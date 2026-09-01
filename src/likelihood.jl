"""
    make_λ0_likelihood(x0, log_p0_nominal, x0_random_coin; multiplicity_thr=0, n_rands=10, device=CPUDevice()) -> DensityFunction

Construct the likelihood of no-light fractions per channel.

We model the fraction of events with no detected light (`λ0`) as follows:

- For each channel, the expected no-light probability `λ0_model` comes from
  the simulation (`log_p0_nominal`) combined with random coincidences
  (`x0_random_coin`) and scaled by per-channel efficiencies
  (the parameters of the model).

- The observed no-light fraction in data is `λ0_data = N0 / N_data`, where `N0`
  is the number of no-light events passing a multiplicity threshold.

- Since `N_data` is large, the binomial distribution `N0 ~ Binomial(N_data,
  λ0_model)` can be approximated by a normal distribution: `λ0_data ~ Normal(μ =
  λ0_model, σ² = λ0_model (1 - λ0_model) / N_data)`.

The likelihood is the sum of log-probabilities across all channels.

# Arguments
- `x0`: observed no-light indicators from data events.
- `log_p0_nominal`: logarithm of the probability to see no light for each
  (event, channel), typically from simulations.
- `x0_random_coin`: observed no-light indicators from random coincidence events.
- `multiplicity_thr`: discard events with multiplicity below this threshold
  (optional, defaults to 0).
- `n_rands`: average forward model results over this amount of random numbers.
- `device`: on which device to run the computation of the forward model. (default
  `CPUDevice()`)

# Examples

Get some data:

```julia
using LegendOpticalFits

runsel = (:p13, :r001)
nev_sim = 10_000
nev_data = 1_000
multiplicity_thr = 6

optmap = load_optical_map("./optmap-p13.lh5", runsel)
log_p0 = log_p0_nominal_ar39(optmap, nev_sim)

x0 = x0_data("l200-p13-r001-ath-tier_evt.lh5", runsel, max_events=nev_data)
x0_rc = x0_data("l200-p13-r001-ant-tier_evt.lh5", runsel, max_events=nev_sim)
```

Build the likelihood (on the CPU by default):

```julia
logl = make_λ0_likelihood(x0, lp0, x0_rc, multiplicity_thr=multiplicity_thr)
```

CUDA via Reactant/XLA:

```julia
using Reactant
Reactant.set_default_backend("cuda")

logl = make_λ0_likelihood(x0, lp0, x0_rc, multiplicity_thr=multiplicity_thr, device=ReactantDevice())
```

# Returns
- A `DensityFunction` object representing the log-likelihood. It can be called
  with a parameter set.
"""
function make_λ0_likelihood(
    x0::Table,
    log_p0_nominal::Table,
    x0_random_coin::Table,
    ;
    multiplicity_thr::Int = 0,
    n_rands::Int = 10,
    device = CPUDevice()
)
    # we choose as channel order the one used in x0
    ϵ_order = columnnames(x0)

    # convert to matrix with the correct order
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵ_order)
    x0_rc, _ = _to_matrix(x0_random_coin, order = ϵ_order)

    # cast x0 to floating point for later computation
    float_x0_rc = eltype(log_p0).(x0_rc)

    # pre-allocate random numbers for forward model evaluation
    n_events, n_channels = size(log_p0)
    rands = rand(n_events, n_channels, n_rands)

    # prepare data
    λ0, N_ev = λ0_data(x0, multiplicity_thr = multiplicity_thr)
    data = [λ0[k] for k in ϵ_order]

    # prepare model for computation on the requested device (CPU or GPU)
    _model(ϵ) = _λ0_model_bulk_ops(ϵ, log_p0, float_x0_rc, rands, multiplicity_thr = multiplicity_thr)
    _model_on_dev = on_device(_model, device, rand(eltype(log_p0), n_channels))

    function _logl(ϵ)
        # ϵ is expected to be a NamedTuple, we just pass the values to
        # the low level routines. make sure the order of the parameters is
        # the correct one
        ϵv = [ϵ[k] for k in ϵ_order]

        # compute the forward model
        model = _model_on_dev(ϵv)

        # and the log-likelihood
        x = data
        μ = model
        # Gaussian approximation of the binomial distribution
        σ = sqrt.((μ .- μ .^ 2) / N_ev)

        return sum(logpdf.(Normal.(μ, σ), x))
    end

    return DensityInterface.logfuncdensity(_logl)
end

export make_λ0_likelihood

"""
    make_λ0_joint_likelihood(datasets; n_rands=10, device=CPUDevice()) -> DensityFunction

Combined no-light-fraction likelihood across several datasets that share a SINGLE
set of per-channel efficiencies `ε`.

Each dataset is scored with the same per-channel, Gaussian-approximated binomial
term as [`make_λ0_likelihood`](@ref), but restricted to its OWN multiplicity
threshold `M ≥ mult_thr`, and the per-dataset log-likelihoods are summed.
Because each dataset contributes its own selected-event count `N` (which sets the
Gaussian width `σ² = μ(1-μ)/N`), datasets with equal `N` carry equal statistical
weight — the basis of a "balanced" design (size each dataset to the same `N`).

This is what allows a joint Ar-39 + calibration-source fit with one shared `ε`:
e.g. Ar-39 (M ≥ 6) plus Cs and Th source selections (higher thresholds), each a
separate entry in `datasets`.

# Arguments
- `datasets`: a vector of NamedTuples, each with fields
    - `x0`      : data no-light indicator `Table` (events × channels);
    - `log_p0`  : simulation `log(p0)` `Table` (events × channels);
    - `x0_rc`   : random-coincidence no-light `Table`;
    - `mult_thr`: multiplicity threshold (keep events with `M ≥ mult_thr`).
  All datasets must share the same channel columns; the channel order (and hence
  the expected parameter names) is taken from `datasets[1].x0`.
- `n_rands`, `device`: as in [`make_λ0_likelihood`](@ref).

# Returns
- A `DensityFunction` of a single `NamedTuple` of per-channel efficiencies.
"""
function make_λ0_joint_likelihood(
    datasets::AbstractVector;
    n_rands::Int = 10,
    device = CPUDevice()
)
    # shared channel order, taken from the first dataset's data table
    ϵ_order = columnnames(datasets[1].x0)

    # pre-build one forward model + one data vector per dataset
    compiled = Any[]
    for ds in datasets
        thr = ds.mult_thr

        log_p0, _ = _to_matrix(ds.log_p0, order = ϵ_order)
        x0_rc, _ = _to_matrix(ds.x0_rc, order = ϵ_order)
        float_x0_rc = eltype(log_p0).(x0_rc)

        n_events, n_channels = size(log_p0)
        rands = rand(eltype(log_p0), n_events, n_channels, n_rands)

        # data no-light fractions above the threshold and the selected-event count N
        λ0, N_ev = λ0_data(ds.x0, multiplicity_thr = thr)
        data = [λ0[k] for k in ϵ_order]

        _model(ϵ) = _λ0_model_bulk_ops(ϵ, log_p0, float_x0_rc, rands, multiplicity_thr = thr)
        _model_on_dev = on_device(_model, device, rand(eltype(log_p0), n_channels))

        @info "joint dataset: nev=$n_events M>=$thr N_data=$N_ev"
        push!(compiled, (_model_on_dev, data, N_ev))
    end

    function _logl(ϵ)
        # ϵ is a NamedTuple; pass values in the shared channel order
        ϵv = [ϵ[k] for k in ϵ_order]

        ll = 0.0
        for (_model_on_dev, data, N_ev) in compiled
            μ = _model_on_dev(ϵv)
            # Gaussian approximation of the binomial
            σ = sqrt.((μ .- μ .^ 2) / N_ev)
            ll += sum(logpdf.(Normal.(μ, σ), data))
        end
        return ll
    end

    return DensityInterface.logfuncdensity(_logl)
end

export make_λ0_joint_likelihood
