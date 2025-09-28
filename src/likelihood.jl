"""
    make_λ0_likelihood(x0, log_p0_nominal, x0_random_coin; multiplicity_thr=0, n_rands=10, smear_factor=0, device=CPUDevice()) -> DensityFunction

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
- `smear_factor`: the width of the likelihood gaussian terms is increased by a
  factor `smear_factor * mean`.
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
    smear_factor::Real = 0,
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
        σ = sqrt.((μ .- μ .^ 2) / N_ev) .+ smear_factor .* μ

        logl = sum(- (x .- μ) .^ 2 ./ (σ .^ 2))

        return logl
    end

    return DensityInterface.logfuncdensity(_logl)
end

export make_λ0_likelihood
