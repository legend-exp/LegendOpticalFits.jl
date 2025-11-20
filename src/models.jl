"""
    λ0_model(efficiencies, log_p0_nominal, x0_random_coin[, multiplicity_thr])

Expected fraction of events in which a SiPM channel sees no light.

Computes the expected fraction by looping over the rows (events) in
`log_p0_nominal`. It samples the expected light/no-light observable by combining
the probability to see no light `p0` with the random coincidence data (at
the same row index). If a multiplicity threshold is specified, events with
lower multiplicity are discarded.

# Arguments
All input data is keyed by detector name (a symbol)

- `efficiencies`: scaling factors for each SiPM channel.
- `log_p0_nominal`: logarithm of the probability to see no light for each
  (event, channel), typically from simulations.
- `x0_random_coin`: presence of light from random coincidences for each
  (event, channel). This is typically coming from a measurement.
- `multiplicity_thr`: discard events with multiplicity below this threshold
  (optional, defaults to 0).
- `n_rands`: average forward model results over this amount of random numbers.

# Returns
- Vector of expectation values for each channel, ordered as the input data
  structures.
"""
function λ0_model(
    efficiencies::NamedTuple,
    log_p0_nominal::Table,
    x0_random_coin::Table,
    ;
    multiplicity_thr::Int = 0,
    n_rands::Int = 10
)::NamedTuple
    # make sure order is consistent with provided scaling_factors
    ϵk = keys(efficiencies);
    ϵv = collect(values(efficiencies))
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵk)
    x0, _ = _to_matrix(x0_random_coin, order = ϵk)

    # get random engine depending on device (CPU/GPU)
    rng = default_device_rng(get_device(log_p0))
    # pre-allocate random numbers for forward model evaluation
    n_events, n_channels = size(log_p0)
    rands = rand(rng, n_events, n_channels, n_rands)

    # call low-level routine
    λ0 = _λ0_model_bulk_ops(ϵv, log_p0, x0, rands, multiplicity_thr = multiplicity_thr)

    # re-label as a NamedTuple keyed by channel symbols
    return NamedTuple{ϵk}(Tuple(λ0))
end

export λ0_model

"""
    _λ0_model_bulk_ops()

Low-level implementation of `λ0_model` using bulk array programming, Reactant /
CUDA compatible.
"""
function _λ0_model_bulk_ops(
    efficiencies::AbstractVector{<:Number}, # Should be Real, but Reactant tracing array elements are subtypes of Number
    log_p0_nominal::AbstractMatrix{<:Number},
    x0_random_coin::AbstractMatrix{<:Number},
    rands::AbstractArray{<:Number,3},
    ;
    multiplicity_thr::Int = 0
)
    T = eltype(efficiencies)

    n_events, n_channels = size(log_p0_nominal)

    # avoid numerical issues
    ϵ = clamp.(efficiencies, 1e-10, 1 - 1e-10)

    # calculate the probability to see no light
    # NOTE: exp is expensive, do it here
    p0 = exp.(log_p0_nominal .* ϵ')

    # draw bernoulli distributed numbers and fold in random coincidences
    drawn = (rands .< p0) .* x0_random_coin

    # compute the event multiplicity
    multiplicity = n_channels .- sum(drawn, dims = 2)

    # and select above a threshold
    weights = one(T) .* (multiplicity .>= multiplicity_thr)

    # calculate expectation for fraction of events with no light in each channel
    return vec(sum(drawn .* weights, dims = (1, 3))) / sum(weights)
end

"""
    _λ0_model_loops()

Low-level implementation of `λ0_model` using for loops.
"""
function _λ0_model_loops(
    efficiencies::AbstractVector{<:AbstractFloat},
    log_p0_nominal::AbstractMatrix{<:AbstractFloat},
    x0_random_coin::AbstractMatrix{Bool},
    ;
    multiplicity_thr::Int = 0
)::Vector{Float64}
    n_events, n_channels = size(log_p0_nominal)

    # avoid numerical issues
    ϵ = clamp.(efficiencies, 1e-10, 1 - 1e-10)

    # no-light probability
    p0 = zeros(n_channels)
    # outcome of bernoulli trials Ber(p0)
    x0 = zeros(Bool, n_channels)
    # number of events with no light, per channel
    N0 = zeros(n_channels)
    # number of events that pass the multiplicity condition
    N = 0

    @inbounds @simd for i in 1:n_events
        # set multiplicity to maximum possible
        multiplicity = n_channels

        @inbounds @simd for k in 1:n_channels
            # compute the no-light probability, folding efficiencies in
            p0[k] = exp(log_p0_nominal[i, k] * ϵ[k])
            # sample the outcome for this event from bernoulli with probability p0
            # and optionally the random coincidences from data
            x0[k] = (rand() < p0[k]) && x0_random_coin[i, k]

            # update the event multiplicity
            multiplicity -= x0[k]
        end

        # finally check if the multiplicity condition is satisfied, if yes use the event
        if multiplicity >= multiplicity_thr
            @inbounds @simd for k in 1:n_channels
                N0[k] += x0[k]
            end
            N += 1
        end
    end

    return N0 ./ N
end

"""
    x0_mock(efficiencies_true, log_p0_nominal, x0_random_coin)

Generate mock x0 data (bernoulli draws folded with random-coincidence mask)
from `efficiencies_true` and `log_p0_nominal` from model.
"""
function x0_mock(
    efficiencies_true::NamedTuple,
    log_p0_nominal::Table,
    x0_random_coin::Table
)
    # make sure order is consistent with provided scaling_factors
    ϵk = keys(efficiencies_true)
    ϵv = collect(values(efficiencies_true))
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵk)
    x0_rc, _ = _to_matrix(x0_random_coin, order = ϵk)

    # get random engine depending on device (CPU/GPU)
    rng = default_device_rng(get_device(log_p0))
    # pre-allocate random numbers for forward model evaluation (2D)
    n_events, n_channels = size(log_p0)
    rands = rand(rng, n_events, n_channels)

    # avoid numerical issues
    ϵ = clamp.(ϵv, 1e-10, 1 - 1e-10)

    # calculate the probability to see no light (duplicate of bulk logic)
    p0 = exp.(log_p0 .* ϵ')

    # draw bernoulli distributed numbers and fold in random coincidences
    drawn = (rands .< p0) .* x0_rc

    # return as a Table 
    return Table(; (ϵk[i] => drawn[:, i] for i in 1:size(drawn, 2))...)
end
