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
    analytical_λ0_curve(log_p0_nominal, x0_random_coin, channel; eps_range=range(0.0,1.0,length=101))

Compute the per-channel no-light probability λ₀ as a function of a single
channel's efficiency. 

This analytical model computes the expectation per event 
λ₀_jk = mean(p0_jk * x0_jk) with p0_jk = exp(log_p0_jk * ε) and x0_jk 
the random coincidence indicators.

# Arguments
- `log_p0_nominal::Table`: table of log p0 values (events × channels).
- `x0_random_coin::Table`: table of random-coincidence indicators (events × channels).
- `channel::Symbol`: Symbol identifying the channel (must be present in
  columnnames(log_p0_nominal)).
- `eps_range`: efficiency values to evaluate (vector or range).

# Returns
- `NamedTuple` with fields :eps (vector of efficiencies) and :λ0
  (vector of no-light probabilities evaluated at those efficiencies).

# Notes
  Can only be used if no multiplicity threshold is applied!
"""
function analytical_λ0_curve(
    log_p0_nominal::Table,
    x0_random_coin::Table,
    channel::Symbol;
    eps_range = range(0.0, 1.0, length = 101)
)
    # derive channel ordering from the provided `log_p0_nominal` table
    order = collect(columnnames(log_p0_nominal))
    log_p0, _ = _to_matrix(log_p0_nominal, order = order)
    x0, _ = _to_matrix(x0_random_coin, order = order)

    # resolve channel index (channel must be a Symbol and present in log_p0_nominal)
    ch_idx = findfirst(==(channel), order)
    if ch_idx === nothing
        throw(ArgumentError("channel $channel not found in efficiencies keys"))
    end

    # columns for the requested channel
    log_p0_col = view(log_p0, :, ch_idx)
    x0_col = view(x0, :, ch_idx)

    eps_vec = collect(eps_range)
    n = length(eps_vec)
    λ0s = Vector{Float64}(undef, n)
    for (i, ε) in enumerate(eps_vec)
        p0 = exp.(log_p0_col .* ε)
        # analytic exact expectation
        λ0s[i] = mean(p0 .* x0_col)
    end

    return (eps = eps_vec, λ0 = λ0s)
end

"""
    analytical_λ0_curve_all(log_p0_nominal, x0_random_coin; eps_range=range(0.0,1.0,length=101))

Compute λ₀(eps) curves for every channel present in `log_p0_nominal`.

Description
- For each channel (columns of `log_p0_nominal`) evaluate the no-light
  probability λ₀ across `eps_range` by calling `compute_λ0_curve`.
- The returned curves are suitable for locating the efficiency at which a
  measured λ₀ would be reproduced (e.g. by interpolation).

Arguments
- see above in description of `compute_λ0_curve`

Returns
- `Dict{Symbol, NamedTuple}` mapping each channel symbol to the NamedTuple
  returned by `compute_λ0_curve` (fields `:eps` and `:λ0`).

Notes
- Channel order is derived from `columnnames(log_p0_nominal)`. Use that
  table to control which channels are computed and their ordering.

"""
function analytical_λ0_curve_all(
    log_p0_nominal::Table,
    x0_random_coin::Table;
    eps_range = range(0.0, 1.0, length = 101)
)
    # derive channel list from the table column names
    ch_list = collect(columnnames(log_p0_nominal))
    results = Dict{Symbol,Any}()
    for ch in ch_list
        results[ch] = analytical_λ0_curve(
            log_p0_nominal,
            x0_random_coin,
            ch;
            eps_range = eps_range
        )
    end

    return results
end
