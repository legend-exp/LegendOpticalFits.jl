"""
    expected_no_light_fraction(scaling_factors, log_p0_nominal, x0_random_coin[, multiplicity_thr])

Expected fraction of events in which a SiPM channel sees no light.

Computes the expected fraction by looping over the rows (events) in
`log_p0_nominal`. It samples the expected light/no-light observable by combining
the probability to see no light `p0` with the random coincidence data (at
the same row index). If a multiplicity threshold is specified, events with
lower multiplicity are discarded.

# Arguments
- `scaling_factors`: scaling factors for each SiPM channel.
- `log_p0_nominal`: logarithm of the probability to see no light for each
  (event, channel), typically from simulations.
- `x0_random_coin`: presence of light from random coincidences for each
  (event, channel). This is typically coming from a measurement.
- `multiplicity_thr`: discard events with multiplicity below this threshold
  (optional, defaults to 0).

# Returns
- Vector of expectation values for each channel, ordered as the input data
  structures.
"""
function expected_no_light_fraction(
    scaling_factors::AbstractVector{<:Float64},
    log_p0_nominal::AbstractMatrix{<:Float64},
    x0_random_coin::AbstractMatrix{Bool},
    ;
    multiplicity_thr::Int = 0
)
    n_events, n_channels = size(log_p0_nominal)

    # avoid numerical issues
    ϵ = clamp.(scaling_factors, 1e-10, 1 - 1e-10)

    # no-light probability
    p0 = zeros(n_channels)
    # outcome of bernoulli trials Ber(p0)
    x0 = zeros(Bool, n_channels)
    # number of events with no light, per channel
    λ0 = zeros(n_channels)
    # number of events that pass the multiplicity condition
    N0 = 0

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
                λ0[k] += x0[k]
            end
            N0 += 1
        end
    end

    return λ0 ./ N0
end

export expected_no_light_fraction


#= version with bulk / array programming operations, FYI
function _version_bulk_ops(
    efficiencies::AbstractVector{<:Float64}, 
    log_p0_nom::Matrix{Float64}, 
    multiplicity_thr::Int,
)
    n_events, n_channels = size(log_p0_nom)

    # avoid numerical issues
    ϵ = clamp!(copy(efficiencies), 1e-10, 1 - 1e-10)

    # calculate the probability to see no light
    # NOTE: exp is expensive, do it here
    p0 = exp.(log_p0_nom .* ϵ')

    # draw bernoulli distributed numbers
    drawn = rand(n_events, n_channels) .< p0

    # compute the event multiplicity
    multiplicity = vec(n_channels .- sum(drawn, dims=2))

    # and select above a threshold
    mask = multiplicity .>= multiplicity_thr

    # calculate expectation for fraction of events with no light in each channel
    return sum(drawn[mask, :], dims=1) / count(mask)
end
=#

function simulate_ar39_log_p0_nominal(optmap, n_events; light_yield = 60)
    # load Ar39 beta spectrum
    dist_ar39 = load_ar39_spectrum("./ar39.csv")

    n_channels = length(collect(keys(optmap)))

    # get n_events beta energies
    sampled_energies = rand(dist_ar39, n_events)
    # get expected number of photons
    mean_ar39_photons = rand.(Poisson.(sampled_energies .* light_yield))

    p0_nom = zeros(n_events, n_channels)
    for event_idx in 1:n_events
        # get valid lar voxel
        point = sample_valid_point(optmap[1])
        # get how many scintillation photons for this event
        n = mean_ar39_photons[event_idx]

        # store map values at point for all channels
        for ch_idx in 1:n_channels
            ξ = optmap[ch_idx][point...]
            # get expected number of detected photons with efficiency 1
            p0_nom[event_idx, ch_idx] = -n * ξ
        end
    end

    return p0_nom
end
