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

# Returns
- Vector of expectation values for each channel, ordered as the input data
  structures.
"""
function λ0_model(
    efficiencies::NamedTuple,
    log_p0_nominal::Table,
    x0_random_coin::Table,
    ;
    multiplicity_thr::Int = 0
)::NamedTuple
    # make sure order is consistent with provided scaling_factors
    ϵk = keys(efficiencies);
    ϵv = collect(values(efficiencies))
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵk)
    x0, _ = _to_matrix(x0_random_coin, order = ϵk)

    # call low-level routine
    λ0 = λ0_model(ϵv, log_p0, x0, multiplicity_thr = multiplicity_thr)

    # re-label as a NamedTuple keyed by channel symbols
    return NamedTuple{ϵk}(Tuple(λ0))
end

"""
    λ0_model(
        efficiencies::AbstractVector{<:AbstractFloat},
        log_p0_nominal::AbstractVector{<:AbstractFloat},
        mask::AbstractMatrix{Bool};
        multiplicity_thr::Int=1)

Low-level version of `λ0_model` working with plain arrays and boolean mask.
"""
function λ0_model(
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

export λ0_model

"""
    λ0_model_bulk_ops(...)

Low-level version of [`λ0_model`](@ref) that uses bulk array programming.
To be used as a CUDA kernel.
"""
function λ0_model_bulk_ops(
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
    log_p0_nominal_ar39(optical_map, n_events; ...)

Logarithm of no-light-probability for simulated Ar-39 events.

Uses the energy distribution of the beta particle emitted in the decay of an
Ar-39 nucleus and the detection probability in liquid argon (through the
optical map) to estimate the probability of seeing no light for a set of random
events. This can be used to compute inputs for [`λ0_model`](@ref).

# Arguments
- `optmap`: liquid argon optical map as loaded by [`load_optical_map`](@ref).
- `n_events`: number of Ar-39 decays to simulate.
- `light_yield`: liquid argon scintillation yield in units of photons / keV.
- `rand_voxel_kwargs...`: optional keyword arguments forwarded to
  [`rand_voxel`](@ref).

# Returns
A table of `log(p0)` for each channel (columns) and event (rows).
"""
function log_p0_nominal_ar39(
    optmap::OpticalMap,
    n_events::Integer
    ;
    light_yield::Integer = 60,
    rand_voxel_kwargs...
)::Table

    # load Ar39 beta spectrum
    dist_ar39 = ar39_beta_energy_dist()

    ch_keys = collect(keys(optmap))

    # get n_events beta energies
    sampled_energies = rand(dist_ar39, n_events)
    # get expected number of photons
    mean_ar39_photons = rand.(Poisson.(sampled_energies .* light_yield))

    # initialize dictionary of zeros for each channel
    p0_nom = Dict(ch => zeros(n_events) for ch in ch_keys)

    for event_idx in 1:n_events
        # get valid lar voxel
        point = rand_voxel(optmap, rand_voxel_kwargs...)

        # get how many scintillation photons for this event
        n = mean_ar39_photons[event_idx]

        # store map values at point for all channels
        for ch in ch_keys
            # map probability at selected voxel
            ξ = getproperty(optmap, ch).weights[point...]
            # get expected number of detected photons with efficiency 1
            p0_nom[ch][event_idx] = -n * ξ
        end
    end

    return Table(p0_nom)
end

export log_p0_nominal_ar39

"""
    ar39_beta_energy_dist() -> MixtureModel{Uniform}

Energy distribution of the beta particle emitted in an Ar-39 nuclear decay.

Return a continuous probability distribution for the beta decay spectrum of
Ar-39. The distribution is constructed from tabulated values from the IAEA
BetaShape database, downloadable at this
[link](https://www-nds.iaea.org/relnsd/v1/data?fields=bin_beta&nuclides=39ar&rad_types=bm).
"""
function ar39_beta_energy_dist()
    csvpath = joinpath(@__DIR__, "..", "data", "ar39-beta-decay-data.csv")

    tbl   = CSV.File(csvpath)
    edges = collect(tbl.bin_en)
    dens  = collect(tbl.dn_de)

    d = diff(edges)
    ΔE = vcat(d, d[end])
    w = dens .* ΔE
    w ./= sum(w)

    comps = [Uniform(edges[i], edges[i] + ΔE[i]) for i in eachindex(edges)]
    return MixtureModel(comps, w)
end

export ar39_beta_energy_dist
