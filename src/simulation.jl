function log_p0_nominal(
    sim_data::Table,
    optmap::OpticalMap,
    ;
    light_yield::Integer = 40,
    rand_voxel_kwargs...
)::Table
    # handy references
    coords = (sim_data.xloc, sim_data.yloc, sim_data.zloc)
    edeps = sim_data.edep

    # output values
    lp0 = Dict{Symbol,Vector{Float64}}()

    # number of scintillation photons
    n = map(edep -> rand.(Poisson.(ustrip.(u"keV", edep) .* light_yield)), edeps)

    for ch in keys(optmap)
        # detection probability
        ξ = detection_prob_vov(optmap[ch], coords...)
        # log p0 for each step
        step_lp0 = map((_ξ, _n) -> -_ξ .* _n, ξ, n)
        # now sum over steps
        lp0[ch] = sum.(step_lp0)
    end

    return Table(lp0)
end

export log_p0_nominal

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
    light_yield::Integer = 40,
    rand_voxel_kwargs...
)::Table

    # load Ar39 beta spectrum
    dist_ar39 = ar39_beta_energy_dist()

    ch_keys = collect(keys(optmap))

    # get n_events beta energies
    sampled_energies = rand(dist_ar39, n_events)
    # sample number of photons
    ar39_photons = rand.(Poisson.(sampled_energies .* light_yield))

    # initialize dictionary of zeros for each channel
    p0_nom = Dict(ch => zeros(n_events) for ch in ch_keys)

    for event_idx in 1:n_events
        # get valid lar voxel
        point = rand_voxel(optmap, rand_voxel_kwargs...)

        # get how many scintillation photons for this event
        n = ar39_photons[event_idx]

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
