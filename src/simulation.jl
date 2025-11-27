"""
    log_p0_nominal(sim_data, optical_map; ...)

Logarithm of no-light-probability for events simulated by remage.

This can be used to compute inputs for [`λ0_model`](@ref).

# Arguments
- `sim_data`: output table (stepping data) for a scintillation detector (the
  same the optical map refers to) from remage. Must contain the fields `xloc`,
  `yloc`, `zloc` and `edep`.
- `optmap`: liquid argon optical map as loaded by [`load_optical_map`](@ref).
- `n_events`: number of Ar-39 decays to simulate.
- `light_yield`: liquid argon scintillation yield.

# Returns
A table of `log(p0)` for each channel (columns) and event (rows).

# Examples
```julia
# load an optical map
optmap = load_optical_map("map.lh5", (:p13, :r001))

# load some remage simulation data
sim_data = lh5open("th228.lh5") do h5
    return h5["stp/lar"][:]
end

# call the function
log_p0_nominal(sim_data, optmap)
```
"""
function log_p0_nominal(
    sim_data::Table,
    optmap::OpticalMap,
    ;
    light_yield::Quantity = 51u"1/keV",
    out_of_bounds_val = nothing,
    map_undefined_val = nothing
)::Table
    # handy references
    coords = (sim_data.xloc, sim_data.yloc, sim_data.zloc)
    edeps = sim_data.edep

    # output values
    lp0 = Dict{Symbol,Vector{Float64}}()

    # number of scintillation photons
    n = map(edep -> rand.(Poisson.(edep .* light_yield)), edeps)

    for ch in keys(optmap)
        # detection probability
        ξ = detection_prob_vov(optmap[ch], coords...; out_of_bounds_val = out_of_bounds_val)
        # throw error if coordinates in voxel where map is not defined (=-1) or set ξ at that coordinate to unser specific value
        if any(x -> x == -1, Iterators.flatten(ξ))
            if map_undefined_val === nothing
                error("coordinate(s) with undefined map")
            else
                ξ = [replace(v, -1 => map_undefined_val) for v in ξ]
            end
        end
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
- `light_yield`: liquid argon scintillation yield.
- `rand_voxel_kwargs...`: optional keyword arguments forwarded to
  [`rand_voxel`](@ref).

# Returns
A table of `log(p0)` for each channel (columns) and event (rows).
"""
function log_p0_nominal_ar39(
    optmap::OpticalMap,
    n_events::Integer
    ;
    light_yield::Quantity = 51u"1/keV",
    rand_voxel_kwargs...
)::Table

    # load Ar39 beta spectrum
    dist_ar39 = ar39_beta_energy_dist()

    # sample beta energies (in keV)
    sampled_energies = rand(dist_ar39, n_events)

    # use the first histogram to get edges for voxel centers
    h = first(values(optmap))
    ex, ey, ez = h.edges

    # helper to compute bin center from edges and index
    center(e, i) = (e[i] + e[i + 1]) / 2

    # build per-event single-step coordinates and energy deposits
    xloc = Vector{Vector{typeof(1.0u"m")}}(undef, n_events)
    yloc = Vector{Vector{typeof(1.0u"m")}}(undef, n_events)
    zloc = Vector{Vector{typeof(1.0u"m")}}(undef, n_events)
    edep = Vector{Vector{typeof(1.0u"keV")}}(undef, n_events)

    for i in 1:n_events
        ix, iy, iz = rand_voxel(optmap, rand_voxel_kwargs...)
        xloc[i] = [center(ex, ix) * u"m"]
        yloc[i] = [center(ey, iy) * u"m"]
        zloc[i] = [center(ez, iz) * u"m"]
        edep[i] = [sampled_energies[i] * u"keV"]
    end

    sim_tbl = Table(xloc = xloc, yloc = yloc, zloc = zloc, edep = edep)

    return log_p0_nominal(sim_tbl, optmap; light_yield = light_yield, rand_voxel_kwargs...)
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
