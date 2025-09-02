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
