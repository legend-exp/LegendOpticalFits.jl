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


function _h5read_revdims_array(fstream::HDF5.File, dset::AbstractString)
    d = fstream[dset]
    A = Array{eltype(d)}(undef, reverse(size(d)))
    permutedims!(A, read(d), ndims(d):-1:1)
    return A
end

function _read_histogram(fstream::LHDataStore, name::AbstractString)
    binning = fstream["$name/binning"]
    isdensity = fstream["$name/isdensity"]
    weights = _h5read_revdims_array(fstream.data_store, "$name/weights")
    return LegendHDF5IO._nt_to_histogram((binning = binning, isdensity = isdensity, weights = weights))
end

function _to_matrix(table::Table; order = nothing)
    _table = table
    if order != nothing
        _table = Table(; (k => getproperty(table, k) for k in order)...)
    end

    return (Tables.matrix(_table), columnnames(_table))
end

function _to_table(matrix::AbstractMatrix, colnames)
    return Table(; (name => matrix[:, j] for (j, name) in enumerate(colnames))...)
end
