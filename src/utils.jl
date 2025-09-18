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
