"""
    OpticalMap ≡ NamedTuple of 3D histograms keyed by channel Symbol

Handy type alias for LEGEND optical maps, i.e. three-dimensional histograms for
each SiPM channel. The field names are the channel symbols (e.g. :S030).
"""
const OpticalMap{names,T} = NamedTuple{names,T} where {names,T<:Tuple{Vararg{Histogram{<:AbstractFloat,3}}}}
export OpticalMap

"""
    load_optical_map(filename, runsel) -> OpticalMap

Load a LEGEND-200 optical map from file.

# Examples
```julia
load_optical_map("./optmap.lh5", (:p13, :r001))
```
"""
function load_optical_map(filename::AbstractString, runsel::RunSelLike)::OpticalMap
    period, run = runsel
    _detname = id -> rawid2detname(CHANNELMAPS[period][run], id)

    lh5open(filename) do file
        names = filter(s -> occursin(r"_\d{7}$", s), keys(file))
        rawids = [parse(Int, match(r"\d+", name).match) for name in names]
        detnames = map(_detname, rawids)
        order = sortperm(string.(detnames))
        kvs = (detnames[i] => _read_histogram(file, "_$(rawids[i])/p_det") for i in order)
        return (; kvs...)
    end
end

export load_optical_map


"""
    rand_voxel(optmap::OpticalMap; xrange = nothing, yrange = nothing, zrange = nothing) -> (ix, iy, iz)

Sample a random valid voxel (bin indices) from an `OpticalMap`.

The function draws random voxel indices `(ix, iy, iz)` within the histogram domain of the optical map.
The histogram of the first channel is used to determine the geometry (all channels share the same dimensions).

# Arguments
- `optmap`: optical map (see [`load_optical_map`](@ref).
- `xrange`, `yrange`, `zrange`: optional `(min,max)` in axis units.  
  If `nothing` (default), the full axis range is used.

# Returns
Tuple `(ix, iy, iz)` of voxel indices.
"""
function rand_voxel(
    optmap::OpticalMap;
    xrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
    yrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
    zrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing
)::Tuple{Int,Int,Int}
    # use the first histogram to get dimensions + edges
    h = first(values(optmap))
    ex, ey, ez = h.edges
    nx, ny, nz = size(h.weights)

    # helper: convert axis range → index range
    function to_index_range(edges, n, r)
        if r === nothing
            return 1, n
        else
            lo, hi = r
            imin = searchsortedfirst(edges, lo)
            imax = searchsortedlast(edges, hi) - 1
            (1 ≤ imin ≤ imax ≤ n) ||
                error("range $r out of bounds for axis with $n bins")
            return imin, imax
        end
    end

    ixmin, ixmax = to_index_range(ex, nx, xrange)
    iymin, iymax = to_index_range(ey, ny, yrange)
    izmin, izmax = to_index_range(ez, nz, zrange)

    for _ in 1:1000
        ix = rand(ixmin:ixmax)
        iy = rand(iymin:iymax)
        iz = rand(izmin:izmax)
        if h.weights[ix, iy, iz] != -1
            return (ix, iy, iz)
        end
    end

    return error("could not find a valid voxel within 1000 trials")
end

export rand_voxel
