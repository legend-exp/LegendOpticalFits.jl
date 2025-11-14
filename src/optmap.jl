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
function load_optical_map(
    filename::AbstractString,
    runsel::RunSelLike;
    exclude_unusable::Bool = false
)::OpticalMap
    period, run = runsel
    chmap = CHANNELMAPS[period][run]
    _detname = id -> rawid2detname(chmap, id)

    lh5open(filename) do file
        names = filter(s -> occursin(r"_\d{7}$", s), keys(file))
        rawids = [parse(Int, match(r"\d+", name).match) for name in names]
        raw_det = [(rawids[i], _detname(rawids[i])) for i in eachindex(rawids)]

        # optionally exclude unusable channels
        if exclude_unusable
            raw_det = filter(pair -> chmap[pair[2]].usable == true, raw_det)
        end

        order = sortperm(string.(last.(raw_det)))
        kvs = ( last(raw_det[i]) => LegendOpticalFits._read_histogram(file, "_$(first(raw_det[i]))/p_det") for i in order )

        return (; kvs...)
    end
end

export load_optical_map

"""
    detection_prob(h, coords...)

Return the bin content of histogram `h` at the given coordinates
`coords` (with units). Coordinates must be inside the histogram
bounds, otherwise an error is thrown.
"""
function detection_prob(h, coords...; out_of_bounds_val=nothing)
    point = ustrip.(u"m", coords)
    idx   = map(searchsortedlast, h.edges, point)

    if any(((i, e),) -> i < 1 || i ≥ length(e), zip(idx, h.edges))
        if out_of_bounds_val === nothing
            error("point $point out of bounds")
        else
            return convert(eltype(h.weights), out_of_bounds_val)
        end
    end

    return h.weights[idx...]
end

"""
    detection_prob_bcast(h, xss, yss, zss)

Apply [`detection_prob_event`] to each triple of vectors
`(xs, ys, zs)` from `(xss, yss, zss)`.

This is useful when `xss`, `yss`, and `zss` are
`VectorOfVectors`, e.g. the coordinates of all hits for many events.
Returns a vector of vectors of bin contents.
"""
function detection_prob_vov(h, xss, yss, zss; out_of_bounds_val=nothing)
    return map((xs, ys, zs) -> detection_prob.(Ref(h), xs, ys, zs; out_of_bounds_val=out_of_bounds_val), xss, yss, zss)
end

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
