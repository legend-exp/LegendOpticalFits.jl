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
    rand_voxel(optmap::OpticalMap; zlims = (20, 180)) -> (x, y, z)

Sample a random valid point from an optical map.

The function draws random `(x, y, z)` coordinates within the histogram domain
of the optical map the histogram of the first channel is used to determine the
geometry (all channels share the same dimensions).

# Arguments
- `optmap`: optical map (see [`load_optical_map`](@ref).
- `zlims`: restrict the sampled z–range. default is `(20, 180)`.
"""
function rand_voxel(
    optmap::OpticalMap
    ;
    zlims::Tuple{<:Integer,<:Integer} = (20, 180)
)::Tuple{Int,Int,Int}
    # histogram for the first channel, as all have the same dimension
    h = first(Tuple(optmap))
    xdim, ydim, zdim = size(h.weights)

    zmin, zmax = zlims
    if zmax > zdim
        error("zmax=$(zmax) exceeds available z-dimension (zdim=$(zdim))")
    end

    for _ in 1:1000
        x = rand(1:xdim)
        y = rand(1:ydim)
        z = rand(zmin:zmax)

        if h.weights[x, y, z] != -1
            return (x, y, z)
        end
    end

    error("could not find a valid voxel with 1000 trials")
    return nothing
end

export rand_voxel
