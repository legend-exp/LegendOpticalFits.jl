"""
    load_optical_map(filename, runsel) -> Dict{String, Histogram}

Load LEGEND-200 optical maps from file.

Returns a mapping of detector names to histograms for each detector.

# Examples
```julia
load_optical_maps("./optmaps.lh5", (:p13, :r001))
```
"""
function load_optical_maps(filename::AbstractString, runsel::RunSelLike)::Dict{Symbol,StatsBase.Histogram}
    period, run = runsel
    _detname = id -> rawid2detname(CHANNELMAPS[period][run], id)

    lh5open(filename) do file
        names = filter(s -> occursin(r"_\d{7}$", s), keys(file))
        rawids = [parse(Int, match(r"\d+", name).match) for name in names]
        return Dict(_detname(id) => _read_histogram(file, "_$id/p_det") for id in rawids)
    end
end

export load_optical_maps


function sample_valid_point(
    optmap::Dict{Symbol,<:StatsBase.Histogram};
    zlims::Tuple{<:Integer,<:Integer} = (20, 180)
)::Tuple{Int,Int,Int}

    # histogram for the first channel, as all have the same dimension
    h = first(values(optmap))
    xdim, ydim, zdim = size(h.weights)

    zmin, zmax = zlims
    if zmax > zdim
        error("zmax=$(zmax) exceeds available z-dimension (zdim=$(zdim))")
    end

    trials = 0
    while true
        trials += 1
        x = rand(1:xdim)
        y = rand(1:ydim)
        z = rand(zmin:zmax)
        hprob = h.weights[x, y, z]
        if hprob != -1
            return (x, y, z)
        end
        # add a break condition to avoid infinite loop
        if trials > 100
            error("too many trials with invalid hprob (-1)")
        end
    end
end

export sample_valid_point
