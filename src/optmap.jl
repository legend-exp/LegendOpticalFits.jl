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
