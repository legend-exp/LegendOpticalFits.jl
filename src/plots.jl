using Plots

"""
    plot_λ0_by_fiber(x0; multiplicity_thr=0, chmap) -> λ0s_plot

Plot per–fiber no-light probabilities λ₀ from an `x0` event table.

# Arguments
- `x0`: a `Table` of booleans as returned by [`x0_data`](@ref), indicating for each event and channel whether a photon was detected.
- `chmap`: channel map (e.g. `CHANNELMAPS[:p13][:r003]`) providing mapping between SiPM channel names and fiber identifiers.
- `multiplicity_thr`: minimum number of channels with light per event required for inclusion in the computation. Defaults to `0`.

"""
function plot_λ0_by_fiber(
    x0::Table,
    runsel::LegendOpticalFits.RunSelLike;
    multiplicity_thr::Int = 0
)
    period, run = runsel
    chmap = LegendOpticalFits.CHANNELMAPS[period][run]

    # compute λ₀ data
    λ0_d, nsel = λ0_data(x0, multiplicity_thr = multiplicity_thr)
    λdict = Dict(pairs(λ0_d))

    # map to fibers and filter to channels present in chmap
    λ0s_plot = Dict(chmap[ch].fiber => λdict[ch] for ch in keys(λdict) if haskey(chmap, ch))

    # sort by fiber symbol (:IB01 … :OB40)
    sorted_keys = sort(collect(keys(λ0s_plot)))
    values = [λ0s_plot[k] for k in sorted_keys]
    xvals = 1:length(sorted_keys)

    # make the scatter plot
    plt = scatter(
        xvals,
        values;
        xticks = (xvals, string.(sorted_keys)),
        title = "multiplicity ≥ $multiplicity_thr",
        ylabel = "no-light probability λ₀",
        legend = false,
        markersize = 6,
        rotation = 90,
        size = (600, 300)
    )

    return plt
end

export plot_λ0_by_fiber
