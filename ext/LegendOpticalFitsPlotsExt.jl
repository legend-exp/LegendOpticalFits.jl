module LegendOpticalFitsPlotsExt

import Plots
using TypedTables: Table
using LegendOpticalFits

"""
    plot_λ0_by_fiber(x0, runsel; multiplicity_thr=0) -> plt

Plot per–fiber no–light probabilities λ₀ from an `x0` event table, using a channel map
looked up from `LegendOpticalFits.CHANNELMAPS` via `runsel`.

# Arguments
- `x0`: a `Table` of booleans as returned by [`x0_data`](@ref), indicating for each event and channel whether a photon was detected.
- `runsel`: a `(period, run)` pair or other `RunSelLike` used to select the channel map from `CHANNELMAPS`.
- `multiplicity_thr`: minimum number of channels with light per event required for inclusion in the computation. Defaults to `0`.

Returns a `Plots.Plot` object.
"""
function plot_λ0_by_fiber(
    x0::Table,
    runsel::LegendOpticalFits.RunSelLike;
    multiplicity_thr::Int = 0
)
    period, run = runsel
    chmap = LegendOpticalFits.CHANNELMAPS[period][run]

    # compute λ₀ from x0
    λ0_d, nsel = LegendOpticalFits.λ0_data(x0, multiplicity_thr = multiplicity_thr)
    λdict = Dict(pairs(λ0_d))

    # map to fibers and filter to channels present in chmap
    λ0s_plot = Dict(chmap[ch].fiber => λdict[ch] for ch in keys(λdict) if haskey(chmap, ch))

    # sort by fiber symbol (:IB01 … :OB40)
    sorted_keys = sort(collect(keys(λ0s_plot)))
    if isempty(sorted_keys)
        plt = Plots.plot(title = "no channels found for run $runsel", size = (600, 300))
        return plt
    end
    values = [λ0s_plot[k] for k in sorted_keys]
    xvals = 1:length(sorted_keys)

    # make the scatter plot
    plt = Plots.scatter(
        xvals,
        values;
        xticks = (xvals, string.(sorted_keys)),
        title = "multiplicity ≥ $multiplicity_thr",
        ylabel = "no-light probability λ₀",
        legend = false,
        markersize = 6,
        xrotation = 90,
        size = (600, 300)
    )

    return plt
end

export plot_λ0_by_fiber

"""
    plot_eff_est_from_curves(runsel, curves, x0; multiplicity_thr=0, xlabel="efficiency", ylabel="no-light probability λ₀") -> (plt, eff_est)

Plot the analytical λ₀ vs efficiency curves (from `curves`, e.g. returned by
`LegendOpticalFits.analytical_λ0_curve_all`) for the channels in the channel map
selected by `runsel`, mark the estimated efficiencies where the horizontal
data λ₀ intersects the curves, and return the plot and the `eff_est` mapping.

# Arguments
- `runsel`: `(period, run)` pair selecting the channel map from `LegendOpticalFits.CHANNELMAPS`.
- `curves`: output of `LegendOpticalFits.analytical_λ0_curve_all` (any NamedTuple/Dict-like mapping channel -> (eps, λ0)).
- `x0`: a `Table` of events (as accepted by `LegendOpticalFits.λ0_data`).

# Returns
- `(plt, eff_est)` where `plt` is a `Plots.Plot` and `eff_est` is a `Dict{Symbol, Union{Float64, Nothing}}` mapping channel -> estimated efficiency (or `nothing` when no estimate was found).
"""
function plot_eff_est_from_curves(
    runsel::LegendOpticalFits.RunSelLike,
    curves::AbstractDict{Symbol,Any},
    x0::Table
)
    period, run = runsel
    chmap = LegendOpticalFits.CHANNELMAPS[period][run]

    # curves to Dict-like for indexing
    curves_dict = Dict(pairs(curves))
    curves_plot = Dict(chmap[ch].fiber => curves_dict[ch] for ch in keys(curves_dict) if haskey(chmap, ch))

    # sort by fiber symbol (:IB01 … :OB40)
    sorted_keys = sort(collect(keys(curves_plot)))
    values = [curves_plot[k] for k in sorted_keys]

    # compute data λ0 and estimates
    λ0_d, nsel = LegendOpticalFits.λ0_data(x0; multiplicity_thr = 0)
    eff_est = LegendOpticalFits.estimate_efficiencies_from_curves(curves, λ0_d)

    plt = Plots.plot(xlabel = "efficiency", ylabel = "no-light probability λ₀")
    seenI = false
    seenO = false

    # build reverse map: fiber_symbol -> channel_symbol
    revmap = Dict{Symbol,Symbol}()
    for ch in keys(chmap)
        revmap[chmap[ch].fiber] = ch
    end

    for (i, k) in enumerate(sorted_keys)
        entry = values[i]
        eps = entry.eps
        λvals = entry.λ0
        s = string(k)
        if startswith(s, "I")
            color = :blue
            lab = seenI ? "" : "IB fibers"
            seenI = true
        else
            color = :red
            lab = seenO ? "" : "OB fibers"
            seenO = true
        end
        Plots.plot!(plt, eps, λvals; color = color, label = lab, alpha = 0.7)

        # find data λ0 for this channel (map fiber -> channel symbol)
        if !haskey(revmap, k)
            continue
        end
        ch_sym = revmap[k]
        if !haskey(λ0_d, ch_sym)
            continue
        end
        target = λ0_d[ch_sym]
        eff = get(eff_est, ch_sym, nothing)

        # mark cross on plot if we have an estimate
        if eff !== nothing
            Plots.plot!(
                plt,
                [eff],
                [target];
                seriestype = :scatter,
                markershape = :x,
                markersize = 5,
                color = color,
                label = false
            )
        end
    end

    return plt, eff_est
end

export plot_λ0_by_fiber, plot_eff_est_from_curves

end # module LegendOpticalFitsPlotsExt
