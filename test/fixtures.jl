using LegendOpticalFits

using Tables
using StatsBase
using DensityInterface

function mock_optmap(runsel)
    # optmap = load_optical_map("map-merged-2cm-p13.lh5", runsel)
    edges = (collect(-0.7:0.02:0.7), collect(-0.7:0.02:0.7), collect(-1.89:0.02:1.89))
    weights = fill(2e-5, (70, 70, 189))
    h = StatsBase.Histogram(edges, weights)
    channels = keys(LegendOpticalFits.CHANNELMAPS[runsel[1]][runsel[2]])
    return (; (ch => h for ch in channels)...)
end

function mock_ar39_data(optmap, n_events; eff = 0.5)
    # some fake x0 Ar39 data
    log_p0 = log_p0_nominal_ar39(optmap, n_events)
    shape = size(Tables.matrix(log_p0))

    # some random efficiencies
    if eff isa Real
        eff = (; (k => eff for k in keys(optmap))...)
    else
        eff = (; (k => v for (k, v) in zip(keys(optmap), eff))...)
    end

    log_p0_m, _ = LegendOpticalFits._to_matrix(log_p0)
    x0_m = rand(shape...) .< exp.(log_p0_m .* collect(values(eff))')

    # some random coincidences
    x0_rc_m = rand(shape...) .< 0.1

    x0 = LegendOpticalFits._to_table(x0_m, keys(eff))
    x0_rc = LegendOpticalFits._to_table(x0_rc_m, keys(eff))

    return x0, x0_rc, log_p0
end
