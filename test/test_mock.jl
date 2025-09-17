using LegendOpticalFits

using Tables
using StatsBase
using DensityInterface
using BenchmarkTools

n_events = 100_000
runsel = (:p13, :r001)

function make_mock_optmap()
    # optmap = load_optical_map("map-merged-2cm-p13.lh5", runsel)
    edges = (collect(-0.7:0.02:0.7), collect(-0.7:0.02:0.7), collect(-1.89:0.02:1.89))
    weights = fill(2e-5, (70, 70, 189))
    h = StatsBase.Histogram(edges, weights)
    channels = keys(LegendOpticalFits.CHANNELMAPS[runsel[1]][runsel[2]])
    return (; (ch => h for ch in channels)...)
end

# some fake x0 Ar39 data
optmap = make_mock_optmap()
log_p0 = log_p0_nominal_ar39(optmap, n_events)
shape = size(Tables.matrix(log_p0))

# some random efficiencies
eff = (; (k => 0.5 for k in keys(optmap))...)

log_p0_m, _ = LegendOpticalFits._to_matrix(log_p0)
x0_m = rand(shape...) .< exp.(log_p0_m .* collect(values(eff))')

# some random coincidences
x0_rc_m = rand(shape...) .< 0.1

x0 = LegendOpticalFits._to_table(x0_m, keys(eff))
x0_rc = LegendOpticalFits._to_table(x0_rc_m, keys(eff))

logl = make_λ0_likelihood(
    x0,
    log_p0,
    x0_rc,
    multiplicity_thr = 8
)

@benchmark logdensityof(logl, eff)
