using LegendOpticalFits

using TypedTables
using Random
using Plots


optmap = mock_optmap((:p13, :r001))
n_channels = length(optmap)
lp0, _ = LegendOpticalFits._to_matrix(log_p0_nominal_ar39(optmap, 500_000))
x0_rc = rand(500_000, n_channels) .< 0.05

ϵ = fill(0.5, n_channels)

plt = plot()
for n_events in [1000, 10_000, 50_000, 100_000, 500_000]
    @info "with $n_events..."
    λ0 = [
        LegendOpticalFits.λ0_model_bulk_ops(ϵ, lp0[1:n_events, :], x0_rc[1:n_events, :], multiplicity_thr = 6)[1]
        for _ in 1:100
    ]
    stephist!(
        plt,
        λ0, bins = 10,
        fill_between = 0, lc = :black,
        ylim = (0, :auto),
        xlabel = "λ0 (first channel)",
        label = "$n_events simulated events"
    )
    display(plt)
end
