using LegendOpticalFits

using TypedTables
using Random
using Plots


optmap = mock_optmap((:p13, :r001))
n_channels = length(optmap)
nev_sim = 10_000

lp0, _ = LegendOpticalFits._to_matrix(log_p0_nominal_ar39(optmap, nev_sim))
x0_rc = trues(nev_sim, n_channels)

ϵ = fill(0.5, n_channels)

plt = plot()
for n_rands in [1, 10, 50, 100]
    @info "with $n_rands..."

    λ0 = Iterators.flatten(
        LegendOpticalFits._λ0_model_bulk_ops(
            ϵ,
            lp0,
            x0_rc,
            rand(nev_sim, n_channels, n_rands),
            multiplicity_thr = 6
        )
        for _ in 1:200
    )
    stephist!(
        plt,
        collect(λ0),
        bins = 20,
        fill_between = 0, lc = :black,
        ylim = (0, :auto),
        xlabel = "λ0_model (all channels, $nev_sim events)",
        label = "$n_rands randoms"
    )
    display(plt)
end


nev_sim = 100_000

x0_rc = trues(nev_sim, n_channels)
rands = rand(nev_sim, n_channels, 10)

ϵ = fill(0.5, n_channels)

plt = plot()
for n_events in [10_000, 50_000, 100_000]
    @info "with $n_events..."

    _x0_rc = x0_rc[1:n_events, :]
    _rands = rands[1:n_events, :, :]

    λ0 = Iterators.flatten(
        LegendOpticalFits._λ0_model_bulk_ops(
            ϵ,
            LegendOpticalFits._to_matrix(log_p0_nominal_ar39(optmap, n_events))[1],
            _x0_rc,
            _rands,
            multiplicity_thr = 6
        )
        for _ in 1:2
    )
    stephist!(
        plt,
        collect(λ0),
        bins = 20,
        fill_between = 0, lc = :black,
        ylim = (0, :auto),
        xlabel = "λ0_model (first channel)",
        label = "$n_events events"
    )
    display(plt)
end
