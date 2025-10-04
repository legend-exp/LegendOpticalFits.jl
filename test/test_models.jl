using Test
using LegendOpticalFits

using TypedTables
using Random


@testset "forward models" begin
    @testset "λ0 model" begin
        n_events = 100_000
        n_channels = 10

        # low level, basic
        ϵ = rand(n_channels)
        lp0 = log.(rand(n_events, n_channels))
        x0_rc = rand(Bool, n_events, n_channels)

        n_events, n_channels = size(lp0)
        rands = rand(n_events, n_channels, 10)

        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)

        @test length(λ0) == n_channels
        @test all(0 .<= λ0 .<= 1)

        # zero efficiencies -> only rc contribute
        ϵ = zeros(n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(isapprox.(λ0, 0.5; atol = 0.01))

        # zero efficiencies and no rc
        x0_rc = trues(n_events, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(λ0 .== 1)

        # zero efficiencies and loads of rc
        x0_rc = falses(n_events, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(λ0 .== 0)

        # efficiency 1, no rc, no light
        ϵ = ones(n_channels)
        lp0 = zeros(n_events, n_channels)
        x0_rc = trues(n_events, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(λ0 .== 1)

        # lights on
        lp0 = fill(-1000, n_events, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(λ0 .== 0)

        # some light
        lp0 = fill(-1, n_events, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(isapprox.(λ0, exp(-1); atol = 0.01))

        # some efficiency
        ϵ = fill(0.5, n_channels)
        λ0 = LegendOpticalFits._λ0_model_bulk_ops(ϵ, lp0, x0_rc, rands)
        @test all(isapprox.(λ0, exp(-0.5); atol = 0.01))

        # high-level
        ϵ = (; (Symbol("S$s") => 0.5 for s in 1:n_channels)...)
        lp0 = Table(; (Symbol("S$s") => fill(-1, n_events) for s in 1:n_channels)...)
        x0_rc = Table(; (Symbol("S$s") => trues(n_events) for s in 1:n_channels)...)

        λ0_hl = λ0_model(ϵ, lp0, x0_rc, n_rands = 10)

        @test λ0_hl isa NamedTuple
        @test all(isapprox.(collect(values(λ0_hl)), λ0; atol = 0.01))
    end
end
