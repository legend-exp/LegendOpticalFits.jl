using Test
using LegendOpticalFits

using TypedTables
using Random
using StatsBase


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

    @testset "log(p0) from Ar-39" begin
        # create a 10x10x10 histogram with random probabilities
        edges = (collect(0.0:1.0:10.0), collect(0.0:1.0:10.0), collect(0.0:1.0:10.0))
        weights = rand(10, 10, 10)

        # set some random voxels to sentinel -1 (invalid)
        for _ in 1:50
            weights[rand(1:10), rand(1:10), rand(1:10)] = -1.0
        end

        h = StatsBase.Histogram(edges, weights)

        # build a NamedTuple optical map with 10 channels
        channels = Tuple(Symbol("S" * lpad(string(i), 3, '0')) for i in 1:10)
        optmap = (; (ch => h for ch in channels)...)

        n_events = 100

        p0 = log_p0_nominal_ar39(optmap, n_events)

        @test p0 isa Table
        @test Set(propertynames(p0)) == Set(channels)
        @test all(length(getproperty(p0, ch)) == n_events for ch in channels)
        @test all(all(getproperty(p0, ch) .<= 0) for ch in channels)
        @test all(all(isfinite.(getproperty(p0, ch))) for ch in channels)

        # with zero light yield, all entries should be exactly zero
        p0_zero = log_p0_nominal_ar39(optmap, n_events; light_yield = 0)
        @test all(all(getproperty(p0_zero, ch) .== 0) for ch in channels)
    end
end
