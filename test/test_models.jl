using Test
using LegendOpticalFits

using TypedTables
using Random
using StatsBase


@testset "forward models" begin
    @testset "λ0 models" begin
        n_events = 1000
        n_channels = 10

        # low level
        scaling_factors = rand(n_channels)
        log_p0_nominal = log.(rand(n_events, n_channels))
        x0_random_coin = rand(Bool, n_events, n_channels)

        λ0 = λ0_model(
            scaling_factors,
            log_p0_nominal,
            x0_random_coin
        )

        @test length(λ0) == n_channels
        @test all(0 .<= λ0 .<= 1)

        # high-level
        scaling_factors = (; (Symbol("S$s") => rand() for s in 1:10)...)
        log_p0_nominal = Table(; (Symbol("S$s") => log.(rand(n_events)) for s in 1:10)...)
        x0_random_coin = Table(; (Symbol("S$s") => rand(Bool, n_events) for s in 1:10)...)

        λ0_hl = λ0_model(
            scaling_factors,
            log_p0_nominal,
            x0_random_coin
        )

        @test λ0_hl isa NamedTuple

    end
end

@testset "log_p0_nominal_ar39" begin
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
