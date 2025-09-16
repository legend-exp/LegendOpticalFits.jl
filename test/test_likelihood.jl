using Test
using LegendOpticalFits

using TypedTables
using Random
using StatsBase
using DensityInterface


@testset "likelihoods" begin
    @testset "λ0" begin
        n_events = 1000
        n_channels = 10

        x0 = Table(; (Symbol("S$s") => rand(Bool, n_events) for s in 1:10)...)
        log_p0_nominal = Table(; (Symbol("S$s") => log.(rand(n_events)) for s in 1:10)...)
        x0_random_coin = Table(; (Symbol("S$s") => rand(Bool, n_events) for s in 1:10)...)

        logf = make_λ0_likelihood(
            x0,
            log_p0_nominal,
            x0_random_coin
        )

        scaling_factors = (; (Symbol("S$s") => rand() for s in 1:10)...)
        @test logdensityof(logf, scaling_factors) isa Real
    end
end
