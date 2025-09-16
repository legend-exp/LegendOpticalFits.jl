using Test
using LegendOpticalFits

using TypedTables


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
        scaling_factors = Dict(Symbol("S$s") => rand() for s in 1:10)
        log_p0_nominal = Table(; (Symbol("S$s") => log.(rand(n_events)) for s in 1:10)...)
        x0_random_coin = Table(; (Symbol("S$s") => rand(Bool, n_events) for s in 1:10)...)

        λ0_hl = λ0_model(
            scaling_factors,
            log_p0_nominal,
            x0_random_coin
        )

        @test λ0_hl isa Dict

    end
end
