using Test
using LegendOpticalFits


@testset "test forward models" begin
    @testset "no-light probability" begin
        n_events = 1000
        n_channels = 10

        scaling_factors = rand(n_channels)
        log_p0_nominal = log.(rand(n_events, n_channels))
        x0_random_coin = rand(Bool, n_events, n_channels)

        fractional_expectations = expected_no_light_fraction(
            scaling_factors,
            log_p0_nominal,
            x0_random_coin
        )

        @test length(fractional_expectations) == n_channels
        @test all(0 .<= fractional_expectations .<= 1)
    end
end
