using Test
using LegendOpticalFits

@testset "utilities" begin
    @testset "ar39 beta energy distribution" begin
        dist = ar39_beta_energy_dist()

        # draw some samples
        samples = rand(dist, 10_000)

        @test all(isfinite, samples)
        @test minimum(samples) ≥ minimum([c.a for c in dist.components])
        @test maximum(samples) ≤ maximum([c.b for c in dist.components])
    end
end
