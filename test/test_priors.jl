using LegendOpticalFits
using Test

using BAT

@testset "priors" begin
    prior = make_efficiencies_prior((:S1, :S2, :S3))
    @test prior isa HierarchicalDistribution
end
