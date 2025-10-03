import Test

Test.@testset verbose=true "Package LegendOpticalFits" begin
    include("test_aqua.jl")
    include("fixtures.jl")
    include("test_data.jl")
    include("test_optmap.jl")
    include("test_simulation.jl")
    include("test_models.jl")
    include("test_priors.jl")
    include("test_likelihood.jl")
    include("test_utils.jl")
    include("test_channelmap.jl")
    include("test_docs.jl")
end # testset
