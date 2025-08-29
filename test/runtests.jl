# This file is a part of LegendOpticalFits.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package LegendOpticalFits" begin
    include("test_aqua.jl")
    include("test_optical_model.jl")
    include("test_docs.jl")
end # testset
