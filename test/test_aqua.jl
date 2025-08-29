# This file is a part of LegendOpticalFits.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import LegendOpticalFits

Test.@testset "Package ambiguities" begin
    Test.@test isempty(Test.detect_ambiguities(LegendOpticalFits))
end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        LegendOpticalFits,
        ambiguities = true
    )
end # testset
