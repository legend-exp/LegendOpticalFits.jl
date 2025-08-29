# This file is a part of LegendOpticalFits.jl, licensed under the MIT License (MIT).

module LegendOpticalFitsTestExt

using Test: @test, @testset
using LegendOpticalFits

function LegendOpticalFits.test_hellofunc(f)
    @testset "test_hellofunc for $(nameof(typeof(f)))" begin
        @test f() == "Hello, World!"
    end
    return nothing
end

end # module LegendOpticalFitsTestExt
