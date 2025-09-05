using Test
using LegendOpticalFits

@testset "channel maps" begin
    @test LegendOpticalFits.rawid2detname(LegendOpticalFits._chmap_p13_r001, 1057602) == :S030
end
