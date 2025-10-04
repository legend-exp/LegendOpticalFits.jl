using Test
using LegendOpticalFits: detection_prob, detection_prob_vov

using LegendTestData
using LegendHDF5IO

@testset "optical map utilities" begin
    @testset "access" begin
        # build a mock optical map
        optmap = mock_optmap((:p13, :r001))

        # load some simulation data
        fn = joinpath(legend_test_data_path(), "data", "remage", "th228-full-optional-v0_13.lh5")
        sim_data = lh5open(fn) do h5
            return h5["stp/scint1"][:]
        end

        coords = (sim_data.xloc, sim_data.yloc, sim_data.zloc)

        # take first hit of first event
        x = only(coords[1][1])
        y = only(coords[2][1])
        z = only(coords[3][1])

        # detection_prob: should return a Float64 bin content
        val = detection_prob(first(optmap), x, y, z)
        @test isa(val, Float64)

        # detection_prob_vov: should return a vector-of-vectors aligned with sim_data
        vals = detection_prob_vov(first(optmap), coords...)
        @test length(vals) == length(sim_data)
        @test all(isa(v, AbstractVector{Float64}) for v in vals)

        # consistency check: first element of vals equals scalar detection_prob
        @test vals[1][1] == detection_prob(first(optmap), x, y, z)
    end
end
