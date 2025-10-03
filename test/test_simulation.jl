using Test
using LegendOpticalFits

using LegendTestData
using StatsBase
using LegendHDF5IO
using TypedTables

@testset "simulation" begin
    @testset "log(p0) from Ar-39" begin
        # create a 10x10x10 histogram with random probabilities
        edges = (collect(0.0:1.0:10.0), collect(0.0:1.0:10.0), collect(0.0:1.0:10.0))
        weights = rand(10, 10, 10)

        # set some random voxels to sentinel -1 (invalid)
        for _ in 1:50
            weights[rand(1:10), rand(1:10), rand(1:10)] = -1.0
        end

        h = StatsBase.Histogram(edges, weights)

        # build a NamedTuple optical map with 10 channels
        channels = Tuple(Symbol("S" * lpad(string(i), 3, '0')) for i in 1:10)
        optmap = (; (ch => h for ch in channels)...)

        n_events = 100

        p0 = log_p0_nominal_ar39(optmap, n_events)

        @test p0 isa Table
        @test Set(propertynames(p0)) == Set(channels)
        @test all(length(getproperty(p0, ch)) == n_events for ch in channels)
        @test all(all(getproperty(p0, ch) .<= 0) for ch in channels)
        @test all(all(isfinite.(getproperty(p0, ch))) for ch in channels)

        # with zero light yield, all entries should be exactly zero
        p0_zero = log_p0_nominal_ar39(optmap, n_events; light_yield = 0)
        @test all(all(getproperty(p0_zero, ch) .== 0) for ch in channels)
    end

    @testset "log(p0) from remage" begin
        # build a mock optical map
        optmap = mock_optmap((:p13, :r001))

        # load some simulation data
        fn = joinpath(legend_test_data_path(), "data", "remage", "th228-full-optional-v0_13.lh5")
        sim_data = lh5open(fn) do h5
            return h5["stp/scint1"][:]
        end

        # call the function
        result = log_p0_nominal(sim_data, optmap; light_yield = 40)

        # basic type/shape checks
        @test result isa Table
        @test length(result) == length(sim_data)  # rows match events

        # columns should match optmap keys
        @test Set(propertynames(result)) == Set(keys(optmap))

        # each column: Vector{Float64} with one value per event
        for ch in propertynames(result)
            vals = getproperty(result, ch)
            @test vals isa Vector{Float64}
            @test length(vals) == length(sim_data)
        end

        # sign check on a sample column (-ξ .* n ≤ 0)
        sample_col = getproperty(result, first(propertynames(result)))
        @test all(v <= 0 for v in sample_col)
    end
end
