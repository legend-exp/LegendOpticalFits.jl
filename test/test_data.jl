using Test
using LegendOpticalFits

using TypedTables
using Random

@testset "data manipulation" begin
    @testset "λ0 from data" begin
        keys = [Symbol("S$s") for s in 1:10]
        x0_random_coin = Table(; (k => rand(Bool, 10_000) for k in keys)...)
        λ0, N = λ0_data(x0_random_coin, multiplicity_thr = 8)

        @test Set(propertynames(λ0)) == Set(keys)
        @test N isa Integer
        @test N > 0
    end
end
