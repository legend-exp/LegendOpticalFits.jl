using Test
using LegendOpticalFits: ustrip_vov

using TypedTables
using Unitful

@testset "utilities" begin
    @testset "ar39 beta energy distribution" begin
        dist = ar39_beta_energy_dist()

        # draw some samples
        samples = rand(dist, 10_000)

        @test all(isfinite, samples)
        @test minimum(samples) ≥ minimum([c.a for c in dist.components])
        @test maximum(samples) ≤ maximum([c.b for c in dist.components])
    end

    @testset "matrix <-> table" begin
        table = Table(; a = [1, 2, 3], b = [4, 5, 6])
        matrix = [[1, 2, 3] [4, 5, 6]]
        colnames = (:a, :b)

        @test LegendOpticalFits._to_matrix(table) == (matrix, colnames)
        @test LegendOpticalFits._to_table(matrix, colnames) == table

        order = (:b, :a)
        @test LegendOpticalFits._to_matrix(table, order = order) == ([[4, 5, 6] [1, 2, 3]], order)
    end

    @testset "ustrip_vov" begin
        v = [[1.0, 2.0]u"keV", [3.0]u"keV"]
        result = ustrip_vov(u"keV", v)
        @test result == [[1.0, 2.0], [3.0]]
    end
end
