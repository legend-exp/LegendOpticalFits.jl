using Test
using LegendOpticalFits

using TypedTables
using Random
using StatsBase
using DensityInterface


@testset "likelihoods" begin
    @testset "λ0" begin
        n_events = 100_000
        runsel = (:p13, :r001)

        optmap = mock_optmap(runsel)
        x0, x0_rc, log_p0 = mock_ar39_data(optmap, n_events, eff = 0.5)

        logl = make_λ0_likelihood(
            x0,
            log_p0,
            x0_rc,
            multiplicity_thr = 6
        )

        eff = (; (k => 0.5 for k in keys(optmap))...)
        @test logdensityof(logl, eff) isa Real
    end
end
