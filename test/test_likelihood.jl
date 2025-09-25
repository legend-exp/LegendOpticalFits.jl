using Test
using LegendOpticalFits

using TypedTables
using Random
using DensityInterface


@testset "likelihoods" begin
    @testset "λ0" begin
        nev_sim = 10_000
        nev_data = 1_000
        n_channels = 50

        # test with asimov dataset
        eff = (; (Symbol("S$s") => 0.5 for s in 1:n_channels)...)
        lp0 = Table(; (Symbol("S$s") => fill(-1, nev_sim) for s in 1:n_channels)...)
        x0_rc = Table(; (Symbol("S$s") => trues(nev_sim) for s in 1:n_channels)...)
        x0 = Table(; (Symbol("S$s") => fill(exp(-0.5), nev_data) for s in 1:n_channels)...)
        # x0 = Table(; (Symbol("S$s") => rand(nev_data) .< exp(-0.5) for s in 1:n_channels)...)

        logl = make_λ0_likelihood(x0, lp0, x0_rc)

        # likelihood at the truth
        ll = logdensityof(logl, eff)
        @test ll isa Real

        # test that we have the max logl
        @test ll > logdensityof(logl, (; eff..., S1 = 0.49))
        @test ll > logdensityof(logl, (; eff..., S1 = 0.51))

        # mock ar39 data
        runsel = (:p13, :r001)
        optmap = mock_optmap(runsel)
        x0, _, lp0 = mock_ar39_data(optmap, nev_sim, eff = 0.5)
        x0 = x0[1:nev_data, :]
        x0_rc = Table(; (k => trues(nev_sim) for k in keys(optmap))...)

        logl = make_λ0_likelihood(
            x0, lp0, x0_rc,
            multiplicity_thr = 6
        )

        eff = (; (k => 0.5 for k in keys(optmap))...)
        @test logdensityof(logl, eff) isa Real
    end
end
