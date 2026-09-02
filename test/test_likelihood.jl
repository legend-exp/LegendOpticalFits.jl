using Test
using LegendOpticalFits

using TypedTables
using Random
using DensityInterface


@testset "likelihoods" begin
    @testset "λ0" begin
        nev_sim = 100_000
        nev_data = 10_000
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

    @testset "multiplicity windows" begin
        Random.seed!(1234)

        nev_sim = 100_000
        nev_data = 20_000

        runsel = (:p13, :r001)
        optmap = mock_optmap(runsel)
        x0, _, lp0 = mock_ar39_data(optmap, nev_sim, eff = 0.5)
        x0 = x0[1:nev_data, :]
        x0_rc = Table(; (k => trues(nev_sim) for k in keys(optmap))...)
        eff = (; (k => 0.5 for k in keys(optmap))...)

        # a window is a strict subset of the corresponding open-ended threshold,
        # and disjoint windows partition it exactly
        _, n_lo = λ0_data(x0, multiplicity_thr = 6, multiplicity_max = 15)
        _, n_hi = λ0_data(x0, multiplicity_thr = 15)
        _, n_all = λ0_data(x0, multiplicity_thr = 6)
        @test n_lo + n_hi == n_all
        @test 0 < n_lo < n_all

        # an empty window is rejected rather than silently giving 0/0
        @test_throws ErrorException λ0_data(x0, multiplicity_thr = 6, multiplicity_max = 6)

        # the forward model honours the upper bound too. A low-multiplicity
        # window has less light, so every channel must see no light MORE often
        # there than in the high-multiplicity window
        λ_lo = λ0_model(eff, lp0, x0_rc, multiplicity_thr = 6, multiplicity_max = 15)
        λ_hi = λ0_model(eff, lp0, x0_rc, multiplicity_thr = 15)
        @test keys(λ_lo) == keys(λ_hi)
        @test all(0 .<= values(λ_lo) .<= 1)
        @test all(values(λ_lo) .> values(λ_hi))

        # slicing one run into two disjoint windows, each its own dataset
        datasets = [
            (x0 = x0, log_p0 = lp0, x0_rc = x0_rc, mult_thr = 6, mult_max = 15),
            (x0 = x0, log_p0 = lp0, x0_rc = x0_rc, mult_thr = 15)
        ]
        logl = make_λ0_joint_likelihood(datasets)

        ll = logdensityof(logl, eff)
        @test ll isa Real
        @test isfinite(ll)

        k1 = first(keys(optmap))
        @test ll > logdensityof(logl, (; eff..., k1 => 0.4))
        @test ll > logdensityof(logl, (; eff..., k1 => 0.6))

        # equal N by construction: the windows have a fixed yield per raw event,
        # so feeding the low window proportionally fewer events makes both
        # datasets carry the same statistical weight
        n_scaled = ceil(Int, nev_data * n_hi / n_lo)
        _, n_lo_scaled = λ0_data(x0[1:n_scaled, :], multiplicity_thr = 6, multiplicity_max = 15)
        @test isapprox(n_lo_scaled, n_hi, rtol = 0.1)
    end
end
