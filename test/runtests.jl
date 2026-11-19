using OrnsteinUhlenbeckJumpProcess
using Test, Dates, Statistics, LinearAlgebra, StableRNGs
using DataFrames, Distributions, MCMCChains, LogExpFunctions, Turing
using CairoMakie
using Aqua

@testset "OrnsteinUhlenbeckJumpProcess.jl" begin

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OrnsteinUhlenbeckJumpProcess)
    end

    @testset "summarize_parameters" begin
        rng = StableRNG(1)
        df = DataFrame(
            Symbol("beta[1]") => randn(rng, 1_000),
            Symbol("beta[2]") => randn(rng, 1_000)
        )
        summary = summarize_parameters(df)

        @test names(summary) == ["Parameter", "Median", "Lower", "Upper", "Index"]
        @test summary.Parameter == ["beta", "beta"]
        @test summary.Index == [1, 2]
        @test summary.Median ≈ [0.0297674, 0.054809] rtol = 1e-5
        @test summary.Lower ≈ [-1.89686, -1.77497] rtol = 1e-5
        @test summary.Upper ≈ [1.93410, 2.03375] rtol = 1e-5
    end

    @testset "posterior_directional_error" begin
        a = [1.0, 2.0, 3.0, 4.0]
        b = [2.0, 3.0, 2.0, 5.0]

        @test posterior_directional_error(a, b) == 0.25
        @test posterior_directional_error([1.0, 2.0], [2.0, 1.0]) == 0.5
        @test posterior_directional_error([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) == 0.0
    end

    @testset "compare_groups" begin
        rng = StableRNG(1)
        df = DataFrame(
            Symbol("group_effect[1]") => randn(rng, 1_000),
            Symbol("group_effect[2]") => randn(rng, 1_000) .+ 1,
            Symbol("group_effect[3]") => randn(rng, 1_000) .+ 2
        )
        comparisons = compare_groups(df)

        @test names(comparisons) ==
              ["Parameter", "IndexA", "IndexB", "PosteriorDirectionalError"]
        @test comparisons.Parameter == ["group_effect", "group_effect", "group_effect"]
        @test comparisons.IndexA == [1, 1, 2]
        @test comparisons.IndexB == [2, 3, 3]
        @test comparisons.PosteriorDirectionalError ≈
              [0.222, 0.07, 0.231] atol = 5e-4
    end

    @testset "hours_elapsed" begin
        x = [
            DateTime(2026, 1, 1, 12, 0, 0),
            DateTime(2026, 1, 1, 15, 0, 0),
            DateTime(2026, 1, 2, 12, 0, 0)
        ]

        @test hours_elapsed(x) == [0.0, 3.0, 21.0]
        @test hours_elapsed([DateTime(2026, 1, 1, 8, 0, 0)]) == [0.0]
        @test hours_elapsed([
            DateTime(2026, 1, 1, 8, 0, 0),
            DateTime(2026, 1, 1, 9, 0, 0, 500)
        ]) ≈ [0.0, 1.000139] rtol = 1e-5
    end

    @testset "logjoint (OrnsteinUhlenbeckModel)" begin
        args = (
            group_indices = [1, 1],
            subject_indices = [1, 1, 2, 2],
            z_obs = [0.5, 0.6, 0.7, 0.4],
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        model = OrnsteinUhlenbeckModel(; args...)

        scenarios = [
            (
                μ_group = [0.5], θ_group = [0.0], γ_group = [0.0],
                κ_μ_group_raw = [0.0], κ_θ_group_raw = [0.0], κ_γ_group_raw = [0.0],
                μ_subject_raw = [0.5, -0.5], θ_subject_raw = [0.2, -0.2],
                γ_subject_raw = [0.1, -0.1], z_raw = [0.3, -0.4, 0.5, -0.6],
                σ_raw = 0.1
            ),
            (
                μ_group = [0.3], θ_group = [0.5], γ_group = [-0.5],
                κ_μ_group_raw = [1.0], κ_θ_group_raw = [-1.0], κ_γ_group_raw = [0.5],
                μ_subject_raw = [1.0, 0.0], θ_subject_raw = [-0.3, 0.4],
                γ_subject_raw = [0.2, -0.6], z_raw = [-0.2, 0.1, -0.3, 0.8],
                σ_raw = -0.5
            ),
            (
                μ_group = [-0.2], θ_group = [-0.36651], γ_group = [-2.30259],
                κ_μ_group_raw = [0.7], κ_θ_group_raw = [-0.2], κ_γ_group_raw = [1.3],
                μ_subject_raw = [-0.8, 0.6], θ_subject_raw = [0.5, 0.1],
                γ_subject_raw = [-0.4, 0.9], z_raw = [1.0, -1.0, 0.2, -0.7],
                σ_raw = 0.9
            )
        ]

        function expected_logjoint(model_args, θ)
            κ(raw) = exp(-2.30259 + 0.5 * raw)

            μ1 = θ.μ_group[1] + κ(θ.κ_μ_group_raw[1]) * θ.μ_subject_raw[1]
            μ2 = θ.μ_group[1] + κ(θ.κ_μ_group_raw[1]) * θ.μ_subject_raw[2]
            θ1 = exp(θ.θ_group[1] + κ(θ.κ_θ_group_raw[1]) * θ.θ_subject_raw[1])
            θ2 = exp(θ.θ_group[1] + κ(θ.κ_θ_group_raw[1]) * θ.θ_subject_raw[2])
            γ1 = exp(θ.γ_group[1] + κ(θ.κ_γ_group_raw[1]) * θ.γ_subject_raw[1])
            γ2 = exp(θ.γ_group[1] + κ(θ.κ_γ_group_raw[1]) * θ.γ_subject_raw[2])

            z1 = μ1 + γ1 * θ.z_raw[1]
            z2 = μ1 + (z1 - μ1) * exp(-θ1 * model_args.Δt[2]) +
                 sqrt(γ1^2 * -expm1(-2 * θ1 * model_args.Δt[2])) * θ.z_raw[2]
            z3 = μ2 + γ2 * θ.z_raw[3]
            z4 = μ2 + (z3 - μ2) * exp(-θ2 * model_args.Δt[4]) +
                 sqrt(γ2^2 * -expm1(-2 * θ2 * model_args.Δt[4])) * θ.z_raw[4]
            z = [z1, z2, z3, z4]
            σ² = exp(-4.60517 + 0.05 * θ.σ_raw)

            raws = vcat(
                θ.κ_μ_group_raw, θ.κ_θ_group_raw, θ.κ_γ_group_raw,
                θ.μ_subject_raw, θ.θ_subject_raw, θ.γ_subject_raw,
                θ.z_raw, θ.σ_raw
            )

            return logpdf(Normal(0.5, 0.2), θ.μ_group[1]) +
                   logpdf(Normal(-0.36651, 0.5), θ.θ_group[1]) +
                   logpdf(Normal(-2.30259, 0.5), θ.γ_group[1]) +
                   sum(logpdf(Normal(0, 1), x) for x in raws) +
                   logpdf(MvNormal(z, σ² * I), model_args.z_obs)
        end

        for θ in scenarios
            @test logjoint(model, θ) ≈ expected_logjoint(args, θ) rtol = 1e-9
        end
    end

    @testset "logjoint (OrnsteinUhlenbeckJumpModel)" begin
        args = (
            group_indices = [1, 1],
            subject_indices = [1, 1, 2, 2],
            z_obs = [0.5, 0.6, 0.7, 0.4],
            Δt = [0.0, 1.0, 0.0, 2.0],
            M⁺ = [0.005, 0.005, 0.003, 0.002],
            M⁻ = [0.001, 0.0, 0.002, 0.001],
            Δτ = [0.0, 0.5, 0.0, 1.0],
            is_transition = [false, true, false, true]
        )
        model = OrnsteinUhlenbeckJumpModel(; args...)
        θ = (
            μ_group = [0.5], θ_group = [0.0], γ_group = [0.0],
            α⁺_group = [2.30259], α⁻_group = [1.60944],
            κ_μ_group_raw = [0.0], κ_θ_group_raw = [0.0], κ_γ_group_raw = [0.0],
            κ_α⁺_group_raw = [0.0], κ_α⁻_group_raw = [0.0],
            μ_subject_raw = [0.5, -0.5], θ_subject_raw = [0.2, -0.2],
            γ_subject_raw = [0.1, -0.1], α⁺_subject_raw = [0.3, -0.3],
            α⁻_subject_raw = [0.2, 0.4],
            z_raw = [0.3, -0.4, 0.5, -0.6], σ_raw = 0.1
        )

        κ0 = exp(-2.30259)
        μ1 = 0.5 + κ0 * 0.5
        μ2 = 0.5 - κ0 * 0.5
        θ1 = exp(κ0 * 0.2)
        θ2 = exp(-κ0 * 0.2)
        γ1 = exp(κ0 * 0.1)
        γ2 = exp(-κ0 * 0.1)
        α⁺1 = exp(2.30259 + κ0 * 0.3)
        α⁺2 = exp(2.30259 - κ0 * 0.3)
        α⁻1 = exp(1.60944 + κ0 * 0.2)
        α⁻2 = exp(1.60944 + κ0 * 0.4)

        z1 = μ1 + γ1 * 0.3
        z2 = μ1 + (z1 - μ1) * exp(-θ1 * 1.0) +
             sqrt(γ1^2 * -expm1(-2 * θ1 * 1.0)) * (-0.4) +
             (α⁺1 * 0.005 - α⁻1 * 0.0) * exp(-θ1 * (1.0 - 0.5))
        z3 = μ2 + γ2 * 0.5
        z4 = μ2 + (z3 - μ2) * exp(-θ2 * 2.0) +
             sqrt(γ2^2 * -expm1(-2 * θ2 * 2.0)) * (-0.6) +
             (α⁺2 * 0.002 - α⁻2 * 0.001) * exp(-θ2 * (2.0 - 1.0))
        z = [z1, z2, z3, z4]
        σ² = exp(-4.60517 + 0.05 * 0.1)
        raws = [
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.5, -0.5, 0.2, -0.2, 0.1, -0.1, 0.3, -0.3, 0.2, 0.4,
            0.3, -0.4, 0.5, -0.6, 0.1
        ]
        expected = logpdf(Normal(0.5, 0.2), 0.5) +
                   logpdf(Normal(-0.36651, 0.5), 0.0) +
                   logpdf(Normal(-2.30259, 0.5), 0.0) +
                   logpdf(Normal(2.30259, 0.5), 2.30259) +
                   logpdf(Normal(2.30259, 0.5), 1.60944) +
                   sum(logpdf(Normal(0, 1), x) for x in raws) +
                   logpdf(MvNormal(z, σ² * I), args.z_obs)

        @test logjoint(model, θ) ≈ expected rtol = 1e-9
    end

    @testset "logjoint (EventProbabilityModel)" begin
        args = (group_indices = [1, 2], subject_indices = [1, 2, 1, 2], e = [1, 0, 0, 1])
        model = EventProbabilityModel(; args...)
        θ = (β_group = [0.5, -0.5], κ_β_group = [0.2, 0.3], β_subject_raw = [0.3, -0.4])
        β_subject = [0.5 + 0.2 * 0.3, -0.5 + 0.3 * (-0.4)]
        β = β_subject[[1, 2, 1, 2]]
        expected = logpdf(Normal(0, 1), 0.5) +
                   logpdf(Normal(0, 1), -0.5) +
                   logpdf(truncated(Normal(0, 1); lower = 0), 0.2) +
                   logpdf(truncated(Normal(0, 1); lower = 0), 0.3) +
                   logpdf(Normal(0, 1), 0.3) +
                   logpdf(Normal(0, 1), -0.4) +
                   sum(logpdf(BernoulliLogit(β[i]), args.e[i]) for i in eachindex(args.e))

        @test logjoint(model, θ) ≈ expected rtol = 1e-12
    end

    @testset "clean_parameters (hand-built chains)" begin
        κ0 = exp(-2.30259)

        # OrnsteinUhlenbeckModel with one group and one subject
        ou_args = (group_indices = [1],)
        ou_chains = Chains(
            [0.5 -0.36651 -2.30259 0.0 0.0 0.0 0.5 0.2 0.3],
            ["μ_group[1]", "θ_group[1]", "γ_group[1]", "κ_μ_group_raw[1]",
             "κ_θ_group_raw[1]", "κ_γ_group_raw[1]", "μ_subject_raw[1]",
             "θ_subject_raw[1]", "γ_subject_raw[1]"]
        )
        ou_cleaned = clean_parameters(OrnsteinUhlenbeckModel, ou_chains, ou_args)

        @test names(ou_cleaned) == [
            "μ_group[1]", "θ_group[1]", "γ_group[1]", "σ_group[1]",
            "κ_μ_group[1]", "κ_θ_group[1]", "κ_γ_group[1]", "μ_subject[1]",
            "θ_subject[1]", "γ_subject[1]", "σ_subject[1]"
        ]
        @test ou_cleaned[1, "μ_group[1]"] == 50.0
        @test ou_cleaned[1, "θ_group[1]"] ≈ log(2) / exp(-0.36651) rtol = 1e-12
        @test ou_cleaned[1, "γ_group[1]"] ≈ 100 * exp(-2.30259) rtol = 1e-12
        @test ou_cleaned[1, "σ_group[1]"] ≈
              100 * exp(-2.30259) * sqrt(2 * exp(-0.36651)) rtol = 1e-12
        @test ou_cleaned[1, "κ_μ_group[1]"] ≈ 100 * κ0 rtol = 1e-12
        @test ou_cleaned[1, "κ_θ_group[1]"] ≈ κ0 rtol = 1e-12
        @test ou_cleaned[1, "κ_γ_group[1]"] ≈ κ0 rtol = 1e-12
        @test ou_cleaned[1, "μ_subject[1]"] ≈ 100 * (0.5 + κ0 * 0.5) rtol = 1e-12
        @test ou_cleaned[1, "θ_subject[1]"] ≈
              log(2) / exp(-0.36651 + κ0 * 0.2) rtol = 1e-12
        @test ou_cleaned[1, "γ_subject[1]"] ≈
              100 * exp(-2.30259 + κ0 * 0.3) rtol = 1e-12
        @test ou_cleaned[1, "σ_subject[1]"] ≈
              100 * exp(-2.30259 + κ0 * 0.3) *
              sqrt(2 * exp(-0.36651 + κ0 * 0.2)) rtol = 1e-12

        # OrnsteinUhlenbeckJumpModel with one group and one subject
        jump_args = (group_indices = [1],)
        jump_chains = Chains(
            [0.5 -0.36651 -2.30259 2.30259 1.60944 0.0 0.0 0.0 0.0 0.0 0.5 0.2 0.3 -0.25 0.4],
            ["μ_group[1]", "θ_group[1]", "γ_group[1]", "α⁺_group[1]", "α⁻_group[1]",
             "κ_μ_group_raw[1]", "κ_θ_group_raw[1]", "κ_γ_group_raw[1]",
             "κ_α⁺_group_raw[1]", "κ_α⁻_group_raw[1]", "μ_subject_raw[1]",
             "θ_subject_raw[1]", "γ_subject_raw[1]", "α⁺_subject_raw[1]",
             "α⁻_subject_raw[1]"]
        )
        jump_cleaned = clean_parameters(OrnsteinUhlenbeckJumpModel, jump_chains, jump_args)

        @test names(jump_cleaned) == [
            "μ_group[1]", "θ_group[1]", "γ_group[1]", "σ_group[1]",
            "α⁺_group[1]", "α⁻_group[1]", "κ_μ_group[1]", "κ_θ_group[1]",
            "κ_γ_group[1]", "κ_α⁺_group[1]", "κ_α⁻_group[1]", "μ_subject[1]",
            "θ_subject[1]", "γ_subject[1]", "σ_subject[1]", "α⁺_subject[1]",
            "α⁻_subject[1]"
        ]
        @test jump_cleaned[1, "α⁺_group[1]"] ≈ exp(2.30259) rtol = 1e-12
        @test jump_cleaned[1, "α⁻_group[1]"] ≈ exp(1.60944) rtol = 1e-12
        @test jump_cleaned[1, "κ_α⁺_group[1]"] ≈ κ0 rtol = 1e-12
        @test jump_cleaned[1, "α⁺_subject[1]"] ≈ exp(2.30259 - κ0 * 0.25) rtol = 1e-12
        @test jump_cleaned[1, "α⁻_subject[1]"] ≈ exp(1.60944 + κ0 * 0.4) rtol = 1e-12
        @test jump_cleaned[1, "μ_subject[1]"] ≈ 100 * (0.5 + κ0 * 0.5) rtol = 1e-12

        # EventProbabilityModel with two groups and two subjects
        event_args = (group_indices = [1, 2],)
        event_chains = Chains(
            [0.5 -0.5 0.2 0.3 0.3 -0.4],
            ["β_group[1]", "β_group[2]", "κ_β_group[1]", "κ_β_group[2]",
             "β_subject_raw[1]", "β_subject_raw[2]"]
        )
        event_cleaned = clean_parameters(EventProbabilityModel, event_chains, event_args)

        @test names(event_cleaned) == [
            "β_group[1]", "β_group[2]", "κ_β_group[1]", "κ_β_group[2]",
            "β_subject[1]", "β_subject[2]"
        ]
        @test event_cleaned[1, "β_group[1]"] ≈ logistic(0.5) rtol = 1e-12
        @test event_cleaned[1, "β_group[2]"] ≈ logistic(-0.5) rtol = 1e-12
        @test event_cleaned[1, "κ_β_group[1]"] == 0.2
        @test event_cleaned[1, "κ_β_group[2]"] == 0.3
        @test event_cleaned[1, "β_subject[1]"] ≈
              logistic(0.5 + 0.2 * 0.3) rtol = 1e-12
        @test event_cleaned[1, "β_subject[2]"] ≈
              logistic(-0.5 + 0.3 * (-0.4)) rtol = 1e-12
    end

    @testset "clean_parameters (prior samples)" begin
        # OrnsteinUhlenbeckModel with two groups and two subjects
        ou_model = OrnsteinUhlenbeckModel(;
            group_indices = [1, 2],
            subject_indices = [1, 2, 1, 2],
            z_obs = missing,
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        θ = rand(StableRNG(1), ou_model)
        ou_data_model = OrnsteinUhlenbeckModel(;
            group_indices = [1, 2],
            subject_indices = [1, 2, 1, 2],
            z_obs = θ.z_obs,
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        ou_prior = sample(StableRNG(1), ou_data_model, Prior(), 50; progress = false)
        ou_cleaned = clean_parameters(
            OrnsteinUhlenbeckModel, ou_prior, (group_indices = [1, 2],)
        )

        @test nrow(ou_cleaned) == 50
        @test names(ou_cleaned) == [
            "μ_group[1]", "μ_group[2]",
            "θ_group[1]", "θ_group[2]",
            "γ_group[1]", "γ_group[2]",
            "σ_group[1]", "σ_group[2]",
            "κ_μ_group[1]", "κ_μ_group[2]",
            "κ_θ_group[1]", "κ_θ_group[2]",
            "κ_γ_group[1]", "κ_γ_group[2]",
            "μ_subject[1]", "μ_subject[2]",
            "θ_subject[1]", "θ_subject[2]",
            "γ_subject[1]", "γ_subject[2]",
            "σ_subject[1]", "σ_subject[2]"
        ]
        ou_raw = DataFrame(ou_prior)
        for i in 1:nrow(ou_cleaned)
            @test ou_cleaned[i, "μ_subject[2]"] ≈
                  100 * (ou_raw[i, "μ_group[2]"] +
                         exp(-2.30259 + 0.5 * ou_raw[i, "κ_μ_group_raw[2]"]) *
                         ou_raw[i, "μ_subject_raw[2]"]) rtol = 1e-9
            @test ou_cleaned[i, "θ_subject[1]"] ≈
                  log(2) / exp(ou_raw[i, "θ_group[1]"] +
                               exp(-2.30259 + 0.5 * ou_raw[i, "κ_θ_group_raw[1]"]) *
                               ou_raw[i, "θ_subject_raw[1]"]) rtol = 1e-9
            @test ou_cleaned[i, "γ_subject[1]"] ≈
                  100 * exp(ou_raw[i, "γ_group[1]"] +
                            exp(-2.30259 + 0.5 * ou_raw[i, "κ_γ_group_raw[1]"]) *
                            ou_raw[i, "γ_subject_raw[1]"]) rtol = 1e-9
            @test ou_cleaned[i, "σ_subject[2]"] ≈
                  100 * exp(ou_raw[i, "γ_group[2]"] +
                            exp(-2.30259 + 0.5 * ou_raw[i, "κ_γ_group_raw[2]"]) *
                            ou_raw[i, "γ_subject_raw[2]"]) *
                  sqrt(2 * exp(ou_raw[i, "θ_group[2]"] +
                               exp(-2.30259 + 0.5 * ou_raw[i, "κ_θ_group_raw[2]"]) *
                               ou_raw[i, "θ_subject_raw[2]"])) rtol = 1e-9
            @test ou_cleaned[i, "κ_μ_group[1]"] ≈
                  100 * exp(-2.30259 + 0.5 * ou_raw[i, "κ_μ_group_raw[1]"]) rtol = 1e-9
        end

        # OrnsteinUhlenbeckJumpModel with one group and two subjects
        jump_model = OrnsteinUhlenbeckJumpModel(;
            group_indices = [1, 1],
            subject_indices = [1, 1, 2, 2],
            z_obs = [0.5, 0.6, 0.7, 0.4],
            Δt = [0.0, 1.0, 0.0, 2.0],
            M⁺ = [0.005, 0.005, 0.003, 0.002],
            M⁻ = [0.001, 0.0, 0.002, 0.001],
            Δτ = [0.0, 0.5, 0.0, 1.0],
            is_transition = [false, true, false, true]
        )
        jump_prior = sample(StableRNG(1), jump_model, Prior(), 20; progress = false)
        jump_cleaned = clean_parameters(
            OrnsteinUhlenbeckJumpModel, jump_prior, (group_indices = [1, 1],)
        )

        @test nrow(jump_cleaned) == 20
        @test names(jump_cleaned) == [
            "μ_group[1]", "θ_group[1]", "γ_group[1]", "σ_group[1]",
            "α⁺_group[1]", "α⁻_group[1]", "κ_μ_group[1]", "κ_θ_group[1]",
            "κ_γ_group[1]", "κ_α⁺_group[1]", "κ_α⁻_group[1]", "μ_subject[1]",
            "μ_subject[2]", "θ_subject[1]", "θ_subject[2]", "γ_subject[1]",
            "γ_subject[2]", "σ_subject[1]", "σ_subject[2]", "α⁺_subject[1]",
            "α⁺_subject[2]", "α⁻_subject[1]", "α⁻_subject[2]"
        ]
        jump_raw = DataFrame(jump_prior)
        for i in 1:nrow(jump_cleaned)
            @test jump_cleaned[i, "α⁺_subject[2]"] ≈
                  exp(jump_raw[i, "α⁺_group[1]"] +
                      exp(-2.30259 + 0.5 * jump_raw[i, "κ_α⁺_group_raw[1]"]) *
                      jump_raw[i, "α⁺_subject_raw[2]"]) rtol = 1e-9
            @test jump_cleaned[i, "θ_subject[1]"] ≈
                  log(2) / exp(jump_raw[i, "θ_group[1]"] +
                               exp(-2.30259 + 0.5 * jump_raw[i, "κ_θ_group_raw[1]"]) *
                               jump_raw[i, "θ_subject_raw[1]"]) rtol = 1e-9
        end

        # EventProbabilityModel with two groups and two subjects
        event_model = EventProbabilityModel(;
            group_indices = [1, 2],
            subject_indices = [1, 2, 1, 2],
            e = [1, 0, 0, 1]
        )
        event_prior = sample(StableRNG(1), event_model, Prior(), 20; progress = false)
        event_cleaned = clean_parameters(
            EventProbabilityModel, event_prior, (group_indices = [1, 2],)
        )

        @test nrow(event_cleaned) == 20
        @test names(event_cleaned) == [
            "β_group[1]", "β_group[2]", "κ_β_group[1]", "κ_β_group[2]",
            "β_subject[1]", "β_subject[2]"
        ]
        event_raw = DataFrame(event_prior)
        for i in 1:nrow(event_cleaned)
            @test event_cleaned[i, "β_subject[2]"] ≈
                  logistic(event_raw[i, "β_group[2]"] +
                           event_raw[i, "κ_β_group[2]"] *
                           event_raw[i, "β_subject_raw[2]"]) rtol = 1e-9
            @test event_cleaned[i, "β_group[1]"] ≈
                  logistic(event_raw[i, "β_group[1]"]) rtol = 1e-9
            @test event_cleaned[i, "κ_β_group[1]"] == event_raw[i, "κ_β_group[1]"]
        end
    end

    @testset "sampling" begin
        ou_model = OrnsteinUhlenbeckModel(;
            group_indices = [1, 1],
            subject_indices = [1, 1, 2, 2],
            z_obs = [0.5, 0.6, 0.7, 0.4],
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        ou_chain = sample(ou_model, NUTS(0.65), 10; progress = false)

        @test size(ou_chain, 1) == 10
        @test Symbol("μ_group[1]") ∈ names(ou_chain)
        @test :σ_raw ∈ names(ou_chain)

        event_model = EventProbabilityModel(;
            group_indices = [1, 2],
            subject_indices = [1, 2, 1, 2],
            e = [1, 0, 0, 1]
        )
        event_chain = sample(event_model, NUTS(0.65), 10; progress = false)

        @test size(event_chain, 1) == 10
        @test Symbol("β_group[1]") ∈ names(event_chain)
    end

    @testset "draw_chains" begin
        chains = Chains(randn(StableRNG(1), 100, 12, 2))
        figure = draw_chains(chains, 5)

        @test figure isa CairoMakie.Figure
        @test first(size(figure.layout)) == 5
        @test first(size(draw_chains(chains).layout)) == 10
    end

    @testset "model evaluation errors" begin
        bad_obs = OrnsteinUhlenbeckModel(;
            group_indices = [1, 1],
            subject_indices = [1, 1, 2, 2],
            z_obs = [0.5, 0.6, 0.7],
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        @test_throws Exception sample(bad_obs, Prior(), 1; progress = false)

        bad_subject = OrnsteinUhlenbeckModel(;
            group_indices = [1, 1],
            subject_indices = [1, 1, 3, 3],
            z_obs = [0.5, 0.6, 0.7, 0.4],
            Δt = [0.0, 1.0, 0.0, 2.0],
            is_transition = [false, true, false, true]
        )
        @test_throws Exception sample(bad_subject, Prior(), 1; progress = false)
    end
end

