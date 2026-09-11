
# versions of OrnsteinUhlenbeckModel and OrnsteinUhlenbeckJumpModel that work with ParetoSmooth

@model function OrnsteinUhlenbeckJumpProcess.OrnsteinUhlenbeckModel(;
        group_indices, subject_indices, z_obs, Δt, is_transition
)
    n_groups = maximum(group_indices)
    n_subjects = length(group_indices)
    n_nodes = length(subject_indices)

    # group-level priors
    μ_group ~ filldist(Normal(0.5, 0.2), n_groups)
    θ_group ~ filldist(Normal(-0.36651, 0.5), n_groups) # log(log(2)); half-time of deviations = 1 hour
    γ_group ~ filldist(Normal(-2.30259, 0.5), n_groups) # log(0.1)

    κ_μ_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_θ_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_γ_group_raw ~ filldist(Normal(0, 1), n_groups)

    κ_μ_group = exp.(-2.30259 .+ 0.5 .* κ_μ_group_raw) # log(0.1)
    κ_θ_group = exp.(-2.30259 .+ 0.5 .* κ_θ_group_raw) # log(0.1)
    κ_γ_group = exp.(-2.30259 .+ 0.5 .* κ_γ_group_raw) # log(0.1)

    # subject-level priors
    μ_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    θ_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    γ_subject_raw ~ filldist(Normal(0, 1), n_subjects)

    @inbounds μ_subject = μ_group[group_indices] .+
                          κ_μ_group[group_indices] .* μ_subject_raw

    @inbounds θ_subject = exp.(θ_group[group_indices] .+
                               κ_θ_group[group_indices] .* θ_subject_raw)

    @inbounds γ_subject = exp.(γ_group[group_indices] .+
                               κ_γ_group[group_indices] .* γ_subject_raw)

    @inbounds μ = μ_subject[subject_indices]
    @inbounds θ = θ_subject[subject_indices]
    @inbounds γ = γ_subject[subject_indices]

    # latent states
    z_raw ~ filldist(Normal(0, 1), n_nodes)

    z = similar(z_raw)

    exp_terms = @. exp(-θ * Δt)
    var_terms = @. γ^2 * -expm1(-2 * θ * Δt)

    # transitions
    @inbounds for i in 1:n_nodes
        if is_transition[i]
            # OU propagation from previous node
            mean = μ[i] + (z[i - 1] - μ[i]) * exp_terms[i]

            z[i] = mean + sqrt(var_terms[i]) * z_raw[i]
        else
            # stationary distribution
            z[i] = μ[i] + γ[i] * z_raw[i]
        end
    end

    # observation noise
    σ_raw ~ Normal(0, 1)
    σ²_obs = exp(-4.60517 + 0.05 * σ_raw) # log(0.01)

    # observation model
    for i in eachindex(z_obs)
        z_obs[i] ~ Normal(z[i], sqrt(σ²_obs))
    end
end

@model function OrnsteinUhlenbeckJumpProcess.OrnsteinUhlenbeckJumpModel(;
        group_indices, subject_indices, z_obs, Δt, M⁺, M⁻, Δτ, is_transition
)
    n_groups = maximum(group_indices)
    n_subjects = length(group_indices)
    n_nodes = length(subject_indices)

    # group-level priors
    μ_group ~ filldist(Normal(0.5, 0.2), n_groups)
    θ_group ~ filldist(Normal(-0.36651, 0.5), n_groups) # log(log(2)); half-time of deviations = 1 hour
    γ_group ~ filldist(Normal(-2.30259, 0.5), n_groups) # log(0.1)
    α⁺_group ~ filldist(Normal(2.30259, 0.5), n_groups) # log(10)
    α⁻_group ~ filldist(Normal(2.30259, 0.5), n_groups) # log(10)

    κ_μ_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_θ_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_γ_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_α⁺_group_raw ~ filldist(Normal(0, 1), n_groups)
    κ_α⁻_group_raw ~ filldist(Normal(0, 1), n_groups)

    κ_μ_group = exp.(-2.30259 .+ 0.5 .* κ_μ_group_raw) # log(0.1)
    κ_θ_group = exp.(-2.30259 .+ 0.5 .* κ_θ_group_raw) # log(0.1)
    κ_γ_group = exp.(-2.30259 .+ 0.5 .* κ_γ_group_raw) # log(0.1)
    κ_α⁺_group = exp.(-2.30259 .+ 0.5 .* κ_α⁺_group_raw) # log(0.1)
    κ_α⁻_group = exp.(-2.30259 .+ 0.5 .* κ_α⁻_group_raw) # log(0.1)

    # subject-level priors
    μ_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    θ_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    γ_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    α⁺_subject_raw ~ filldist(Normal(0, 1), n_subjects)
    α⁻_subject_raw ~ filldist(Normal(0, 1), n_subjects)

    @inbounds μ_subject = μ_group[group_indices] .+
                          κ_μ_group[group_indices] .* μ_subject_raw

    @inbounds θ_subject = exp.(θ_group[group_indices] .+
                               κ_θ_group[group_indices] .* θ_subject_raw)

    @inbounds γ_subject = exp.(γ_group[group_indices] .+
                               κ_γ_group[group_indices] .* γ_subject_raw)

    @inbounds α⁺_subject = exp.(α⁺_group[group_indices] .+
                                κ_α⁺_group[group_indices] .* α⁺_subject_raw)

    @inbounds α⁻_subject = exp.(α⁻_group[group_indices] .+
                                κ_α⁻_group[group_indices] .* α⁻_subject_raw)

    @inbounds μ = μ_subject[subject_indices]
    @inbounds θ = θ_subject[subject_indices]
    @inbounds γ = γ_subject[subject_indices]
    @inbounds α⁺ = α⁺_subject[subject_indices]
    @inbounds α⁻ = α⁻_subject[subject_indices]

    # latent states
    z_raw ~ filldist(Normal(0, 1), n_nodes)

    z = similar(z_raw)

    exp_terms = @. exp(-θ * Δt)
    var_terms = @. γ^2 * -expm1(-2 * θ * Δt)

    impulses = @. (α⁺ * M⁺ - α⁻ * M⁻) * exp(-θ * (Δt - Δτ))

    # transitions
    @inbounds for i in 1:n_nodes
        if is_transition[i]
            # OU propagation from previous node
            mean = μ[i] + (z[i - 1] - μ[i]) * exp_terms[i]

            z[i] = mean + sqrt(var_terms[i]) * z_raw[i] + impulses[i]
        else
            # stationary distribution
            z[i] = μ[i] + γ[i] * z_raw[i]
        end
    end

    # observation noise
    σ_raw ~ Normal(0, 1)
    σ²_obs = exp(-4.60517 + 0.05 * σ_raw) # log(0.01)

    # observation model
    for i in eachindex(z_obs)
        z_obs[i] ~ Normal(z[i], sqrt(σ²_obs))
    end
end