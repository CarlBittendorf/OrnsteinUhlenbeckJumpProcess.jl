
####################################################################################################
# EXPORTS
####################################################################################################

export OrnsteinUhlenbeckJumpModel

####################################################################################################
# IMPLEMENTATIONS
####################################################################################################

"""
    OrnsteinUhlenbeckJumpModel(;
        group_indices,
        subject_indices,
        z_obs,
        Δt,
        M⁺,
        M⁻,
        Δτ,
        is_transition,
        kwargs...
    )

Define a hierarchical Ornstein-Uhlenbeck model with positive and negative
impulsive jumps.

The model estimates group-level and subject-level equilibrium means, mean
reversion rates, stationary standard deviations, and jump-response
coefficients. Latent states follow Ornstein-Uhlenbeck transitions augmented by
exponentially decaying impulses.

# Arguments
- `group_indices`: Group index for each subject.
- `subject_indices`: Subject index associated with each observation node.
- `z_obs`: Vector of observed values.
- `Δt`: Time elapsed since the preceding node.
- `M⁺`: Magnitude of the positive impulse associated with each node.
- `M⁻`: Magnitude of the negative impulse associated with each node.
- `Δτ`: Time from the preceding node to the impulse.
- `is_transition`: Boolean vector indicating whether each node follows from
  the preceding node. A value of `false` initializes the node from its
  stationary distribution.
- `kwargs...`: Additional keyword arguments reserved for model extensions.

# Returns
- A Turing model representing the hierarchical Ornstein-Uhlenbeck jump
  process.

# Details
- `μ_subject` is the subject-specific equilibrium mean.
- `θ_subject` is the positive subject-specific mean reversion rate.
- `γ_subject` is the positive subject-specific stationary standard deviation.
- `α⁺_subject` and `α⁻_subject` are positive response coefficients for positive
  and negative impulses.
- Subject-level parameters use non-centered hierarchical parameterizations.
- For transition nodes, the Ornstein-Uhlenbeck conditional mean is
  `μ + (z_previous - μ) * exp(-θ * Δt)`.
- The conditional transition variance is
  `γ² * (1 - exp(-2θ * Δt))`.
- The impulse contribution is
  `(α⁺ * M⁺ - α⁻ * M⁻) * exp(-θ * (Δt - Δτ))`.
- Nodes that do not represent transitions are initialized as
  `Normal(μ, γ)`.
- Observations are modeled jointly with Gaussian observation noise.

# Example
```julia
julia> model = OrnsteinUhlenbeckJumpModel(;
           group_indices = group_indices,
           subject_indices = subject_indices,
           z_obs = z_obs,
           Δt = Δt,
           M⁺ = M⁺,
           M⁻ = M⁻,
           Δτ = Δτ,
           is_transition = is_transition
       )
```
"""
@model function OrnsteinUhlenbeckJumpModel(;
        group_indices, subject_indices, z_obs, Δt, M⁺, M⁻, Δτ, is_transition, kwargs...
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
    z_obs ~ MvNormal(z, σ²_obs * I)
end

"""
    clean_parameters(::typeof(OrnsteinUhlenbeckJumpModel), chains, args)

Clean posterior samples from `OrnsteinUhlenbeckJumpModel`.

The method reconstructs subject-level Ornstein-Uhlenbeck and jump-response
parameters, derives diffusion scales, and converts parameters to interpretable
units.

# Arguments
- `model::typeof(OrnsteinUhlenbeckJumpModel)`: The model function used for
  dispatch.
- `chains`: Posterior samples from `OrnsteinUhlenbeckJumpModel`.
- `args`: Model arguments containing `group_indices`, with one group index per
  subject.

# Returns
- `DataFrame`: Posterior samples containing cleaned group-level and
  subject-level Ornstein-Uhlenbeck and jump-response parameters.

# Details
- Reconstructs subject-level `μ`, `θ`, `γ`, `α⁺`, and `α⁻` using the model's
  non-centered hierarchical parameterization.
- Converts `μ` to percentage points.
- Converts log mean reversion rates `θ` to half-lives using
  `log(2) / exp(θ)`.
- Converts log stationary standard deviations `γ` to percentage points.
- Exponentiates `α⁺` and `α⁻` to obtain positive jump-response coefficients.
- Derives diffusion scales as `100 * exp(γ) * sqrt(2 * exp(θ))`.
- Transforms raw group-level variation parameters `κ` to their positive model
  scales, with `κ_μ_group` expressed in percentage points.
- Retains the cleaned `μ`, `θ`, `γ`, `σ`, `α⁺`, `α⁻`, and `κ` parameters for
  all groups and subjects.
"""
function clean_parameters(::typeof(OrnsteinUhlenbeckJumpModel), chains, args)
    n_groups = maximum(args.group_indices)
    n_subjects = length(args.group_indices)

    param(name, i) = Symbol("$name[$i]")
    κ(raw) = exp(-2.30259 + 0.5 * raw)

    subject_rules = mapreduce(vcat, 1:n_subjects) do s
        g = args.group_indices[s]

        [
            [
                param("μ_group", g),
                param("κ_μ_group_raw", g),
                param("μ_subject_raw", s)
            ] => ByRow((μ, κ_raw, z) -> μ + κ(κ_raw) * z) => param("μ_subject", s),
            [
                param("θ_group", g),
                param("κ_θ_group_raw", g),
                param("θ_subject_raw", s)
            ] => ByRow((θ, κ_raw, z) -> θ + κ(κ_raw) * z) => param("θ_subject", s),
            [
                param("γ_group", g),
                param("κ_γ_group_raw", g),
                param("γ_subject_raw", s)
            ] => ByRow((γ, κ_raw, z) -> γ + κ(κ_raw) * z) => param("γ_subject", s),
            [
                param("α⁺_group", g),
                param("κ_α⁺_group_raw", g),
                param("α⁺_subject_raw", s)
            ] => ByRow((α⁺, κ_raw, z) -> α⁺ + κ(κ_raw) * z) => param("α⁺_subject", s),
            [
                param("α⁻_group", g),
                param("κ_α⁻_group_raw", g),
                param("α⁻_subject_raw", s)
            ] => ByRow((α⁻, κ_raw, z) -> α⁻ + κ(κ_raw) * z) => param("α⁻_subject", s)
        ]
    end

    function σ_rules(level, n)
        [[param("θ_$level", i), param("γ_$level", i)] =>
             ByRow((θ, γ) -> 100 * exp(γ) * sqrt(2 * exp(θ))) =>
                 param("σ_$level", i)
         for i in 1:n]
    end

    function κ_rules(name, scale = identity)
        [param("κ_$(name)_group_raw", g) =>
             ByRow(raw -> scale(κ(raw))) =>
                 param("κ_$(name)_group", g)
         for g in 1:n_groups]
    end

    columns(name, level, n) = [param("$(name)_$(level)", i) for i in 1:n]

    μ_group = columns("μ", "group", n_groups)
    θ_group = columns("θ", "group", n_groups)
    γ_group = columns("γ", "group", n_groups)
    σ_group = columns("σ", "group", n_groups)
    α⁺_group = columns("α⁺", "group", n_groups)
    α⁻_group = columns("α⁻", "group", n_groups)

    κ_μ_group = columns("κ_μ", "group", n_groups)
    κ_θ_group = columns("κ_θ", "group", n_groups)
    κ_γ_group = columns("κ_γ", "group", n_groups)
    κ_α⁺_group = columns("κ_α⁺", "group", n_groups)
    κ_α⁻_group = columns("κ_α⁻", "group", n_groups)

    μ_subject = columns("μ", "subject", n_subjects)
    θ_subject = columns("θ", "subject", n_subjects)
    γ_subject = columns("γ", "subject", n_subjects)
    σ_subject = columns("σ", "subject", n_subjects)
    α⁺_subject = columns("α⁺", "subject", n_subjects)
    α⁻_subject = columns("α⁻", "subject", n_subjects)

    parameter_columns = vcat(
        μ_group,
        θ_group,
        γ_group,
        σ_group,
        α⁺_group,
        α⁻_group,
        κ_μ_group,
        κ_θ_group,
        κ_γ_group,
        κ_α⁺_group,
        κ_α⁻_group,
        μ_subject,
        θ_subject,
        γ_subject,
        σ_subject,
        α⁺_subject,
        α⁻_subject
    )

    @chain chains begin
        # convert Chains object to DataFrame, where the parameters are columns and the iterations are rows
        DataFrame

        # construct subject-level parameters on their model scales
        transform(subject_rules...)

        # compute diffusion scales before transforming θ and γ
        transform(
            σ_rules("group", n_groups)...,
            σ_rules("subject", n_subjects)...
        )

        # transform group-level between-subject standard deviation
        transform(
            κ_rules("μ", x -> 100 * x)...,
            κ_rules("θ")...,
            κ_rules("γ")...,
            κ_rules("α⁺")...,
            κ_rules("α⁻")...
        )

        # convert group- and subject-level parameters to interpretable scales
        # μ: percentage points
        # θ: half-life
        # γ: percentage points
        # α⁺: percentage points
        # α⁻: percentage points
        transform(
            μ_group .=> ByRow(μ -> 100 * μ),
            μ_subject .=> ByRow(μ -> 100 * μ),
            θ_group .=> ByRow(θ -> log(2) / exp(θ)),
            θ_subject .=> ByRow(θ -> log(2) / exp(θ)),
            γ_group .=> ByRow(γ -> 100 * exp(γ)),
            γ_subject .=> ByRow(γ -> 100 * exp(γ)),
            α⁺_group .=> ByRow(α⁺ -> exp(α⁺)),
            α⁺_subject .=> ByRow(α⁺ -> exp(α⁺)),
            α⁻_group .=> ByRow(α⁻ -> exp(α⁻)),
            α⁻_subject .=> ByRow(α⁻ -> exp(α⁻));
            renamecols = false
        )

        # retain only interpretable parameters
        select(parameter_columns)
    end
end