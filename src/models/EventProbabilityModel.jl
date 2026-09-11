
####################################################################################################
# EXPORTS
####################################################################################################

export EventProbabilityModel

####################################################################################################
# IMPLEMENTATIONS
####################################################################################################

"""
    EventProbabilityModel(;
        group_indices,
        subject_indices,
        e,
        kwargs...
    )

Define a hierarchical logistic model for binary event probabilities.

The model estimates group-level log-odds and subject-level deviations, then
models each observed event using a Bernoulli distribution with a logit
parameterization.

# Arguments
- `group_indices`: Group index for each subject.
- `subject_indices`: Subject index associated with each observation.
- `e`: Vector of binary event indicators.
- `kwargs...`: Additional keyword arguments reserved for model extensions.

# Returns
- A Turing model representing the hierarchical event-probability model.

# Details
- `β_group` represents the group-specific event log-odds.
- `κ_β_group` represents the positive group-specific scale of subject-level
  variation.
- `β_subject` uses a non-centered hierarchical parameterization.
- The event probability for observation `i` is
  `logistic(β_subject[subject_indices[i]])`.
- Each event is modeled independently using `BernoulliLogit`.

# Example
```julia
julia> model = EventProbabilityModel(;
           group_indices = group_indices,
           subject_indices = subject_indices,
           e = e
       )
```
"""
@model function EventProbabilityModel(; group_indices, subject_indices, e, kwargs...)
    n_groups = maximum(group_indices)
    n_subjects = length(group_indices)

    # group-level priors
    β_group ~ filldist(Normal(0, 1), n_groups)
    κ_β_group ~ filldist(truncated(Normal(0, 1); lower = 0), n_groups)

    # subject-level priors
    β_subject_raw ~ filldist(Normal(0, 1), n_subjects)

    @inbounds β_subject = β_group[group_indices] .+
                          κ_β_group[group_indices] .* β_subject_raw

    @inbounds β = β_subject[subject_indices]

    # observation model
    @inbounds for i in eachindex(e)
        e[i] ~ BernoulliLogit(β[i])
    end
end

"""
    clean_parameters(::typeof(EventProbabilityModel), chains, args)

Clean posterior samples from `EventProbabilityModel`.

The method reconstructs subject-level event log-odds and converts group-level
and subject-level effects to event probabilities.

# Arguments
- `model::typeof(EventProbabilityModel)`: The model function used for
  dispatch.
- `chains`: Posterior samples from `EventProbabilityModel`.
- `args`: Model arguments containing `group_indices`, with one group index per
  subject.

# Returns
- `DataFrame`: Posterior samples containing group-level event probabilities,
  group-level between-subject scales, and subject-level event probabilities.

# Details
- Reconstructs subject-level `β` from the corresponding group-level log-odds,
  between-subject scale, and standardized subject-level effect.
- Converts group-level and subject-level `β` values from log-odds to
  probabilities using the logistic function.
- Leaves `κ_β_group` on its positive log-odds scale.
- Retains `β_group`, `κ_β_group`, and `β_subject` for all groups and subjects.
"""
function clean_parameters(::typeof(EventProbabilityModel), chains, args)
    n_groups = maximum(args.group_indices)
    n_subjects = length(args.group_indices)

    param(name, i) = Symbol("$name[$i]")

    subject_rules = mapreduce(vcat, 1:n_subjects) do s
        g = args.group_indices[s]

        [
            [
            param("β_group", g),
            param("κ_β_group", g),
            param("β_subject_raw", s)
        ] => ByRow((β, κ, z) -> β + κ * z) => param("β_subject", s)
        ]
    end

    columns(name, level, n) = [param("$(name)_$(level)", i) for i in 1:n]

    β_group = columns("β", "group", n_groups)

    κ_β_group = columns("κ_β", "group", n_groups)

    β_subject = columns("β", "subject", n_subjects)

    parameter_columns = vcat(
        β_group,
        κ_β_group,
        β_subject
    )

    @chain chains begin
        # convert Chains object to DataFrame, where the parameters are columns and the iterations are rows
        DataFrame

        # construct subject-level parameters on their model scales
        transform(subject_rules...)

        # convert group- and subject-level parameters to interpretable scales
        # β: probability
        transform(
            β_group .=> ByRow(β -> logistic(β)),
            β_subject .=> ByRow(β -> logistic(β));
            renamecols = false
        )

        # retain only interpretable parameters
        select(parameter_columns)
    end
end