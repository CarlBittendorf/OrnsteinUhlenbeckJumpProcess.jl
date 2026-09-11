
function prepare_args(df::DataFrame)
    df_model = @chain df begin
        sort([:Group, :Participant])

        groupby(:Participant)
        transform(
            groupindices => :Subject;
            renamecols = false
        )
    end

    groups = @chain df_model begin
        groupby(:Group)
        combine(groupindices => :GroupIndex)

        getproperty(:Group)
    end

    group_indices = @chain df_model begin
        groupby(:Group)
        transform(groupindices => :GroupIndex)

        groupby(:Subject)
        combine(:GroupIndex => first; renamecols = false)

        getproperty(:GroupIndex)
    end

    subject_indices = df_model.Subject
    Δt = df_model.HoursSinceObservation
    M⁺ = df_model.PositiveEventIntensity
    M⁻ = df_model.NegativeEventIntensity
    Δτ = df_model.HoursSinceEvent

    first_indices = findall(diff([0; subject_indices]) .!= 0)
    is_transition = trues(length(Δt))
    is_transition[first_indices] .= false

    return (; groups, group_indices, subject_indices, Δt, M⁺, M⁻, Δτ, is_transition)
end

function prepare_args(df::DataFrame, variable::Symbol)
    z_obs = @chain df begin
        sort([:Group, :Participant])

        groupby(:Participant)
        transform(
            groupindices => :Subject;
            renamecols = false
        )

        getproperty(variable)
    end

    return merge(prepare_args(df), (; z_obs))
end

function fit_model(
        f::Function, filename::AbstractString, args;
        iterations = 2500, nchains = 4
)
    rng = StableRNG(1)
    model = f(; args...)

    chains = sample(
        rng,
        model,
        NUTS(; adtype = AutoReverseDiff(; compile = true)),
        MCMCThreads(),
        iterations,
        nchains
    )

    mkpath(dirname(filename))
    jldsave(filename; chains)

    @info "Saved model to " * filename
end