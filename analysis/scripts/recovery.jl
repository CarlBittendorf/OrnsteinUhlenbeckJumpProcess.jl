include("../startup.jl")
include("simulation.jl")

####################################################################################################
# PARAMETER RECOVERY
####################################################################################################

const descriptions = [
    "μ_group" => "Attractor",
    "θ_group" => "Half-Life of Deviations",
    "γ_group" => "Variability",
    "α⁺_group" => "Reactivity (Positive Events)",
    "α⁻_group" => "Reactivity (Negative Events)"
]

df_params = DataFrame()

for j in eachindex(MODELS)
    folder = joinpath("models", "simulation", string(MODELS[j]))

    for model in MODELS
        df_model = DataFrame()

        for i in eachindex(args)
            filename = joinpath(folder, string(model) * " ($i).jld2")

            df_iteration = @chain filename begin
                load("chains")
                DataFrame

                # use only group-level parameters
                select(Cols(x -> contains(x, "group")))

                # calculate the median for each parameter
                combine(All() .=> median; renamecols = false)

                # stack the DataFrame so that each parameter and its median correspond to one row
                stack
                rename([:Parameter, :RecoveredValue])

                # separate parameter names from indices
                transform(
                    :Parameter => ByRow(x -> split(x, ['[', ']'])[1]);
                    renamecols = false
                )

                transform(
                    :Parameter => ByRow(string),
                    All() => ((x...) -> string(MODELS[j])) => :ModelSimulation,
                    All() => ((x...) -> string(model)) => :ModelFit,
                    All() => ((x...) -> i) => :Iteration;
                    renamecols = false
                )
            end

            df_model = vcat(df_model, df_iteration)
        end

        df_combined = @chain params[:, j] begin
            DataFrame

            # use only group-level parameters
            select(Cols(x -> contains(x, "group")))

            # turn the entries of group-level parameters into numbers
            transform(
                All() .=> ByRow(first),
                eachindex => :Iteration;
                renamecols = false
            )

            # stack the DataFrame so that each parameter and its median correspond to one row
            stack
            rename(
                :variable => :Parameter,
                :value => :TrueValue
            )

            rightjoin(df_model; on = [:Parameter, :Iteration])

            subset(:Parameter => ByRow(!startswith("κ")))

            transform(
                [:ModelSimulation, :ModelFit] .=>
                    ByRow(x -> x == "OrnsteinUhlenbeckModel" ? "OU" : "OU-J"),
                :Parameter => (x -> replace(x, descriptions...));
                renamecols = false
            )

            # transform the parameters to an interpretable scale
            transform(
                [
                [:Parameter, :TrueValue],
                [:Parameter, :RecoveredValue]
            ] .=>
                ByRow(
                    (
                    p,
                    x
                ) -> p == "Attractor" ? 100 * x :
                     p == "Half-Life of Deviations" ? log(2) / exp(x) :
                     p == "Variability" ? 100 * exp(x) :
                     exp(x)
                ) .=> [:TrueValue, :RecoveredValue]
            )
        end

        df_params = vcat(df_params, df_combined)
    end
end

# compare the recovered parameters to their true values using scatter plots
figure = @chain df_params begin
    subset([:ModelSimulation, :ModelFit] .=> ByRow(isequal("OU-J")))
    dropmissing

    draw(
        data(_) *
        mapping(
            :TrueValue => "True Value", :RecoveredValue => "Recovered Value";
            layout = :Parameter => presorted
        );
        facet = (; linkxaxes = :none, linkyaxes = :none),
        axis = (; width = 250, height = 250)
    )
end

save(
    joinpath("figures", "paper", "Parameter Recovery.png"), figure;
    px_per_unit = 3
)

# calculate the Pearson correlation for each parameter
df_within = @chain df_params begin
    subset([:ModelSimulation, :ModelFit] => ByRow(isequal))
    rename(:ModelSimulation => :Model)

    groupby([:Model, :Parameter])
    combine([:TrueValue, :RecoveredValue] => cor => :Correlation)
end

CSV.write(joinpath("data", "Parameter Recovery Within.csv"), df_within)

@chain df_within begin
    subset(:Model => ByRow(isequal("OU-J")))
    transform(
        :Correlation => ByRow(x -> replace(string(round(x; digits = 2)), "0." => "."));
        renamecols = false
    )
    select(:Parameter, :Correlation)

    simple_table
    print_latex_table
end

group_params = last.(descriptions)

# create a correlation matrix between recovered parameters
df_between = @chain df_params begin
    subset([:ModelSimulation, :ModelFit] .=> ByRow(isequal("OU-J")))

    unstack(:Iteration, :Parameter, :RecoveredValue)

    combine(
        All() => ((x...) -> first.(combinations(group_params, 2))) => :FirstParameter,
        All() => ((x...) -> last.(combinations(group_params, 2))) => :SecondParameter,
        group_params =>
            ((x...) -> [cor(x[i], x[j])
                        for (i, j) in combinations(eachindex(group_params), 2)]) =>
                :Correlation
    )

    groupby(:SecondParameter)
    combine(
        ([:FirstParameter, :Correlation] =>
             ((x, c) -> any(x .== p) ? only(c[x .== p]) : missing) => p
    for p in group_params)...
    )

    rename(:SecondParameter => :Parameter)
end

CSV.write("data/Parameter Recovery Between.csv", df_between)

@chain df_between begin
    transform(
        Not(:Parameter) .=> ByRow(x -> ismissing(x) ? "" :
                   replace(string(round(x; digits = 2)), "0." => "."));
        renamecols = false
    )
    select(1:5)

    simple_table
    print_latex_table
end

####################################################################################################
# MODEL RECOVERY
####################################################################################################

df_psisloo = CSV.read(joinpath("data", "Model Recovery PSIS-LOO.csv"), DataFrame)

df_confusion = @chain df_psisloo begin
    subset(
        :Statistic => ByRow(isequal("cv_elpd")),
        :Column => ByRow(isequal("total"))
    )

    transform(
        [:ModelSimulation, :ModelFit] .=>
            ByRow(x -> x == "OrnsteinUhlenbeckModel" ? "OU" : "OU-J");
        renamecols = false
    )

    groupby([:ModelSimulation, :Iteration])
    combine(
        [:ModelFit, :Value] =>
        ((m, d) -> first(m[d .== maximum(filter(!isnan, d))])) => :BestModel
    )

    groupby(:ModelSimulation)
    combine(
        (:BestModel => (m -> count(m .== model)) => model for model in ["OU", "OU-J"])...
    )
end

CSV.write(joinpath("data", "Model Recovery Confusion Matrix.csv"), df_confusion)

@chain df_confusion begin
    transform(
        ["OU", "OU-J"] .=>
            ByRow(x -> replace(string(round(x / 100; digits = 2)), "0." => "."));
        renamecols = false
    )
    rename("ModelSimulation" => "")

    simple_table
    print_latex_table
end