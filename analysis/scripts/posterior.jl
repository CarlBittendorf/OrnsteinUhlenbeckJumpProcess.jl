include("../startup.jl")

using AlgebraOfGraphics: density

folder = joinpath("models", "real")

df_data = d"Dataset"

model = OrnsteinUhlenbeckJumpModel
variable = :Valence

args = prepare_args(df_data, variable)
predictive_model = model(; merge(args, (; z_obs = missing))...)

chains = @chain begin
    make_filename(folder, model, variable)
    load("chains")
end

Random.seed!(1)

predictive_samples = predict(
    predictive_model,
    chains[rand(1:2500, 10), Colon(), rand(1:4, 10)]
)

z_pred = @chain predictive_samples begin
    DataFrame
    select(Not([:iteration, :chain]))
    Matrix
    transpose
end

# calculate RMSE
sqrt(mean(abs2, 100 .* (z_pred .- args.z_obs)))

df_model = @chain df_data begin
    sort([:Group, :Participant])

    groupby(:Participant)
    transform(
        groupindices => :Subject,
        eachindex => :Prompts;
        renamecols = false
    )
end

df_figure = vcat(
    transform(
        df_model,
        :Valence => (x -> "Original") => :Data,
        :Valence => (x -> 100 * x);
        renamecols = false
    ),
    (transform(
         df_model,
         :Valence => (x -> "Predicted") => :Data,
         :Valence => (x -> i) => :Iteration,
         :Valence => (x -> 100 * z_pred[:, i]);
         renamecols = false
     ) for i in 1:5)...;
    cols = :union
)

figure = @chain df_figure begin
    subset(:Participant => ByRow(x -> x in [181080, 187005, 182167, 184013]))
    transform(
        :Participant => (x -> replace(x,
            181080 => "BPD",
            184013 => "BPD-REM",
            187005 => "ANX",
            182167 => "HC"
        ));
        renamecols = false
    )

    draw(
        mapping(
            :Prompts, :Valence;
            row = :Participant => sorter("BPD", "BPD-REM", "ANX", "HC")
        ) * visual(Lines) *
        (
            data(subset(_, :Data => ByRow(isequal("Original")))) *
            mapping(; group = :Iteration, color = :Data) +
            data(subset(_, :Data => ByRow(isequal("Predicted")))) *
            mapping(; col = :Iteration, color = :Data)
        );
        axis = (; width = 250, height = 250)
    )
end

save(
    joinpath("figures", "paper", "Posterior Predictive Checks.png"), figure;
    px_per_unit = 3
)