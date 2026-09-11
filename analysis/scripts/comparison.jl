include("../startup.jl")

folder = joinpath("models", "real")

df_data = d"Dataset"

args = prepare_args(df_data, :Valence)

####################################################################################################
# MEDIANS AND HPD INTERVALS
####################################################################################################

df_params = @chain begin
    vcat(
        load_parameters(OrnsteinUhlenbeckModel, folder, :Valence, args),
        load_parameters(OrnsteinUhlenbeckJumpModel, folder, :Valence, args)
    )

    transform(
        :Model => ByRow(x -> x == "OrnsteinUhlenbeckModel" ? "OU" : "OU-J");
        renamecols = false
    )
end

figure = @chain df_params begin
    transform(
        [:Lower, :Median] => ByRow((x, m) -> m - x) => :Lower,
        [:Upper, :Median] => ByRow((x, m) -> x - m) => :Upper,
        [:Model, :Parameter, :Index] =>
            ByRow((
                m,
                p,
                i
            ) -> contains(p, "Reactivity") ? convert(Float64, i) :
                 m == "OU" ? i - 0.1 : i + 0.1) => :X;
        renamecols = false
    )

    draw(
        data(_) *
        mapping(;
            layout = :Parameter => presorted,
            color = :Model => presorted
        ) *
        (
            mapping(:X, :Median, :Lower, :Upper) * visual(Errorbars; whiskerwidth = 5) +
            mapping(:X, :Median)
        );
        facet = (; linkxaxes = :none, linkyaxes = :none),
        axis = (;
            width = 250,
            height = 250,
            limits = ((0.75, 4.25), nothing),
            xlabel = "Group",
            ylabel = "",
            xtickformat = x -> unique(_.Group)
        )
    )
end

@chain df_params begin
    subset(:Model => ByRow(isequal("OU")))
    select(:Parameter, :Group, :Median, :Lower, :Upper)

    simple_table
    print_latex_table
end

@chain df_params begin
    subset(:Model => ByRow(isequal("OU-J")))
    select(:Parameter, :Group, :Median, :Lower, :Upper)

    simple_table
    print_latex_table
end

####################################################################################################
# POSTERIOR DIRECTIONAL ERRORS
####################################################################################################

df_comparisons = @chain begin
    vcat(
        load_comparisons(OrnsteinUhlenbeckModel, folder, :Valence, args),
        load_comparisons(OrnsteinUhlenbeckJumpModel, folder, :Valence, args)
    )

    select(:Parameter, :Comparison, :PosteriorDirectionalError, :Model)
    transform(
        :Model => ByRow(x -> x == "OrnsteinUhlenbeckModel" ? "OU" : "OU-J"),
        :PosteriorDirectionalError =>
            ByRow(x -> x < 0.001 ? "<0.001" : string(round(x; digits = 3)));
        renamecols = false
    )
    transform(
        :PosteriorDirectionalError => ByRow(x -> replace(x, "0." => "."));
        renamecols = false
    )

    sort(:Parameter; by = x -> findfirst(isequal(x), last.(DESCRIPTIONS)))
    sort(
        :Comparison;
        by = x -> findfirst(
            isequal(x),
            [
                "BPD vs BPD-REM", "BPD vs ANX", "BPD vs HC",
                "BPD-REM vs ANX", "BPD-REM vs HC", "ANX vs HC"
            ]
        )
    )
end

@chain df_comparisons begin
    unstack(:Comparison, :PosteriorDirectionalError)
    dropmissing

    CSV.write(joinpath("data", "Group Comparisons.csv"), _)
end

@chain df_comparisons begin
    subset(:Model => ByRow(isequal("OU-J")))
    select(Not(:Model))

    unstack(:Parameter, :PosteriorDirectionalError)

    simple_table
    print_latex_table
end

####################################################################################################
# RELATIVE PARAMETER CHANGES
####################################################################################################

@chain df_params begin
    subset(:Parameter => ByRow(!startswith("Reactivity")))

    groupby([:Parameter, :Group])
    combine(
        :Median =>
        (x -> round((x[2] - x[1]) / x[1] * 100; digits = 2)) => :RelativeChange
    )

    rename(:RelativeChange => "Relative Change [%]")

    simple_table
    print_latex_table
end

####################################################################################################
# EVENTS
####################################################################################################

df_figure = @chain begin
    vcat(
        load_parameters(EventProbabilityModel, folder, "PositiveEvents", args),
        load_parameters(EventProbabilityModel, folder, "NegativeEvents", args)
    )

    transform(
        :Variable =>
        (x -> replace(x,
            "PositiveEvents" => "Probability of a Positive Event",
            "NegativeEvents" => "Probability of a Negative Event"
        )) => :Model
    )
    transform(
        [:Lower, :Median] => ByRow((x, m) -> m - x) => :Lower,
        [:Upper, :Median] => ByRow((x, m) -> x - m) => :Upper
    )
end

draw!(
    figure.figure[3, 1:2],
    data(df_figure) *
    mapping(; layout = :Model => presorted) *
    (
        mapping(:Index, :Median, :Lower, :Upper) * visual(Errorbars; whiskerwidth = 5) +
        mapping(:Index, :Median)
    );
    facet = (; linkxaxes = :none, linkyaxes = :none),
    axis = (;
        width = 250,
        height = 250,
        limits = ((0.75, 4.25), nothing),
        xlabel = "Group",
        ylabel = "",
        xtickformat = x -> unique(df_figure.Group)
    )
)

resize!(figure.figure, 1000, 1025)

Label(figure.figure[1:2, 1:4, TopLeft()], "a", font = :bold, fontsize = 24, halign = :left)
Label(figure.figure[3, 1:2, TopLeft()], "b", font = :bold, fontsize = 24, halign = :left)

save(joinpath("figures", "paper", "Group-Level Parameters.png"), figure; px_per_unit = 3)

@chain begin
    vcat(
        load_comparisons(EventProbabilityModel, folder, "PositiveEvents", args),
        load_comparisons(EventProbabilityModel, folder, "NegativeEvents", args)
    )

    transform(
        :Variable =>
            (x -> replace(x,
                "PositiveEvents" => "Probability of a Positive Event",
                "NegativeEvents" => "Probability of a Negative Event"
            )) => :Parameter,
        :PosteriorDirectionalError =>
            ByRow(x -> x < 0.001 ? "<0.001" : string(round(x; digits = 3)));
        renamecols = false
    )
    transform(
        :PosteriorDirectionalError => ByRow(x -> replace(x, "0." => "."));
        renamecols = false
    )

    select(:Parameter, :Comparison, :PosteriorDirectionalError)

    sort(
        :Comparison;
        by = x -> findfirst(
            isequal(x),
            [
                "BPD vs BPD-REM", "BPD vs ANX", "BPD vs HC",
                "BPD-REM vs ANX", "BPD-REM vs HC", "ANX vs HC"
            ]
        )
    )

    unstack(:Parameter, :PosteriorDirectionalError)

    simple_table
    print_latex_table
end