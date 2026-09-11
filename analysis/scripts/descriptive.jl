include("../startup.jl")

####################################################################################################
# SAMPLE CHARACTERISTICS
####################################################################################################

df_summary = @chain d"Dataset" begin
    groupby(:Participant)
    combine(
        :Group => first,
        :FormStart => (x -> 100 * length(x) / 48) => :Compliance;
        renamecols = false
    )

    # add :Age and :Medication columns
    leftjoin(d"BPD Affect & Self-Esteem Questionnaires"; on = :Participant)

    transform(:Medication => ByRow(x -> !ismissing(x) && x); renamecols = false)

    sort(:Group; by = x -> findfirst(isequal(x), ["BPD", "BPD-REM", "ANX", "HC"]))
end

@chain df_summary begin
    table_one(
        [
            :Age => extrema_analysis => "Age [years]",
            :Medication => boolean_analysis => "Medication",
            :Compliance => mean_analysis => "Compliance [%]"
        ];
        groupby = :Group => "Group",
        show_n = true,
        sort = false
    )
    print_latex_table
end