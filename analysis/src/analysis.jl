
function load_parameters(model, folder, variable, args)
    @chain begin
        make_filename(folder, model, variable)
        load("chains")

        clean_parameters(model, _, args)
        summarize_parameters

        subset(
            :Parameter => ByRow(!endswith("subject")),
            :Parameter => ByRow(!startswith("γ")),
            :Parameter => ByRow(!startswith("κ"))
        )
        transform(
            All() => ((x...) -> string(model)) => :Model,
            All() => ((x...) -> string(variable)) => :Variable,
            :Parameter => (x -> replace(x, DESCRIPTIONS...)),
            :Index => ByRow(x -> args.groups[x]) => :Group,
            :Index => (x -> replace(x, 1 => 3, 2 => 1, 3 => 2));
            renamecols = false
        )
        sort(:Index)
    end
end

function load_comparisons(model, folder, variable, args)
    @chain begin
        make_filename(folder, model, variable)
        load("chains")

        clean_parameters(model, _, args)
        compare_groups

        subset(
            :Parameter => ByRow(!startswith("γ")),
            :Parameter => ByRow(!startswith("κ"))
        )
        transform(
            All() => ((x...) -> string(model)) => :Model,
            All() => ((x...) -> string(variable)) => :Variable,
            :Parameter => (x -> replace(x, DESCRIPTIONS...)),
            [:IndexA, :IndexB] .=> ByRow(x -> args.groups[x]) .=> [:GroupA, :GroupB];
            renamecols = false
        )
        transform([:GroupA, :GroupB] => ByRow((a, b) -> a * " vs " * b) => :Comparison)
        transform(
            :Comparison =>
                (x -> replace(x,
                    "ANX vs BPD" => "BPD vs ANX",
                    "ANX vs BPD-REM" => "BPD-REM vs ANX"
                ));
            renamecols = false
        )
    end
end