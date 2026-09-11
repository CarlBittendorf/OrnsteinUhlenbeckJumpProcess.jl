
function make_filename(folder, model, variable)
    joinpath(folder, string(model), string(variable) * ".jld2")
end

enumerate_days(x) = Dates.value.(Date.(x) .- Date(minimum(x))) .+ 1

function hat(x)
    if length(x) == 1
        return x * '̂'
    else
        return join([first(x), '̂', x[collect(eachindex(x))[2:end]]...])
    end
end

mean_analysis(column) = (Concat(mean(column), " (", std(column), ")") => "Mean (SD)",)

function extrema_analysis(column)
    (
        Concat(mean(column), " (", std(column), ")") => "Mean (SD)",
        Concat("[", minimum(column), ", ", maximum(column), "]") => "[Min, Max]"
    )
end

function boolean_analysis(column)
    (Concat(count(column), " (", 100 * mean(column), "%)") => "n (%)",)
end

function print_latex_table(table)
    io = IOBuffer()

    show(io, MIME"text/latex"(), table)

    @chain io begin
        seekstart
        read(String)
        replace(raw"\textbackslash{}" => "\\", raw"\$" => raw"$", raw"\_" => "_")
        println
    end
end