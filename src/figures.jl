
####################################################################################################
# EXPORTS
####################################################################################################

export draw_chains

####################################################################################################
# IMPLEMENTATIONS
####################################################################################################

"""
    draw_chains(chains, indices)
    draw_chains(chains, n::Int = 10)

Create trace and density plots for selected parameters of an MCMC chain.

For each parameter, the function displays a trace plot across iterations and a
density plot, with separate colors for each chain.

# Arguments
- `chains`: A chain object containing parameter names and sampled values.
- `indices`: The indices of the parameters to plot.
- `n::Int`: The number of parameters to plot, starting from the first. Defaults
  to `10`.

# Returns
- `Figure`: A figure with one row per parameter, containing a trace plot and a
  density plot.

# Details
- Parameter names are obtained from `chains.name_map.parameters`.
- The first column shows sampled values across iterations.
- The second column shows the estimated density of the sampled values.
- Colors distinguish the individual chains.

# Example
```julia
julia> figure = draw_chains(chains, 10:20)

julia> figure = draw_chains(chains, 5)
```
"""
function draw_chains(chains, indices)
    params = string.(chains.name_map.parameters)[indices]
    iterations = size(chains.value, 1)
    n = size(chains.value, 3)

    figure = Figure(; size = (800, 250 * length(params)))

    for (i, (index, label)) in enumerate(zip(indices, params))
        df = DataFrame(
            :X => repeat(collect(1:iterations), n),
            :Y => vec(chains.value[:, index, :]),
            :Chain => repeat(string.(collect(1:n)); inner = iterations)
        )

        draw!(
            figure[i, 1],
            data(df) * mapping(:X => "", :Y => ""; color = :Chain) * visual(Lines);
            axis = (; title = label)
        )
        draw!(
            figure[i, 2],
            data(df) * mapping(:Y => ""; color = :Chain) * AlgebraOfGraphics.density();
            axis = (; title = label, ylabel = "")
        )
    end

    return figure
end

draw_chains(chains, n::Int = 10) = draw_chains(chains, 1:n)