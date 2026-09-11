
####################################################################################################
# EXPORTS
####################################################################################################

export clean_parameters, summarize_parameters, posterior_directional_error, compare_groups

####################################################################################################
# IMPLEMENTATIONS
####################################################################################################

"""
    clean_parameters(model, chains, args)

Convert posterior samples into a data frame of interpretable model parameters.

Implementations reconstruct subject-level parameters from their hierarchical
parameterization, transform parameters to meaningful scales, and retain only
the columns relevant to interpretation and reporting.

# Arguments
- `model`: The model function used to generate `chains`.
- `chains`: Posterior samples, typically represented as an `MCMCChains.Chains`
  object.
- `args`: Model arguments containing the indexing information required to
  reconstruct hierarchical parameters.

# Returns
- `DataFrame`: Posterior samples with one row per iteration and one column per
  cleaned group-level or subject-level parameter.

# Details
- Dispatches on `model` to apply model-specific transformations.
- Reconstructs subject-level parameters from group-level parameters,
  between-subject scales, and standardized subject-level effects.
- Removes raw and latent parameters that are not intended for direct
  interpretation.
"""
function clean_parameters end

"""
    summarize_parameters(df)

Summarize posterior samples using medians and highest posterior density
intervals.

The function computes the median and highest posterior density interval for
each parameter, then separates indexed parameter names into parameter and
index columns.

# Arguments
- `df::AbstractDataFrame`: A data frame whose columns contain posterior samples
  for individual parameters. Indexed parameter names should use the form
  `"parameter[index]"`.

# Returns
- `DataFrame`: A summary with the columns `:Parameter`, `:Index`, `:Median`,
  `:Lower`, and `:Upper`.

# Details
- Converts `df` to a `Chains` object to calculate highest posterior density
  intervals with `hpd`.
- Calculates the median of each parameter column.
- Removes parameters for which no complete summary is available.
- Parses parameter names such as `"beta[2]"` into `"beta"` and the integer
  index `2`.

# Example
```julia
julia> df = DataFrame(
           Symbol("beta[1]") => randn(1_000),
           Symbol("beta[2]") => randn(1_000)
       );

julia> summarize_parameters(df)
2×5 DataFrame
 Row │ Parameter  Index  Median     Lower      Upper
     │ String     Int64  Float64    Float64    Float64
─────┼──────────────────────────────────────────────────
   1 │ beta           1  0.012      -1.95       1.94
   2 │ beta           2  0.008      -1.91       1.97
```
"""
function summarize_parameters(df)
    df_hpd = @chain df begin
        # convert DataFrame back to Chains
        Chains(Array(_), names(_))

        # calculate the highest posterior density interval for each parameter
        hpd

        DataFrame

        rename([:Parameter, :Lower, :Upper])
        transform(:Parameter => ByRow(string); renamecols = false)
    end

    @chain df begin
        # calculate the median for each parameter
        combine(All() .=> median; renamecols = false)

        # stack the DataFrame so that each parameter and its median correspond to one row
        stack
        rename([:Parameter, :Median])

        # add highest posterior density intervals
        leftjoin(df_hpd; on = :Parameter)
        dropmissing

        # separate parameter names and indices into two columns
        transform(
            :Parameter =>
            ByRow(x -> split(x, ['[', ']'])[1:2]) => [:Parameter, :Index]
        )
        transform(
            :Parameter => ByRow(string),
            :Index => ByRow(x -> parse(Int, x));
            renamecols = false
        )
    end
end

"""
    posterior_directional_error(a, b)

Estimate the posterior probability that the inferred direction of the
difference between two quantities is incorrect.

The function calculates the posterior proportions for `a < b` and `a > b`,
then returns the smaller of the two.

# Arguments
- `a`: Posterior samples for the first quantity.
- `b`: Posterior samples for the second quantity. Must be broadcast-compatible
  with `a`.

# Returns
- A value between `0` and `0.5` representing the posterior directional error.
  Values near `0` indicate strong evidence for one direction, while values near
  `0.5` indicate high directional uncertainty.

# Details
- Comparisons are performed elementwise using broadcasting.
- Samples where `a == b` contribute to neither directional probability.

# Example
```julia
julia> a = [1.0, 2.0, 3.0, 4.0];

julia> b = [2.0, 3.0, 2.0, 5.0];

julia> posterior_directional_error(a, b)
0.25
```
"""
posterior_directional_error(a, b) = min(mean(a .< b), mean(a .> b))

"""
    compare_groups(df)

Compare group-level posterior distributions using posterior directional errors.

The function identifies group-level parameters, extracts the posterior samples
for each group, and calculates the posterior directional error for every
pairwise group comparison within each parameter.

# Arguments
- `df::AbstractDataFrame`: A data frame whose columns contain posterior samples.
  Group-level columns must include `"group"` in their names and use the naming
  convention `"parameter[index]"`.

# Returns
- `DataFrame`: A pairwise comparison table with the columns `:Parameter`,
  `:IndexA`, `:IndexB`, and `:PosteriorDirectionalError`.

# Details
- Selects only columns whose names contain `"group"`.
- Parses names such as `"group_effect[2]"` into a parameter name and group
  index.
- Groups samples by parameter and generates all unique pairs of groups.
- Uses `posterior_directional_error` to quantify the directional uncertainty
  for each pair.
- `IndexA` and `IndexB` refer to the positions of the groups within each
  parameter, not necessarily the parsed group indices.

# Example
```julia
julia> df = DataFrame(
           Symbol("group_effect[1]") => randn(1_000),
           Symbol("group_effect[2]") => randn(1_000) .+ 1,
           Symbol("group_effect[3]") => randn(1_000) .+ 2
       );

julia> compare_groups(df)
3×4 DataFrame
 Row │ Parameter     IndexA  IndexB  PosteriorDirectionalError
     │ String        Int64   Int64   Float64
─────┼──────────────────────────────────────────────────────────
   1 │ group_effect       1       2                      0.24
   2 │ group_effect       1       3                      0.08
   3 │ group_effect       2       3                      0.23
```
"""
function compare_groups(df)
    @chain df begin
        # use only group-level parameters
        select(Cols(x -> contains(x, "group")))

        # collapse all data points into a single row
        combine(All() .=> (x -> [x]); renamecols = false)

        # transpose the DataFrame so that each parameter and its associated chains correspond to one row
        DataFrame(:Parameter => names(_), :Chains => first.(eachcol(_)))

        # separate parameter names and indices into two columns
        transform(
            :Parameter =>
            ByRow(x -> split(x, ['[', ']'])[1:2]) => [:Parameter, :Index]
        )
        transform(
            :Parameter => ByRow(string),
            :Index => ByRow(x -> parse(Int, x));
            renamecols = false
        )

        # calculate posterior directional errors for each combination of two groups
        groupby(:Parameter)
        combine(
            :Chains => (x -> collect(combinations(eachindex(x), 2))) => [:IndexA, :IndexB],
            :Chains =>
                (x -> [posterior_directional_error(x[a], x[b])
                       for (a, b) in combinations(eachindex(x), 2)]) =>
                    :PosteriorDirectionalError
        )
    end
end