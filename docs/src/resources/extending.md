# Extending the package

The package is organized around a small set of conventions, so additional
stochastic-process models can be added without reorganizing the existing code.
This page describes the pattern shared by all models; the [Models](@ref) pages
and the [API reference](@ref) show the resulting structure.

## Model conventions

New models follow the structure of the existing ones:

1. **Define the model as an exported Turing `@model` function.** The function
   name is the model name (for example `OrnsteinUhlenbeckModel`), data are
   passed as keyword arguments, and a trailing `kwargs...` keeps calls
   forward-compatible when new arguments are added to existing models.
2. **Use non-centered parameterizations.** Group-level parameters,
   between-subject scales, and standardized subject-level effects keep the
   posterior geometry friendly for gradient-based samplers.
3. **Follow the shared data conventions.** Observations live on the unit
   interval, time differences are measured in hours, and event magnitudes are
   scaled to $[0, 0.01]$; see [Models](@ref).
4. **Provide a `clean_parameters` method.** Dispatch on the model *function*,
   `clean_parameters(::typeof(NewModel), chains, args)`, and return a
   `DataFrame` with one column per interpretable parameter. Reconstruct
   subject-level parameters from their hierarchical parameterization and
   transform parameters to meaningful scales.
5. **Document and test.** Document the likelihood and priors on a page under
   `docs/src/models/`, add the page to `pages` in `docs/make.jl`, and add
   `logjoint`, reconstruction, and sampling tests to `test/runtests.jl`.

## Skeleton

```julia
export MyProcessModel

@model function MyProcessModel(; group_indices, subject_indices, z_obs, kwargs...)
    # group-level priors, between-subject scales, standardized
    # subject-level effects, and the observation model
    ...
end

function clean_parameters(::typeof(MyProcessModel), chains, args)
    param(name, i) = Symbol("$name[$i]")
    ...
    @chain chains begin
        # convert the chain to a DataFrame and reconstruct the
        # interpretable parameters
        DataFrame
        transform(subject_rules...)
        select(parameter_columns)
    end
end
```

The `param`, `columns`, `subject_rules`, and `κ_rules` helpers in the existing
`clean_parameters` methods show how indexed parameter names are constructed
and how subject-level parameters are reconstructed; see
[`OrnsteinUhlenbeckModel`](@ref) for a complete example.

## Fitting at scale

The repository also ships an auxiliary `analysis/` subpackage with scripts for
simulation and fitting at scale, including compiled reverse-mode automatic
differentiation. It is a good reference for running large models efficiently.
