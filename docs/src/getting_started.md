# Getting started

`OrnsteinUhlenbeckJumpProcess` provides hierarchical Bayesian models for
longitudinal data and event probabilities, built with
[Turing.jl](https://turinglang.org). Fitting a model produces posterior samples
for group-level and subject-level parameters, which the package's analysis
functions convert into interpretable quantities such as equilibrium means,
deviation half-lives, and event reactivity.

This page walks through the core workflow with synthetic data: simulating
observations, choosing a model, sampling with NUTS, and cleaning the posterior
samples.

## Installation

The package is not yet registered in the Julia General registry. Install it
directly from GitHub:

```julia
using Pkg

Pkg.add(url = "https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl")
```

If you are working from a local clone, the package can also be developed
directly from its path:

```julia
pkg> develop /path/to/OrnsteinUhlenbeckJumpProcess
```

## Data conventions

All models operate on **nodes**: each element of the data vectors corresponds
to one observation node. Nodes that follow their predecessor in time form a
time series, while the first node of each series is drawn from the stationary
distribution. Three conventions apply to every model:

- Observations live on the **unit interval** (0, 1), regardless of the
  original rating scale.
- Time differences are measured in **hours** and are computed from the
  preceding node of the same subject.
- Event magnitudes (jump model) are scaled to $[0, 0.01]$, so that the
  reactivity parameters can be read as percentage points of affect change.

See [Models](@ref) for the full conventions and the likelihood of each model.

## A first synthetic example

We simulate two subjects, each observed on an 8-hour grid over one day. The
`hours_elapsed` utility converts timestamps into the required `Δt` vector:

```@example getting-started
using OrnsteinUhlenbeckJumpProcess, Dates

times = DateTime("2026-01-01T00:00:00"):Hour(8):DateTime("2026-01-01T16:00:00")
Δt = repeat(hours_elapsed(collect(times)), 2)
is_transition = repeat([false, true, true], 2)
```

Next, we draw observations from a data-generating model with the default
priors. Setting `z_obs = missing` marks the observations as unknown, so that
`rand` samples them along with all other random variables:

```@example getting-started
using StableRNGs

data = (;
    group_indices = [1, 1],
    subject_indices = [1, 1, 1, 2, 2, 2],
    Δt = Δt,
    is_transition = is_transition
)
generative = OrnsteinUhlenbeckModel(; data..., z_obs = missing)
z_obs = rand(StableRNG(1), generative).z_obs
```

The resulting values are the synthetic observations. We fix them in a model
that we can fit:

```@example getting-started
model = OrnsteinUhlenbeckModel(; data..., z_obs = z_obs);
```

## Choose a model

The package provides three models. Choose the one that matches the data:

| Model | Use it for |
| --- | --- |
| [`OrnsteinUhlenbeckModel`](@ref) | Longitudinal observations with continuous mean-reverting dynamics |
| [`OrnsteinUhlenbeckJumpModel`](@ref) | Longitudinal observations with measured positive and negative events |
| [`EventProbabilityModel`](@ref) | Binary event indicators, without a latent time series |

Each model is an ordinary Turing model, so the data-generating trick above
works for all of them: pass the observations as `missing` and call `rand` on
the model to sample from the joint prior.

## Fitting with NUTS

Any Turing-compatible sampler works. The default NUTS sampler is a good
starting point for small models:

```@example getting-started
using Turing

chains = sample(StableRNG(1), model, NUTS(), 200; progress = false);
size(chains)
```

For larger models, consider compiling the reverse-mode automatic
differentiation and running multiple chains in parallel; see the
[Analysis workflow](@ref) and [Jump model workflow](@ref) examples and the
README for the recommended setup.

## Cleaning the posterior samples

The chain contains raw, non-centered parameters. `clean_parameters`
reconstructs the interpretable subject-level parameters:

```@example getting-started
cleaned = clean_parameters(OrnsteinUhlenbeckModel, chains, data);
names(cleaned)
```

Finally, summarize the cleaned parameters with medians and 95% highest
posterior density intervals:

```@example getting-started
summarize_parameters(cleaned[:, [
    "μ_group[1]", "θ_group[1]", "γ_group[1]", "σ_group[1]",
    "μ_subject[1]", "θ_subject[1]", "γ_subject[1]", "σ_subject[1]",
    "μ_subject[2]", "θ_subject[2]", "γ_subject[2]", "σ_subject[2]"
]])
```

## Common configuration choices

The workflow usually changes in a few predictable ways:

- Increase the number of samples, or use `MCMCThreads()` with several parallel
  chains, to stabilize the posterior summaries.
- Compile the automatic differentiation (`AutoReverseDiff(; compile = true)`)
  when the number of nodes grows.
- Pass the event data (`M⁺`, `M⁻`, `Δτ`) to
  [`OrnsteinUhlenbeckJumpModel`](@ref) to include measured events.
- Keep the time unit in hours, so that the reported half-lives are in hours.

## Common mistakes

A few patterns usually lead to confusion:

- Passing observations outside the unit interval, which breaks the prior
  scales and the interpretation of the cleaned parameters.
- Forgetting that `Δt` must be `0` at the first node of each series and is
  measured from the preceding node of the same subject.
- Marking the first node of a series with `is_transition = true`, which
  references a node that has no predecessor.
- Reading the raw chain columns (`*_raw`, `z_raw`) directly instead of using
  `clean_parameters`.
- Mixing up the scales: `μ` and `γ` are reported on a 0-100 percentage scale,
  `θ` as a half-life in hours, and `α⁺`, `α⁻` as percentage points per unit of
  impulse magnitude.

## Next steps

The rest of the documentation expands on the same ideas with a bit more
detail:

- [Models](@ref): data conventions, likelihoods, and priors of each model.
- [Analysis workflow](@ref): group comparisons and chain visualization.
- [Jump model workflow](@ref): fitting the model with event jumps.
- [Extending the package](@ref): adding new stochastic-process models.
- [API reference](@ref): complete reference of the public interface.
