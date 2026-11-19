```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "OrnsteinUhlenbeckJumpProcess.jl"
  text:
  tagline: Hierarchical Ornstein-Uhlenbeck process with event-based jumps for modeling affect dynamics
  actions:
    - theme: brand
      text: Getting started
      link: /getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl
---
```

# Welcome to OrnsteinUhlenbeckJumpProcess!

Hierarchical Bayesian models for longitudinal data and event probabilities,
built with [Turing.jl](https://turinglang.org). The package provides:

- **OrnsteinUhlenbeckModel**: a hierarchical Ornstein-Uhlenbeck process for
  longitudinal observations.
- **OrnsteinUhlenbeckJumpModel**: the Ornstein-Uhlenbeck process extended with
  positive and negative impulsive jumps.
- **EventProbabilityModel**: a hierarchical logistic model for binary events.
- Utilities for preparing data, cleaning and summarizing posterior samples,
  comparing groups, and visualizing MCMC chains.

## Why an Ornstein-Uhlenbeck process with jumps?

Ambulatory assessment studies often quantify affective instability from
intensive longitudinal data without separating responses to measured events
from the continuous background dynamics. The
[`OrnsteinUhlenbeckJumpModel`](@ref) addresses this by decomposing affective
dynamics into four interpretable components:

- an **attractor**: the equilibrium state to which affect reverts;
- the **half-life of deviations**: how quickly affect returns to the attractor;
- **volatility**: how much affect fluctuates around the attractor;
- **event reactivity**: separate positive and negative response coefficients
  for measured events.

The process is formulated in continuous time with a closed-form transition
distribution, so irregularly spaced observations require no resampling or
interpolation, and events that occur *between* two measurements still enter the
likelihood through their exponentially decaying impulse. Because event
responses are modeled separately from the continuous dynamics, the estimated
half-life and volatility reflect background affect dynamics that are not
confounded with reactions to measured events.

## Installation

The package is not yet registered in the Julia General registry. Install it
directly from GitHub:

```julia
using Pkg

Pkg.add(url = "https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl")
```

## Quick start

The models operate on **nodes**: each element of the data vectors corresponds
to one observation node, time differences are measured in **hours**, and
observations live on the **unit interval**. A quick tour of the core workflow:

```@example index
using OrnsteinUhlenbeckJumpProcess, Dates, StableRNGs

data = (;
    group_indices = [1, 1],
    subject_indices = [1, 1, 1, 2, 2, 2],
    Δt = repeat(hours_elapsed(collect(
        DateTime("2026-01-01T00:00:00"):Hour(8):DateTime("2026-01-01T16:00:00")
    )), 2),
    is_transition = repeat([false, true, true], 2)
)

# draw synthetic observations from the prior
generative = OrnsteinUhlenbeckModel(; data..., z_obs = missing)
z_obs = rand(StableRNG(1), generative).z_obs

# fix the observations in a model that can be fitted
model = OrnsteinUhlenbeckModel(; data..., z_obs = z_obs)

round.(z_obs; digits = 3)
```

## Documentation

- [Getting started](@ref): install the package and fit a first model.
- [Models](@ref): the data conventions, likelihoods, and priors of each model.
- [Analysis workflow](@ref) and [Jump model workflow](@ref): end-to-end
  examples that fit models and analyze the posterior samples.
- [Extending the package](@ref): add new stochastic-process models following
  the package conventions.
- [API reference](@ref): complete reference of the public interface.

## Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research
Foundation) – GRK2739/2 – Project Nr. 447089431 – Research Training Group:
KD²School – Designing Biosignal-Adaptive Systems for Decision-Making Processes

```@meta
CurrentModule = OrnsteinUhlenbeckJumpProcess
```

