# Context Matters: Extending the Ornstein–Uhlenbeck Process of Affect Dynamics with Event‑Based Jumps

[![Build Status](https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia implementation of a hierarchical Ornstein-Uhlenbeck process with event-based jumps for modeling affect dynamics.

## Abstract

Affective instability is a transdiagnostic feature of several mental disorders, yet ambulatory assessment studies often quantify it without distinguishing responses to measured events from continuous background dynamics. We extend the Ornstein-Uhlenbeck (OU) process with an event-related jump component.

The resulting OU-Jump (OU-J) model decomposes affective dynamics into an attractor, the half-life of deviations, volatility, and positive- and negative-event reactivity. Its closed-form transition distribution accommodates irregularly spaced observations and events occurring between measurements.

Parameter- and model-recovery simulations showed strong recovery of group-level parameters and high sensitivity to jump dynamics. In an application to intensive longitudinal data from 314 participants, incorporating measured events resulted in shorter estimated half-lives and lower volatility than the standard OU model.

## Models

The package provides three hierarchical Bayesian models:

- `OrnsteinUhlenbeckModel`: continuous mean-reverting dynamics
- `OrnsteinUhlenbeckJumpModel`: OU dynamics with positive and negative event effects
- `EventProbabilityModel`: hierarchical probabilities of binary events

## Installation

Follow the instructions on [https://julialang.org/install](https://julialang.org/install) to download and install Julia (if you have not already).

You can then install this package with the following command:

```julia
using Pkg

Pkg.add("https://github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl")
```

## Usage

### Data Preparation

Before fitting a model, rescale all affect ratings to the unit interval, regardless of their original measurement scale:

```julia
z_obs = (affect .- affect_min) ./ (affect_max - affect_min)
```

Here, affect_min and affect_max are the theoretical endpoints of the original response scale. For example, ratings recorded from 0 to 100 can be converted with:

```julia
z_obs = affect ./ 100
```

Positive and negative event magnitudes should likewise be represented on a common scale and rescaled to [0, 0.01]:

```julia
M⁺ = 0.01 .* positive_event_magnitude ./ maximum_event_magnitude
M⁻ = 0.01 .* negative_event_magnitude ./ maximum_event_magnitude
```

Under this scaling, an event with the maximum possible intensity has magnitude 0.01. Because affect is modeled on the unit interval and subsequently reported on a 0 to 100 scale, the event-reactivity parameters α⁺ and α⁻ can be interpreted as the expected affect change, in percentage points, immediately following a maximally intense event.

Use theoretical scale endpoints rather than sample-specific minima and maxima wherever possible. This preserves comparability across participants, groups, and datasets.

### Model Fitting

Construct aligned vectors for the observations, hierarchical indices, elapsed times, and event information. See the analysis directory, especially analysis/src/fitting, for an example of how to generate these vectors from a DataFrame.

```julia
using OrnsteinUhlenbeckJumpProcess, Turing, ReverseDiff, StableRNGs

args = (
    group_indices = group_indices,
    subject_indices = subject_indices,
    z_obs = z_obs,
    Δt = Δt,
    M⁺ = M⁺,
    M⁻ = M⁻,
    Δτ = Δτ,
    is_transition = is_transition
)

rng = StableRNG(1)
model = OrnsteinUhlenbeckJumpModel(; args...)

chains = sample(
    rng,
    model,
    NUTS(; adtype = AutoReverseDiff(; compile = true)),
    MCMCThreads(),
    2500,
    4
)
```

The time unit used for Δt and Δτ should be hours and determines the unit of the estimated half-lives.

## Posterior processing

Convert samples to interpretable parameter scales:

```julia
parameters = clean_parameters(OrnsteinUhlenbeckJumpModel, chains, args)
```

Summarize posterior medians and highest posterior density intervals:

```julia
summary = summarize_parameters(parameters)
```

Compare group-level posterior distributions:

```julia
comparisons = compare_groups(parameters)
```

Inspect trace and density plots:

```julia
figure = draw_chains(chains, 10)
```

## Citation

If you use this package or model, please cite the associated manuscript.
Complete citation information will be added upon publication.

## Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – GRK2739/2 – Project Nr. 447089431 – Research Training Group: KD²School – Designing Biosignal-Adaptive Systems for Decision-Making Processes