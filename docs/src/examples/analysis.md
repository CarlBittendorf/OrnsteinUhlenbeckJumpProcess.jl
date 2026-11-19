# Analysis workflow

This page walks through the full analysis workflow: fitting a model, cleaning
the posterior samples into interpretable parameters, summarizing them,
comparing groups, and visualizing the chains.

## Fitting a model

We fit the Ornstein-Uhlenbeck model to synthetic data from two subjects nested
in two groups, using Turing's NUTS sampler:

```@example analysis
using OrnsteinUhlenbeckJumpProcess, Turing, StableRNGs
using OrnsteinUhlenbeckJumpProcess.AlgebraOfGraphics: set_aog_theme!

set_aog_theme!()

data = (;
    group_indices = [1, 2],
    subject_indices = [1, 1, 1, 2, 2, 2],
    Δt = repeat([0.0, 8.0, 16.0], 2),
    is_transition = repeat([false, true, true], 2)
)
generative = OrnsteinUhlenbeckModel(; data..., z_obs = missing)
z_obs = rand(StableRNG(1), generative).z_obs;

model = OrnsteinUhlenbeckModel(; data..., z_obs = z_obs);
chains = sample(StableRNG(1), model, NUTS(), 200; progress = false);
```

For real data, the vectors would come from the measurements instead; see the
README for the recommended data preparation and sampler configuration at
scale.

## Cleaning posterior samples

The chain contains raw, non-centered parameters. `clean_parameters`
reconstructs the interpretable parameters and drops the latent and raw
variables:

```@example analysis
cleaned = clean_parameters(OrnsteinUhlenbeckModel, chains, data);
names(cleaned)
```

## Summarizing parameters

`summarize_parameters` computes medians and 95% highest posterior density
intervals for every cleaned parameter:

```@example analysis
summarize_parameters(cleaned[:, [
    "μ_group[1]", "μ_group[2]", "θ_group[1]", "θ_group[2]",
    "γ_group[1]", "γ_group[2]", "σ_group[1]", "σ_group[2]"
]])
```

## Comparing groups

`compare_groups` quantifies how clearly separated the group-level posteriors
are. For every pair of groups within a parameter, it reports the posterior
directional error:

```@example analysis
compare_groups(cleaned)
```

The **posterior directional error** of two posterior samples is the smaller of
the posterior probabilities of each direction, so values near `0` indicate
strong evidence for one direction and values near `0.5` indicate that the data
cannot distinguish the two groups:

```@example analysis
posterior_directional_error(
    cleaned[:, "μ_subject[1]"], cleaned[:, "μ_subject[2]"]
)
```

## Visualizing chains

`draw_chains` plots trace plots and densities for selected parameters of a
chain. Indices refer to the positions of the parameters in the chain:

```@example analysis
draw_chains(chains, 1:4)
```

## Next steps

- [Jump model workflow](@ref): the same workflow for the model with event
  jumps.
- [Models](@ref): the interpretation of the cleaned parameter scales.
- [Extending the package](@ref): the auxiliary `analysis/` subpackage for
  simulation and fitting at scale.
- [API reference](@ref): the full function reference.
