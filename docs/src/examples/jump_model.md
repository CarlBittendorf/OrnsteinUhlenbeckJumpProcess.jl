# Jump model workflow

This page fits [`OrnsteinUhlenbeckJumpModel`](@ref) to synthetic data with
measured events and analyzes the event-reactivity parameters. The workflow is
identical to [Analysis workflow](@ref); here we focus on the event-related
parts.

## Simulating event-related jumps

Each node carries the magnitude of the positive or negative event associated
with it, scaled to [0, 0.01], and the time between the preceding node and the
event. We draw observations from a data-generating model with the default
priors:

```@example jump-model
using OrnsteinUhlenbeckJumpProcess, StableRNGs, Turing
using OrnsteinUhlenbeckJumpProcess.AlgebraOfGraphics: set_aog_theme!

set_aog_theme!()

data = (;
    group_indices = [1, 2],
    subject_indices = [1, 1, 1, 2, 2, 2],
    Δt = repeat([0.0, 8.0, 16.0], 2),
    M⁺ = [0.0, 0.01, 0.0, 0.0, 0.0, 0.005],
    M⁻ = [0.0, 0.0, 0.002, 0.0, 0.0, 0.0],
    Δτ = [0.0, 4.0, 2.0, 0.0, 0.0, 6.0],
    is_transition = repeat([false, true, true], 2)
)
generative = OrnsteinUhlenbeckJumpModel(; data..., z_obs = missing)
z_obs = rand(StableRNG(1), generative).z_obs;

model = OrnsteinUhlenbeckJumpModel(; data..., z_obs = z_obs);
```

Zeros in `M⁺` and `M⁻` mark nodes without an event of that kind, so their
magnitudes do not influence the likelihood.

## Fitting the model

```@example jump-model
chains = sample(StableRNG(1), model, NUTS(), 200; progress = false);
size(chains)
```

## Cleaning the posterior samples

```@example jump-model
cleaned = clean_parameters(OrnsteinUhlenbeckJumpModel, chains, data);
names(cleaned)
```

The cleaned columns include the event-reactivity parameters. Summarize the
group-level reactivity:

```@example jump-model
summarize_parameters(cleaned[:, [
    "α⁺_group[1]", "α⁻_group[1]", "α⁺_group[2]", "α⁻_group[2]"
]])
```

Under the event scaling, `α⁺` and `α⁻` can be read as the expected affect
change, in percentage points, immediately following a maximally intense event.

## Comparing groups

```@example jump-model
compare_groups(cleaned[:, [
    "α⁺_group[1]", "α⁺_group[2]", "α⁻_group[1]", "α⁻_group[2]"
]])
```

## Visualizing chains

```@example jump-model
draw_chains(chains, 1:4)
```

## Next steps

- [The Ornstein–Uhlenbeck jump model](@ref): the likelihood and priors.
- [Analysis workflow](@ref): the same workflow for the plain model.
- [API reference](@ref): the full function reference.
