# Models

All models are hierarchical Bayesian models built with
[Turing.jl](https://turinglang.org) and share the same node-based data
conventions. They differ in the process they describe:

| Model | Process |
| --- | --- |
| [`OrnsteinUhlenbeckModel`](@ref) | Continuous mean-reverting dynamics for longitudinal observations |
| [`OrnsteinUhlenbeckJumpModel`](@ref) | Ornstein-Uhlenbeck dynamics with positive and negative event jumps |
| [`EventProbabilityModel`](@ref) | Hierarchical probabilities of binary events |

Follow the links for the likelihood and priors of each model. The remainder of
this page describes the conventions shared by all of them.

## Data conventions

All models operate on **nodes**: each element of the data vectors corresponds
to one observation node. Nodes belonging to the same subject form a time
series; the first node of each series is drawn from the stationary
distribution, and every later node follows an Ornstein-Uhlenbeck transition
from its predecessor.

| Argument | Meaning |
| --- | --- |
| `subject_indices` | One subject index per node (length = number of nodes). |
| `group_indices` | One group index per subject (length = number of subjects). |
| `Δt` | Hours elapsed since the preceding node of the same subject; `0` at series starts. Compute with `hours_elapsed`. |
| `is_transition` | `is_transition[i] == false` initializes node `i` from its stationary distribution. |
| `z_obs` | Observed values on the unit interval (0, 1). |
| `M⁺`, `M⁻` | Positive and negative impulse magnitudes in $[0, 0.01]$ (jump model only). |
| `Δτ` | Time from the preceding node to the impulse (jump model only). |
| `e` | Binary event indicator per observation (event model only). |

Additional keyword arguments are accepted and ignored (`kwargs...`), which
keeps model calls forward-compatible.

Throughout the model pages, $s$ indexes subjects, $g$ indexes groups, and node
$i$ belongs to subject $s_i$, which belongs to group $g_{s_i}$.

## Simulating data

Every model can generate data: pass the observations as `missing` and call
`rand` on the model to sample all random variables, including the observations
(see [Getting started](@ref) for a complete example).

## Interpreting cleaned parameters

`clean_parameters` (see [Analysis workflow](@ref)) converts posterior samples
into interpretable units:

- `μ_...` columns: equilibrium means on the 0-100 percentage scale.
- `θ_...` columns: deviation half-lives in hours, $\log 2 / \theta$.
- `γ_...` columns: stationary standard deviations on the percentage scale.
- `σ_...` columns: diffusion scale of the SDE $dz = \theta (\mu - z)\, dt + \sigma\, dW$, i.e. $\gamma \sqrt{2 \theta}$ on the percentage scale.
- `κ_μ_...` columns: percentage-point shifts per standardized subject effect.
- `κ_θ_...`, `κ_γ_...`, `κ_α⁺_...`, `κ_α⁻_...` columns: multiplicative scales of subject-level deviations.
- `α⁺_...`, `α⁻_...` columns: response coefficients per unit of impulse magnitude $M^\pm$.
- Event model: `β_group_...` and `β_subject_...` columns are event probabilities on the (0, 1) scale; `κ_β_...` columns are multiplicative scales of subject-level log-odds deviations.
