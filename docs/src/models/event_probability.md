# The event probability model

A hierarchical logistic model for binary events, exported as
[`EventProbabilityModel`](@ref). Each subject has a log-odds $\beta_s$ nested
within its group. Unlike the other models, there is no latent time series:
every observation is an independent Bernoulli trial conditional on its
subject's log-odds. The indexing conventions are shared with the other models;
see [Models](@ref).

```@example event-probability
using OrnsteinUhlenbeckJumpProcess

model = EventProbabilityModel(;
    group_indices = [1, 2],
    subject_indices = [1, 2, 1, 2],
    e = [1, 0, 0, 1]
);
```

## Likelihood

$$\beta_s = \beta_{g_s} + \kappa_{\beta, g_s}\, \beta^{\mathrm{raw}}_s, \qquad p_i = \mathrm{logistic}(\beta_{s_i}), \qquad e_i \sim \mathrm{Bernoulli}(p_i).$$

## Priors

| Prior | Interpretation |
| --- | --- |
| $\beta_g \sim \mathcal{N}(0,\ 1)$ | group log-odds |
| $\kappa_{\beta, g} \sim \mathrm{truncated}(\mathcal{N}(0,\ 1);\ \mathrm{lower} = 0)$ | between-subject scale; non-negative by construction |
| $\beta^{\mathrm{raw}}_s \sim \mathcal{N}(0,\ 1)$ | standardized subject-level effect |

## Interpreting cleaned parameters

`clean_parameters` reconstructs the subject-level log-odds and converts the
log-odds to probabilities:

- `β_group_...` and `β_subject_...` columns: event probabilities on the (0, 1)
  scale.
- `κ_β_...` columns: multiplicative scales of subject-level log-odds
  deviations, kept on the log-odds scale.
