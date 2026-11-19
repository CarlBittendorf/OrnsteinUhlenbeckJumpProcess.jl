# The Ornstein–Uhlenbeck model

A hierarchical Ornstein-Uhlenbeck model for longitudinal observations,
exported as [`OrnsteinUhlenbeckModel`](@ref). The model estimates group-level
and subject-level equilibrium means, mean-reversion rates, and stationary
standard deviations. The data conventions are shared with all other models;
see [Models](@ref).

```@example ornstein-uhlenbeck
using OrnsteinUhlenbeckJumpProcess

model = OrnsteinUhlenbeckModel(;
    group_indices = [1, 1],
    subject_indices = [1, 1, 2, 2],
    z_obs = [0.5, 0.6, 0.7, 0.4],
    Δt = [0.0, 1.0, 0.0, 2.0],
    is_transition = [false, true, false, true]
);
```

## Likelihood

For transition nodes, the latent state follows an Ornstein-Uhlenbeck
transition from its predecessor:

$$z_i \mid z_{i-1} \sim \mathcal{N}\left(\mu_{s_i} + (z_{i-1} - \mu_{s_i})\, e^{-\theta_{s_i} \Delta t_i},\ \gamma_{s_i}^2 \left(1 - e^{-2 \theta_{s_i} \Delta t_i}\right)\right),$$

and the first node of each series is drawn from the stationary distribution:

$$z_i \sim \mathcal{N}(\mu_{s_i},\ \gamma_{s_i}^2).$$

Observations are noisy versions of the latent states:

$$z^{\mathrm{obs}}_i \sim \mathcal{N}(z_i,\ \sigma^2_{\mathrm{obs}}), \qquad \sigma^2_{\mathrm{obs}} = e^{-4.60517 + 0.05\, \sigma_{\mathrm{raw}}} \approx 0.01.$$

## Hierarchy and priors

| Prior | Interpretation |
| --- | --- |
| $\mu_g \sim \mathcal{N}(0.5,\ 0.2)$ | group equilibrium mean |
| $\log \theta_g \sim \mathcal{N}(\log\log 2,\ 0.5)$ | group mean-reversion rate; the prior median gives a deviation half-life of one hour |
| $\log \gamma_g \sim \mathcal{N}(\log 0.1,\ 0.5)$ | group stationary standard deviation |
| $\kappa^{\mathrm{raw}} \sim \mathcal{N}(0,\ 1)$, $\kappa = e^{-2.30259 + 0.5\, \kappa^{\mathrm{raw}}}$ | between-subject scale; prior median $0.1$ |
| $\mu^{\mathrm{raw}}_s,\ \theta^{\mathrm{raw}}_s,\ \gamma^{\mathrm{raw}}_s \sim \mathcal{N}(0,\ 1)$ | standardized subject-level effects |
| $\sigma_{\mathrm{raw}} \sim \mathcal{N}(0,\ 1)$ | observation-noise modifier |

Subject-level parameters are reconstructed as

$$\mu_s = \mu_{g_s} + \kappa_{\mu, g_s}\, \mu^{\mathrm{raw}}_s, \qquad \theta_s = e^{\log \theta_{g_s} + \kappa_{\theta, g_s}\, \theta^{\mathrm{raw}}_s}, \qquad \gamma_s = e^{\log \gamma_{g_s} + \kappa_{\gamma, g_s}\, \gamma^{\mathrm{raw}}_s}.$$

## Implementation notes

The implementation is non-centered: the transition noise enters through a
standardized `z_raw` vector and the subject-level parameters are
reconstructed from standardized effects, which improves the geometry of the
posterior for gradient-based samplers.
