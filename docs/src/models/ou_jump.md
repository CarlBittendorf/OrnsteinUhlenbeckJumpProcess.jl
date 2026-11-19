# The Ornstein–Uhlenbeck jump model

A hierarchical Ornstein-Uhlenbeck model with positive and negative impulsive
jumps, exported as [`OrnsteinUhlenbeckJumpModel`](@ref). The model estimates
the same parameters as the plain Ornstein-Uhlenbeck model, plus positive and
negative response coefficients for measured events. The data conventions are
shared with all other models; see [Models](@ref).

```@example ou-jump
using OrnsteinUhlenbeckJumpProcess

model = OrnsteinUhlenbeckJumpModel(;
    group_indices = [1, 1],
    subject_indices = [1, 1, 2, 2],
    z_obs = [0.5, 0.6, 0.7, 0.4],
    Δt = [0.0, 1.0, 0.0, 2.0],
    M⁺ = [0.005, 0.005, 0.003, 0.002],
    M⁻ = [0.001, 0.0, 0.002, 0.001],
    Δτ = [0.0, 0.5, 0.0, 1.0],
    is_transition = [false, true, false, true]
);
```

## Likelihood

For transition nodes, the latent state follows an Ornstein-Uhlenbeck
transition whose conditional mean gains an exponentially decaying impulse:

$$z_i \mid z_{i-1} \sim \mathcal{N}\left(\mu_{s_i} + (z_{i-1} - \mu_{s_i})\, e^{-\theta_{s_i} \Delta t_i} + \left(\alpha^+_{s_i} M^+_i - \alpha^-_{s_i} M^-_i\right) e^{-\theta_{s_i} (\Delta t_i - \Delta \tau_i)},\ \gamma_{s_i}^2 \left(1 - e^{-2 \theta_{s_i} \Delta t_i}\right)\right),$$

where $M^+_i, M^-_i \in [0, 0.01]$ are the measured impulse magnitudes and
$\Delta \tau_i$ is the time between the preceding node and the impulse. The
first node of each series is drawn from the stationary distribution, and the
observation model is identical to the Ornstein-Uhlenbeck model.

The impulse acts on the transition mean, so its effect decays with the same
rate as the Ornstein-Uhlenbeck memory.

## Event effects

The response coefficients are positive and shared across subjects through an
additional log-normal hierarchy:

$$\log \alpha^+_g, \log \alpha^-_g \sim \mathcal{N}(\log 10,\ 0.5), \qquad \alpha^\pm_s = e^{\log \alpha^\pm_{g_s} + \kappa_{\alpha^\pm, g_s}\, \alpha^{\pm, \mathrm{raw}}_s}.$$

The between-subject scales $\kappa_{\alpha^+, g}$ and $\kappa_{\alpha^-, g}$
follow the same $\kappa = e^{-2.30259 + 0.5\, \kappa^{\mathrm{raw}}}$
parameterization as the other scales, with prior median $0.1$.

## Priors

All other priors, the stationary initialization, and the observation model
are identical to the Ornstein-Uhlenbeck model; see its
[priors](@ref "The Ornstein–Uhlenbeck model") for the full table.
