from collections import defaultdict
from functools import partial
from torch import tensor

import statistics
import pyro
from pyro import condition, do, sample
from pyro.optim import SGD
import torch.distributions.constraints as constraints

from pyro.distributions import Normal, Delta
from pyro.infer import EmpiricalMarginal, Importance, SVI, Trace_ELBO


def g(a, b):
    return a / (a + b)


class GF_SCM():
    def __init__(self, rates, totals, spike_width):
        self.rates = rates
        self.totals = totals

        def f_SOS(EGFR, IGFR, N):
            p = g(
                rates['SOS_activation_by_EGFR']*EGFR +
                rates['SOS_activation_by_IGFR']*IGFR,
                rates['SOS_deactivation']
            )
            mu = totals['SOS'] * p
            sigma = mu * (1. - p)
            SOS = N * sigma + mu
            return SOS

        def f_Ras(SOS, N):
            p = g(
                rates['Ras_activation_by_SOS']*SOS,
                rates['Ras_deactivation']
            )
            mu = totals['Ras'] * p
            sigma = mu * (1. - p)
            Ras = N * sigma + mu
            return Ras

        def f_PI3K(EGFR, IGFR, Ras, N):
            p = g(
                rates['PI3K_activation_by_EGFR'] * EGFR +
                rates['PI3K_activation_by_IGFR'] * IGFR +
                rates['PI3K_activation_by_Ras'] * Ras,
                rates['PI3K_deactivation']
            )
            mu = totals['PI3K'] * p
            sigma = mu * (1. - p)
            PI3K = N * sigma + mu
            return PI3K

        def f_AKT(PI3K, N):
            p = g(
                rates['AKT_activation_by_PI3K'] * PI3K,
                rates['AKT_deactivation']
            )
            mu = totals['AKT'] * p
            sigma = mu * (1. - p)
            AKT = N * sigma + mu
            return AKT

        def f_Raf(Ras, AKT, N):
            p = g(
                rates['Raf_activation_by_Ras'] * Ras,
                rates['Raf_deactivation_by_phosphotase'] +
                rates['Raf_deactivation_by_AKT'] * AKT
            )
            mu = totals['Raf'] * p
            sigma = mu * (1. - p)
            Raf = N * sigma + mu
            return Raf

        def f_Mek(Raf, N):
            p = g(
                rates['Mek_activation_by_Raf'] * Raf,
                rates['Mek_deactivation']
            )
            mu = totals['Mek'] * p
            sigma = mu * (1. - p)
            Mek = N * sigma + mu
            return Mek

        def f_Erk(Mek, N):
            p = g(
                rates['Erk_activation_by_Mek'] * Mek,
                rates['Erk_deactivation']
            )
            mu = totals['Erk'] * p
            sigma = mu * (1. - p)
            Erk = N * sigma + mu
            return Erk

        def model(noise):
            N_SOS = sample('N_SOS', noise['N_SOS'])
            N_Ras = sample('N_Ras', noise['N_Ras'])
            N_PI3K = sample('N_PI3K', noise['N_PI3K'])
            N_AKT = sample('N_AKT', noise['N_AKT'])
            N_Raf = sample('N_Raf', noise['N_Raf'])
            N_Mek = sample('N_Mek', noise['N_Mek'])
            N_Erk = sample('N_Erk', noise['N_Erk'])

            EGFR = tensor(37.)
            IGFR = tensor(5.)
            SOS = sample('SOS', Delta(f_SOS(EGFR, IGFR, N_SOS)))
            Ras = sample('Ras', Delta(f_Ras(SOS, N_Ras)))
            PI3K = sample('PI3K', Delta(f_PI3K(EGFR, IGFR, Ras, N_PI3K)))
            AKT = sample('AKT', Delta(f_AKT(PI3K, N_AKT)))
            Raf = sample('Raf', Delta(f_Raf(Ras, AKT, N_Raf)))
            Mek = sample('Mek', Delta(f_Mek(Raf, N_Mek)))
            Erk = sample('Erk', Delta(f_Erk(Mek, N_Erk)))
            noise_samples = N_SOS, N_Ras, N_PI3K, N_AKT, N_Raf, N_Mek, N_Erk
            proteins = SOS, Ras, PI3K, AKT, Raf, Mek, Erk
            return proteins, noise_samples

        Spike = partial(Normal, scale=tensor(spike_width))

        def noisy_model(noise):
            N_SOS = sample('N_SOS', noise['N_SOS'])
            N_Ras = sample('N_Ras', noise['N_Ras'])
            N_PI3K = sample('N_PI3K', noise['N_PI3K'])
            N_AKT = sample('N_AKT', noise['N_AKT'])
            N_Raf = sample('N_Raf', noise['N_Raf'])
            N_Mek = sample('N_Mek', noise['N_Mek'])
            N_Erk = sample('N_Erk', noise['N_Erk'])

            EGFR = tensor(37.)
            IGFR = tensor(5.)
            SOS = sample('SOS', Spike(f_SOS(EGFR, IGFR, N_SOS)))
            Ras = sample('Ras', Spike(f_Ras(SOS, N_Ras)))
            PI3K = sample('PI3K', Spike(f_PI3K(EGFR, IGFR, Ras, N_PI3K)))
            AKT = sample('AKT', Spike(f_AKT(PI3K, N_AKT)))
            Raf = sample('Raf', Spike(f_Raf(Ras, AKT, N_Raf)))
            Mek = sample('Mek', Spike(f_Mek(Raf, N_Mek)))
            Erk = sample('Erk', Spike(f_Erk(Mek, N_Erk)))
            noise_samples = N_SOS, N_Ras, N_PI3K, N_AKT, N_Raf, N_Mek, N_Erk
            proteins = SOS, Ras, PI3K, AKT, Raf, Mek, Erk
            return proteins, noise_samples

        self.model = model
        self.noisy_model = noisy_model

    def infer(self, model, noise):
        return Importance(model, num_samples=1000).run(noise)

    def update_noise_svi(self, observed_steady_state, initial_noise):
        def guide(noise):
            noise_terms = list(noise.keys())
            mu_constraints = constraints.interval(-3., 3.)
            sigma_constraints = constraints.interval(.0001, 3)
            mu = {
                k: pyro.param(
                    '{}_mu'.format(k),
                    tensor(0.),
                    constraint=mu_constraints
                ) for k in noise_terms
            }
            sigma = {
                k: pyro.param(
                    '{}_sigma'.format(k),
                    tensor(1.),
                    constraint=sigma_constraints
                ) for k in noise_terms
            }
            for noise in noise_terms:
                sample(noise, Normal(mu[noise], sigma[noise]))

        observation_model = condition(self.noisy_model, observed_steady_state)
        pyro.clear_param_store()
        svi = SVI(
            model=observation_model,
            guide=guide,
            optim=SGD({"lr": 0.001, "momentum": 0.1}),
            loss=Trace_ELBO()
        )

        losses = []
        num_steps = 1000
        samples = defaultdict(list)
        for t in range(num_steps):
            losses.append(svi.step(initial_noise))
            for noise in initial_noise.keys():
                mu = '{}_mu'.format(noise)
                sigma = '{}_sigma'.format(noise)
                samples[mu].append(pyro.param(mu).item())
                samples[sigma].append(pyro.param(sigma).item())
        means = {k: statistics.mean(v) for k, v in samples.items()}
        updated_noise = {
            'N_SOS': Normal(means['N_SOS_mu'], means['N_SOS_sigma']),
            'N_Ras': Normal(means['N_Ras_mu'], means['N_Ras_sigma']),
            'N_PI3K': Normal(means['N_PI3K_mu'], means['N_PI3K_sigma']),
            'N_AKT': Normal(means['N_AKT_mu'], means['N_AKT_sigma']),
            'N_Raf': Normal(means['N_Raf_mu'], means['N_Raf_sigma']),
            'N_Mek': Normal(means['N_Mek_mu'], means['N_Mek_sigma']),
            'N_Erk': Normal(means['N_Erk_mu'], means['N_Erk_sigma'])
        }

        return updated_noise, losses

    def update_noise_importance(self, observed_steady_state, initial_noise):
        observation_model = condition(self.noisy_model, observed_steady_state)
        posterior = self.infer(observation_model, initial_noise)
        updated_noise = {
            k: EmpiricalMarginal(posterior, sites=k)
            for k in initial_noise.keys()
        }
        return updated_noise

def scm_ras_erk_counterfactual(
    rates,
    totals,
    observation,
    ras_intervention,
    spike_width=1.0,
    svi=True
):
    gf_scm = GF_SCM(rates, totals, spike_width)
    noise = {
        'N_SOS': Normal(0., 1.),
        'N_Ras': Normal(0., 1.),
        'N_PI3K': Normal(0., 1.),
        'N_AKT': Normal(0., 1.),
        'N_Raf': Normal(0., 1.),
        'N_Mek': Normal(0., 1.),
        'N_Erk': Normal(0., 1.)
    }
    if svi:
        updated_noise, _ = gf_scm.update_noise_svi(observation, noise)
    else:
        updated_noise = gf_scm.update_noise_importance(observation, noise)
    counterfactual_model = do(gf_scm.model, ras_intervention)
    cf_posterior = gf_scm.infer(counterfactual_model, updated_noise)
    cf_erk_marginal = EmpiricalMarginal(cf_posterior, sites='Erk')

    scm_causal_effect_samples = [
        observation['Erk'] - float(cf_erk_marginal.sample())
        for _ in range(500)
    ]
    return scm_causal_effect_samples
