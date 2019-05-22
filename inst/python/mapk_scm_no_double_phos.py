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


def g1(a):
    return a / (a + tensor(1.))


def g2(a):
    return a**2 / (a**2 + a + tensor(1.))


class MAPK_SCM():
    def __init__(self, rates, totals, spike_width):
        self.rates = rates
        self.totals = totals

        def f_Raf(N):
            E1 = tensor(1.)
            weight = self.rates['raf_activate']/rates['raf_deactivate']
            p = g1(E1 * weight)
            mu = totals['Raf'] * p
            sigma = mu * (1. - p)
            Raf = N * sigma + mu
            return Raf

        def f_Mek(Raf, N):
            weight = rates['mek_activate']/rates['mek_deactivate']
            p = g1(Raf * weight)
            mu = totals['Mek'] * p
            sigma = mu * (1. - p)
            Mek = N * sigma + mu
            return Mek

        def f_Erk(Mek, N):
            weight = rates['erk_activate']/rates['erk_deactivate']
            p = g1(Mek * weight)
            mu = totals['Erk'] * p
            sigma = mu * (1. - p)
            Erk = N * sigma + mu
            return Erk

        def model(noise):
            N_Raf = sample('N_Raf', noise['N_Raf'])
            N_Mek = sample('N_Mek', noise['N_Mek'])
            N_Erk = sample('N_Erk', noise['N_Erk'])
            Raf = sample('Raf', Delta(f_Raf(N_Raf)))
            Mek = sample('Mek', Delta(f_Mek(Raf, N_Mek)))
            Erk = sample('Erk', Delta(f_Erk(Mek, N_Erk)))
            return Raf, Mek, Erk, N_Raf, N_Mek, N_Erk

        Spike = partial(Normal, scale=tensor(spike_width))

        def noisy_model(noise):
            N_Raf = sample('N_Raf', noise['N_Raf'])
            N_Mek = sample('N_Mek', noise['N_Mek'])
            N_Erk = sample('N_Erk', noise['N_Erk'])
            Raf = sample('Raf', Spike(f_Raf(N_Raf)))
            Mek = sample('Mek', Spike(f_Mek(Raf, N_Mek)))
            Erk = sample('Erk', Spike(f_Erk(Mek, N_Erk)))
            return Raf, Mek, Erk

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


def scm_erk_counterfactual(
    rates,
    totals,
    observation,
    raf_intervention,
    spike_width=1.0,
    svi=True
):
    mapk_scm = MAPK_SCM(rates, totals, spike_width)
    noise = {
        'N_Raf': Normal(0., 1.),
        'N_Mek': Normal(0., 1.),
        'N_Erk': Normal(0., 1.)
    }
    if svi:
        updated_noise, _ = mapk_scm.update_noise_svi(observation, noise)
    else:
        updated_noise = mapk_scm.update_noise_importance(observation, noise)
    counterfactual_model = do(mapk_scm.model, raf_intervention)
    cf_posterior = mapk_scm.infer(counterfactual_model, updated_noise)
    cf_erk_marginal = EmpiricalMarginal(cf_posterior, sites='Erk')

    scm_causal_effect_samples = [
        observation['Erk'] - float(cf_erk_marginal.sample())
        for _ in range(500)
    ]
    return scm_causal_effect_samples
