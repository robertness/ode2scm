from mapk_scm import scm_erk_counterfactual


def main():
    rates = {
        'raf_activate': 0.1,
        'raf_deactivate': 0.1,
        'mek_activate': 0.1,
        'mek_deactivate': 2,
        'erk_activate': 0.1,
        'erk_deactivate': 1
    }

    totals = {
        'Raf': 100,
        'Mek': 100,
        'Erk': 100
    }

    observation = {
        'Raf': 49.99773,
        'Mek': 64.10115,
        'Erk': 84.7213
    }

    raf_intervention = {
        'Raf': 24.96818
    }

    out = scm_erk_counterfactual(
        rates,
        totals,
        observation,
        raf_intervention,
        spike_width=1.0,
        svi=True
    )
    print(out)


main()
