r"""
Mediation analysis: does mLOY mediate the effect of smoking on cancer risk?

The paper uses the unified approach of Lange, Vansteelandt & Bekaert (2012,
Am J Epidemiol) to estimate the natural direct effect (NDE) and the natural
indirect effect (NIE).

Causal structure under test:
    smoking (X)  -->  mLOY / mLRR (M)  -->  solid tumour (Y)
                 \-------------------------------/
                          direct effect (NDE)

The implementation below uses an equivalent regression/counterfactual
formulation for a rare time-to-event outcome, so the resulting ratios
approximate effects on the hazard scale. Paper's finding: the NIE is
essentially null (mLOY is NOT an important mediator) while the direct effect
of smoking remains strong.

Note: this is a simplified but conceptually aligned version of Lange et al.
For publication, use a fully validated implementation (e.g. the SAS macro or
the `CMAverse` R package).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def natural_effects(
    df: pd.DataFrame,
    exposure: str = "smoke_current",
    mediator: str = "mlrr_sd",
    outcome: str = "event",
    confounders: tuple[str, ...] = ("bmi", "alcohol_heavy", "edu_degree", "race_nonwhite"),
) -> dict:
    """Estimate the NDE and NIE (Lange-style product-of-coefficients approach
    on the log-hazard/log-risk scale for a rare outcome).

    Returns hazard ratios for the direct and indirect effects.
    """
    conf = " + ".join(confounders)

    # Mediator model: M ~ X + C  (linear, since mLRR is continuous)
    m_formula = f"{mediator} ~ {exposure} + {conf}"
    m_model = smf.ols(m_formula, data=df).fit()
    alpha = m_model.params[exposure]  # effect of X -> M

    # Outcome model: Y ~ X + M + C  (logistic as a log-risk proxy for a rare outcome)
    y_formula = f"{outcome} ~ {exposure} + {mediator} + {conf}"
    y_model = smf.glm(y_formula, data=df, family=sm.families.Binomial()).fit()
    beta_x = y_model.params[exposure]   # direct effect X -> Y | M
    beta_m = y_model.params[mediator]   # effect M -> Y | X

    nde_log = beta_x                    # natural direct effect (log scale)
    nie_log = alpha * beta_m            # natural indirect effect (product method)

    # Standard error of the NIE via the delta method
    var_alpha = m_model.bse[exposure] ** 2
    var_beta_m = y_model.bse[mediator] ** 2
    se_nie = np.sqrt(
        (beta_m ** 2) * var_alpha + (alpha ** 2) * var_beta_m
    )
    z = nie_log / se_nie if se_nie > 0 else np.nan

    from scipy.stats import norm
    p_nie = 2 * (1 - norm.cdf(abs(z))) if not np.isnan(z) else np.nan

    prop_mediated = nie_log / (nde_log + nie_log) if (nde_log + nie_log) != 0 else np.nan

    return {
        "exposure": exposure,
        "mediator": mediator,
        "NDE_HR": float(np.exp(nde_log)),
        "NIE_HR": float(np.exp(nie_log)),
        "NIE_ci_low": float(np.exp(nie_log - 1.96 * se_nie)),
        "NIE_ci_high": float(np.exp(nie_log + 1.96 * se_nie)),
        "NIE_p": float(p_nie),
        "proportion_mediated": float(prop_mediated),
    }


def format_mediation(res: dict) -> str:
    return (
        f"[{res['exposure']} -> {res['mediator']} -> outcome]\n"
        f"  DIRECT effect (NDE):        HR={res['NDE_HR']:.3f}\n"
        f"  INDIRECT effect (NIE):      HR={res['NIE_HR']:.3f} "
        f"(95% CI {res['NIE_ci_low']:.3f}-{res['NIE_ci_high']:.3f}); "
        f"P={res['NIE_p']:.2g}\n"
        f"  Proportion mediated:        {100*res['proportion_mediated']:.1f}%"
    )
