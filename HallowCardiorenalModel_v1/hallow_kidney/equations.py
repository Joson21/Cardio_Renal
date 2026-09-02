"""Small reusable mathematical relationships used by the kidney model."""

from __future__ import annotations

from math import exp, log10

from .parameters import ModelParameters


def safe_exp(value: float) -> float:
    """Evaluate an exponential without floating-point overflow."""

    return exp(max(-700.0, min(700.0, value)))


def positive(value: float, floor: float) -> float:
    """Return a positive value for logarithms and divisions."""

    return max(floor, value)


def fmol_per_mL_to_pg_per_mL(
    concentration_fmol_per_mL: float,
    molecular_weight_g_per_mol: float,
) -> float:
    """Convert fmol/mL to pg/mL for a peptide of known molecular weight."""

    return concentration_fmol_per_mL * molecular_weight_g_per_mol * 1.0e-3


def at1_vascular_effect(
    at1_bound_pg_per_mL: float,
    a: float,
    b: float,
    c: float,
    floor: float,
) -> float:
    x = positive(at1_bound_pg_per_mL, floor)
    return a + b * x - c / x


def myogenic_multiplier(p_gh_mmHg: float, p: ModelParameters) -> float:
    """Corrected equation containing the missing leading +1."""

    return 1.0 + p.c_gp_autoreg * (p_gh_mmHg / p.p_gh_nom_mmHg - 1.0)


def tgf_multiplier(phi_md_sodium_mEq_per_min: float) -> float:
    exponent = (phi_md_sodium_mEq_per_min - 3.859) / -0.9617
    return 0.3408 + 3.449 / (3.88 + safe_exp(exponent))


def raw_md_renin_feedback(phi_md_sodium_mEq_per_min: float, p: ModelParameters) -> float:
    exponent = (
        phi_md_sodium_mEq_per_min - p.offset_md_renin_mEq_per_min
    ) / p.x_md_renin_mEq_per_min
    return p.a_md_renin + p.b_md_renin / (p.c_md_renin + safe_exp(exponent))


def raw_rsna_renin_feedback(rsna: float) -> float:
    return 1.89 - 2.056 / (1.358 + safe_exp(rsna - 0.8667))


def at1_renin_feedback(at1_bound_pg_per_mL: float, p: ModelParameters) -> float:
    ratio = positive(at1_bound_pg_per_mL, p.minimum_positive) / p.at1_feedback_eq_pg_per_mL
    exponent = p.a_at1_renin - p.b_at1_renin * log10(ratio)
    return 10.0**exponent


def rihp_multiplier(delta_rihp_mmHg: float, sensitivity: float, p: ModelParameters) -> float:
    """Provisional pressure-natriuresis multiplier supplied by the user."""

    exponent = delta_rihp_mmHg / p.rihp_pressure_scale_mmHg
    logistic = 1.0 / (1.0 + safe_exp(exponent))
    return 1.0 + sensitivity * (logistic - 0.5)
