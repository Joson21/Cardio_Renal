"""Equilibrium calculations and residual reporting."""

from __future__ import annotations

from math import log

import numpy as np
from scipy.optimize import least_squares

from .model import KidneyModel
from .parameters import ModelParameters
from .states import STATE_NAMES, initial_state


def isolated_raas_equilibrium(p: ModelParameters) -> dict[str, float]:
    """Solve the systemic RAAS cascade for unit renin feedback multipliers."""

    prc = p.n_rs_pg_per_mL_per_min / (log(2.0) / p.h_renin_min)
    pra = prc * p.x_prc_pra_fmol_per_h_per_pg / 60.0
    agt = (p.k_agt_fmol_per_mL_per_min - pra) / (log(2.0) / p.h_agt_min)
    angi = pra / (
        p.c_ace_per_min
        + p.c_chym_per_min
        + p.c_nep_per_min
        + log(2.0) / p.h_angi_min
    )
    angii = (
        (p.c_ace_per_min + p.c_chym_per_min) * angi
        / (
            p.c_ace2_per_min
            + p.c_angii_to_angiv_per_min
            + p.c_at1_per_min
            + p.c_at2_per_min
            + log(2.0) / p.h_angii_min
        )
    )
    ang17 = (
        p.c_nep_per_min * angi + p.c_ace2_per_min * angii
    ) / (log(2.0) / p.h_ang17_min)
    angiv = p.c_angii_to_angiv_per_min * angii / (
        log(2.0) / p.h_angiv_min
    )
    at1 = p.c_at1_per_min * angii / (log(2.0) / p.h_at1_min)
    at2 = p.c_at2_per_min * angii / (log(2.0) / p.h_at2_min)
    return {
        "prc_pg_per_mL": prc,
        "pra_fmol_per_mL_per_min": pra,
        "agt_fmol_per_mL": agt,
        "angi_fmol_per_mL": angi,
        "angii_fmol_per_mL": angii,
        "ang17_fmol_per_mL": ang17,
        "angiv_fmol_per_mL": angiv,
        "at1_bound_fmol_per_mL": at1,
        "at2_bound_fmol_per_mL": at2,
    }


def derivative_report(model: KidneyModel, y: np.ndarray | None = None) -> dict[str, float]:
    """Return every state derivative at the reference time."""

    state = initial_state(model.p) if y is None else np.asarray(y, dtype=float)
    derivative = model.rhs(0.0, state)
    return {name: float(value) for name, value in zip(STATE_NAMES, derivative)}


def refine_full_steady_state(
    model: KidneyModel,
    guess: np.ndarray | None = None,
    *,
    max_function_evaluations: int = 4000,
) -> tuple[np.ndarray, dict[str, float]]:
    """Numerically refine a full coupled equilibrium for constant inputs.

    Positive physiological states are represented logarithmically during the
    optimization. The two adaptation integrals remain unconstrained. This
    function is deliberately separate from ``initial_state``: the agreed RAAS
    initialization is source/branch-derived, while a coupled steady state also
    depends on the selected external inputs and extension switches.
    """

    y_guess = np.asarray(
        guess if guess is not None else initial_state(model.p), dtype=float
    )
    if np.any(y_guess[:12] <= 0.0):
        raise ValueError("The first 12 state guesses must be positive")

    z0 = np.concatenate([np.log(y_guess[:12]), y_guess[12:]])
    derivative_scales = np.array(
        [
            0.126,
            0.0015,
            1.0,
            100.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.03,
            0.15,
            0.1,
            0.1,
        ]
    )

    def decode(z: np.ndarray) -> np.ndarray:
        return np.concatenate([np.exp(z[:12]), z[12:]])

    def residual(z: np.ndarray) -> np.ndarray:
        return model.rhs(0.0, decode(z)) / derivative_scales

    solution = least_squares(
        residual,
        z0,
        max_nfev=max_function_evaluations,
        xtol=1.0e-11,
        ftol=1.0e-11,
        gtol=1.0e-11,
    )
    steady = decode(solution.x)
    report = {
        "success": float(solution.success),
        "cost": float(solution.cost),
        "optimality": float(solution.optimality),
        "max_abs_scaled_residual": float(np.max(np.abs(residual(solution.x)))),
        "function_evaluations": float(solution.nfev),
    }
    return steady, report
