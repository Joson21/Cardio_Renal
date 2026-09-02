"""Automated mathematical and baseline-consistency checks."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np

from .model import KidneyModel
from .states import STATE_INDEX, initial_state
from .steady_state import isolated_raas_equilibrium


@dataclass(frozen=True)
class ValidationCheck:
    name: str
    passed: bool
    value: float
    tolerance: float | None
    explanation: str


def run_validation(model: KidneyModel) -> list[ValidationCheck]:
    """Run checks that should hold for every supplied configuration."""

    y0 = initial_state(model.p)
    out = model.algebraic_outputs(0.0, y0)
    derivative = model.rhs(0.0, y0)
    equilibrium = isolated_raas_equilibrium(model.p)

    raas_state_names = (
        "prc_pg_per_mL",
        "agt_fmol_per_mL",
        "angi_fmol_per_mL",
        "angii_fmol_per_mL",
        "ang17_fmol_per_mL",
        "angiv_fmol_per_mL",
        "at1_bound_fmol_per_mL",
        "at2_bound_fmol_per_mL",
    )
    raas_initial_error = max(
        abs(y0[STATE_INDEX[name]] - equilibrium[name]) for name in raas_state_names
    )
    raas_derivative_indices = [STATE_INDEX[name] for name in raas_state_names]
    raas_derivative_max = float(
        np.max(np.abs(derivative[raas_derivative_indices]))
    )
    finite_outputs = all(isfinite(value) for value in out.values())
    extension_neutrality = max(
        abs(out["p_b_cvp_mmHg"]),
        abs(out["gamma_rihp_pt"] - 1.0),
        abs(out["gamma_rihp_dt"] - 1.0),
        abs(out["gamma_rihp_cd"] - 1.0),
    )

    return [
        ValidationCheck(
            "corrected_afferent_resistance",
            model.p.r_aa_0_mmHg_min_per_L_per_nephron == 6.0e7,
            model.p.r_aa_0_mmHg_min_per_L_per_nephron,
            None,
            "Single-nephron nominal resistance must include 2e6 parallel nephrons.",
        ),
        ValidationCheck(
            "corrected_efferent_resistance",
            model.p.r_ea_0_mmHg_min_per_L_per_nephron == 1.0e8,
            model.p.r_ea_0_mmHg_min_per_L_per_nephron,
            None,
            "Single-nephron nominal resistance must include 2e6 parallel nephrons.",
        ),
        ValidationCheck(
            "isolated_raas_initial_condition",
            raas_initial_error < 1.0e-9,
            raas_initial_error,
            1.0e-9,
            "Agreed peptide initial conditions must equal the unit-feedback equilibrium.",
        ),
        ValidationCheck(
            "raas_initial_derivative",
            raas_derivative_max < 1.0e-8,
            raas_derivative_max,
            1.0e-8,
            "In baseline-normalized mode, the RAAS branch starts at equilibrium.",
        ),
        ValidationCheck(
            "renal_algebraic_residual",
            abs(out["p_gh_residual_mmHg"]) < 1.0e-7,
            abs(out["p_gh_residual_mmHg"]),
            1.0e-7,
            "The implicit renal hemodynamic loop must be solved consistently.",
        ),
        ValidationCheck(
            "baseline_extension_neutrality",
            extension_neutrality < 1.0e-12,
            extension_neutrality,
            1.0e-12,
            "CVP and RIHP additions must be exactly neutral at their references.",
        ),
        ValidationCheck(
            "finite_baseline_outputs",
            finite_outputs,
            1.0 if finite_outputs else 0.0,
            None,
            "Every algebraic output must be finite.",
        ),
        ValidationCheck(
            "positive_gfr",
            out["gfr_mL_per_min"] > 0.0,
            out["gfr_mL_per_min"],
            None,
            "The reference state should have positive filtration.",
        ),
    ]


def assert_validation(model: KidneyModel) -> list[ValidationCheck]:
    checks = run_validation(model)
    failed = [check for check in checks if not check.passed]
    if failed:
        details = "; ".join(
            f"{check.name}={check.value:.6g}" for check in failed
        )
        raise AssertionError(f"Model validation failed: {details}")
    return checks
