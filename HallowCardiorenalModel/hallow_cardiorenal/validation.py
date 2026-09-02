"""Automated equation, coupling, and equilibrium checks."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np

from hallow_kidney.states import initial_state as kidney_initial_state
from hallow_kidney.steady_state import isolated_raas_equilibrium

from .model import CoupledCardiorenalModel
from .states import STATE_INDEX, initial_state


@dataclass(frozen=True)
class ValidationCheck:
    name: str
    passed: bool
    value: float
    tolerance: float | None
    explanation: str


def run_validation(
    model: CoupledCardiorenalModel,
    *,
    equilibrium_state: np.ndarray | None = None,
) -> list[ValidationCheck]:
    y0 = initial_state(model.p)
    output = model.algebraic_outputs(0.0, y0)
    derivative = model.rhs(0.0, y0)
    kidney_y0 = kidney_initial_state(model.p.kidney)
    raas = isolated_raas_equilibrium(model.p.kidney)
    raas_names = (
        "prc_pg_per_mL",
        "agt_fmol_per_mL",
        "angi_fmol_per_mL",
        "angii_fmol_per_mL",
        "ang17_fmol_per_mL",
        "angiv_fmol_per_mL",
        "at1_bound_fmol_per_mL",
        "at2_bound_fmol_per_mL",
    )
    raas_error = max(
        abs(kidney_y0[2 + index] - raas[name])
        for index, name in enumerate(raas_names)
    )
    finite = all(isfinite(value) for value in output.values())
    checks = [
        ValidationCheck(
            "cardiovascular_algebraic_residual",
            abs(output["cardiovascular_algebraic_residual_L_per_min"]) < 1.0e-8,
            abs(output["cardiovascular_algebraic_residual_L_per_min"]),
            1.0e-8,
            "Cardiac output must equal venous return.",
        ),
        ValidationCheck(
            "cvp_bowman_reference_neutrality",
            abs(output["p_b_cvp_mmHg"]) < 1.0e-12,
            abs(output["p_b_cvp_mmHg"]),
            1.0e-12,
            "The CVP novelty must remain neutral at Pra=Pra_ref.",
        ),
        ValidationCheck(
            "raas_initial_equilibrium",
            raas_error < 1.0e-9,
            raas_error,
            1.0e-9,
            "The agreed systemic RAAS initial values must match the analytical unit-feedback equilibrium.",
        ),
        ValidationCheck(
            "finite_coupled_outputs",
            finite,
            1.0 if finite else 0.0,
            None,
            "Every coupled algebraic output must be finite.",
        ),
    ]
    if model.p.cardiovascular.closure_mode == "baseline_centered":
        checks[1:1] = [
            ValidationCheck(
                "centered_reference_map",
                abs(output["p_ma_mmHg"] - 100.0) < 1.0e-9,
                output["p_ma_mmHg"],
                1.0e-9,
                "The recommended closure must preserve MAP=100 mmHg at the nominal reference.",
            ),
            ValidationCheck(
                "centered_reference_pra",
                abs(output["p_ra_mmHg"]) < 1.0e-9,
                output["p_ra_mmHg"],
                1.0e-9,
                "The recommended closure must preserve Pra=0 mmHg at the nominal reference.",
            ),
            ValidationCheck(
                "centered_reference_cardiac_output",
                abs(output["cardiac_output_L_per_min"] - 5.0) < 1.0e-9,
                output["cardiac_output_L_per_min"],
                1.0e-9,
                "The recommended closure must preserve cardiac output=5 L/min.",
            ),
            ValidationCheck(
                "vascularity_reference_equilibrium",
                abs(derivative[STATE_INDEX["vascularity"]]) < 1.0e-12,
                abs(derivative[STATE_INDEX["vascularity"]]),
                1.0e-12,
                "Baseline-centered vascular formation and destruction must balance at the reference.",
            ),
        ]
    if equilibrium_state is not None:
        equilibrium_derivative = model.rhs(0.0, equilibrium_state)
        max_derivative = float(np.max(np.abs(equilibrium_derivative)))
        equilibrium_output = model.algebraic_outputs(0.0, equilibrium_state)
        sodium_balance = abs(
            equilibrium_output["sodium_intake_mEq_per_min"]
            - equilibrium_output["phi_urine_sodium_mEq_per_min"]
        )
        water_balance = abs(
            equilibrium_output["phi_water_in_L_per_min"]
            - equilibrium_output["phi_urine_L_per_min"]
        )
        checks.extend(
            [
                ValidationCheck(
                    "coupled_equilibrium_derivative",
                    max_derivative < 1.0e-8,
                    max_derivative,
                    1.0e-8,
                    "Every state derivative must vanish at the supplied coupled equilibrium.",
                ),
                ValidationCheck(
                    "coupled_sodium_balance",
                    sodium_balance < 1.0e-9,
                    sodium_balance,
                    1.0e-9,
                    "Sodium intake and urinary sodium output must balance at equilibrium.",
                ),
                ValidationCheck(
                    "coupled_water_balance",
                    water_balance < 1.0e-10,
                    water_balance,
                    1.0e-10,
                    "Water intake and urine output must balance at equilibrium.",
                ),
            ]
        )
    return checks


def assert_validation(
    model: CoupledCardiorenalModel,
    *,
    equilibrium_state: np.ndarray | None = None,
) -> list[ValidationCheck]:
    checks = run_validation(model, equilibrium_state=equilibrium_state)
    failed = [check for check in checks if not check.passed]
    if failed:
        details = "; ".join(
            f"{check.name}={check.value:.6g}" for check in failed
        )
        raise AssertionError(f"Coupled validation failed: {details}")
    return checks
