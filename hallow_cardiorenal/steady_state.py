"""Coupled-equilibrium calculation and derivative reporting."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .model import CoupledCardiorenalModel
from .simulation import simulate_coupled
from .states import STATE_NAMES, initial_state


POSITIVE_STATE_INDICES = tuple(range(12)) + (14,)
FREE_STATE_INDICES = (12, 13)


def derivative_report(
    model: CoupledCardiorenalModel,
    y: np.ndarray | None = None,
) -> dict[str, float]:
    state = initial_state(model.p) if y is None else np.asarray(y, dtype=float)
    derivative = model.rhs(0.0, state)
    return {name: float(value) for name, value in zip(STATE_NAMES, derivative)}


def refine_coupled_equilibrium(
    model: CoupledCardiorenalModel,
    guess: np.ndarray,
    *,
    max_function_evaluations: int = 6000,
) -> tuple[np.ndarray, dict[str, float]]:
    """Refine a positive, physically constrained coupled equilibrium."""

    y_guess = np.asarray(guess, dtype=float)
    if y_guess.shape != (len(STATE_NAMES),):
        raise ValueError(f"Expected {len(STATE_NAMES)} state values")
    if np.any(y_guess[list(POSITIVE_STATE_INDICES)] <= 0.0):
        raise ValueError("All concentration, amount, volume, and vascularity guesses must be positive")

    positive_values = np.log(y_guess[list(POSITIVE_STATE_INDICES)])
    free_values = y_guess[list(FREE_STATE_INDICES)]
    z0 = np.concatenate([positive_values, free_values])
    derivative_scales = np.array(
        [
            0.126,
            0.001,
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
            1.0e-5,
        ],
        dtype=float,
    )

    def decode(z: np.ndarray) -> np.ndarray:
        y = np.empty(len(STATE_NAMES), dtype=float)
        y[list(POSITIVE_STATE_INDICES)] = np.exp(z[: len(POSITIVE_STATE_INDICES)])
        y[list(FREE_STATE_INDICES)] = z[len(POSITIVE_STATE_INDICES) :]
        return y

    def residual(z: np.ndarray) -> np.ndarray:
        return model.rhs(0.0, decode(z)) / derivative_scales

    solution = least_squares(
        residual,
        z0,
        max_nfev=max_function_evaluations,
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
    )
    equilibrium = decode(solution.x)
    scaled = residual(solution.x)
    report = {
        "success": float(solution.success),
        "cost": float(solution.cost),
        "optimality": float(solution.optimality),
        "max_abs_scaled_residual": float(np.max(np.abs(scaled))),
        "max_abs_unscaled_derivative": float(
            np.max(np.abs(model.rhs(0.0, equilibrium)))
        ),
        "function_evaluations": float(solution.nfev),
    }
    return equilibrium, report


def calculate_coupled_equilibrium(
    model: CoupledCardiorenalModel,
    *,
    guess: np.ndarray | None = None,
    settle_days: float = 300.0,
) -> tuple[np.ndarray, dict[str, float]]:
    """Settle the model and then refine all coupled derivatives to zero."""

    starting_state = (
        np.asarray(guess, dtype=float)
        if guess is not None
        else initial_state(model.p)
    )
    settled = simulate_coupled(
        model,
        days=settle_days,
        y0=starting_state,
        samples_per_day=1,
        method="BDF",
        relative_tolerance=1.0e-9,
        absolute_tolerance=1.0e-11,
    )
    if not settled.success:
        raise RuntimeError(settled.message)
    equilibrium, report = refine_coupled_equilibrium(
        model,
        settled.states[:, -1],
    )
    report["settle_days"] = float(settle_days)
    return equilibrium, report
