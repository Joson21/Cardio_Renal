"""Replaceable model inputs and experiment protocol helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


InputFunction = Callable[[float], float]


def constant_input(value: float) -> InputFunction:
    """Return an input function that always produces ``value``."""

    return lambda _time_min: float(value)


def step_input(
    baseline: float,
    changed: float,
    start_min: float,
    end_min: float | None = None,
) -> InputFunction:
    """Return a constant-step or finite-duration pulse input."""

    def profile(time_min: float) -> float:
        after_start = time_min >= start_min
        before_end = end_min is None or time_min < end_min
        return float(changed if after_start and before_end else baseline)

    return profile


def ramp_input(
    baseline: float,
    changed: float,
    start_min: float,
    duration_min: float,
) -> InputFunction:
    """Return a linear ramp that remains at the final value afterward."""

    if duration_min <= 0.0:
        raise ValueError("duration_min must be positive")

    def profile(time_min: float) -> float:
        if time_min <= start_min:
            return float(baseline)
        fraction = min(1.0, (time_min - start_min) / duration_min)
        return float(baseline + fraction * (changed - baseline))

    return profile


@dataclass(frozen=True)
class ModelInputs:
    """External signals supplied to the standalone kidney model.

    ``p_peritubular_mmHg`` is intentionally an input because the supplied
    equations define only its effect on RIHP, not an equation that predicts the
    absolute peritubular pressure. It can later be replaced by a coupled model.
    """

    p_ma_mmHg: InputFunction
    p_ra_mmHg: InputFunction
    sodium_intake_mEq_per_min: InputFunction
    p_peritubular_mmHg: InputFunction

    @classmethod
    def baseline(cls) -> "ModelInputs":
        return cls(
            p_ma_mmHg=constant_input(100.0),
            p_ra_mmHg=constant_input(0.0),
            sodium_intake_mEq_per_min=constant_input(0.126),
            p_peritubular_mmHg=constant_input(0.0),
        )
