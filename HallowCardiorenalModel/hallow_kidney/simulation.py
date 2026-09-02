"""Numerical integration and result export helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .model import KidneyModel
from .states import STATE_NAMES, initial_state


@dataclass
class SimulationResult:
    """State trajectories plus recalculated algebraic outputs."""

    time_min: np.ndarray
    states: np.ndarray
    outputs: dict[str, np.ndarray]
    success: bool
    message: str

    @property
    def time_days(self) -> np.ndarray:
        return self.time_min / 1440.0

    def state(self, name: str) -> np.ndarray:
        index = STATE_NAMES.index(name)
        return self.states[index]

    def save_csv(self, path: str | Path) -> Path:
        """Save all states and outputs as a plain numeric CSV file."""

        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_names = sorted(self.outputs)
        columns = [self.time_min, self.time_days]
        columns.extend(self.state(name) for name in STATE_NAMES)
        columns.extend(self.outputs[name] for name in output_names)
        header = ["time_min", "time_days", *STATE_NAMES, *output_names]
        np.savetxt(
            output_path,
            np.column_stack(columns),
            delimiter=",",
            header=",".join(header),
            comments="",
        )
        return output_path


def time_grid(days: float, samples_per_day: int = 48) -> np.ndarray:
    """Create an evenly spaced output grid expressed internally in minutes."""

    if days <= 0.0:
        raise ValueError("days must be positive")
    points = max(2, int(days * samples_per_day) + 1)
    return np.linspace(0.0, days * 1440.0, points)


def simulate(
    model: KidneyModel,
    *,
    days: float,
    y0: np.ndarray | None = None,
    samples_per_day: int = 48,
    method: str = "BDF",
    relative_tolerance: float = 1.0e-7,
    absolute_tolerance: float = 1.0e-9,
) -> SimulationResult:
    """Integrate the model with a stiff ODE solver."""

    initial = np.asarray(y0 if y0 is not None else initial_state(model.p), dtype=float)
    evaluation_times = time_grid(days, samples_per_day)
    solution = solve_ivp(
        model.rhs,
        (evaluation_times[0], evaluation_times[-1]),
        initial,
        method=method,
        t_eval=evaluation_times,
        rtol=relative_tolerance,
        atol=absolute_tolerance,
    )
    if not solution.success:
        return SimulationResult(
            time_min=solution.t,
            states=solution.y,
            outputs={},
            success=False,
            message=solution.message,
        )

    snapshots = [
        model.algebraic_outputs(time, solution.y[:, index])
        for index, time in enumerate(solution.t)
    ]
    output_names = snapshots[0].keys()
    outputs = {
        name: np.array([snapshot[name] for snapshot in snapshots], dtype=float)
        for name in output_names
    }
    return SimulationResult(
        time_min=solution.t,
        states=solution.y,
        outputs=outputs,
        success=True,
        message=solution.message,
    )
