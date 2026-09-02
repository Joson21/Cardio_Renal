"""Bidirectionally coupled long-term cardiovascular-kidney model."""

from __future__ import annotations

import numpy as np

from hallow_kidney import KidneyModel, ModelInputs
from hallow_kidney.inputs import constant_input

from .cardiovascular import CardiovascularSubsystem
from .inputs import CardiorenalInputs
from .parameters import CardiorenalParameters
from .states import STATE_NAMES, split_state


class CoupledCardiorenalModel:
    """Close MAP and Pra around the existing kidney/RAAS module."""

    def __init__(
        self,
        parameters: CardiorenalParameters,
        inputs: CardiorenalInputs | None = None,
    ) -> None:
        self.p = parameters
        self.inputs = inputs or CardiorenalInputs.baseline()
        kidney_inputs = ModelInputs(
            p_ma_mmHg=constant_input(parameters.cardiovascular.p_ma_ref_mmHg),
            p_ra_mmHg=constant_input(parameters.cardiovascular.p_ra_ref_mmHg),
            sodium_intake_mEq_per_min=self.inputs.sodium_intake_mEq_per_min,
            p_peritubular_mmHg=self.inputs.p_peritubular_mmHg,
        )
        self.kidney = KidneyModel(parameters.kidney, kidney_inputs)
        self.cardiovascular = CardiovascularSubsystem(
            parameters.cardiovascular,
            self.kidney,
        )

    def _coupling_values(
        self,
        kidney_state: np.ndarray,
        vascularity: float,
    ) -> dict[str, float]:
        return self.cardiovascular.solve(
            v_ecf_L=float(kidney_state[1]),
            vascularity=vascularity,
            baro_integral_min=float(kidney_state[12]),
        )

    def algebraic_outputs(self, time_min: float, y: np.ndarray) -> dict[str, float]:
        kidney_state, vascularity = split_state(y)
        cardio = self._coupling_values(kidney_state, vascularity)
        kidney = self.kidney.algebraic_outputs(
            time_min,
            kidney_state,
            input_overrides={
                "p_ma_mmHg": cardio["p_ma_mmHg"],
                "p_ra_mmHg": cardio["p_ra_mmHg"],
            },
        )
        # Kidney values take precedence for shared autonomic names because its
        # evaluation uses the solved right-atrial pressure.  The cardiovascular
        # algebraic loop uses a Pra=0 placeholder for epsilon_aum-only work.
        outputs = {**cardio, **kidney}
        return {name: float(value) for name, value in outputs.items()}

    def rhs(self, time_min: float, y: np.ndarray) -> np.ndarray:
        kidney_state, vascularity = split_state(y)
        cardio = self._coupling_values(kidney_state, vascularity)
        kidney_derivative = self.kidney.rhs(
            time_min,
            kidney_state,
            input_overrides={
                "p_ma_mmHg": cardio["p_ma_mmHg"],
                "p_ra_mmHg": cardio["p_ra_mmHg"],
            },
        )
        derivative = np.concatenate(
            [kidney_derivative, np.array([cardio["d_vascularity_per_min"]])]
        )
        if derivative.shape != (len(STATE_NAMES),):
            raise RuntimeError("Coupled derivative ordering does not match STATE_NAMES")
        return derivative
