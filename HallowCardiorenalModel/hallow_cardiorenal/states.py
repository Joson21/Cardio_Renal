"""State ordering and initial conditions for the coupled model."""

from __future__ import annotations

import numpy as np

from hallow_kidney.states import STATE_NAMES as KIDNEY_STATE_NAMES
from hallow_kidney.states import initial_state as kidney_initial_state

from .parameters import CardiorenalParameters


CARDIOVASCULAR_STATE_NAMES = ("vascularity",)
STATE_NAMES = (*KIDNEY_STATE_NAMES, *CARDIOVASCULAR_STATE_NAMES)
STATE_INDEX = {name: index for index, name in enumerate(STATE_NAMES)}
KIDNEY_STATE_COUNT = len(KIDNEY_STATE_NAMES)


def initial_state(parameters: CardiorenalParameters | None = None) -> np.ndarray:
    """Return the nominal paper/agreed state before whole-system settling."""

    kidney_parameters = parameters.kidney if parameters is not None else None
    return np.concatenate(
        [
            kidney_initial_state(kidney_parameters),
            np.array([1.0], dtype=float),
        ]
    )


def cvp_reference_equilibrium_state() -> np.ndarray:
    """Return the verified equilibrium of the default CVP-coupled profile.

    This state was obtained by 300 days of stiff settling followed by a
    constrained nonlinear residual solve. It belongs specifically to
    ``cvp_extended_cardiorenal_parameters()`` with baseline external inputs.
    Recalculate it whenever equations or parameters change.
    """

    return np.array(
        [
            2078.06772219359,
            14.5538978853956,
            31.2084493147996,
            471996.215368979,
            13.6488910514231,
            7.20193521776991,
            23.2983885566133,
            2.03474676043154,
            24.5208630874133,
            8.10435305431458,
            1.39843091679175,
            0.601991304698793,
            976.073779120935,
            -8.9707788694985,
            1.08021901986482,
        ],
        dtype=float,
    )


def corrected_hallow_reference_equilibrium_state() -> np.ndarray:
    """Verified equilibrium for the corrected configuration with CVP off."""

    return np.array(
        [
            2078.14079361099,
            14.5544071921238,
            31.2095064306459,
            471995.285060547,
            13.6493533768296,
            7.20217916708525,
            23.2991777369884,
            2.03481568288976,
            24.5216936763565,
            8.1046275709992,
            1.39796705524511,
            0.601991304698793,
            974.990854776849,
            -8.96107189235108,
            1.08012715685293,
        ],
        dtype=float,
    )


def split_state(y: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray(y, dtype=float)
    if values.shape != (len(STATE_NAMES),):
        raise ValueError(
            f"Expected {len(STATE_NAMES)} coupled states, got shape {values.shape}"
        )
    return values[:KIDNEY_STATE_COUNT], float(values[-1])
