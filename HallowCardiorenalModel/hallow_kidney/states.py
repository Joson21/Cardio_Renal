"""ODE state ordering and agreed initial conditions."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .parameters import ModelParameters


STATE_NAMES = (
    "m_sodium_mEq",
    "v_ecf_L",
    "prc_pg_per_mL",
    "agt_fmol_per_mL",
    "angi_fmol_per_mL",
    "angii_fmol_per_mL",
    "ang17_fmol_per_mL",
    "angiv_fmol_per_mL",
    "at1_bound_fmol_per_mL",
    "at2_bound_fmol_per_mL",
    "n_aldo",
    "n_adh",
    "baro_integral_min",
    "ra_integral_min",
)

STATE_INDEX = {name: index for index, name in enumerate(STATE_NAMES)}


def initial_state(parameters: ModelParameters | None = None) -> np.ndarray:
    """Return the agreed baseline state vector.

    The peptide values solve the isolated Hallow systemic RAAS equilibrium for
    N_rs = 97 pg/mL/h and nu_MD = nu_RSNA = nu_AT1 = 1.
    """

    del parameters  # Reserved so future profiles can provide different states.
    return np.array(
        [
            2160.0,             # total extracellular sodium, mEq
            15.0,               # extracellular fluid volume, L
            27.988283793246,    # PRC, pg/mL
            474830.102634539,   # AGT, fmol/mL (= 474.8301026 pmol/mL)
            12.240564481657,    # AngI, fmol/mL
            6.458821606363,     # AngII, fmol/mL
            20.894402803236,    # Ang(1-7), fmol/mL
            1.824796522374,     # AngIV, fmol/mL
            21.990739367508,    # AT1-bound AngII, fmol/mL
            7.268125723160,     # AT2-bound AngII, fmol/mL
            1.0,                # normalized aldosterone
            1.0,                # normalized ADH
            0.0,                # corrected baroreceptor adaptation integral
            0.0,                # right-atrial-pressure adaptation integral
        ],
        dtype=float,
    )


def unpack_state(y: np.ndarray) -> dict[str, float]:
    """Convert the solver array into a readable name-to-value dictionary."""

    return {name: float(y[index]) for index, name in enumerate(STATE_NAMES)}


def pack_state(values: Mapping[str, float]) -> np.ndarray:
    """Create a solver array from a mapping containing every state name."""

    missing = [name for name in STATE_NAMES if name not in values]
    if missing:
        raise KeyError(f"Missing state values: {missing}")
    return np.array([values[name] for name in STATE_NAMES], dtype=float)
