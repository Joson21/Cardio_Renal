"""Karaaslan/Hallow long-term cardiac and systemic-vascular closure."""

from __future__ import annotations

from math import exp

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from hallow_kidney.model import KidneyModel

from .parameters import CardiovascularParameters


class CardiovascularSubsystem:
    """Evaluate Eqs. 32-43 (Karaaslan) / A18-A26 (Hallow).

    The source calls this a heart/vasculature block, but it is a long-term
    lumped closure. Cardiac output, right-atrial pressure, and MAP are
    algebraic variables; vascularity is its sole additional ODE state.
    """

    def __init__(
        self,
        parameters: CardiovascularParameters,
        kidney_model: KidneyModel,
    ) -> None:
        self.p = parameters
        self.kidney_model = kidney_model
        if self.p.closure_mode not in {"baseline_centered", "published_literal"}:
            raise ValueError(
                "closure_mode must be 'baseline_centered' or 'published_literal'"
            )

    def blood_volume_L(self, v_ecf_L: float) -> float:
        p = self.p
        exponent = -(v_ecf_L - p.vb_midpoint_L) * p.vb_slope_per_L
        return p.vb_offset_L + p.vb_amplitude_L / (1.0 + exp(exponent))

    def _values_at_pressure(
        self,
        p_ma_mmHg: float,
        v_ecf_L: float,
        vascularity: float,
        baro_integral_min: float,
    ) -> dict[str, float]:
        p = self.p
        if vascularity <= 0.0:
            raise ValueError("vascularity must remain positive")

        # Pra does not enter epsilon_aum, so a placeholder is sufficient while
        # solving the cardiovascular algebraic loop.
        autonomic = self.kidney_model.autonomic_outputs(
            p_ma_mmHg,
            0.0,
            baro_integral_min,
        )
        epsilon_aum = autonomic["epsilon_aum"]
        blood_volume = self.blood_volume_L(v_ecf_L)
        arterial_resistance = (
            p.arterial_resistance_constant_mmHg_min_per_L
            / vascularity
            * epsilon_aum
        )
        total_peripheral_resistance = (
            arterial_resistance + p.basic_venous_resistance_mmHg_min_per_L
        )
        cardiac_output = p_ma_mmHg / total_peripheral_resistance

        p_ra_raw = p.pra_scale_mmHg * exp(
            p.pra_flow_coefficient_min_per_L * cardiac_output
        )
        if p.closure_mode == "baseline_centered":
            p_ra = p_ra_raw - p.p_ra_raw_reference_mmHg + p.p_ra_ref_mmHg
            pmf_shift = p.mean_filling_pressure_reference_shift_mmHg
        else:
            p_ra = p_ra_raw
            pmf_shift = 0.0

        mean_filling_pressure = (
            p.pmf_slope_mmHg_per_L * blood_volume
            - p.pmf_offset_mmHg
        ) * epsilon_aum + pmf_shift
        venous_return_resistance = (
            8.0 * p.basic_venous_resistance_mmHg_min_per_L
            + arterial_resistance
        ) / p.venous_return_denominator
        venous_return = (
            mean_filling_pressure - p_ra
        ) / venous_return_resistance
        residual = cardiac_output - venous_return
        return {
            **autonomic,
            "blood_volume_L": blood_volume,
            "mean_filling_pressure_mmHg": mean_filling_pressure,
            "arterial_resistance_mmHg_min_per_L": arterial_resistance,
            "basic_venous_resistance_mmHg_min_per_L": (
                p.basic_venous_resistance_mmHg_min_per_L
            ),
            "total_peripheral_resistance_mmHg_min_per_L": (
                total_peripheral_resistance
            ),
            "venous_return_resistance_mmHg_min_per_L": venous_return_resistance,
            "cardiac_output_L_per_min": cardiac_output,
            "venous_return_L_per_min": venous_return,
            "p_ra_raw_mmHg": p_ra_raw,
            "p_ra_mmHg": p_ra,
            "p_ma_mmHg": p_ma_mmHg,
            "cardiovascular_algebraic_residual_L_per_min": residual,
        }

    def solve(
        self,
        v_ecf_L: float,
        vascularity: float,
        baro_integral_min: float,
    ) -> dict[str, float]:
        """Solve the instantaneous closed cardiac/vascular algebraic loop."""

        p = self.p

        def residual(p_ma: float) -> float:
            return self._values_at_pressure(
                p_ma,
                v_ecf_L,
                vascularity,
                baro_integral_min,
            )["cardiovascular_algebraic_residual_L_per_min"]

        root: float | None = None
        lower_residual = residual(p.p_ma_min_mmHg)
        upper_residual = residual(p.p_ma_max_mmHg)
        if lower_residual == 0.0:
            root = p.p_ma_min_mmHg
        elif lower_residual * upper_residual < 0.0:
            # The physiological cases normally bracket a unique solution over
            # the full configured interval. Trying this first avoids an
            # unnecessary dense scan during every stiff-solver evaluation.
            root = float(
                brentq(
                    residual,
                    p.p_ma_min_mmHg,
                    p.p_ma_max_mmHg,
                    xtol=p.algebraic_tolerance,
                )
            )
        else:
            # Retain a grid search for unusual parameter sets whose residual
            # can cross zero more than once inside the interval.
            grid = np.linspace(p.p_ma_min_mmHg, p.p_ma_max_mmHg, 180)
            values = np.array([residual(value) for value in grid])
            for left, right, f_left, f_right in zip(
                grid[:-1], grid[1:], values[:-1], values[1:]
            ):
                if f_left == 0.0:
                    root = float(left)
                    break
                if f_left * f_right < 0.0:
                    root = float(
                        brentq(
                            residual,
                            float(left),
                            float(right),
                            xtol=p.algebraic_tolerance,
                        )
                    )
                    break

        if root is None:
            optimum = minimize_scalar(
                lambda value: abs(residual(float(value))),
                bounds=(p.p_ma_min_mmHg, p.p_ma_max_mmHg),
                method="bounded",
            )
            root = float(optimum.x)
            if abs(residual(root)) > 1.0e-7:
                raise RuntimeError(
                    "The cardiovascular algebraic loop has no acceptable "
                    f"solution for V_ECF={v_ecf_L:.6g} L, "
                    f"vascularity={vascularity:.6g}"
                )

        output = self._values_at_pressure(
            root,
            v_ecf_L,
            vascularity,
            baro_integral_min,
        )
        cardiac_output = output["cardiac_output_L_per_min"]
        if p.closure_mode == "baseline_centered":
            vascularity_formation = (
                p.vascularity_destruction_per_min
                * p.vascularity_ref
                * exp(
                    -p.vascularity_flow_coefficient_min_per_L
                    * (cardiac_output - p.cardiac_output_ref_L_per_min)
                )
            )
        else:
            vascularity_formation = (
                p.vascularity_formation_prefactor_per_min
                * exp(
                    -p.vascularity_flow_coefficient_min_per_L
                    * cardiac_output
                )
            )
        vascularity_destruction = (
            p.vascularity_destruction_per_min * vascularity
        )
        output.update(
            {
                "vascularity": float(vascularity),
                "vascularity_formation_per_min": vascularity_formation,
                "vascularity_destruction_per_min": vascularity_destruction,
                "d_vascularity_per_min": (
                    vascularity_formation - vascularity_destruction
                ),
            }
        )
        return {name: float(value) for name, value in output.items()}
