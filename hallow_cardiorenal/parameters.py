"""Parameters and named configurations for the coupled model."""

from __future__ import annotations

from dataclasses import dataclass
from math import exp

from hallow_kidney.parameters import (
    ModelParameters,
    corrected_hallow_parameters,
    cvp_extended_parameters,
    hallow_table3_calibrated_parameters,
    rihp_demo_parameters,
)


@dataclass(frozen=True)
class CardiovascularParameters:
    """Karaaslan/Hallow long-term heart and systemic-vascular parameters.

    ``baseline_centered`` is the recommended executable reconstruction. It
    retains the published curve shapes but resolves the paper's conflict
    between the exponential right-atrial-pressure equation and the tabulated
    reference ``P_ra = 0 mmHg``. ``published_literal`` evaluates Eqs. A18-A26
    without those reference corrections.
    """

    closure_mode: str = "baseline_centered"

    # Blood-volume relationship, Karaaslan Eq. 32 / Hallow Eq. A18.
    vb_offset_L: float = 4.560227
    vb_amplitude_L: float = 2.431217
    vb_midpoint_L: float = 18.11278
    vb_slope_per_L: float = 0.47437

    # Mean filling pressure, Eq. 33 / A19.
    pmf_slope_mmHg_per_L: float = 7.436
    pmf_offset_mmHg: float = 30.18

    # Frank-Starling/right-atrial relationship, Eq. 35 / A22.
    pra_scale_mmHg: float = 0.2787
    pra_flow_coefficient_min_per_L: float = 0.2281

    # Vascularity, arterial resistance, and venous resistance.
    vascularity_formation_prefactor_per_min: float = 11.312 / 100000.0
    vascularity_flow_coefficient_min_per_L: float = 0.4799
    vascularity_destruction_per_min: float = 0.00001
    arterial_resistance_constant_mmHg_min_per_L: float = 16.6
    basic_venous_resistance_mmHg_min_per_L: float = 3.4
    venous_return_denominator: float = 31.0

    # Agreed reference used by the baseline-centered reconstruction.
    v_ecf_ref_L: float = 15.0
    vascularity_ref: float = 1.0
    cardiac_output_ref_L_per_min: float = 5.0
    p_ma_ref_mmHg: float = 100.0
    p_ra_ref_mmHg: float = 0.0

    # Algebraic solver controls.
    p_ma_min_mmHg: float = 20.0
    p_ma_max_mmHg: float = 250.0
    algebraic_tolerance: float = 1.0e-10

    @property
    def p_ra_raw_reference_mmHg(self) -> float:
        return self.pra_scale_mmHg * exp(
            self.pra_flow_coefficient_min_per_L
            * self.cardiac_output_ref_L_per_min
        )

    @property
    def blood_volume_reference_L(self) -> float:
        exponent = -(
            self.v_ecf_ref_L - self.vb_midpoint_L
        ) * self.vb_slope_per_L
        return self.vb_offset_L + self.vb_amplitude_L / (1.0 + exp(exponent))

    @property
    def venous_return_resistance_reference_mmHg_min_per_L(self) -> float:
        arterial_resistance_reference = (
            self.arterial_resistance_constant_mmHg_min_per_L
            / self.vascularity_ref
        )
        return (
            8.0 * self.basic_venous_resistance_mmHg_min_per_L
            + arterial_resistance_reference
        ) / self.venous_return_denominator

    @property
    def mean_filling_pressure_reference_shift_mmHg(self) -> float:
        raw_reference = (
            self.pmf_slope_mmHg_per_L * self.blood_volume_reference_L
            - self.pmf_offset_mmHg
        )
        required_reference = (
            self.cardiac_output_ref_L_per_min
            * self.venous_return_resistance_reference_mmHg_min_per_L
            + self.p_ra_ref_mmHg
        )
        return required_reference - raw_reference


@dataclass(frozen=True)
class CardiorenalParameters:
    """Complete, explicit parameter bundle for one coupled configuration."""

    configuration_name: str
    kidney: ModelParameters
    cardiovascular: CardiovascularParameters


def corrected_hallow_cardiorenal_parameters() -> CardiorenalParameters:
    """Corrected Hallow kidney plus baseline-centered cardiovascular closure."""

    return CardiorenalParameters(
        configuration_name="corrected_hallow_cardiorenal",
        kidney=corrected_hallow_parameters(),
        cardiovascular=CardiovascularParameters(),
    )


def cvp_extended_cardiorenal_parameters() -> CardiorenalParameters:
    """Recommended coupled configuration with the fixed CVP-Bowman novelty."""

    return CardiorenalParameters(
        configuration_name="cvp_extended_cardiorenal",
        kidney=cvp_extended_parameters(),
        cardiovascular=CardiovascularParameters(),
    )


def hallow_table3_calibrated_cardiorenal_parameters() -> CardiorenalParameters:
    """Healthy coupled profile jointly calibrated to Hallow Table 3.

    The source-reproduction profiles above are retained unchanged.  This
    profile uses the paper's MAP and cardiac-output targets as the exact
    cardiovascular reference and the positive pressure predicted by its
    printed right-atrial-pressure equation at that cardiac output.
    """

    return CardiorenalParameters(
        configuration_name="hallow_table3_calibrated_cardiorenal",
        kidney=hallow_table3_calibrated_parameters(),
        cardiovascular=CardiovascularParameters(
            v_ecf_ref_L=15.0,
            vascularity_ref=1.3053901358986109,
            cardiac_output_ref_L_per_min=5.15,
            p_ma_ref_mmHg=83.0,
            p_ra_ref_mmHg=0.9022132710207813,
        ),
    )


def rihp_demo_cardiorenal_parameters() -> CardiorenalParameters:
    """Coupled CVP model plus the provisional RIHP transport multipliers."""

    return CardiorenalParameters(
        configuration_name="cvp_and_rihp_cardiorenal_demo",
        kidney=rihp_demo_parameters(),
        cardiovascular=CardiovascularParameters(),
    )


def published_literal_cardiorenal_parameters() -> CardiorenalParameters:
    """Literal cardiovascular equations for source-comparison diagnostics."""

    return CardiorenalParameters(
        configuration_name="published_literal_cardiovascular_comparison",
        kidney=corrected_hallow_parameters(),
        cardiovascular=CardiovascularParameters(closure_mode="published_literal"),
    )
