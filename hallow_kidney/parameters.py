"""Parameter definitions for the corrected and extended model configurations.

All rate constants in this file use the model's internal time unit: minutes.
The source papers often report RAAS rates in h^-1; those values are converted
here by dividing by 60. Parameter comments identify corrections or derived
values that differ from a literal reading of the published tables.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import log


@dataclass(frozen=True)
class ModelParameters:
    """Flat, explicit parameter collection.

    A flat dataclass is intentionally used so a new Python user can search for
    a symbol and change it without navigating several nested objects.
    """

    # Configuration switches -------------------------------------------------
    configuration_name: str = "corrected_hallow"
    enable_cvp_bowman: bool = False
    enable_rihp: bool = False
    renin_feedback_mode: str = "baseline_normalized"

    # Optional physical reference point for normalized renin feedbacks.  When
    # all three values are supplied, the published feedback curves are divided
    # by their values at this point instead of at ``initial_state()``.  This is
    # useful for a separately calibrated baseline while leaving the original
    # reconstruction unchanged.
    renin_md_reference_mEq_per_min: float | None = None
    renin_rsna_reference: float | None = None
    renin_at1_reference_pg_per_mL: float | None = None

    # Reference inputs --------------------------------------------------------
    p_ma_ref_mmHg: float = 100.0
    p_ra_ref_mmHg: float = 0.0
    p_peritubular_ref_mmHg: float = 0.0
    sodium_intake_ref_mEq_per_min: float = 0.126

    # Renal hemodynamics -------------------------------------------------------
    n_nephrons: float = 2.0e6
    r_preglom_0_mmHg_min_per_L: float = 15.0
    # Corrected factor-of-ten values: 30 and 50 mmHg min/L multiplied by 2e6.
    r_aa_0_mmHg_min_per_L_per_nephron: float = 6.0e7
    r_ea_0_mmHg_min_per_L_per_nephron: float = 1.0e8
    kf_nL_per_min_per_mmHg_per_nephron: float = 3.905
    p_b_tubular_mmHg: float = 18.0
    p_go_mmHg: float = 28.0
    p_gh_nom_mmHg: float = 60.0
    c_gp_autoreg: float = 5.0

    # AT1-bound AngII vascular effects ---------------------------------------
    a_at1_preglom: float = 0.8
    b_at1_preglom: float = 0.055
    c_at1_preglom: float = 0.185
    a_at1_aa: float = 0.8
    b_at1_aa: float = 0.055
    c_at1_aa: float = 0.185
    a_at1_ea: float = 0.925
    b_at1_ea: float = 0.05
    c_at1_ea: float = 0.17

    # Renin feedback ----------------------------------------------------------
    n_rs_pg_per_mL_per_min: float = 97.0 / 60.0
    a_md_renin: float = 0.2262
    b_md_renin: float = 28.04
    c_md_renin: float = 2.0
    offset_md_renin_mEq_per_min: float = 5.0
    x_md_renin_mEq_per_min: float = 1.15
    a_at1_renin: float = 0.0102
    b_at1_renin: float = 0.95
    # Derived so the agreed AT1 equilibrium produces nu_AT1 = 1 exactly.
    at1_feedback_eq_pg_per_mL: float = 22.4444725
    h_renin_min: float = 12.0
    x_prc_pra_fmol_per_h_per_pg: float = 61.0

    # Systemic RAAS kinetics (all converted to min^-1) ------------------------
    k_agt_fmol_per_mL_per_min: float = 34620.0 / 60.0
    h_agt_min: float = 10.0 * 60.0
    c_ace_per_min: float = 54.1 / 60.0
    c_chym_per_min: float = 1.1 / 60.0
    c_nep_per_min: float = 1.1 / 60.0
    c_ace2_per_min: float = 2.4 / 60.0
    c_angii_to_angiv_per_min: float = 23.5 / 60.0
    c_at1_per_min: float = 11.8 / 60.0
    c_at2_per_min: float = 3.9 / 60.0
    h_angi_min: float = 0.5
    h_angii_min: float = 0.66
    h_ang17_min: float = 30.0
    h_angiv_min: float = 0.5
    h_at1_min: float = 12.0
    h_at2_min: float = 12.0
    angii_molecular_weight_g_per_mol: float = 1046.18

    # Tubular sodium handling -------------------------------------------------
    n_eta_pt: float = 0.8
    n_eta_dt: float = 0.5
    n_eta_cd: float = 0.93
    a_aldo_dt: float = 0.24
    b_aldo_dt: float = 0.3
    a_aldo_cd: float = 0.74
    b_aldo_cd: float = 0.06

    # Aldosterone, ADH, and autonomic adaptation -----------------------------
    c_k_mEq_per_L: float = 5.0
    c_aldo_ref_ng_per_L: float = 85.0
    c_adh_ref_mU_per_L: float = 4.0
    t_aldo_min: float = 30.0
    t_adh_min: float = 6.0
    k_baro_per_min: float = 0.0000667
    k_ra_adapt_per_min: float = 0.0007
    ra_adh_gain_per_mmHg: float = 0.2
    urine_flow_floor_L_per_min: float = 0.0003
    tubular_water_reab_base_L_per_min: float = 0.0251
    n_rsna: float = 1.0

    # CVP novelty -------------------------------------------------------------
    p_b_cvp_slope: float = 0.1468

    # Optional RIHP example ---------------------------------------------------
    rihp_pressure_scale_mmHg: float = 1.0
    s_pn_pt: float = 3.0
    s_pn_dt: float = 3.0
    s_pn_cd: float = 3.0

    # Numerical safeguards ----------------------------------------------------
    minimum_positive: float = 1.0e-12
    algebraic_tolerance_mmHg: float = 1.0e-9

    @property
    def renin_clearance_per_min(self) -> float:
        return log(2.0) / self.h_renin_min


def corrected_hallow_parameters() -> ModelParameters:
    """Corrected Hallow reproduction, without the new CVP or RIHP pathways."""

    return ModelParameters()


def cvp_extended_parameters() -> ModelParameters:
    """Corrected Hallow model with the fixed CVP-to-Bowman pressure novelty."""

    return replace(
        corrected_hallow_parameters(),
        configuration_name="cvp_extended",
        enable_cvp_bowman=True,
    )


def hallow_table3_calibrated_parameters() -> ModelParameters:
    """Healthy baseline calibrated jointly to Hallow Table 3.

    This is intentionally separate from :func:`corrected_hallow_parameters`.
    The values are the result of the constrained baseline calibration described
    in ``CALIBRATION_REPORT.md``; source-paper values remain available through
    the original constructors.
    """

    return replace(
        cvp_extended_parameters(),
        configuration_name="hallow_table3_calibrated",
        p_ma_ref_mmHg=83.0,
        p_ra_ref_mmHg=0.9022132710207813,
        n_rs_pg_per_mL_per_min=0.3314228322349465,
        at1_feedback_eq_pg_per_mL=15.8,
        renin_md_reference_mEq_per_min=5.6611005854660235,
        renin_rsna_reference=1.3224841838351185,
        renin_at1_reference_pg_per_mL=4.716373441108469,
        r_preglom_0_mmHg_min_per_L=5.917079624209335,
        r_aa_0_mmHg_min_per_L_per_nephron=2.3668318496837342e7,
        r_ea_0_mmHg_min_per_L_per_nephron=1.1592660660155734e8,
        n_eta_cd=1.0346820276676627,
        c_aldo_ref_ng_per_L=87.62718001985094,
        c_adh_ref_mU_per_L=5.210847998533489,
        tubular_water_reab_base_L_per_min=0.019385600350148517,
    )


def rihp_demo_parameters() -> ModelParameters:
    """CVP model plus the provisional segment-specific RIHP multipliers."""

    return replace(
        cvp_extended_parameters(),
        configuration_name="cvp_and_rihp_demo",
        enable_rihp=True,
    )


def with_published_renin_feedback(parameters: ModelParameters) -> ModelParameters:
    """Return a copy using the literal, unnormalized Hallow feedback equations.

    The published functions do not evaluate to one at the agreed initial state.
    This option is retained for equation-level reproduction and comparison.
    """

    return replace(parameters, renin_feedback_mode="published_dynamic")
