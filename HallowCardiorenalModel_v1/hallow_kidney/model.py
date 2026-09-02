"""Coupled algebraic and differential equations for the kidney model."""

from __future__ import annotations

from math import exp, log, log10
from typing import Any

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .equations import (
    at1_renin_feedback,
    at1_vascular_effect,
    fmol_per_mL_to_pg_per_mL,
    myogenic_multiplier,
    positive,
    raw_md_renin_feedback,
    raw_rsna_renin_feedback,
    rihp_multiplier,
    safe_exp,
    tgf_multiplier,
)
from .inputs import ModelInputs
from .parameters import ModelParameters
from .states import initial_state, unpack_state


class KidneyModel:
    """Executable corrected Hallow model.

    Mean arterial pressure, right atrial pressure, sodium intake, and
    peritubular pressure are inputs. The state vector contains only quantities
    integrated in time. All other physiological quantities are recomputed from
    the current state by :meth:`algebraic_outputs`.
    """

    def __init__(
        self,
        parameters: ModelParameters,
        inputs: ModelInputs | None = None,
    ) -> None:
        self.p = parameters
        self.inputs = inputs or ModelInputs.baseline()
        allowed_modes = {"baseline_normalized", "published_dynamic", "fixed_unit"}
        if self.p.renin_feedback_mode not in allowed_modes:
            raise ValueError(
                f"Unknown renin_feedback_mode={self.p.renin_feedback_mode!r}; "
                f"choose one of {sorted(allowed_modes)}"
            )

        # The user-selected peptide initial conditions assume all three renin
        # feedback multipliers are one. The normalized mode preserves the shape
        # of each published feedback while dividing by its value at this agreed
        # reference condition.
        explicit_feedback_reference = (
            self.p.renin_md_reference_mEq_per_min,
            self.p.renin_rsna_reference,
            self.p.renin_at1_reference_pg_per_mL,
        )
        if any(value is not None for value in explicit_feedback_reference):
            if not all(value is not None for value in explicit_feedback_reference):
                raise ValueError(
                    "renin feedback reference values must be supplied together"
                )
            md_reference, rsna_reference, at1_reference = (
                float(value) for value in explicit_feedback_reference
            )
        else:
            reference = self.algebraic_outputs(
                0.0,
                initial_state(self.p),
                include_renin_feedback=False,
            )
            md_reference = reference["phi_md_sodium_mEq_per_min"]
            rsna_reference = reference["rsna"]
            at1_reference = reference["at1_bound_pg_per_mL"]
        self._raw_v_md_reference = raw_md_renin_feedback(md_reference, self.p)
        self._raw_v_rsna_reference = raw_rsna_renin_feedback(rsna_reference)
        self._raw_v_at1_reference = at1_renin_feedback(at1_reference, self.p)

    @property
    def feedback_reference_values(self) -> dict[str, float]:
        """Raw published feedback values used by baseline normalization."""

        return {
            "raw_v_md_reference": self._raw_v_md_reference,
            "raw_v_rsna_reference": self._raw_v_rsna_reference,
            "raw_v_at1_reference": self._raw_v_at1_reference,
        }

    def _input_values(
        self,
        time_min: float,
        overrides: dict[str, float] | None = None,
    ) -> dict[str, float]:
        values = {
            "p_ma_mmHg": float(self.inputs.p_ma_mmHg(time_min)),
            "p_ra_mmHg": float(self.inputs.p_ra_mmHg(time_min)),
            "sodium_intake_mEq_per_min": float(
                self.inputs.sodium_intake_mEq_per_min(time_min)
            ),
            "p_peritubular_mmHg": float(
                self.inputs.p_peritubular_mmHg(time_min)
            ),
        }
        if overrides:
            unknown = set(overrides) - set(values)
            if unknown:
                raise KeyError(f"Unknown kidney input overrides: {sorted(unknown)}")
            values.update({name: float(value) for name, value in overrides.items()})
        return values

    def autonomic_outputs(
        self,
        p_ma_mmHg: float,
        p_ra_mmHg: float,
        baro_integral_min: float,
    ) -> dict[str, float]:
        # Corrected normalization: a_auto is exactly one at MAP = 100 mmHg.
        a_auto = exp(-0.011 * (p_ma_mmHg - self.p.p_ma_ref_mmHg))
        a_chemo = 0.25 * a_auto

        # Corrected weighted adaptation: the weighted baroreceptor contribution
        # adapts toward 0.75, not toward 1.
        a_baro = 0.75 * a_auto - self.p.k_baro_per_min * baro_integral_min
        epsilon_aum = a_chemo + a_baro

        alpha_map = 0.5 + 1.1 / (
            1.0 + safe_exp((p_ma_mmHg - 100.0) / 15.0)
        )
        alpha_rap = 1.0 - 0.008 * p_ra_mmHg
        rsna = self.p.n_rsna * alpha_map * alpha_rap
        beta_rsna = 1.5 * (rsna - 1.0) + 1.0
        return {
            "a_auto": a_auto,
            "a_chemo": a_chemo,
            "a_baro": a_baro,
            "epsilon_aum": epsilon_aum,
            "epsilon_aum_adh": max(1.0, epsilon_aum),
            "alpha_map": alpha_map,
            "alpha_rap": alpha_rap,
            "rsna": rsna,
            "beta_rsna": beta_rsna,
        }

    def algebraic_outputs(
        self,
        time_min: float,
        y: np.ndarray,
        *,
        include_renin_feedback: bool = True,
        input_overrides: dict[str, float] | None = None,
    ) -> dict[str, float]:
        """Evaluate every non-state equation for one time and state vector."""

        p = self.p
        state = unpack_state(np.asarray(y, dtype=float))
        inputs = self._input_values(time_min, input_overrides)
        p_ma = inputs["p_ma_mmHg"]
        p_ra = inputs["p_ra_mmHg"]

        if state["v_ecf_L"] <= 0.0:
            raise ValueError("V_ECF must remain positive")

        c_sodium = state["m_sodium_mEq"] / state["v_ecf_L"]
        c_aldo = p.c_aldo_ref_ng_per_L * state["n_aldo"]
        c_adh = p.c_adh_ref_mU_per_L * state["n_adh"]
        at1_bound_pg_per_mL = fmol_per_mL_to_pg_per_mL(
            state["at1_bound_fmol_per_mL"],
            p.angii_molecular_weight_g_per_mol,
        )

        autonomic = self.autonomic_outputs(
            p_ma,
            p_ra,
            state["baro_integral_min"],
        )

        p_b_cvp = (
            p.p_b_cvp_slope * (p_ra - p.p_ra_ref_mmHg)
            if p.enable_cvp_bowman
            else 0.0
        )
        p_b_total = p.p_b_tubular_mmHg + p_b_cvp

        delta_rihp = (
            inputs["p_peritubular_mmHg"] - p.p_peritubular_ref_mmHg
        )
        if p.enable_rihp:
            gamma_rihp_pt = rihp_multiplier(delta_rihp, p.s_pn_pt, p)
            gamma_rihp_dt = rihp_multiplier(delta_rihp, p.s_pn_dt, p)
            gamma_rihp_cd = rihp_multiplier(delta_rihp, p.s_pn_cd, p)
        else:
            gamma_rihp_pt = gamma_rihp_dt = gamma_rihp_cd = 1.0

        psi_at1_preglom = at1_vascular_effect(
            at1_bound_pg_per_mL,
            p.a_at1_preglom,
            p.b_at1_preglom,
            p.c_at1_preglom,
            p.minimum_positive,
        )
        psi_at1_aa = at1_vascular_effect(
            at1_bound_pg_per_mL,
            p.a_at1_aa,
            p.b_at1_aa,
            p.c_at1_aa,
            p.minimum_positive,
        )
        psi_at1_ea = at1_vascular_effect(
            at1_bound_pg_per_mL,
            p.a_at1_ea,
            p.b_at1_ea,
            p.c_at1_ea,
            p.minimum_positive,
        )

        def renal_values(p_gh_mmHg: float) -> dict[str, float]:
            sigma_myo = myogenic_multiplier(p_gh_mmHg, p)
            net_filtration_pressure = p_gh_mmHg - p_b_total - p.p_go_mmHg
            sngfr_nL_per_min = max(
                0.0,
                p.kf_nL_per_min_per_mmHg_per_nephron
                * net_filtration_pressure,
            )
            gfr_L_per_min = sngfr_nL_per_min * p.n_nephrons * 1.0e-9
            phi_filtered_sodium = gfr_L_per_min * c_sodium

            gamma_filtered = 0.8 + 0.3 / (
                1.0
                + safe_exp(
                    1.0 + (phi_filtered_sodium - 14.0) / 138.0
                )
            )
            gamma_at1_tubule = 0.95 + 0.12 / (
                1.0
                + safe_exp(
                    2.6
                    - 1.8
                    * log10(
                        positive(at1_bound_pg_per_mL, p.minimum_positive)
                    )
                )
            )
            gamma_rsna = 0.5 + 0.7 / (
                1.0 + safe_exp((1.0 - autonomic["rsna"]) / 2.18)
            )
            eta_pt = (
                p.n_eta_pt
                * gamma_filtered
                * gamma_at1_tubule
                * gamma_rsna
                * gamma_rihp_pt
            )
            phi_pt_reab = phi_filtered_sodium * eta_pt
            phi_md_sodium = phi_filtered_sodium - phi_pt_reab
            sigma_tgf = tgf_multiplier(phi_md_sodium)

            r_preglom = (
                p.r_preglom_0_mmHg_min_per_L
                * autonomic["beta_rsna"]
                * psi_at1_preglom
                * sigma_myo
            )
            r_aa = (
                p.r_aa_0_mmHg_min_per_L_per_nephron
                * autonomic["beta_rsna"]
                * sigma_tgf
                * psi_at1_aa
                * sigma_myo
            )
            r_ea = p.r_ea_0_mmHg_min_per_L_per_nephron * psi_at1_ea
            r_renal = r_preglom + (r_aa + r_ea) / p.n_nephrons
            phi_renal_blood = p_ma / r_renal
            implied_p_gh = p_ma - phi_renal_blood * (
                r_preglom + r_aa / p.n_nephrons
            )
            return {
                "p_gh_mmHg": p_gh_mmHg,
                "implied_p_gh_mmHg": implied_p_gh,
                "p_gh_residual_mmHg": p_gh_mmHg - implied_p_gh,
                "sigma_myo": sigma_myo,
                "net_filtration_pressure_mmHg": net_filtration_pressure,
                "sngfr_nL_per_min": sngfr_nL_per_min,
                "gfr_L_per_min": gfr_L_per_min,
                "gfr_mL_per_min": gfr_L_per_min * 1000.0,
                "phi_filtered_sodium_mEq_per_min": phi_filtered_sodium,
                "gamma_filtered_sodium": gamma_filtered,
                "gamma_at1_tubule": gamma_at1_tubule,
                "gamma_rsna_tubule": gamma_rsna,
                "eta_pt": eta_pt,
                "phi_pt_reab_mEq_per_min": phi_pt_reab,
                "phi_md_sodium_mEq_per_min": phi_md_sodium,
                "sigma_tgf": sigma_tgf,
                "r_preglom_mmHg_min_per_L": r_preglom,
                "r_aa_single_mmHg_min_per_L": r_aa,
                "r_ea_single_mmHg_min_per_L": r_ea,
                "r_aa_equivalent_mmHg_min_per_L": r_aa / p.n_nephrons,
                "r_ea_equivalent_mmHg_min_per_L": r_ea / p.n_nephrons,
                "r_renal_mmHg_min_per_L": r_renal,
                "renal_blood_flow_L_per_min": phi_renal_blood,
            }

        # A positive myogenic multiplier requires Pgh > 48 mmHg for the
        # selected gain. The upper bound is just below arterial pressure.
        myogenic_zero = p.p_gh_nom_mmHg * (1.0 - 1.0 / p.c_gp_autoreg)
        lower = max(
            myogenic_zero + 1.0e-8,
            p_b_total + p.p_go_mmHg + 1.0e-8,
        )
        upper = p_ma - 1.0e-8
        if upper <= lower:
            raise ValueError(
                "No physical glomerular-pressure interval exists for the "
                f"current inputs: lower={lower:.6g}, upper={upper:.6g} mmHg"
            )

        residual_lower = renal_values(lower)["p_gh_residual_mmHg"]
        residual_upper = renal_values(upper)["p_gh_residual_mmHg"]
        if residual_lower * residual_upper <= 0.0:
            p_gh = brentq(
                lambda value: renal_values(value)["p_gh_residual_mmHg"],
                lower,
                upper,
                xtol=p.algebraic_tolerance_mmHg,
            )
        else:
            optimum = minimize_scalar(
                lambda value: abs(
                    renal_values(value)["p_gh_residual_mmHg"]
                ),
                bounds=(lower, upper),
                method="bounded",
            )
            p_gh = float(optimum.x)
            if abs(renal_values(p_gh)["p_gh_residual_mmHg"]) > 1.0e-5:
                raise RuntimeError(
                    "The renal algebraic loop has no acceptable solution at "
                    f"t={time_min:.6g} min"
                )

        renal = renal_values(p_gh)

        # Distal tubule and collecting duct sodium handling -----------------
        psi_aldo_dt = p.a_aldo_dt * positive(c_aldo, p.minimum_positive) ** p.b_aldo_dt
        eta_dt = p.n_eta_dt * psi_aldo_dt * gamma_rihp_dt
        phi_dt_reab = renal["phi_md_sodium_mEq_per_min"] * eta_dt
        phi_dt_sodium = renal["phi_md_sodium_mEq_per_min"] - phi_dt_reab

        lambda_dt = 0.82 + 0.39 / (
            1.0 + safe_exp((phi_dt_sodium - 1.6) / 2.0)
        )
        c_anp_hat = 7.427 - 6.554 / (
            1.0 + safe_exp((p_ra - 3.762) / 1.0)
        )
        lambda_anp = -0.1 * c_anp_hat + 1.1199
        lambda_aldo_cd = p.a_aldo_cd * positive(
            c_aldo, p.minimum_positive
        ) ** p.b_aldo_cd
        eta_cd = (
            p.n_eta_cd
            * lambda_dt
            * lambda_anp
            * lambda_aldo_cd
            * gamma_rihp_cd
        )
        phi_cd_reab = phi_dt_sodium * eta_cd
        phi_urine_sodium = phi_dt_sodium - phi_cd_reab

        # Water handling ------------------------------------------------------
        mu_aldo = 0.17 + 0.94 / (
            1.0
            + safe_exp(
                (0.48 - 1.2 * log10(positive(c_aldo, p.minimum_positive)))
                / 0.88
            )
        )
        mu_adh = 0.37 + 0.8 / (
            1.0
            + safe_exp(
                0.6 - 3.7 * log10(positive(c_adh, p.minimum_positive))
            )
        )
        phi_tubular_water_reab = (
            p.tubular_water_reab_base_L_per_min
            - 0.0011 / (mu_aldo * mu_adh)
            + 0.8 * renal["gfr_L_per_min"]
        )
        phi_urine = max(
            p.urine_flow_floor_L_per_min,
            renal["gfr_L_per_min"] - phi_tubular_water_reab,
        )
        phi_water_in = 0.0081 / (
            1.0
            + 1.822 * positive(c_adh, p.minimum_positive) ** -1.607
        ) - 0.00531

        # Aldosterone and ADH stimuli ----------------------------------------
        xi_k_sodium = max(
            0.0,
            (p.c_k_mEq_per_L / positive(c_sodium, p.minimum_positive))
            / 0.003525
            - 9.0,
        )
        xi_map = (
            69.03 * exp(-0.0425 * p_ma)
            if p_ma <= 100.0
            else 1.0
        )
        xi_at1 = 0.4 + 2.4 / (
            1.0
            + safe_exp(
                (
                    2.82
                    - 1.5
                    * log10(
                        positive(at1_bound_pg_per_mL, p.minimum_positive)
                    )
                )
                / 0.8
            )
        )
        n_aldo_stimulus = xi_k_sodium * xi_map * xi_at1

        delta_ra = (
            p.ra_adh_gain_per_mmHg * (p_ra - p.p_ra_ref_mmHg)
            - p.k_ra_adapt_per_min * state["ra_integral_min"]
        )
        sodium_adh_stimulus = max(0.0, c_sodium - 141.0)
        autonomic_adh_stimulus = autonomic["epsilon_aum_adh"] - 1.0
        n_adh_stimulus = (
            sodium_adh_stimulus + autonomic_adh_stimulus - delta_ra
        ) / 3.0

        # Renin feedbacks -----------------------------------------------------
        raw_v_md = raw_md_renin_feedback(
            renal["phi_md_sodium_mEq_per_min"], p
        )
        raw_v_rsna = raw_rsna_renin_feedback(autonomic["rsna"])
        raw_v_at1 = at1_renin_feedback(at1_bound_pg_per_mL, p)
        if not include_renin_feedback or p.renin_feedback_mode == "fixed_unit":
            v_md = v_rsna = v_at1 = 1.0
        elif p.renin_feedback_mode == "published_dynamic":
            v_md, v_rsna, v_at1 = raw_v_md, raw_v_rsna, raw_v_at1
        else:
            v_md = raw_v_md / self._raw_v_md_reference
            v_rsna = raw_v_rsna / self._raw_v_rsna_reference
            v_at1 = raw_v_at1 / self._raw_v_at1_reference

        renin_secretion = p.n_rs_pg_per_mL_per_min * v_md * v_rsna * v_at1
        pra_fmol_per_mL_per_min = (
            state["prc_pg_per_mL"] * p.x_prc_pra_fmol_per_h_per_pg / 60.0
        )

        outputs: dict[str, Any] = {
            **inputs,
            **autonomic,
            **renal,
            "c_sodium_mEq_per_L": c_sodium,
            "c_aldo_ng_per_L": c_aldo,
            "c_adh_mU_per_L": c_adh,
            "at1_bound_pg_per_mL": at1_bound_pg_per_mL,
            "p_b_cvp_mmHg": p_b_cvp,
            "p_b_total_mmHg": p_b_total,
            "delta_rihp_mmHg": delta_rihp,
            "gamma_rihp_pt": gamma_rihp_pt,
            "gamma_rihp_dt": gamma_rihp_dt,
            "gamma_rihp_cd": gamma_rihp_cd,
            "psi_at1_preglom": psi_at1_preglom,
            "psi_at1_aa": psi_at1_aa,
            "psi_at1_ea": psi_at1_ea,
            "psi_aldo_dt": psi_aldo_dt,
            "eta_dt": eta_dt,
            "phi_dt_reab_mEq_per_min": phi_dt_reab,
            "phi_dt_sodium_mEq_per_min": phi_dt_sodium,
            "lambda_dt": lambda_dt,
            "c_anp_hat": c_anp_hat,
            "lambda_anp": lambda_anp,
            "lambda_aldo_cd": lambda_aldo_cd,
            "eta_cd": eta_cd,
            "phi_cd_reab_mEq_per_min": phi_cd_reab,
            "phi_urine_sodium_mEq_per_min": phi_urine_sodium,
            "mu_aldo": mu_aldo,
            "mu_adh": mu_adh,
            "phi_tubular_water_reab_L_per_min": phi_tubular_water_reab,
            "phi_urine_L_per_min": phi_urine,
            "phi_water_in_L_per_min": phi_water_in,
            "xi_k_sodium": xi_k_sodium,
            "xi_map": xi_map,
            "xi_at1": xi_at1,
            "n_aldo_stimulus": n_aldo_stimulus,
            "delta_ra": delta_ra,
            "sodium_adh_stimulus": sodium_adh_stimulus,
            "autonomic_adh_stimulus": autonomic_adh_stimulus,
            "n_adh_stimulus": n_adh_stimulus,
            "raw_v_md": raw_v_md,
            "raw_v_rsna": raw_v_rsna,
            "raw_v_at1": raw_v_at1,
            "v_md": v_md,
            "v_rsna": v_rsna,
            "v_at1": v_at1,
            "renin_secretion_pg_per_mL_per_min": renin_secretion,
            "pra_fmol_per_mL_per_min": pra_fmol_per_mL_per_min,
        }
        return {name: float(value) for name, value in outputs.items()}

    def rhs(
        self,
        time_min: float,
        y: np.ndarray,
        *,
        input_overrides: dict[str, float] | None = None,
    ) -> np.ndarray:
        """Return the ODE derivatives in the ordering from ``STATE_NAMES``."""

        p = self.p
        state = unpack_state(np.asarray(y, dtype=float))
        out = self.algebraic_outputs(
            time_min,
            y,
            input_overrides=input_overrides,
        )

        d_m_sodium = (
            out["sodium_intake_mEq_per_min"]
            - out["phi_urine_sodium_mEq_per_min"]
        )
        d_v_ecf = out["phi_water_in_L_per_min"] - out["phi_urine_L_per_min"]
        d_prc = (
            out["renin_secretion_pg_per_mL_per_min"]
            - p.renin_clearance_per_min * state["prc_pg_per_mL"]
        )
        d_agt = (
            p.k_agt_fmol_per_mL_per_min
            - out["pra_fmol_per_mL_per_min"]
            - log(2.0) / p.h_agt_min * state["agt_fmol_per_mL"]
        )
        d_angi = (
            out["pra_fmol_per_mL_per_min"]
            - (p.c_ace_per_min + p.c_chym_per_min + p.c_nep_per_min)
            * state["angi_fmol_per_mL"]
            - log(2.0) / p.h_angi_min * state["angi_fmol_per_mL"]
        )
        d_angii = (
            (p.c_ace_per_min + p.c_chym_per_min)
            * state["angi_fmol_per_mL"]
            - (
                p.c_ace2_per_min
                + p.c_angii_to_angiv_per_min
                + p.c_at1_per_min
                + p.c_at2_per_min
            )
            * state["angii_fmol_per_mL"]
            - log(2.0) / p.h_angii_min * state["angii_fmol_per_mL"]
        )
        d_ang17 = (
            p.c_nep_per_min * state["angi_fmol_per_mL"]
            + p.c_ace2_per_min * state["angii_fmol_per_mL"]
            - log(2.0) / p.h_ang17_min * state["ang17_fmol_per_mL"]
        )
        d_angiv = (
            p.c_angii_to_angiv_per_min * state["angii_fmol_per_mL"]
            - log(2.0) / p.h_angiv_min * state["angiv_fmol_per_mL"]
        )
        d_at1 = (
            p.c_at1_per_min * state["angii_fmol_per_mL"]
            - log(2.0) / p.h_at1_min * state["at1_bound_fmol_per_mL"]
        )
        d_at2 = (
            p.c_at2_per_min * state["angii_fmol_per_mL"]
            - log(2.0) / p.h_at2_min * state["at2_bound_fmol_per_mL"]
        )
        d_n_aldo = (out["n_aldo_stimulus"] - state["n_aldo"]) / p.t_aldo_min
        d_n_adh = (out["n_adh_stimulus"] - state["n_adh"]) / p.t_adh_min
        d_baro_integral = out["a_baro"] - 0.75
        d_ra_integral = out["delta_ra"]

        return np.array(
            [
                d_m_sodium,
                d_v_ecf,
                d_prc,
                d_agt,
                d_angi,
                d_angii,
                d_ang17,
                d_angiv,
                d_at1,
                d_at2,
                d_n_aldo,
                d_n_adh,
                d_baro_integral,
                d_ra_integral,
            ],
            dtype=float,
        )
