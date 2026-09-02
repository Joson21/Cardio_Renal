import unittest
from dataclasses import replace

import numpy as np

from hallow_kidney import (
    KidneyModel,
    ModelInputs,
    corrected_hallow_parameters,
    cvp_extended_parameters,
    rihp_demo_parameters,
    simulate,
)
from hallow_kidney.equations import myogenic_multiplier, rihp_multiplier
from hallow_kidney.inputs import constant_input, step_input
from hallow_kidney.states import initial_state
from hallow_kidney.steady_state import isolated_raas_equilibrium
from hallow_kidney.validation import run_validation


class KidneyModelTests(unittest.TestCase):
    def test_corrected_resistances(self) -> None:
        p = corrected_hallow_parameters()
        self.assertEqual(p.r_aa_0_mmHg_min_per_L_per_nephron, 6.0e7)
        self.assertEqual(p.r_ea_0_mmHg_min_per_L_per_nephron, 1.0e8)

    def test_corrected_myogenic_baseline(self) -> None:
        p = corrected_hallow_parameters()
        self.assertEqual(myogenic_multiplier(60.0, p), 1.0)

    def test_agreed_raas_initial_conditions_are_equilibrium(self) -> None:
        p = corrected_hallow_parameters()
        equilibrium = isolated_raas_equilibrium(p)
        y0 = initial_state(p)
        expected = [
            equilibrium["prc_pg_per_mL"],
            equilibrium["agt_fmol_per_mL"],
            equilibrium["angi_fmol_per_mL"],
            equilibrium["angii_fmol_per_mL"],
            equilibrium["ang17_fmol_per_mL"],
            equilibrium["angiv_fmol_per_mL"],
            equilibrium["at1_bound_fmol_per_mL"],
            equilibrium["at2_bound_fmol_per_mL"],
        ]
        np.testing.assert_allclose(y0[2:10], expected, rtol=0.0, atol=1.0e-9)

    def test_cvp_extension_is_reference_centered(self) -> None:
        inputs = ModelInputs(
            constant_input(100.0),
            step_input(0.0, 5.0, 1.0),
            constant_input(0.126),
            constant_input(0.0),
        )
        model = KidneyModel(cvp_extended_parameters(), inputs)
        y0 = initial_state(model.p)
        self.assertEqual(model.algebraic_outputs(0.0, y0)["p_b_cvp_mmHg"], 0.0)
        self.assertAlmostEqual(
            model.algebraic_outputs(2.0, y0)["p_b_cvp_mmHg"],
            0.734,
            places=12,
        )

    def test_cvp_bowman_pathway_reduces_instantaneous_gfr(self) -> None:
        inputs = ModelInputs(
            constant_input(100.0),
            constant_input(5.0),
            constant_input(0.126),
            constant_input(0.0),
        )
        reproduction = KidneyModel(corrected_hallow_parameters(), inputs)
        extension = KidneyModel(cvp_extended_parameters(), inputs)
        y0 = initial_state(extension.p)
        without_bowman_path = reproduction.algebraic_outputs(0.0, y0)
        with_bowman_path = extension.algebraic_outputs(0.0, y0)
        self.assertLess(
            with_bowman_path["gfr_mL_per_min"],
            without_bowman_path["gfr_mL_per_min"],
        )

    def test_rihp_multiplier_is_baseline_neutral(self) -> None:
        p = rihp_demo_parameters()
        self.assertEqual(rihp_multiplier(0.0, 3.0, p), 1.0)
        self.assertAlmostEqual(
            rihp_multiplier(1.0, 3.0, p),
            0.306824264110,
            places=12,
        )

    def test_autonomic_weighting_and_adh_local_clamp_at_reference(self) -> None:
        model = KidneyModel(cvp_extended_parameters())
        output = model.algebraic_outputs(0.0, initial_state(model.p))
        self.assertEqual(output["a_chemo"], 0.25)
        self.assertEqual(output["a_baro"], 0.75)
        self.assertEqual(output["epsilon_aum"], 1.0)
        self.assertEqual(output["epsilon_aum_adh"], 1.0)

    def test_doubled_sodium_input_only_changes_protocol(self) -> None:
        profile = step_input(0.126, 0.252, start_min=10.0)
        self.assertEqual(profile(9.0), 0.126)
        self.assertEqual(profile(10.0), 0.252)

    def test_validation_suite_passes(self) -> None:
        model = KidneyModel(cvp_extended_parameters())
        self.assertTrue(all(check.passed for check in run_validation(model)))

    def test_published_feedback_option_is_not_silently_normalized(self) -> None:
        p = replace(
            corrected_hallow_parameters(),
            renin_feedback_mode="published_dynamic",
        )
        model = KidneyModel(p)
        out = model.algebraic_outputs(0.0, initial_state(p))
        self.assertNotEqual(out["v_md"], 1.0)
        self.assertNotEqual(out["v_rsna"], 1.0)

    def test_short_stiff_simulation_runs(self) -> None:
        model = KidneyModel(cvp_extended_parameters())
        result = simulate(model, days=0.02, samples_per_day=10)
        self.assertTrue(result.success)
        self.assertTrue(np.all(np.isfinite(result.states)))


if __name__ == "__main__":
    unittest.main()
