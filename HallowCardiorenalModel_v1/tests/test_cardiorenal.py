import unittest

import numpy as np

from hallow_cardiorenal import (
    CoupledCardiorenalModel,
    corrected_hallow_cardiorenal_parameters,
    cvp_extended_cardiorenal_parameters,
    cvp_reference_equilibrium_state,
    hallow_table3_calibrated_cardiorenal_parameters,
    hallow_table3_reference_equilibrium_state,
    initial_state,
    published_literal_cardiorenal_parameters,
    simulate_coupled,
)
from hallow_cardiorenal.states import STATE_INDEX
from hallow_cardiorenal.validation import run_validation


class CardiorenalModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.model = CoupledCardiorenalModel(
            cvp_extended_cardiorenal_parameters()
        )

    def test_centered_nominal_cardiovascular_reference(self) -> None:
        output = self.model.algebraic_outputs(0.0, initial_state(self.model.p))
        self.assertAlmostEqual(output["p_ma_mmHg"], 100.0, places=9)
        self.assertAlmostEqual(output["p_ra_mmHg"], 0.0, places=9)
        self.assertAlmostEqual(
            output["cardiac_output_L_per_min"], 5.0, places=9
        )
        self.assertAlmostEqual(
            output["cardiovascular_algebraic_residual_L_per_min"],
            0.0,
            places=9,
        )

    def test_coupling_passes_computed_pressures_to_kidney(self) -> None:
        output = self.model.algebraic_outputs(0.0, initial_state(self.model.p))
        self.assertAlmostEqual(
            output["p_b_cvp_mmHg"],
            0.1468 * output["p_ra_mmHg"],
            places=12,
        )

    def test_more_ecf_volume_raises_map_pra_and_cardiac_output(self) -> None:
        lower = initial_state(self.model.p)
        upper = lower.copy()
        lower[STATE_INDEX["v_ecf_L"]] = 14.0
        upper[STATE_INDEX["v_ecf_L"]] = 16.0
        low_output = self.model.algebraic_outputs(0.0, lower)
        high_output = self.model.algebraic_outputs(0.0, upper)
        for name in (
            "blood_volume_L",
            "p_ma_mmHg",
            "p_ra_mmHg",
            "cardiac_output_L_per_min",
        ):
            self.assertGreater(high_output[name], low_output[name])

    def test_published_literal_mode_exposes_source_reference_conflict(self) -> None:
        literal = CoupledCardiorenalModel(
            published_literal_cardiorenal_parameters()
        )
        output = literal.algebraic_outputs(0.0, initial_state(literal.p))
        self.assertAlmostEqual(
            output["p_ra_mmHg"], output["p_ra_raw_mmHg"], places=12
        )
        self.assertGreater(output["p_ra_mmHg"], 0.0)
        self.assertNotAlmostEqual(output["p_ma_mmHg"], 100.0, places=3)

    def test_cvp_reference_state_is_a_full_equilibrium(self) -> None:
        state = cvp_reference_equilibrium_state()
        derivative = self.model.rhs(0.0, state)
        self.assertLess(np.max(np.abs(derivative)), 1.0e-8)
        self.assertTrue(
            all(
                check.passed
                for check in run_validation(
                    self.model,
                    equilibrium_state=state,
                )
            )
        )

    def test_cvp_off_and_on_are_distinct_named_configurations(self) -> None:
        off = corrected_hallow_cardiorenal_parameters()
        on = cvp_extended_cardiorenal_parameters()
        self.assertFalse(off.kidney.enable_cvp_bowman)
        self.assertTrue(on.kidney.enable_cvp_bowman)
        self.assertNotEqual(off.configuration_name, on.configuration_name)

    def test_short_stiff_coupled_simulation_runs(self) -> None:
        result = simulate_coupled(
            self.model,
            days=0.02,
            y0=cvp_reference_equilibrium_state(),
            samples_per_day=10,
        )
        self.assertTrue(result.success)
        self.assertTrue(np.all(np.isfinite(result.states)))
        self.assertTrue(np.all(np.isfinite(result.outputs["p_ma_mmHg"])))

    def test_hallow_table3_calibrated_equilibrium(self) -> None:
        model = CoupledCardiorenalModel(
            hallow_table3_calibrated_cardiorenal_parameters()
        )
        state = hallow_table3_reference_equilibrium_state()
        output = model.algebraic_outputs(0.0, state)
        self.assertLess(np.max(np.abs(model.rhs(0.0, state))), 1.0e-10)
        targets = {
            "p_ma_mmHg": 83.0,
            "p_ra_mmHg": 0.9022132710207813,
            "cardiac_output_L_per_min": 5.15,
            "gfr_mL_per_min": 99.0,
            "renal_blood_flow_L_per_min": 0.9,
            "c_sodium_mEq_per_L": 143.3,
            "phi_urine_sodium_mEq_per_min": 0.126,
            "p_b_cvp_mmHg": 0.0,
            "p_b_total_mmHg": 18.0,
            "c_aldo_ng_per_L": 100.0,
        }
        for name, target in targets.items():
            self.assertAlmostEqual(output[name], target, places=9)
        self.assertTrue(
            all(
                check.passed
                for check in run_validation(model, equilibrium_state=state)
            )
        )
        self.assertAlmostEqual(
            output["phi_urine_L_per_min"] * 1440.0,
            2.1,
            places=9,
        )
        self.assertAlmostEqual(
            output["pra_fmol_per_mL_per_min"] * 60.0,
            350.0,
            places=9,
        )

    def test_coupled_output_rsna_uses_solved_pra(self) -> None:
        model = CoupledCardiorenalModel(
            hallow_table3_calibrated_cardiorenal_parameters()
        )
        state = hallow_table3_reference_equilibrium_state()
        output = model.algebraic_outputs(0.0, state)
        kidney_state = state[:-1]
        expected = model.kidney.algebraic_outputs(
            0.0,
            kidney_state,
            input_overrides={
                "p_ma_mmHg": output["p_ma_mmHg"],
                "p_ra_mmHg": output["p_ra_mmHg"],
            },
        )
        self.assertAlmostEqual(output["rsna"], expected["rsna"], places=12)


if __name__ == "__main__":
    unittest.main()
