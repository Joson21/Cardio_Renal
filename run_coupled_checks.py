"""Run fast mathematical checks for the coupled cardiorenal model."""

from hallow_cardiorenal import (
    CoupledCardiorenalModel,
    hallow_table3_calibrated_cardiorenal_parameters,
    hallow_table3_reference_equilibrium_state,
)
from hallow_cardiorenal.steady_state import derivative_report
from hallow_cardiorenal.validation import assert_validation


def main() -> None:
    model = CoupledCardiorenalModel(
        hallow_table3_calibrated_cardiorenal_parameters()
    )
    equilibrium = hallow_table3_reference_equilibrium_state()
    checks = assert_validation(model, equilibrium_state=equilibrium)
    print("All coupled checks passed:\n")
    for check in checks:
        print(f"  PASS  {check.name}: {check.value:.10g}")

    output = model.algebraic_outputs(0.0, equilibrium)
    print("\nCalibrated coupled equilibrium outputs:")
    for name, unit in (
        ("p_ma_mmHg", "mmHg"),
        ("p_ra_mmHg", "mmHg"),
        ("cardiac_output_L_per_min", "L/min"),
        ("gfr_mL_per_min", "mL/min"),
        ("renal_blood_flow_L_per_min", "L/min"),
        ("phi_urine_sodium_mEq_per_min", "mEq/min"),
    ):
        print(f"  {name:38s} {output[name]: .10g} {unit}")

    print("\nEquilibrium state derivatives:")
    for name, value in derivative_report(model, equilibrium).items():
        print(f"  {name:38s} {value: .6e}")


if __name__ == "__main__":
    main()
