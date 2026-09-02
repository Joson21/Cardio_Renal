"""Demonstrate the fixed CVP-to-Bowman pressure extension."""

from pathlib import Path

from hallow_kidney import KidneyModel, ModelInputs, cvp_extended_parameters, simulate
from hallow_kidney.inputs import constant_input, step_input
from hallow_kidney.plotting import plot_standard_summary


def main() -> None:
    start_min = 0.1 * 1440.0
    inputs = ModelInputs(
        p_ma_mmHg=constant_input(100.0),
        p_ra_mmHg=step_input(0.0, 5.0, start_min=start_min),
        sodium_intake_mEq_per_min=constant_input(0.126),
        p_peritubular_mmHg=constant_input(0.0),
    )
    model = KidneyModel(cvp_extended_parameters(), inputs)
    result = simulate(model, days=0.5, samples_per_day=96)
    if not result.success:
        raise RuntimeError(result.message)

    directory = Path("results")
    result.save_csv(directory / "cvp_step.csv")
    plot_standard_summary(
        result,
        directory / "cvp_step.png",
        title="Right atrial pressure step: 0 to 5 mmHg",
    )
    final_pb_cvp = result.outputs["p_b_cvp_mmHg"][-1]
    print(f"Expected P_B,CVP = 0.1468 x 5 = 0.734 mmHg")
    print(f"Simulated P_B,CVP = {final_pb_cvp:.3f} mmHg")


if __name__ == "__main__":
    main()
