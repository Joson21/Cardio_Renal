"""Worked provisional RIHP extension inside the coupled model."""

from dataclasses import replace
from pathlib import Path

from hallow_cardiorenal import (
    CardiorenalInputs,
    CoupledCardiorenalModel,
    hallow_table3_calibrated_cardiorenal_parameters,
    hallow_table3_reference_equilibrium_state,
    plot_coupled_rihp_summary,
    simulate_coupled,
)
from hallow_kidney.inputs import constant_input, step_input


def main() -> None:
    inputs = CardiorenalInputs(
        sodium_intake_mEq_per_min=constant_input(0.126),
        p_peritubular_mmHg=step_input(
            baseline=0.0,
            changed=0.1,
            start_min=1.0 * 1440.0,
        ),
    )
    base_parameters = hallow_table3_calibrated_cardiorenal_parameters()
    parameters = replace(
        base_parameters,
        configuration_name="hallow_table3_calibrated_rihp_demo",
        kidney=replace(
            base_parameters.kidney,
            configuration_name="hallow_table3_calibrated_rihp_demo",
            enable_rihp=True,
        ),
    )
    model = CoupledCardiorenalModel(parameters, inputs)
    result = simulate_coupled(
        model,
        days=4.0,
        y0=hallow_table3_reference_equilibrium_state(),
        samples_per_day=24,
    )
    if not result.success:
        raise RuntimeError(result.message)
    results = Path("results/cardiorenal")
    csv_path = result.save_csv(results / "coupled_rihp_demo.csv")
    figure_path = plot_coupled_rihp_summary(
        result,
        results / "coupled_rihp_demo.png",
        title="Provisional RIHP transport extension (+0.1 mmHg at day 1)",
    )
    print(f"Saved {csv_path}")
    print(f"Saved {figure_path}")
    print(
        "Final PT/DT/CD RIHP multipliers: "
        f"{result.outputs['gamma_rihp_pt'][-1]:.6f}, "
        f"{result.outputs['gamma_rihp_dt'][-1]:.6f}, "
        f"{result.outputs['gamma_rihp_cd'][-1]:.6f}"
    )
    print("S_P-N,i=3 is provisional; this example is not a calibrated prediction.")


if __name__ == "__main__":
    main()
