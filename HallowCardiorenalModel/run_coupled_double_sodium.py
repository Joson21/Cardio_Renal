"""Compare baseline sodium intake with a twofold step in the closed model."""

from pathlib import Path

from hallow_cardiorenal import (
    CardiorenalInputs,
    CoupledCardiorenalModel,
    cvp_extended_cardiorenal_parameters,
    cvp_reference_equilibrium_state,
    plot_coupled_comparison,
    simulate_coupled,
)
from hallow_kidney.inputs import constant_input, step_input


def main() -> None:
    start_day = 2.0
    control_inputs = CardiorenalInputs.baseline()
    doubled_inputs = CardiorenalInputs(
        sodium_intake_mEq_per_min=step_input(
            baseline=0.126,
            changed=0.252,
            start_min=start_day * 1440.0,
        ),
        p_peritubular_mmHg=constant_input(0.0),
    )
    parameters = cvp_extended_cardiorenal_parameters()
    equilibrium = cvp_reference_equilibrium_state()
    control = simulate_coupled(
        CoupledCardiorenalModel(parameters, control_inputs),
        days=30.0,
        y0=equilibrium,
        samples_per_day=12,
    )
    doubled = simulate_coupled(
        CoupledCardiorenalModel(parameters, doubled_inputs),
        days=30.0,
        y0=equilibrium,
        samples_per_day=12,
    )
    if not control.success or not doubled.success:
        raise RuntimeError(control.message if not control.success else doubled.message)

    results = Path("results/cardiorenal")
    control_path = control.save_csv(results / "coupled_sodium_control.csv")
    doubled_path = doubled.save_csv(results / "coupled_sodium_doubled.csv")
    figure_path = plot_coupled_comparison(
        control,
        doubled,
        results / "coupled_sodium_doubled_comparison.png",
        title="Coupled response to doubled sodium intake at day 2",
        changed_label="Sodium intake doubled",
    )
    print(f"Saved {control_path}")
    print(f"Saved {doubled_path}")
    print(f"Saved {figure_path}")
    for name, unit in (
        ("p_ma_mmHg", "mmHg"),
        ("p_ra_mmHg", "mmHg"),
        ("cardiac_output_L_per_min", "L/min"),
        ("gfr_mL_per_min", "mL/min"),
    ):
        print(
            f"Final {name}: control={control.outputs[name][-1]:.6g}, "
            f"doubled={doubled.outputs[name][-1]:.6g} {unit}"
        )


if __name__ == "__main__":
    main()
