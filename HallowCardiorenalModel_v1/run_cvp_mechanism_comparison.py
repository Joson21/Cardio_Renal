"""Compare the same sodium challenge with the CVP-Bowman path off and on."""

from dataclasses import replace
from pathlib import Path

from hallow_cardiorenal import (
    CardiorenalInputs,
    CoupledCardiorenalModel,
    hallow_table3_calibrated_cardiorenal_parameters,
    hallow_table3_reference_equilibrium_state,
    plot_coupled_comparison,
    simulate_coupled,
)
from hallow_kidney.inputs import constant_input, step_input


def main() -> None:
    inputs = CardiorenalInputs(
        sodium_intake_mEq_per_min=step_input(
            baseline=0.126,
            changed=0.252,
            start_min=2.0 * 1440.0,
        ),
        p_peritubular_mmHg=constant_input(0.0),
    )
    with_parameters = hallow_table3_calibrated_cardiorenal_parameters()
    without_parameters = replace(
        with_parameters,
        configuration_name="hallow_table3_calibrated_cvp_off",
        kidney=replace(
            with_parameters.kidney,
            configuration_name="hallow_table3_calibrated_cvp_off",
            enable_cvp_bowman=False,
        ),
    )
    equilibrium = hallow_table3_reference_equilibrium_state()
    without_cvp = simulate_coupled(
        CoupledCardiorenalModel(without_parameters, inputs),
        days=30.0,
        y0=equilibrium,
        samples_per_day=12,
    )
    with_cvp = simulate_coupled(
        CoupledCardiorenalModel(with_parameters, inputs),
        days=30.0,
        y0=equilibrium,
        samples_per_day=12,
    )
    if not without_cvp.success or not with_cvp.success:
        raise RuntimeError(
            without_cvp.message if not without_cvp.success else with_cvp.message
        )
    results = Path("results/cardiorenal")
    off_path = without_cvp.save_csv(results / "sodium_step_cvp_bowman_off.csv")
    on_path = with_cvp.save_csv(results / "sodium_step_cvp_bowman_on.csv")
    figure_path = plot_coupled_comparison(
        without_cvp,
        with_cvp,
        results / "cvp_bowman_mechanism_comparison.png",
        title="Doubled sodium intake: CVP-Bowman pathway off versus on",
        changed_label="CVP-Bowman enabled",
        control_label="CVP-Bowman disabled",
    )
    print(f"Saved {off_path}")
    print(f"Saved {on_path}")
    print(f"Saved {figure_path}")
    print(
        "Final GFR without/with CVP-Bowman: "
        f"{without_cvp.outputs['gfr_mL_per_min'][-1]:.6f} / "
        f"{with_cvp.outputs['gfr_mL_per_min'][-1]:.6f} mL/min"
    )


if __name__ == "__main__":
    main()
