"""Compare the same sodium challenge with the CVP-Bowman path off and on."""

from pathlib import Path

from hallow_cardiorenal import (
    CardiorenalInputs,
    CoupledCardiorenalModel,
    corrected_hallow_cardiorenal_parameters,
    corrected_hallow_reference_equilibrium_state,
    cvp_extended_cardiorenal_parameters,
    cvp_reference_equilibrium_state,
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
    without_cvp = simulate_coupled(
        CoupledCardiorenalModel(
            corrected_hallow_cardiorenal_parameters(), inputs
        ),
        days=30.0,
        y0=corrected_hallow_reference_equilibrium_state(),
        samples_per_day=12,
    )
    with_cvp = simulate_coupled(
        CoupledCardiorenalModel(cvp_extended_cardiorenal_parameters(), inputs),
        days=30.0,
        y0=cvp_reference_equilibrium_state(),
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
