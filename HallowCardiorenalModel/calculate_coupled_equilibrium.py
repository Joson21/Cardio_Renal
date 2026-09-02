"""Recalculate a coupled equilibrium after changing equations or parameters."""

from hallow_cardiorenal import (
    CoupledCardiorenalModel,
    cvp_extended_cardiorenal_parameters,
    cvp_reference_equilibrium_state,
)
from hallow_cardiorenal.states import STATE_NAMES
from hallow_cardiorenal.steady_state import calculate_coupled_equilibrium


def main() -> None:
    model = CoupledCardiorenalModel(cvp_extended_cardiorenal_parameters())
    state, report = calculate_coupled_equilibrium(
        model,
        # The supplied reference makes recalculation fast. Replace this with
        # initial_state(model.p) after a major structural change if necessary.
        guess=cvp_reference_equilibrium_state(),
        settle_days=30.0,
    )
    print("Equilibrium solver report:")
    for name, value in report.items():
        print(f"  {name:34s} {value:.12g}")
    print("\nState vector to copy into hallow_cardiorenal/states.py:")
    for name, value in zip(STATE_NAMES, state):
        print(f"  {name:38s} {value:.15g}")
    print("\nKey outputs:")
    output = model.algebraic_outputs(0.0, state)
    for name in (
        "p_ma_mmHg",
        "p_ra_mmHg",
        "cardiac_output_L_per_min",
        "gfr_mL_per_min",
        "renal_blood_flow_L_per_min",
    ):
        print(f"  {name:38s} {output[name]:.12g}")


if __name__ == "__main__":
    main()
