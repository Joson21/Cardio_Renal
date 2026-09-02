"""Compare BDF and Radau on the 30-day doubled-sodium protocol."""

from pathlib import Path

import numpy as np

from hallow_cardiorenal import (
    CardiorenalInputs,
    CoupledCardiorenalModel,
    cvp_extended_cardiorenal_parameters,
    cvp_reference_equilibrium_state,
    simulate_coupled,
)
from hallow_kidney.inputs import constant_input, step_input


def main() -> None:
    inputs = CardiorenalInputs(
        sodium_intake_mEq_per_min=step_input(
            0.126,
            0.252,
            start_min=2.0 * 1440.0,
        ),
        p_peritubular_mmHg=constant_input(0.0),
    )
    model = CoupledCardiorenalModel(
        cvp_extended_cardiorenal_parameters(), inputs
    )
    common = dict(
        days=30.0,
        y0=cvp_reference_equilibrium_state(),
        samples_per_day=8,
        relative_tolerance=1.0e-8,
        absolute_tolerance=1.0e-10,
    )
    bdf = simulate_coupled(model, method="BDF", **common)
    radau = simulate_coupled(model, method="Radau", **common)
    if not bdf.success or not radau.success:
        raise RuntimeError(bdf.message if not bdf.success else radau.message)

    output_names = (
        "p_ma_mmHg",
        "p_ra_mmHg",
        "cardiac_output_L_per_min",
        "gfr_mL_per_min",
        "renal_blood_flow_L_per_min",
        "phi_urine_sodium_mEq_per_min",
        "p_b_cvp_mmHg",
    )
    rows = []
    print("Maximum BDF-Radau output differences:\n")
    for name in output_names:
        absolute = float(np.max(np.abs(bdf.outputs[name] - radau.outputs[name])))
        scale = max(1.0e-12, float(np.max(np.abs(radau.outputs[name]))))
        relative = absolute / scale
        rows.append((name, absolute, relative))
        print(f"  {name:38s} absolute={absolute:.8g} relative={relative:.8g}")

    results = Path("results/cardiorenal")
    results.mkdir(parents=True, exist_ok=True)
    path = results / "bdf_radau_comparison.csv"
    np.savetxt(
        path,
        np.array(rows, dtype=object),
        delimiter=",",
        fmt="%s",
        header="output,max_absolute_difference,max_relative_difference",
        comments="",
    )
    print(f"\nSaved {path}")


if __name__ == "__main__":
    main()
