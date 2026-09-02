"""Evaluate the literal and reference-centered cardiovascular equations."""

from pathlib import Path

import numpy as np

from hallow_cardiorenal import (
    CoupledCardiorenalModel,
    corrected_hallow_cardiorenal_parameters,
    initial_state,
    published_literal_cardiorenal_parameters,
)


def main() -> None:
    models = (
        (
            "baseline_centered",
            CoupledCardiorenalModel(corrected_hallow_cardiorenal_parameters()),
        ),
        (
            "published_literal",
            CoupledCardiorenalModel(published_literal_cardiorenal_parameters()),
        ),
    )
    names = (
        "p_ma_mmHg",
        "p_ra_mmHg",
        "p_ra_raw_mmHg",
        "cardiac_output_L_per_min",
        "mean_filling_pressure_mmHg",
        "total_peripheral_resistance_mmHg_min_per_L",
        "d_vascularity_per_min",
        "gfr_mL_per_min",
    )
    rows: list[list[object]] = []
    for mode, model in models:
        output = model.algebraic_outputs(0.0, initial_state(model.p))
        rows.append([mode, *(output[name] for name in names)])
        print(f"\n{mode}")
        for name in names:
            print(f"  {name:52s} {output[name]: .10g}")

    results = Path("results/cardiorenal")
    results.mkdir(parents=True, exist_ok=True)
    path = results / "published_vs_centered_cardiovascular.csv"
    array = np.array([[str(value) for value in row] for row in rows])
    np.savetxt(
        path,
        array,
        delimiter=",",
        fmt="%s",
        header=",".join(("closure_mode", *names)),
        comments="",
    )
    print(f"\nSaved {path}")


if __name__ == "__main__":
    main()
