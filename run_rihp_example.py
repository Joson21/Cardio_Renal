"""Worked RIHP/peritubular-pressure extension example."""

from pathlib import Path

from hallow_kidney import KidneyModel, ModelInputs, rihp_demo_parameters, simulate
from hallow_kidney.inputs import constant_input, step_input
from hallow_kidney.plotting import plot_rihp_summary


def main() -> None:
    inputs = ModelInputs(
        p_ma_mmHg=constant_input(100.0),
        p_ra_mmHg=constant_input(0.0),
        sodium_intake_mEq_per_min=constant_input(0.126),
        # A small +0.1 mmHg step is used because S_P-N=3 is still provisional.
        p_peritubular_mmHg=step_input(
            0.0,
            0.1,
            start_min=0.1 * 1440.0,
        ),
    )
    model = KidneyModel(rihp_demo_parameters(), inputs)
    result = simulate(model, days=0.5, samples_per_day=96)
    if not result.success:
        raise RuntimeError(result.message)

    directory = Path("results")
    result.save_csv(directory / "rihp_step.csv")
    plot_rihp_summary(
        result,
        directory / "rihp_step.png",
        title="Provisional RIHP example: peritubular pressure +0.1 mmHg",
    )
    gamma = result.outputs["gamma_rihp_pt"][-1]
    print(f"Final RIHP multiplier (S_P-N = 3): {gamma:.6f}")
    if gamma <= 0.0:
        print("WARNING: the provisional unbounded multiplier is non-positive.")


if __name__ == "__main__":
    main()
