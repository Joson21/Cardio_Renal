"""Run a short standalone baseline diagnostic simulation."""

from pathlib import Path

from hallow_kidney import KidneyModel, cvp_extended_parameters, simulate
from hallow_kidney.plotting import plot_standard_summary
from hallow_kidney.validation import assert_validation


def main() -> None:
    model = KidneyModel(cvp_extended_parameters())
    assert_validation(model)

    # The standalone model treats MAP and Pra as external inputs. Until it is
    # coupled to the cardiovascular model, use this run as a short diagnostic,
    # not as a claim of a stable long-term whole-body equilibrium.
    result = simulate(model, days=0.5, samples_per_day=48)
    if not result.success:
        raise RuntimeError(result.message)

    results_directory = Path("results")
    csv_path = result.save_csv(results_directory / "baseline.csv")
    figure_path = plot_standard_summary(
        result,
        results_directory / "baseline.png",
        title="Corrected Hallow + CVP model: baseline diagnostic",
    )
    print(f"Saved {csv_path}")
    print(f"Saved {figure_path}")
    print(f"Initial GFR: {result.outputs['gfr_mL_per_min'][0]:.3f} mL/min")
    print(f"Final GFR:   {result.outputs['gfr_mL_per_min'][-1]:.3f} mL/min")
    print(
        "Initial sodium balance: "
        f"{result.outputs['sodium_intake_mEq_per_min'][0] - result.outputs['phi_urine_sodium_mEq_per_min'][0]:.6f} mEq/min"
    )


if __name__ == "__main__":
    main()
