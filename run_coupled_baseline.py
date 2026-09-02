"""Run the verified coupled reference state without a perturbation."""

from pathlib import Path

from hallow_cardiorenal import (
    CoupledCardiorenalModel,
    hallow_table3_calibrated_cardiorenal_parameters,
    hallow_table3_reference_equilibrium_state,
    plot_coupled_summary,
    simulate_coupled,
)


def main() -> None:
    model = CoupledCardiorenalModel(
        hallow_table3_calibrated_cardiorenal_parameters()
    )
    result = simulate_coupled(
        model,
        days=10.0,
        y0=hallow_table3_reference_equilibrium_state(),
        samples_per_day=12,
    )
    if not result.success:
        raise RuntimeError(result.message)

    results = Path("results/cardiorenal")
    csv_path = result.save_csv(results / "coupled_baseline.csv")
    figure_path = plot_coupled_summary(
        result,
        results / "coupled_baseline.png",
        title="Coupled Hallow/Karaaslan calibrated healthy baseline",
    )
    print(f"Saved {csv_path}")
    print(f"Saved {figure_path}")
    print(
        "MAP, Pra, CO, GFR, total Bowman pressure: "
        f"{result.outputs['p_ma_mmHg'][-1]:.3f} mmHg, "
        f"{result.outputs['p_ra_mmHg'][-1]:.4f} mmHg, "
        f"{result.outputs['cardiac_output_L_per_min'][-1]:.3f} L/min, "
        f"{result.outputs['gfr_mL_per_min'][-1]:.3f} mL/min, "
        f"{result.outputs['p_b_total_mmHg'][-1]:.3f} mmHg"
    )


if __name__ == "__main__":
    main()
