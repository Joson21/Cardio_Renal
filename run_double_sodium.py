"""Compare reference sodium intake with a twofold step in sodium intake."""

from pathlib import Path

from hallow_kidney import KidneyModel, ModelInputs, cvp_extended_parameters, simulate
from hallow_kidney.inputs import constant_input, step_input
from hallow_kidney.plotting import plot_two_scenarios


def make_inputs(*, double_sodium: bool) -> ModelInputs:
    sodium = (
        step_input(0.126, 0.252, start_min=0.1 * 1440.0)
        if double_sodium
        else constant_input(0.126)
    )
    return ModelInputs(
        p_ma_mmHg=constant_input(100.0),
        p_ra_mmHg=constant_input(0.0),
        sodium_intake_mEq_per_min=sodium,
        p_peritubular_mmHg=constant_input(0.0),
    )


def main() -> None:
    parameters = cvp_extended_parameters()
    control = simulate(
        KidneyModel(parameters, make_inputs(double_sodium=False)),
        days=0.5,
        samples_per_day=96,
    )
    doubled = simulate(
        KidneyModel(parameters, make_inputs(double_sodium=True)),
        days=0.5,
        samples_per_day=96,
    )
    if not control.success or not doubled.success:
        raise RuntimeError(f"Control: {control.message}; doubled: {doubled.message}")

    directory = Path("results")
    control.save_csv(directory / "sodium_control.csv")
    doubled.save_csv(directory / "sodium_doubled.csv")
    plot_two_scenarios(
        control,
        doubled,
        directory / "sodium_doubled_comparison.png",
        title="Paired experiment: sodium intake doubled",
        changed_label="Sodium intake doubled",
    )
    print("Sodium input changes from 0.126 to 0.252 mEq/min at day 0.1.")
    print(
        "Final difference in total sodium: "
        f"{doubled.state('m_sodium_mEq')[-1] - control.state('m_sodium_mEq')[-1]:.3f} mEq"
    )


if __name__ == "__main__":
    main()
