"""Run the fast mathematical consistency checks."""

from hallow_kidney import KidneyModel, cvp_extended_parameters
from hallow_kidney.steady_state import derivative_report
from hallow_kidney.validation import assert_validation


def main() -> None:
    model = KidneyModel(cvp_extended_parameters())
    checks = assert_validation(model)
    print("All required checks passed:\n")
    for check in checks:
        print(f"  PASS  {check.name}: {check.value:.8g}")

    print("\nInitial full-model derivatives (not all are expected to be zero):")
    for name, value in derivative_report(model).items():
        print(f"  {name:32s} {value: .8e}")


if __name__ == "__main__":
    main()
