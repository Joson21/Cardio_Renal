"""Standard figures for coupled cardiovascular-kidney simulations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from .simulation import CoupledSimulationResult


def _make_flat_trace_readable(axis: plt.Axes, values) -> None:
    """Avoid misleading scientific-offset axes for equilibrium traces."""

    lower = float(min(values))
    upper = float(max(values))
    center = 0.5 * (lower + upper)
    if upper - lower < max(abs(center) * 1.0e-6, 1.0e-10):
        padding = max(abs(center) * 0.02, 0.01)
        axis.set_ylim(center - padding, center + padding)
    axis.ticklabel_format(axis="y", style="plain", useOffset=False)


def plot_coupled_summary(
    result: CoupledSimulationResult,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Plot the principal cardiovascular, renal, and volume outputs."""

    if not result.success:
        raise RuntimeError(result.message)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        ("p_ma_mmHg", "MAP (mmHg)"),
        ("p_ra_mmHg", r"$P_{ra}$ (mmHg)"),
        ("cardiac_output_L_per_min", "Cardiac output (L/min)"),
        ("gfr_mL_per_min", "GFR (mL/min)"),
        ("renal_blood_flow_L_per_min", "Renal blood flow (L/min)"),
        ("phi_urine_sodium_mEq_per_min", "Urinary Na (mEq/min)"),
        ("blood_volume_L", "Blood volume (L)"),
        ("p_b_total_mmHg", r"Total $P_B$ (mmHg)"),
    )
    fig, axes = plt.subplots(4, 2, figsize=(11, 12), sharex=True)
    for axis, (name, label) in zip(axes.flat, panels):
        values = result.outputs[name]
        axis.plot(result.time_days, values, linewidth=2)
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
        _make_flat_trace_readable(axis, values)
    axes[-1, 0].set_xlabel("Time (days)")
    axes[-1, 1].set_xlabel("Time (days)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_coupled_comparison(
    control: CoupledSimulationResult,
    changed: CoupledSimulationResult,
    path: str | Path,
    *,
    title: str,
    changed_label: str,
    control_label: str = "Control",
) -> Path:
    """Compare a coupled perturbation with a matched control run."""

    if not control.success or not changed.success:
        raise RuntimeError("Both simulations must finish successfully")
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        ("p_ma_mmHg", "MAP (mmHg)"),
        ("p_ra_mmHg", r"$P_{ra}$ (mmHg)"),
        ("cardiac_output_L_per_min", "Cardiac output (L/min)"),
        ("gfr_mL_per_min", "GFR (mL/min)"),
        ("phi_urine_sodium_mEq_per_min", "Urinary Na (mEq/min)"),
        ("__state_v_ecf_L", "ECF volume (L)"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for axis, (name, label) in zip(axes.flat, panels):
        control_values = (
            control.state("v_ecf_L")
            if name == "__state_v_ecf_L"
            else control.outputs[name]
        )
        changed_values = (
            changed.state("v_ecf_L")
            if name == "__state_v_ecf_L"
            else changed.outputs[name]
        )
        axis.plot(control.time_days, control_values, label=control_label)
        axis.plot(changed.time_days, changed_values, label=changed_label)
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    axes[0, 0].legend()
    axes[-1, 0].set_xlabel("Time (days)")
    axes[-1, 1].set_xlabel("Time (days)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_coupled_rihp_summary(
    result: CoupledSimulationResult,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Show the provisional RIHP mechanism inside the coupled model."""

    if not result.success:
        raise RuntimeError(result.message)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    days = result.time_days
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)

    axes[0, 0].plot(days, result.outputs["delta_rihp_mmHg"], linewidth=2)
    axes[0, 0].set_ylabel(r"$\Delta$RIHP (mmHg)")

    for name, label in (
        ("gamma_rihp_pt", "PT"),
        ("gamma_rihp_dt", "DT"),
        ("gamma_rihp_cd", "CD"),
    ):
        axes[0, 1].plot(days, result.outputs[name], label=label)
    axes[0, 1].set_ylabel("RIHP multiplier")
    axes[0, 1].legend()

    panels = (
        (axes[1, 0], "p_ma_mmHg", "MAP (mmHg)"),
        (axes[1, 1], "p_ra_mmHg", r"$P_{ra}$ (mmHg)"),
        (
            axes[2, 0],
            "phi_urine_sodium_mEq_per_min",
            "Urinary Na (mEq/min)",
        ),
        (axes[2, 1], "gfr_mL_per_min", "GFR (mL/min)"),
    )
    for axis, name, label in panels:
        axis.plot(days, result.outputs[name], linewidth=2)
        axis.set_ylabel(label)
    axes[-1, 0].set_xlabel("Time (days)")
    axes[-1, 1].set_xlabel("Time (days)")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path
