"""Standard plots used by the example scripts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from .simulation import SimulationResult


def plot_standard_summary(
    result: SimulationResult,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Create a compact six-panel physiological summary figure."""

    if not result.success:
        raise RuntimeError(result.message)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    days = result.time_days

    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    panels = (
        ("gfr_mL_per_min", "GFR (mL/min)"),
        ("renal_blood_flow_L_per_min", "Renal blood flow (L/min)"),
        ("phi_urine_sodium_mEq_per_min", "Urinary Na (mEq/min)"),
        ("c_sodium_mEq_per_L", "Plasma Na (mEq/L)"),
        ("at1_bound_pg_per_mL", "AT1-bound AngII (pg/mL)"),
        ("p_b_total_mmHg", "Bowman pressure (mmHg)"),
    )
    for axis, (name, label) in zip(axes.flat, panels):
        axis.plot(days, result.outputs[name], linewidth=2)
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("Time (days)")
    axes[-1, 1].set_xlabel("Time (days)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_two_scenarios(
    control: SimulationResult,
    changed: SimulationResult,
    path: str | Path,
    *,
    title: str,
    changed_label: str,
) -> Path:
    """Compare a perturbation against a simultaneously simulated control."""

    if not control.success or not changed.success:
        raise RuntimeError("Both simulations must finish successfully")
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        ("gfr_mL_per_min", "GFR (mL/min)"),
        ("phi_urine_sodium_mEq_per_min", "Urinary Na (mEq/min)"),
        ("c_sodium_mEq_per_L", "Plasma Na (mEq/L)"),
        ("at1_bound_pg_per_mL", "AT1-bound AngII (pg/mL)"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for axis, (name, label) in zip(axes.flat, panels):
        axis.plot(control.time_days, control.outputs[name], label="Control")
        axis.plot(changed.time_days, changed.outputs[name], label=changed_label)
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


def plot_rihp_summary(
    result: SimulationResult,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Plot the RIHP driver, segmental effects, and main downstream outputs."""

    if not result.success:
        raise RuntimeError(result.message)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    days = result.time_days
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

    axes[0, 0].plot(days, result.outputs["delta_rihp_mmHg"], linewidth=2)
    axes[0, 0].set_ylabel("Delta RIHP (mmHg)")

    for name, label in (
        ("gamma_rihp_pt", "PT"),
        ("gamma_rihp_dt", "DT"),
        ("gamma_rihp_cd", "CD"),
    ):
        axes[0, 1].plot(days, result.outputs[name], label=label)
    axes[0, 1].set_ylabel("RIHP multiplier")
    axes[0, 1].legend()

    axes[1, 0].plot(
        days,
        result.outputs["phi_urine_sodium_mEq_per_min"],
        linewidth=2,
    )
    axes[1, 0].set_ylabel("Urinary Na (mEq/min)")
    axes[1, 0].set_xlabel("Time (days)")

    axes[1, 1].plot(days, result.outputs["gfr_mL_per_min"], linewidth=2)
    axes[1, 1].set_ylabel("GFR (mL/min)")
    axes[1, 1].set_xlabel("Time (days)")

    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path
