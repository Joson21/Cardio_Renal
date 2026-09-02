"""External experiment inputs remaining after cardiovascular coupling."""

from __future__ import annotations

from dataclasses import dataclass

from hallow_kidney.inputs import InputFunction, constant_input


@dataclass(frozen=True)
class CardiorenalInputs:
    """Signals prescribed from outside the closed long-term model.

    MAP and right-atrial pressure are deliberately absent: the cardiovascular
    closure now calculates them from volume, autonomic activity, and
    vascularity.
    """

    sodium_intake_mEq_per_min: InputFunction
    p_peritubular_mmHg: InputFunction

    @classmethod
    def baseline(cls) -> "CardiorenalInputs":
        return cls(
            sodium_intake_mEq_per_min=constant_input(0.126),
            p_peritubular_mmHg=constant_input(0.0),
        )
