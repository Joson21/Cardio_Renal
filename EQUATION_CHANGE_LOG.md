# Equation and Assumption Change Log

## 1. Why this log exists

The code combines printed source equations, corrections required for an
executable reconstruction, a baseline reconciliation, a fixed novel mechanism,
and a provisional worked extension. These categories must not be conflated.

Status labels:

- **Source:** implemented from Karaaslan/Hallow/Lo.
- **Correction:** chosen resolution of a publication/translation problem.
- **Reconciliation:** source curve retained but shifted/normalized to preserve
  a stated nominal reference.
- **Novelty:** investigator-defined mechanism intended for the research model.
- **Provisional:** illustrative extension requiring calibration.

## 2. Kidney and RAAS decisions retained from version 0.1

| Item | Implemented choice | Status |
|---|---|---|
| Nominal afferent resistance | `6.0e7`, from `30 x 2e6` | Correction |
| Nominal efferent resistance | `1.0e8`, from `50 x 2e6` | Correction |
| Myogenic multiplier | add missing leading `1 +` so reference multiplier is one | Correction |
| Renin concentration unit/value | PRC reference `27.9883 pg/mL`, not `28 ng/mL` | Correction |
| Renin secretion | `97 pg/mL/h` | Project decision |
| RAAS peptide initial values | solve exact unit-feedback Hallow equilibrium | Correction/project decision |
| AngII loss terms | ACE2, AngIV, AT1, AT2, and half-life clearance all negative | Correction |
| Renin feedback default | normalize printed dynamic functions to one at agreed nominal state | Reconciliation |
| Baroreceptor adaptation | weighted state adapts to `0.75`; integrand is `a_baro-0.75` | Correction |
| ADH autonomic lower bound | apply `max(1,epsilon_aum)` only within ADH block | Correction |
| Right-atrial ADH adaptation | implement its integral as an ODE state | Source implementation |

### 2.1 Baroreceptor equation implemented

\[
a_{baro}=0.75a_{auto}-k_b I_b,
\qquad
\dot I_b=a_{baro}-0.75.
\]

This is the expanded form of a `3/4` factor applied to the entire braced
expression. If the paper is written as

\[
a_{baro}=\frac34\left[a_{auto}-k\int(a_{baro}-1)dt\right],
\]

the outside factor multiplies the integral term as well, giving an effective
integral coefficient `3k/4`. The implemented reference-centered integral is the
equivalent weighted-state formulation selected during reconstruction.

## 3. Cardiovascular connection added in version 0.2

### 3.1 Source equations added

The following Karaaslan/Hallow long-term equations are now implemented:

- ECF volume to blood volume;
- blood volume/autonomic tone to mean filling pressure;
- venous return;
- cardiac output equals venous return;
- exponential cardiac-output/right-atrial-pressure curve;
- arterial, venous-return, and total peripheral resistance;
- MAP equals cardiac output times total peripheral resistance;
- vascularity formation/destruction dynamics.

`vascularity` is appended as state 15. MAP, `P_ra`, cardiac output, blood volume,
and all resistances are algebraic.

### 3.2 Right-atrial-pressure inconsistency

Printed:

\[
P_{ra}^{raw}=0.2787\exp(0.2281CO).
\]

At `CO=5`, the printed equation yields `0.871866 mmHg`, while the source nominal
table states `P_ra=0`. Two modes are therefore provided.

- `published_literal`: no change to the printed equation.
- `baseline_centered`:

  \[
  P_{ra}=0.2787[\exp(0.2281CO)-\exp(0.2281\cdot5)].
  \]

The latter is a **reconciliation**, not a claim that the centered equation was
printed in the paper.

### 3.3 Mean-filling-pressure shift

The blood-volume and mean-filling-pressure equations give `7.091371 mmHg` at
`V_ECF=15`, whereas `CO=5`, `P_ra=0`, and the printed venous-return resistance
require `7.064516 mmHg`. The baseline-centered mode adds
`-0.0268552364 mmHg` after autonomic scaling. Status: **reconciliation**.

### 3.4 Vascularity formation normalization

The printed formation equation gives `1.0267147e-5 /min` at `CO=5`; destruction
is `1e-5 /min` at vascularity one. The baseline-centered mode retains the
published flow slope and normalizes formation to equal destruction at the
stated nominal reference. Status: **reconciliation**.

### 3.5 Cardiovascular algebraic solution

The implementation solves `CO(MAP)-VR(MAP)=0` at every ODE evaluation. This is
an implementation requirement, not a physiological equation change.

## 4. CVP-to-Bowman-pressure novelty

Implemented only when `enable_cvp_bowman=True`:

\[
P_{B,CVP}=0.1468(P_{ra}-P_{ra,ref}),
\qquad P_{ra,ref}=0\ \mathrm{mmHg},
\]

\[
P_{B,total}=18+P_{B,CVP}\quad\mathrm{mmHg}.
\]

Status: **Novelty**. The slope is fixed by the investigator's experimental
digitization and is not changed during the current implementation.

Reference behavior:

- at `P_ra=0`, `P_B,CVP=0`;
- at `P_ra=5 mmHg`, `P_B,CVP=0.734 mmHg`.

## 5. Renal perfusion-pressure decision

The current model retains

\[
\Phi_{rb}=P_{ma}/R_r.
\]

It does **not** subtract `P_ra` or renal venous pressure. Status: **faithful
reproduction boundary for the current novelty**.

A future configuration may add

\[
\Phi_{rb}=(P_{ma}-P_{rv})/R_r,
\qquad P_{rv}\approx P_{ra},
\]

but must be separately named, tested, and calibrated. This would represent a
vascular perfusion mechanism distinct from elevated Bowman/interstitial
pressure.

## 6. RIHP worked extension

When `enable_rihp=True`:

\[
\mathrm{RIHP}-\mathrm{RIHP}_{ref}
=P_{peritubular}-P_{peritubular,ref},
\]

\[
\gamma_{rihp,i}=1+S_{P-N,i}
\left[
\frac{1}{1+\exp((\mathrm{RIHP}-\mathrm{RIHP}_{ref})/(1\ \mathrm{mmHg}))}
-0.5
\right].
\]

Current example values are `S_P-N,PT=S_P-N,DT=S_P-N,CD=3`. The multipliers are
applied to proximal, distal, and collecting-duct sodium reabsorption.

Status: **Provisional**. The reference is neutral (`gamma=1`), but segment
sensitivities and the pressure scale require data-based calibration.

## 7. Explicit non-features

The following are not present and must not be inferred from the project name:

- no four heart chambers;
- no valves or cardiac-cycle timing;
- no ventricular elastance or pressure-volume loops;
- no explicit contractility or ejection fraction;
- no pulmonary circulation;
- no renal venous pressure subtraction;
- no dynamic RIHP state;
- no diabetes, nephron-loss, drug-PK, or explicit heart-failure module unless
  added in a future named configuration.

## 8. Initial-condition changes in version 0.2

The nominal kidney/RAAS state from version 0.1 remains available and now has
`vascularity=1` appended. Two full coupled equilibria are additionally stored:

- corrected CVP-off equilibrium;
- recommended CVP-enabled equilibrium.

These equilibria were found by stiff settling followed by constrained
nonlinear residual refinement. They are computed consequences of the current
model, not source-paper parameter values.

## 9. Calibration status

No parameters were silently adjusted to force the full equilibrium to Hallow's
reported healthy baseline outputs. The reference-centered cardiovascular
changes only resolve internal nominal conflicts. The resulting full equilibrium
has lower GFR and renal blood flow than the paper's healthy targets and reaches
the implemented urine-flow floor. Status: **mathematically verified,
physiologically not fully recalibrated**.

## 10. Version summary

| Version | Main content |
|---|---|
| 0.1.0 | corrected standalone kidney/RAAS, CVP-Bowman novelty, provisional RIHP example |
| 0.2.0 | bidirectional Karaaslan/Hallow cardiovascular coupling, centered/literal modes, 15-state equilibria, coupled examples and validation |
