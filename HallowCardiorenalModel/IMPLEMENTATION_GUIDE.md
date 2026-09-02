# Implementation and Model-Regeneration Guide

## 1. Purpose

Use this document to rebuild the coupled model from source papers and the
latest project equation PDF. It records decisions that are not uniquely
recoverable from the publications, defines the software boundary, and gives
acceptance tests for a regenerated implementation.

The target is the Karaaslan/Hallow **long-term cardiorenal system** with the
Hallow systemic RAAS extension, the agreed equation corrections, and the fixed
CVP-to-Bowman-pressure novelty. It is not a beat-to-beat cardiac model.

## 2. Required source bundle

Provide all of the following when regenerating the code:

1. Hallow et al. 2014 main paper and appendix/supplement.
2. Karaaslan et al. 2005 paper.
3. Lo et al. systemic RAAS paper used by Hallow.
4. The newest combined kidney-equations PDF containing highlighted project
   corrections and extensions.
5. The experimental data or digitization record behind the fixed `0.1468`
   `P_B,CVP` versus `P_ra` slope.
6. This guide.
7. `CARDIOVASCULAR_EQUATIONS_AND_COUPLING.md`.
8. `PARAMETERS_AND_INITIAL_CONDITIONS.md`.
9. `EQUATION_CHANGE_LOG.md`.
10. `VALIDATION_REPORT.md`.

Do not regenerate from the equation PDF alone. Prose-only bounds, units,
initialization decisions, and publication conflicts are necessary for an
executable reconstruction.

## 3. Source hierarchy

When sources disagree, use this order:

1. Explicit project decisions in `EQUATION_CHANGE_LOG.md`.
2. The newest project equation PDF for highlighted corrections and novelty.
3. Hallow for modified renal hemodynamics, systemic RAAS coupling, AT1 effects,
   and Hallow-refitted tubular/aldosterone relationships.
4. Karaaslan for inherited sodium, water, autonomic, ADH, ANP, and
   cardiovascular equations.
5. Lo for systemic RAAS topology and kinetic interpretation.
6. Physiology only when the above do not define an executable choice.

Never silently select one conflicting value. Record the alternatives, chosen
value, units, and reason in the equation-change log.

## 4. Non-negotiable project decisions

The following choices define this implementation and must survive code
regeneration unless the investigator explicitly revises them.

| Topic | Required choice |
|---|---|
| Afferent nominal resistance | `6.0e7 mmHg min/L` per single-nephron resistance term |
| Efferent nominal resistance | `1.0e8 mmHg min/L` per single-nephron resistance term |
| Myogenic multiplier | include the leading `1 +` |
| Renin secretion | `N_rs=97 pg/mL/h = 97/60 pg/mL/min` |
| Reference PRC | `97/3.466 = 27.988283793246 pg/mL` |
| RAAS initialization | solve the exact Hallow rate equations at unit renin feedback; do not copy Lo table values blindly |
| Default renin feedback | normalize each published dynamic curve to one at the agreed nominal reference |
| Baroreceptor adaptation | `a_baro=0.75*a_auto-k_b*I_b`, `dI_b/dt=a_baro-0.75` |
| ADH autonomic clamp | apply `max(1,epsilon_aum)` only inside the ADH secretion block |
| CVP-Bowman novelty | `P_B,CVP=0.1468(P_ra-P_ra,ref)` with `P_ra,ref=0` |
| Renal blood flow | retain `P_ma/R_renal` in the current novelty configuration |
| Cardiovascular default | use the documented reference-centered closure |
| Literal source option | retain a separately named `published_literal` mode |
| RIHP | keep disabled by default; `S_P-N,i=3` is a provisional example |

## 5. Model boundary

### 5.1 External inputs after coupling

The closed long-term model has two experiment inputs:

| Input | Baseline | Meaning |
|---|---:|---|
| sodium intake | `0.126 mEq/min` | external dietary sodium protocol |
| peritubular-pressure signal | `0 mmHg` relative to reference | optional provisional RIHP driver |

`P_ma` and `P_ra` are no longer external inputs in the coupled model. They are
algebraic solutions of the cardiovascular closure. The standalone
`hallow_kidney` module still accepts pressure functions so it can be tested in
isolation.

### 5.2 State vector

The coupled state vector contains 15 values in this exact order:

1. total extracellular sodium;
2. extracellular fluid volume;
3. PRC;
4. AGT;
5. AngI;
6. AngII;
7. Ang(1-7);
8. AngIV;
9. AT1-bound AngII;
10. AT2-bound AngII;
11. normalized aldosterone;
12. normalized ADH;
13. baroreceptor adaptation integral;
14. right-atrial ADH adaptation integral;
15. vascularity.

All pressures, flows, resistances, hormone stimuli, and multipliers are
algebraic outputs. Do not assign ODE initial conditions to algebraic quantities
such as MAP, `P_ra`, cardiac output, GFR, or `P_B,total`.

## 6. Unit convention

Use minutes as the internal time unit.

### 6.1 Rate conversion

\[
k_{/min}=k_{/h}/60.
\]

For a half-life `h` in minutes:

\[
k_{clear}=\ln(2)/h.
\]

Do not apply both a supplied first-order rate and a half-life clearance for the
same process unless the equation explicitly contains both.

### 6.2 Filtration

`K_f` is in `nL/min/mmHg/nephron`:

\[
GFR_{L/min}=K_f(P_{gh}-P_{B,total}-P_{go})N_{nephron}10^{-9}.
\]

Filtered sodium is

\[
\Phi_{filsod}=GFR_{L/min}C_{sodium,mEq/L}.
\]

### 6.3 AT1-bound AngII conversion

The RAAS states use `fmol/mL`; downstream empirical equations use a mass
concentration. With AngII molecular weight `1046.18 g/mol`:

\[
[AT1]_{pg/mL}=[AT1]_{fmol/mL}\times1046.18\times10^{-3}.
\]

### 6.4 Exponential arguments

Every exponential argument must be dimensionless. The RIHP example therefore
uses `delta_RIHP/(1 mmHg)`, even though its numerical value is unchanged when
pressures are stored in mmHg.

## 7. Regeneration sequence

### Step 1: create immutable parameter bundles

Create a flat kidney parameter dataclass and a cardiovascular dataclass. Build
named configuration constructors instead of using global variables. Every
simulation should retain a configuration name in its output metadata or code
path.

Required constructors:

- corrected kidney/cardiorenal, CVP off;
- recommended CVP extension;
- provisional CVP+RIHP example; and
- published-literal cardiovascular comparison.

### Step 2: implement pure helper equations

Implement and unit-test small equations independently:

- concentration conversion;
- corrected myogenic multiplier;
- TGF multiplier;
- three renin-feedback functions;
- AT1 vascular effects;
- RIHP multiplier;
- constant, step, and ramp input protocols.

### Step 3: implement the kidney algebra

Given current states, MAP, `P_ra`, sodium intake, and optional peritubular
pressure:

1. calculate sodium, aldosterone, ADH, and AT1 mass concentrations;
2. calculate autonomic and RSNA outputs;
3. calculate `P_B,CVP` and `P_B,total`;
4. calculate optional RIHP multipliers;
5. solve the renal hemodynamic/TGF/myogenic algebraic loop;
6. calculate distal and collecting-duct sodium handling;
7. calculate water handling;
8. calculate aldosterone, ADH, and renin stimuli;
9. return a dictionary containing every diagnostic output.

The renal loop must solve

\[
f(P_{gh})=P_{gh}-P_{gh,implied}=0
\]

at every ODE evaluation. A one-pass calculation using a stale glomerular
pressure is not equivalent.

### Step 4: implement kidney derivatives

Use algebraic outputs to construct one derivative for each of the first 14
states. In the AngII equation, ACE2 conversion, AngIV conversion, AT1 binding,
AT2 binding, and half-life clearance are all losses from free AngII. Their
product equations receive the corresponding positive formation terms.

Implement adaptation integrals as states:

\[
\dot I_b=a_{baro}-0.75,
\qquad
a_{baro}=0.75a_{auto}-k_bI_b,
\]

\[
\dot I_{ra}=\delta_{ra},
\qquad
\delta_{ra}=0.2(P_{ra}-P_{ra,ref})-k_{ra}I_{ra}.
\]

### Step 5: implement the cardiovascular algebra

Implement the equations in
`CARDIOVASCULAR_EQUATIONS_AND_COUPLING.md`. For a candidate MAP, calculate
autonomic tone, resistance, cardiac output, `P_ra`, mean filling pressure, and
venous return. Solve

\[
CO(MAP)-VR(MAP)=0.
\]

The reference-centered mode must give `MAP=100`, `P_ra=0`, and `CO=5` at
`V_ECF=15`, `vascularity=1`, and zero adaptation integrals.

### Step 6: close the two modules

During every full-model right-hand-side evaluation:

1. solve cardiovascular outputs from `V_ECF`, vascularity, and baroreceptor
   adaptation;
2. pass calculated MAP and `P_ra` to the kidney;
3. append the vascularity derivative to the 14 kidney derivatives.

Do not duplicate a pressure state. MAP and `P_ra` remain algebraic.

### Step 7: use a stiff ODE solver

Use BDF by default and retain Radau as an independent comparison. Inside each
solver callback, solve both scalar algebraic loops to their configured
tolerances. Export both ODE states and algebraic outputs on a regular reporting
grid.

### Step 8: calculate equilibria explicitly

Keep two concepts separate:

- nominal paper/agreed state: readable source reference, not a whole-system
  equilibrium;
- full coupled equilibrium: all 15 state derivatives vanish under one exact
  configuration and input set.

To calculate a new equilibrium:

1. integrate the stiff model for a long settling period;
2. transform positive states to logarithms;
3. refine all scaled derivatives with constrained nonlinear least squares;
4. verify sodium and water balance independently;
5. store the state with its exact configuration name;
6. invalidate and recalculate it after any structural or parameter change.

The delivered helper is `calculate_coupled_equilibrium.py`.

## 8. RAAS initialization procedure

Use

\[
N_{rs}=97\ \mathrm{pg/mL/h},
\qquad
PRC^*=97/3.466=27.988283793246\ \mathrm{pg/mL},
\]

with `v_MD=v_RSNA=v_AT1=1`. Solve the steady cascade using the exact Hallow
rates and half-lives. Substitute the resulting concentrations into every RAAS
ODE and verify residuals numerically.

The Lo table values are source references, not automatically an equilibrium of
the Hallow-modified parameter set. Only use them directly if they satisfy the
implemented equations and units.

## 9. Required automated tests

A regenerated version is not accepted until all of these pass:

1. nominal afferent and efferent resistances equal `6e7` and `1e8`;
2. the corrected myogenic multiplier equals one at `P_gh=60 mmHg`;
3. the analytical unit-feedback RAAS equilibrium matches the nominal peptide
   initial values;
4. `P_B,CVP=0` at `P_ra=P_ra,ref` and equals `0.734 mmHg` at `P_ra=5 mmHg`;
5. the corrected baroreceptor contributions are `0.75` and `0.25` at reference;
6. the ADH clamp is local, not global;
7. the cardiovascular residual `CO-VR` is below `1e-8 L/min`;
8. reference-centered nominal MAP, `P_ra`, and CO equal `100`, `0`, and `5`;
9. vascularity formation equals destruction at the centered nominal reference;
10. raising `V_ECF` at fixed states raises blood volume, MAP, `P_ra`, and CO;
11. the published-literal mode exposes, rather than hides, the source baseline
    mismatch;
12. the supplied full equilibrium has every derivative below `1e-8` in its
    native units;
13. sodium intake equals urinary sodium and water intake equals urine flow at
    that equilibrium;
14. a short BDF simulation remains finite;
15. BDF and Radau agree for a 30-day doubled-sodium protocol within documented
    numerical tolerances.

Run:

```bash
python run_coupled_checks.py
python -m unittest discover -s tests -v
```

## 10. Reproduction versus calibration

Passing mathematical checks does not mean physiological calibration is
complete. Use three levels of verification:

1. **Equation verification:** signs, units, source structure, algebraic
   residuals, and state ordering.
2. **Reference verification:** published nominal values and equilibrium
   consistency.
3. **Predictive calibration/validation:** fit only identified parameters to
   independent baseline and perturbation data, then test on held-out data.

The current full coupled equilibrium is mathematically valid but does not match
all of Hallow's reported healthy baseline outputs. Do not force GFR, renal blood
flow, MAP, and cardiac output simultaneously by arbitrary independent edits.
Define calibration targets, plausible parameter bounds, an objective function,
and identifiability checks first.

## 11. Rules for future structural extensions

For each new mechanism:

1. state the physiological hypothesis in one sentence;
2. write the dimensional equation and define every reference value;
3. decide whether it is an input, parameter, algebraic output, or ODE state;
4. ensure it is neutral at the reference unless a baseline shift is intended;
5. put it behind a named configuration switch;
6. write a direct unit test;
7. compare enabled versus disabled versions from their own equilibria;
8. perform sensitivity analysis;
9. recalibrate only parameters affected by the new mechanism unless residuals
   demonstrate a broader need;
10. update all five guides and the validation report.

Adding an algebraic RIHP relation does not automatically require a new state.
Adding delayed RIHP dynamics does. Adding renal venous pressure to renal blood
flow is a separate pathway from the Bowman-pressure term. Adding a
four-chamber heart requires a new subsystem, states, parameters, and data.

## 12. Packaging acceptance

Before delivering a regenerated project:

1. install it in a clean virtual environment;
2. run all tests and coupled scripts;
3. open every generated PNG and inspect labels, units, and axis scaling;
4. verify CSV headers are unique and contain time in minutes and days;
5. confirm the archive contains source code, tests, guides, and representative
   results but no virtual environment or cache folders;
6. record package version and SHA-256 checksum.
