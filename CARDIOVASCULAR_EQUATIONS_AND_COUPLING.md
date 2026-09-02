# Cardiovascular Equations and Cardiorenal Coupling

## 1. Scope

The connected cardiovascular subsystem is the long-term closure inherited from
Karaaslan and printed in the Hallow appendix as Eqs. A18-A26. It converts
extracellular volume, autonomic tone, and vascularity into blood volume, mean
filling pressure, venous return, cardiac output, right-atrial pressure, total
peripheral resistance, and mean arterial pressure.

Only vascularity is a new differential state. Cardiac output, right-atrial
pressure, and MAP are instantaneous algebraic variables. Consequently, this
block should be described as a **lumped long-term cardiovascular closure**, not
as a four-chamber cardiac mechanics model.

## 2. Published cardiovascular equations

The code uses litres, minutes, and mmHg internally.

### 2.1 Blood volume

\[
V_b=4.560227+
\frac{2.431217}{1+\exp[-0.47437(V_{ECF}-18.11278)]}
\quad\mathrm{L}.
\]

At the stated reference `V_ECF=15 L`, this gives

\[
V_b^*=5.012287704\ \mathrm{L}.
\]

### 2.2 Mean filling pressure

The printed form is

\[
P_{mf}^{raw}=(7.436V_b-30.18)\epsilon_{aum}.
\]

At `V_ECF=15 L` and `epsilon_aum=1`, it gives
`7.091371365 mmHg`.

### 2.3 Venous return and cardiac output

\[
\Phi_{vr}=\frac{P_{mf}-P_{ra}}{R_{vr}},
\qquad
\Phi_{co}=\Phi_{vr}.
\]

### 2.4 Right-atrial pressure relation

The printed equation is

\[
P_{ra}^{raw}=0.2787\exp(0.2281\Phi_{co})\quad\mathrm{mmHg}.
\]

At the paper's nominal cardiac output of `5 L/min`, this evaluates to

\[
P_{ra}^{raw,*}=0.8718661675\ \mathrm{mmHg},
\]

not the `0 mmHg` nominal right-atrial pressure stated in the Karaaslan baseline
table. This is a source-level inconsistency; it cannot be removed by more
accurate numerical solution.

### 2.5 Vascularity

The printed dynamics are

\[
\frac{d\mathrm{vas}}{dt}=\mathrm{vas}_f-\mathrm{vas}_d,
\]

\[
\mathrm{vas}_f=11.312\times10^{-5}\exp(-0.4799\Phi_{co}),
\qquad
\mathrm{vas}_d=10^{-5}\mathrm{vas}.
\]

At `CO=5 L/min` and `vas=1`, the printed formation rate is
`1.026714717e-5 /min`, while destruction is `1.0e-5 /min`; hence the printed
nominal vascularity is also not an exact equilibrium.

### 2.6 Resistances and pressure

\[
R_a=\frac{16.6}{\mathrm{vas}}\epsilon_{aum},
\]

\[
R_{vr}=\frac{8R_{bv}+R_a}{31},
\qquad R_{bv}=3.4\ \mathrm{mmHg\,min/L},
\]

\[
R_{tp}=R_a+R_{bv},
\qquad
P_{ma}=\Phi_{co}R_{tp}.
\]

At `vas=epsilon_aum=1`, `R_tp=20 mmHg min/L`; therefore `CO=5 L/min`
produces exactly `MAP=100 mmHg`.

## 3. Recommended reference-centered reconstruction

The default mode, `baseline_centered`, preserves the shapes and slopes of the
published functions while making the stated nominal triplet executable:

\[
P_{ma}^*=100\ \mathrm{mmHg},\qquad
P_{ra}^*=0\ \mathrm{mmHg},\qquad
\Phi_{co}^*=5\ \mathrm{L/min}.
\]

### 3.1 Centered right-atrial pressure

\[
P_{ra}=0.2787\left[
\exp(0.2281\Phi_{co})-\exp(0.2281\Phi_{co}^*)
\right]+P_{ra}^*.
\]

This is a vertical translation. Its derivative with respect to cardiac output
is unchanged. Values slightly below zero are allowed when cardiac output is
below `5 L/min`; small negative right-atrial pressures are mathematically
consistent with this centered empirical curve and are not automatically
clamped.

### 3.2 Mean-filling-pressure reconciliation

At the nominal reference,

\[
R_{vr}^*=\frac{8(3.4)+16.6}{31}=1.412903226
\ \mathrm{mmHg\,min/L}.
\]

The filling pressure required for `CO=5` and `P_ra=0` is

\[
P_{mf}^{required}=5R_{vr}^*=7.064516129\ \mathrm{mmHg}.
\]

The code therefore uses the constant shift

\[
\Delta P_{mf}=7.064516129-7.091371365
=-0.0268552364\ \mathrm{mmHg},
\]

and evaluates

\[
P_{mf}=(7.436V_b-30.18)\epsilon_{aum}+\Delta P_{mf}.
\]

The change is less than `0.03 mmHg` and exists only to close the stated
reference exactly.

### 3.3 Reference-normalized vascularity formation

The recommended mode uses

\[
\mathrm{vas}_f=k_d\mathrm{vas}^*
\exp[-0.4799(\Phi_{co}-\Phi_{co}^*)],
\]

where `k_d=1e-5 /min`, `vas*=1`, and `CO*=5 L/min`. Thus formation and
destruction are exactly equal at the nominal reference without changing the
published flow sensitivity `0.4799 min/L`.

## 4. Published-literal comparison mode

`published_literal_cardiorenal_parameters()` evaluates the equations as
printed, without any of the three reference reconciliations above. At the
nominal state, the closed algebraic solution is approximately:

| Quantity | Published literal | Reference-centered |
|---|---:|---:|
| MAP | 95.8013 mmHg | 100.0000 mmHg |
| `P_ra` | 0.797515 mmHg | 0.000000 mmHg |
| Cardiac output | 4.60923 L/min | 5.00000 L/min |
| Mean filling pressure | 7.42658 mmHg | 7.06452 mmHg |
| Total peripheral resistance | 20.78467 mmHg min/L | 20.00000 mmHg min/L |
| `d(vascularity)/dt` | 2.38494e-6 /min | approximately 0 /min |

Use literal mode for equation auditing, not as the default perturbation model.

## 5. Cardiovascular algebraic loop

For a given `V_ECF`, vascularity, and baroreceptor adaptation state, MAP is not
known in advance. The code solves one scalar closure residual.

For each candidate MAP:

1. calculate `a_auto`, `a_baro`, `a_chemo`, and `epsilon_aum`;
2. calculate arterial and total peripheral resistance;
3. calculate `CO = MAP/R_tp`;
4. calculate `P_ra(CO)`;
5. calculate blood volume and mean filling pressure;
6. calculate venous-return resistance and venous return; and
7. solve

\[
f(P_{ma})=\frac{P_{ma}}{R_{tp}}-
\frac{P_{mf}-P_{ra}}{R_{vr}}=0.
\]

The accepted root therefore satisfies cardiac output equals venous return.
Brent's method is used when a sign-changing interval exists; a bounded residual
minimization is retained only as a diagnostic fallback.

## 6. Bidirectional coupling to the kidney

The coupling layer performs this sequence during every ODE evaluation:

1. read kidney `V_ECF` and cardiovascular `vascularity`;
2. solve the cardiovascular algebraic loop;
3. pass calculated `P_ma` and `P_ra` to the kidney algebraic equations;
4. calculate kidney/RAAS/autonomic/ADH outputs;
5. calculate all 14 kidney-state derivatives;
6. append the vascularity derivative as state 15.

Important feedback paths include:

- `V_ECF -> V_b -> P_mf -> CO/P_ra/MAP`;
- `MAP -> autonomic multiplier -> systemic resistance`;
- `MAP/P_ra -> RSNA, ADH, ANP, aldosterone and renal hemodynamics`;
- `P_ra -> P_B,CVP -> filtration pressure -> GFR` when novelty is enabled;
- renal sodium/water excretion -> total sodium and `V_ECF`.

## 7. CVP-Bowman mechanism and renal perfusion pressure

The novelty is

\[
P_{B,\mathrm{CVP}}=0.1468(P_{ra}-0),
\qquad
P_{B,total}=18+P_{B,\mathrm{CVP}}.
\]

The current renal blood-flow equation remains

\[
\Phi_{rb}=\frac{P_{ma}}{R_r}.
\]

If renal venous pressure is later approximated by right-atrial pressure, a
separate structural configuration could use

\[
\Phi_{rb}=\frac{P_{ma}-P_{rv}}{R_r},
\qquad P_{rv}\approx P_{ra}.
\]

That is not mathematical double counting: reduced vascular perfusion pressure
and increased interstitial/Bowman pressure are distinct physical pathways.
They can, however, be empirically overlapping if the fitted `0.1468` slope was
derived from data in which both mechanisms were already present. Adding both
therefore requires data-based identification and should be tested as a
separately named hypothesis.

## 8. What this block cannot simulate

The implemented closure has no ventricular elastance, contractility state,
preload-dependent stroke volume, valve dynamics, pulmonary circulation,
left/right chamber volumes, heart rate, or explicit heart-failure parameter.
An experiment involving ejection fraction, beat-to-beat waveforms, or pump
failure requires a new cardiac module and new calibration data; it cannot be
obtained by merely changing one of the current constants.
