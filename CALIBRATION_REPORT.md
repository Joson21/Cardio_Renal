# Hallow Table 3 Joint Baseline Calibration

## 1. Purpose and scope

This calibration creates a healthy, internally balanced baseline for the
coupled Hallow/Karaaslan reconstruction.  It addresses the previous equilibrium
with `P_ra=-0.0314 mmHg`, `GFR=69.67 mL/min`, renal blood flow `0.499 L/min`,
and water balance at the hard urine-flow floor.

The original equation-reproduction profiles were not overwritten.  They remain
available through `corrected_hallow_cardiorenal_parameters()` and
`cvp_extended_cardiorenal_parameters()`.  The fitted profile is explicitly
named `hallow_table3_calibrated_cardiorenal_parameters()` and its matching
state is `hallow_table3_reference_equilibrium_state()`.

The primary source is Hallow et al., *A model-based approach to investigating
the pathophysiological mechanisms of hypertension and response to
antihypertensive therapies: extending the Guyton model*, 2014,
[doi:10.1152/ajpregu.00039.2013](https://doi.org/10.1152/ajpregu.00039.2013).
The long-term cardiovascular, autonomic, sodium, and water equations originate
from Karaaslan et al., *Long-Term Mathematical Model Involving Renal
Sympathetic Nerve Activity, Arterial Pressure, and Sodium Excretion*, 2005,
[doi:10.1007/s10439-005-5976-4](https://doi.org/10.1007/s10439-005-5976-4).

## 2. Calibration targets

The following observable Table 3 values were used as exact targets:

| Quantity | Table 3 target | Calibrated result |
|---|---:|---:|
| MAP | 83 mmHg | 83.000000 mmHg |
| GFR | 99 mL/min | 99.000000 mL/min |
| urinary sodium | 0.126 mEq/min | 0.126000 mEq/min |
| ECF volume | 15 L | 15.000000 L |
| daily urine flow | 2.1 L/day | 2.100000 L/day |
| sodium concentration | 143.3 mEq/L | 143.300000 mEq/L |
| cardiac output | 5.15 L/min | 5.150000 L/min |
| renal blood flow | 0.9 L/min | 0.900000 L/min |
| aldosterone | 100 pg/mL | 100.000000 pg/mL |
| PRA | 350 fmol/mL/h | 350.000000 fmol/mL/h |

The calibrated baseline additionally has:

| Derived quantity | Result |
|---|---:|
| right-atrial pressure | +0.902213 mmHg |
| CVP-induced Bowman increment | 0 mmHg |
| total Bowman pressure | 18.000000 mmHg |
| glomerular pressure | 58.676056 mmHg |
| net filtration pressure | 12.676056 mmHg |
| single-nephron GFR | 49.5 nL/min |
| renal vascular resistance | 92.222222 mmHg min/L |
| total peripheral resistance | 16.116505 mmHg min/L |
| maximum absolute state derivative | below `7e-16` in native units |

## 3. Rounded source values that cannot all be exact simultaneously

Several printed Table 3 values conflict with the equations when treated as
unrounded exact numbers:

1. `MAP = CO * TPR`: `83/5.15 = 16.1165`, not the printed TPR value 17.
2. `RBF = MAP/R_renal`: `83/0.9 = 92.2222`, not the printed renal resistance
   value 86.
3. `99 mL/min / 0.9 L/min = 0.11`, not the printed filtration fraction 0.19.
   A filtration fraction near 0.19 can be obtained using renal plasma flow,
   but this reduced model has no hematocrit or renal-plasma-flow state.

MAP, cardiac output, GFR, and renal blood flow were therefore treated as the
primary targets.  TPR and renal resistance were calculated from model
identities.  The filtration-fraction value was not independently fitted.

The source also contains two important parameter inconsistencies:

- Table 1 prints single-nephron afferent/efferent resistances of `6e6` and
  `1e7 mmHg min/L`.  The accompanying text says the typical whole-kidney
  resistances `30` and `50 mmHg min/L` are multiplied by `2e6` nephrons,
  which gives `6e7` and `1e8`.  The source-reproduction profile uses the
  dimensionally consistent latter pair, and the calibration starts from it.
- Table 3 PRA is `350 fmol/mL/h`; with the printed conversion factor
  `X_PRC-PRA=61`, this implies PRC `5.7377 pg/mL`.  Table 1 separately prints
  PRC(0) `15 pg/mL`, while the adjacent text discusses a typical value of 28.
  The observable Table 3 PRA was prioritized.
- Hallow Table S1 prints nominal renin secretion as `97 ng/mL/h`, even though
  PRC and the downstream equilibrium are expressed on a `pg/mL` scale.  The
  pre-calibration reconstruction interpreted the numerical value as
  `97 pg/mL/h`.  The calibrated value is independently derived from the Table
  3 PRA target, the conversion factor, and the 12-minute half-life, so it does
  not depend on that ambiguous printed unit.

## 4. How the constrained joint calibration was performed

### 4.1 Cardiovascular reference

The printed right-atrial equation was retained:

`P_ra = 0.2787 * exp(0.2281 * CO)`.

At `CO=5.15 L/min`, this gives `P_ra=0.9022132710 mmHg`.  Using that value as
both the cardiovascular and kidney reference makes the pressure absolute and
positive while keeping the CVP/Bowman increment exactly zero at baseline.

The required total peripheral resistance is fixed by `MAP/CO`:

`TPR = 83/5.15 = 16.11650485 mmHg min/L`.

The source arterial-resistance constant `K_bar=16.6` and basic venous
resistance `R_bv=3.4` were retained.  Vascularity was therefore set by

`vascularity_ref = K_bar / (TPR - R_bv) = 1.3053901359`.

The mean-filling-pressure shift is then derived by the existing cardiovascular
closure so that venous return equals cardiac output at the same reference.

### 4.2 RAAS and PRA

The target PRC follows directly from the unchanged conversion factor:

`PRC = PRA/X_PRC-PRA = 350/61 = 5.737704918 pg/mL`.

With a 12-minute renin half-life and feedback multipliers normalized to one at
the calibrated operating point, the fitted nominal renin secretion is

`N_rs = PRC*ln(2)/12 = 0.3314228322 pg/mL/min`

or `19.88536993 pg/mL/h`.

All RAAS production, conversion, binding, and clearance rates were left
unchanged.  Those kinetics produce baseline AT1-bound AngII of
`4.716373441 pg/mL`.  The shapes of the published macula-densa, RSNA, and AT1
renin-feedback curves were retained and normalized at the calibrated physical
reference (`MD sodium=5.661100585 mEq/min`, `RSNA=1.322484184`, and
`AT1-bound AngII=4.716373441 pg/mL`).

### 4.3 Renal blood flow and GFR

The source values `K_f=3.905 nL/min/mmHg/nephron`, `N=2e6`, total Bowman
pressure `18 mmHg`, and oncotic pressure `28 mmHg` were retained.  Therefore,
the glomerular pressure needed for `GFR=99 mL/min` is determined rather than
fitted:

`P_gh = 18 + 28 + 99/(3.905*2) = 58.67605634 mmHg`.

At `MAP=83` and `RBF=0.9`, the required upstream pressure-drop resistance is
`(83-58.67605634)/0.9 = 27.02660407 mmHg min/L`; the remaining efferent
resistance is `65.19561815 mmHg min/L`.  The preglomerular and afferent nominal
resistances were scaled together, preserving their source ratio, while the
efferent resistance was fitted separately.  This is the minimum two-group
resistance calibration that identifies both total renal flow and glomerular
pressure without changing `K_f`.

### 4.4 Sodium, aldosterone, ADH, and water balance

Once GFR, plasma sodium, PRA, and pressures were fixed, collecting-duct sodium
reabsorption was adjusted to make urinary sodium equal intake.  Although the
fitted nominal factor is slightly above one, the complete collecting-duct
reabsorption fraction at baseline is `0.957384`, so the physical fractional
output remains below one.

The aldosterone concentration scale was adjusted slightly so that the
steady-state aldosterone stimulus produces `100 pg/mL`.  The ADH concentration
scale was adjusted so the steady-state sodium/osmotic stimulus produces the
water intake associated with `2.1 L/day`.

Finally, the source tubular-water-reabsorption intercept `0.0251 L/min` was
reduced.  With the source value, calculated tubular reabsorption exceeds GFR at
this target and urine flow is forced onto the `0.0003 L/min` numerical floor.
The calibrated intercept makes urine flow `2.1 L/day` with the floor inactive.

## 5. Parameters changed from the source/reproduction profile

The `Source/reproduction` column below is the value used by the executable
profile before this calibration.  Where the papers disagree, the distinction
is important: Karaaslan's nominal cardiovascular state is `MAP=100 mmHg`,
`P_ra=0 mmHg`, and `CO=5 L/min`, whereas Hallow Table 3 reports `MAP=83 mmHg`
and `CO=5.15 L/min`.  Karaaslan gives `n_eta-CD=0.93`, aldosterone `85 ng/L`,
ADH `4 mU/L`, and a water-reabsorption intercept of `0.025 L/min`; Hallow
retains the first three relationships but prints `0.0251 L/min` in Eq. A11.
The present calibration uses the Hallow Table 3 phenotype as its target.

### Cardiovascular and pressure references

| Parameter | Source/reproduction | Calibrated | Change | Reason |
|---|---:|---:|---:|---|
| kidney and cardiovascular `p_ma_ref_mmHg` | 100 | 83 | -17.0% | Match Table 3 MAP and center autonomic adaptation there. |
| `cardiac_output_ref_L_per_min` | 5.0 | 5.15 | +3.0% | Match Table 3 cardiac output. |
| kidney and cardiovascular `p_ra_ref_mmHg` | 0 | 0.902213 | derived | Use the positive pressure from printed Eq. A22 and make the CVP increment neutral there. |
| `vascularity_ref` | 1.0 | 1.305390 | +30.54% | Match the MAP/CO-derived TPR without changing `K_bar` or venous resistance. |
| derived mean-filling-pressure shift | -0.026855 mmHg | +0.442132 mmHg | derived | Close venous return at the calibrated MAP, CO, and `P_ra`. |

### Renal, RAAS, and balance parameters

| Parameter | Source/reproduction | Calibrated | Change | Reason |
|---|---:|---:|---:|---|
| `N_rs`, pg/mL/h | 97 | 19.885370 | -79.50% | Match Table 3 PRA=350 with unchanged PRA/PRC conversion and half-life. |
| `AT1_bound_ANGII_EQ`, pg/mL | 22.444473 | 15.8 | -29.60% | Restore Hallow Table S1's printed value; this is a reconstruction correction, not a departure from the paper. Normalized feedback remains one at the calibrated operating point. |
| `R_preglom_0`, mmHg min/L | 15 | 5.917080 | -60.55% | Reduce upstream pressure loss while matching RBF and glomerular pressure. |
| `r_aa_0`, mmHg min/L per nephron | `6.0e7` | `2.36683185e7` | -60.55% | Same upstream scale as preglomerular resistance; preserves their ratio. |
| `r_ea_0`, mmHg min/L per nephron | `1.0e8` | `1.15926607e8` | +15.93% | Maintain the glomerular pressure needed for GFR=99 while total RBF is 0.9. |
| `n_eta_cd` | 0.93 | 1.034682 | +11.26% | Make urinary sodium equal 0.126 mEq/min at the target filtration load. |
| aldosterone concentration scale | 85 | 87.627180 pg/mL | +3.09% | Produce Table 3 aldosterone=100 at the steady-state stimulus. |
| ADH concentration scale | 4.0 | 5.210848 mU/L | +30.27% | Produce the ADH concentration required for 2.1 L/day water intake. |
| tubular-water-reabsorption intercept | 0.0251 | 0.0193856 L/min | -22.77% | Match 2.1 L/day urine flow with the numerical floor inactive. |

The calibrated initial sodium amount is `2149.5 mEq` instead of `2160 mEq`,
because `143.3 mEq/L * 15 L = 2149.5 mEq`.  RAAS, aldosterone, ADH, and
vascularity state values were also replaced by their calibrated equilibrium
values; these are state initializations, not additional parameter fits.

## 6. Source values deliberately left unchanged

- nephron count `2e6`;
- single-nephron filtration coefficient `K_f=3.905`;
- tubular Bowman pressure `18 mmHg` and oncotic pressure `28 mmHg`;
- CVP-to-Bowman slope `0.1468`;
- sodium intake `0.126 mEq/min`;
- basic venous resistance `3.4` and arterial-resistance constant `16.6`;
- printed blood-volume, mean-filling-pressure, and right-atrial-pressure curve
  coefficients;
- all RAAS cascade kinetic rate constants and peptide half-lives;
- PRA/PRC conversion factor `61`;
- proximal and distal baseline sodium-reabsorption factors.

## 7. Implementation corrections found during calibration

Two coding issues became visible when the reference point was generalized:

1. The reference venous-return resistance assumed vascularity was always one.
   It now uses `K_bar/vascularity_ref`, which is identical for the old profiles
   and correct for the calibrated one.
2. The coupled output dictionary previously allowed cardiovascular placeholder
   autonomic values calculated at `P_ra=0` to overwrite kidney values calculated
   at solved `P_ra`.  Kidney autonomic diagnostics now take precedence.  This
   changes reported RSNA-related outputs but not the already-solved pressure
   loop.

The baseline plot now displays total Bowman pressure rather than the small
CVP-induced increment, preventing a negative increment from being mistaken for
negative physical Bowman pressure.

## 8. Interpretation and limitations

This is an exact baseline calibration, not an independent validation.  A single
steady state cannot uniquely identify every underlying physiological
parameter.  Identifiability was improved by changing only two renal resistance
groups, keeping the preglomerular/afferent ratio fixed, and retaining `K_f`.

The fitted right-atrial pressure follows the source equation and is positive,
but it should not be interpreted as a universal clinical CVP target.  The
calibrated CVP/Bowman slope remains unvalidated against independent pressure
data.  Drug, sodium-loading, heart-failure, and other perturbation responses
must be tested against data not used for this baseline fit before the profile is
used predictively.
