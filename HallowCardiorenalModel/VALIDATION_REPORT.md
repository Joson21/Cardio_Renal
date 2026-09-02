# Validation Report

## 1. Status

The delivered version passes equation-level, coupling, equilibrium, stiff
integration, and cross-solver checks. These results establish mathematical
consistency of the implementation. They do not establish that the current
parameter set is fully recalibrated to all healthy-human outputs reported by
Hallow.

Validation environment used for the recorded run:

| Component | Version |
|---|---:|
| Python | 3.12.13 |
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| Matplotlib | 3.10.8 |

## 2. Commands

```bash
python run_coupled_checks.py
python -m unittest discover -s tests -v
python run_coupled_baseline.py
python run_coupled_double_sodium.py
python run_cvp_mechanism_comparison.py
python run_coupled_rihp_example.py
python run_published_cardiovascular_comparison.py
python run_solver_comparison.py
```

All scripts completed successfully. The unit-test suite contained 18 passing
tests: 11 standalone kidney tests and 7 cardiorenal tests.

## 3. Direct equation and interface checks

| Check | Result | Acceptance |
|---|---:|---:|
| cardiovascular `CO-VR` residual at nominal state | about `2.7e-15 L/min` | `<1e-8` |
| reference-centered MAP | 100.000000 mmHg | within `1e-9` |
| reference-centered `P_ra` | about `1.1e-16 mmHg` | within `1e-9` of zero |
| reference-centered cardiac output | 5.000000 L/min | within `1e-9` |
| centered vascularity derivative | about `8.5e-21 /min` | `<1e-12` |
| `P_B,CVP` at `P_ra=0` | 0 mmHg | exact reference neutrality |
| `P_B,CVP` at `P_ra=5` | 0.734 mmHg | exact slope check |
| RAAS nominal-versus-analytical equilibrium max error | `2.91e-10` | `<1e-9` |
| outputs finite | yes | all finite |

Additional interface tests confirm that raising `V_ECF` from 14 to 16 L at the
nominal remaining states raises blood volume, MAP, `P_ra`, and cardiac output.
This demonstrates the cardiovascular-to-kidney link is not one-way.

## 4. Source-reference inconsistency exposed by literal mode

At the nominal paper/agreed state:

| Output | Published literal | Reference-centered |
|---|---:|---:|
| MAP, mmHg | 95.80126489 | 100.00000000 |
| `P_ra`, mmHg | 0.79751500 | approximately 0 |
| raw `P_ra`, mmHg | 0.79751500 | 0.87186617 before centering |
| cardiac output, L/min | 4.60922714 | 5.00000000 |
| mean filling pressure, mmHg | 7.42657535 | 7.06451613 |
| TPR, mmHg min/L | 20.78466997 | 20.00000000 |
| vascularity derivative, 1/min | 2.38494e-6 | approximately 0 |
| GFR, mL/min | 75.61860548 | 83.41725476 |

The literal option is therefore useful for auditing but is not used as the
recommended reference closure.

## 5. Full coupled equilibrium verification

For `cvp_extended_cardiorenal_parameters()` with baseline sodium and
peritubular-pressure inputs:

| Metric | Result |
|---|---:|
| maximum absolute native state derivative | `4.55e-13` |
| sodium intake minus urinary sodium | `7.0e-15 mEq/min` |
| water intake minus urine output | `1.79e-18 L/min` |
| cardiovascular `CO-VR` residual | less than `1e-13 L/min` |

Principal equilibrium outputs:

| Output | Result |
|---|---:|
| MAP | 92.43248596 mmHg |
| `P_ra` | -0.03139773 mmHg |
| cardiac output | 4.83920855 L/min |
| blood volume | 4.93951912 L |
| total peripheral resistance | 19.10074447 mmHg min/L |
| GFR | 69.67300933 mL/min |
| renal blood flow | 0.49900516 L/min |
| plasma sodium | 142.78427254 mEq/L |
| urinary sodium | 0.12600000 mEq/min |
| daily urine volume | 0.432 L/day |
| aldosterone | 118.8666 ng/L (numerically pg/mL) |
| PRA from implemented PRC relation | 1903.715 fmol/mL/h |
| `P_B,CVP` | -0.00460919 mmHg |

The baseline script remains flat for ten days at numerical precision when
started from this state.

## 6. Comparison with Hallow's reported baseline simulation

Hallow Table 3 reports a calibrated baseline simulation. The current
reconstruction has deliberately not been tuned to match it.

| Quantity | Hallow Table 3 value | Current full equilibrium |
|---|---:|---:|
| MAP, mmHg | 83 | 92.43 |
| GFR, mL/min | 99 | 69.67 |
| urinary sodium, mEq/min | 0.126 | 0.126 |
| ECF volume, L | 15 | 14.554 |
| daily urine, L/day | 2.1 | 0.432 |
| sodium concentration, mEq/L | 143.3 | 142.784 |
| cardiac output, L/min | 5.15 | 4.839 |
| renal blood flow, L/min | 0.9 | 0.499 |
| renal vascular resistance, mmHg min/L | 86 | 185.23 |
| total peripheral resistance, mmHg min/L | 17 | 19.10 |
| aldosterone, pg/mL | 100 | 118.87 |
| PRA, fmol/mL/h | 350 | 1903.72 |

Interpretation:

- the implementation achieves exact sodium balance and a stable closed-loop
  pressure, but renal vascular resistance is high and GFR/RBF are low;
- water balance occurs at the hard urine-flow floor, producing an implausibly
  low daily urine output relative to the Hallow target;
- the user-selected renin branch and dynamic whole-system equilibrium produce
  a much higher PRA than Hallow's reported baseline;
- therefore the model is suitable for equation reconstruction, sensitivity
  work, and mechanism-development, but needs a defined calibration exercise
  before quantitative healthy-human prediction.

## 7. Doubled-sodium coupled experiment

Protocol: sodium intake steps from `0.126` to `0.252 mEq/min` at day 2. Both
intervention and control begin at the recommended full equilibrium.

| Day | MAP, mmHg | `P_ra`, mmHg | CO, L/min | GFR, mL/min | RBF, L/min | urinary Na, mEq/min | `V_ECF`, L | `P_B,CVP`, mmHg |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 92.4325 | -0.03140 | 4.8392 | 69.6730 | 0.49901 | 0.12600 | 14.5539 | -0.00461 |
| 2 | 92.4325 | -0.03140 | 4.8392 | 69.6730 | 0.49901 | 0.12600 | 14.5539 | -0.00461 |
| 5 | 98.3935 | 0.07750 | 5.3734 | 80.1134 | 0.53387 | 0.25506 | 15.4329 | 0.01138 |
| 10 | 98.4190 | 0.05024 | 5.2456 | 80.1790 | 0.53405 | 0.25489 | 15.2579 | 0.00737 |
| 30 | 98.3678 | 0.00749 | 5.0375 | 80.0982 | 0.53366 | 0.25252 | 14.9826 | 0.00110 |

The `P_ra` response is small in this long-term empirical closure. Accordingly,
the fixed CVP-Bowman pressure contribution is also small in this particular
protocol.

## 8. CVP-Bowman mechanism isolation

The doubled-sodium protocol was repeated with the CVP-Bowman pathway disabled
and enabled. Each configuration started from its own equilibrium.

At day 30:

| Output | CVP-Bowman off | CVP-Bowman on |
|---|---:|---:|
| MAP, mmHg | 98.365795 | 98.367848 |
| `P_ra`, mmHg | 0.007507 | 0.007492 |
| cardiac output, L/min | 5.037588 | 5.037512 |
| GFR, mL/min | 80.101626 | 80.098236 |
| renal blood flow, L/min | 0.533652 | 0.533655 |
| `P_B,CVP`, mmHg | 0 | 0.00109983 |

The sign is correct: elevated `P_ra` raises Bowman pressure and slightly lowers
GFR. The numerical size is modest because the driver is only about `0.0075
mmHg` at day 30.

## 9. RIHP worked-example check

A `+0.1 mmHg` peritubular-pressure step at day 1 with `S_P-N,i=3` produces:

\[
\gamma_{rihp,PT}=\gamma_{rihp,DT}=\gamma_{rihp,CD}=0.925062.
\]

Urinary sodium rises immediately because all three segmental reabsorption
fractions decrease. The subsequent MAP, `P_ra`, and GFR changes are model
responses to the sodium/water loss. This is a code-path demonstration, not a
calibrated experiment.

## 10. BDF versus Radau

Both solvers used `rtol=1e-8`, `atol=1e-10`, a common 30-day output grid, and the
doubled-sodium protocol.

| Output | Maximum absolute difference | Maximum relative difference |
|---|---:|---:|
| MAP | 2.63625e-6 mmHg | 2.67810e-8 |
| `P_ra` | 5.16390e-8 mmHg | 6.49885e-7 |
| cardiac output | 2.61910e-7 L/min | 4.86607e-8 |
| GFR | 4.30313e-6 mL/min | 5.36505e-8 |
| renal blood flow | 1.31947e-8 L/min | 2.47003e-8 |
| urinary sodium | 5.79709e-8 mEq/min | 2.26553e-7 |
| `P_B,CVP` | 7.58061e-9 mmHg | 6.49885e-7 |

These differences are negligible relative to the model's physiological scale
and support numerical robustness of the reported protocol.

## 11. Remaining validation work before predictive use

1. Reconcile the current user-selected renin branch with Hallow's baseline PRA
   and PRC targets or explicitly define a different target population.
2. Calibrate renal vascular/hemodynamic parameters jointly to GFR, RBF, MAP,
   filtration fraction, and nephron-number assumptions.
3. Revisit the water-handling floor and fit daily urine output without using a
   hard-bound equilibrium as the target mechanism.
4. Validate the `0.1468` relationship against independent CVP/interstitial or
   Bowman-pressure data, including uncertainty in the digitized slope.
5. If renal venous perfusion pressure is added, constrain it with simultaneous
   RBF and GFR data to separate vascular and interstitial pathways.
6. Calibrate RIHP segment sensitivities only after segmental or otherwise
   identifiable sodium-transport data are available.
7. Validate long-term responses against independent sodium-loading, volume,
   CVP, MAP, GFR, and hormone time courses.

Until those steps are completed, figures should be labeled as simulations of
the current reconstructed/extended parameterization, not as validated clinical
predictions.
