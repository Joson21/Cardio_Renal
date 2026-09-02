# Parameters and Initial Conditions

## 1. Configuration to use

For the calibrated healthy baseline and new coupled experiments, use:

```python
parameters = hallow_table3_calibrated_cardiorenal_parameters()
y0 = hallow_table3_reference_equilibrium_state()
```

This selects the jointly calibrated cardiovascular, renal-hemodynamic, RAAS,
sodium, and water baseline. It retains the fixed CVP-Bowman slope and leaves
RIHP off. The fitted values, derivations, and reasons for every departure from
the source-paper reconstruction are documented in `CALIBRATION_REPORT.md`.

For the pre-calibration source reconstruction with the CVP-Bowman mechanism,
use `cvp_extended_cardiorenal_parameters()` with
`cvp_reference_equilibrium_state()`. For source comparison without the novelty,
use
`corrected_hallow_cardiorenal_parameters()` with
`corrected_hallow_reference_equilibrium_state()`.

### Calibrated overrides

The detailed tables below retain the source-reconstruction values for audit.
The calibrated profile overrides only the following quantities:

| Code name | Source reconstruction | Calibrated value |
|---|---:|---:|
| `P_ma_ref`, mmHg | 100 | 83 |
| `P_ra_ref`, mmHg | 0 | 0.9022132710 |
| `cardiac_output_ref`, L/min | 5 | 5.15 |
| `vascularity_ref` | 1 | 1.3053901359 |
| `N_rs`, pg/mL/min | 1.6166666667 | 0.3314228322 |
| `AT1_bound_ANGII_EQ`, pg/mL | 22.444473 | 15.8 |
| `R_preglom,0`, mmHg min/L | 15 | 5.9170796242 |
| `r_aa,0` | 6.0e7 | 2.3668318497e7 |
| `r_ea,0` | 1.0e8 | 1.1592660660e8 |
| `n_eta_cd` | 0.93 | 1.0346820277 |
| aldosterone scaling reference, ng/L | 85 | 87.6271800199 |
| ADH scaling reference, mU/L | 4 | 5.2108479985 |
| water-reabsorption intercept, L/min | 0.0251 | 0.01938560035 |

The calibrated initial total sodium is `2149.5 mEq`; the RAAS, hormone,
adaptation, and vascularity states are also set to their joint equilibrium.

## 2. External baseline inputs

| Symbol/code name | Value | Unit |
|---|---:|---|
| sodium intake | 0.126 | mEq/min |
| peritubular-pressure difference input | 0 | mmHg relative to reference |

MAP and `P_ra` are algebraic outputs in the coupled model, not external inputs.

## 3. Cardiovascular parameters

All values are in `CardiovascularParameters`.

| Code name | Value | Unit/meaning |
|---|---:|---|
| `closure_mode` | `baseline_centered` | recommended executable closure |
| `vb_offset_L` | 4.560227 | L |
| `vb_amplitude_L` | 2.431217 | L |
| `vb_midpoint_L` | 18.11278 | L ECF |
| `vb_slope_per_L` | 0.47437 | 1/L |
| `pmf_slope_mmHg_per_L` | 7.436 | mmHg/L |
| `pmf_offset_mmHg` | 30.18 | mmHg |
| `pra_scale_mmHg` | 0.2787 | mmHg |
| `pra_flow_coefficient_min_per_L` | 0.2281 | min/L |
| literal vascularity formation prefactor | 1.1312e-4 | 1/min |
| vascularity flow coefficient | 0.4799 | min/L |
| vascularity destruction constant | 1.0e-5 | 1/min |
| arterial resistance constant | 16.6 | mmHg min/L |
| basic venous resistance `R_bv` | 3.4 | mmHg min/L |
| venous-return denominator | 31 | dimensionless |
| `V_ECF_ref` | 15 | L |
| `vascularity_ref` | 1 | dimensionless |
| `cardiac_output_ref` | 5 | L/min |
| `P_ma_ref` | 100 | mmHg |
| `P_ra_ref` | 0 | mmHg |
| MAP algebraic search range | 20 to 250 | mmHg |
| algebraic root tolerance | 1.0e-10 | MAP coordinate tolerance |

Derived reference values, calculated rather than independently fitted:

| Derived value | Value |
|---|---:|
| blood volume at `V_ECF=15` | 5.0122877038 L |
| raw `P_ra` at `CO=5` | 0.8718661675 mmHg |
| venous-return resistance at reference | 1.4129032258 mmHg min/L |
| raw mean filling pressure at reference | 7.0913713654 mmHg |
| required mean filling pressure at reference | 7.0645161290 mmHg |
| centered mean-filling-pressure shift | -0.0268552364 mmHg |

## 4. Renal hemodynamic parameters

| Parameter | Value | Unit/comment |
|---|---:|---|
| functioning nephrons | 2.0e6 | nephrons |
| preglomerular resistance | 15.0 | mmHg min/L |
| nominal afferent resistance | 6.0e7 | single-nephron-equivalent mmHg min/L term |
| nominal efferent resistance | 1.0e8 | single-nephron-equivalent mmHg min/L term |
| filtration coefficient `K_f` | 3.905 | nL/min/mmHg/nephron |
| tubular Bowman pressure | 18.0 | mmHg |
| glomerular oncotic pressure | 28.0 | mmHg |
| nominal glomerular pressure | 60.0 | mmHg |
| myogenic gain | 5.0 | dimensionless |
| glomerular algebraic tolerance | 1.0e-9 | mmHg |

The afferent/efferent values are the corrected factor-of-ten choices: `30` and
`50 mmHg min/L` multiplied by `2e6` nephrons.

## 5. AT1 vascular-effect parameters

Each empirical multiplier has form `A + B*x - C/x`, where `x` is AT1-bound
AngII in pg/mL.

| Site | A | B | C |
|---|---:|---:|---:|
| preglomerular | 0.8 | 0.055 | 0.185 |
| afferent | 0.8 | 0.055 | 0.185 |
| efferent | 0.925 | 0.05 | 0.17 |

## 6. Renin and systemic RAAS parameters

### 6.1 Agreed renin branch

| Quantity | Value | Unit |
|---|---:|---|
| nominal renin secretion `N_rs` | 97 | pg/mL/h |
| code value `N_rs` | 1.6166666667 | pg/mL/min |
| PRC half-life | 12 | min |
| renin clearance | 0.05776226505 | 1/min |
| reference PRC | 27.988283793246 | pg/mL |
| PRC-to-PRA factor | 61 | fmol/mL/h per pg/mL PRC |
| reference feedbacks | `v_MD=v_RSNA=v_AT1=1` | dimensionless |

### 6.2 RAAS production/conversion rates

| Process | Source-scale value | Code value | Internal unit |
|---|---:|---:|---|
| AGT production | 34620 | 577.0 | fmol/mL/min |
| ACE | 54.1 /h | 0.9016666667 | 1/min |
| chymase | 1.1 /h | 0.0183333333 | 1/min |
| NEP | 1.1 /h | 0.0183333333 | 1/min |
| ACE2 | 2.4 /h | 0.04 | 1/min |
| AngII to AngIV | 23.5 /h | 0.3916666667 | 1/min |
| AT1 binding | 11.8 /h | 0.1966666667 | 1/min |
| AT2 binding | 3.9 /h | 0.065 | 1/min |

### 6.3 RAAS half-lives

| Species | Half-life |
|---|---:|
| AGT | 600 min |
| AngI | 0.5 min |
| AngII | 0.66 min |
| Ang(1-7) | 30 min |
| AngIV | 0.5 min |
| AT1-bound AngII | 12 min |
| AT2-bound AngII | 12 min |

AngII molecular weight is `1046.18 g/mol` for conversion from fmol/mL to
pg/mL.

## 7. Tubular, hormone, and autonomic parameters

| Parameter | Value | Meaning |
|---|---:|---|
| nominal proximal fraction `n_eta_pt` | 0.8 | dimensionless |
| nominal distal fraction `n_eta_dt` | 0.5 | dimensionless |
| nominal collecting-duct fraction `n_eta_cd` | 0.93 | dimensionless |
| distal aldosterone A/B | 0.24 / 0.3 | empirical |
| collecting-duct aldosterone A/B | 0.74 / 0.06 | empirical |
| potassium reference | 5.0 | mEq/L |
| aldosterone reference | 85.0 | ng/L |
| ADH reference | 4.0 | mU/L |
| aldosterone time constant | 30 | min |
| ADH time constant | 6 | min |
| baroreceptor adaptation rate | 0.0000667 | 1/min |
| right-atrial ADH adaptation rate | 0.0007 | 1/min |
| right-atrial ADH gain | 0.2 | 1/mmHg |
| urine-flow floor | 0.0003 | L/min |
| nominal RSNA multiplier | 1 | dimensionless |

## 8. Novel and provisional parameters

| Parameter | Value | Status |
|---|---:|---|
| `P_B,CVP` slope | 0.1468 | fixed current novelty |
| `P_ra,ref` | 0 mmHg | required reference centering |
| RIHP pressure scale | 1 mmHg | provisional normalization |
| `S_P-N,PT` | 3 | provisional |
| `S_P-N,DT` | 3 | provisional |
| `S_P-N,CD` | 3 | provisional |

The `0.1468` slope maps mmHg of `P_ra` change to mmHg of Bowman-pressure
change, so it is numerically dimensionless (`mmHg/mmHg`).

## 9. Nominal paper/agreed initial state

`initial_state(parameters)` is useful for source-reference checks. It is not a
full cardiorenal equilibrium.

| State | Initial value | Unit |
|---|---:|---|
| total sodium | 2160.0 | mEq |
| `V_ECF` | 15.0 | L |
| PRC | 27.988283793246 | pg/mL |
| AGT | 474830.102634539 | fmol/mL |
| AngI | 12.240564481657 | fmol/mL |
| AngII | 6.458821606363 | fmol/mL |
| Ang(1-7) | 20.894402803236 | fmol/mL |
| AngIV | 1.824796522374 | fmol/mL |
| AT1-bound AngII | 21.990739367508 | fmol/mL |
| AT2-bound AngII | 7.268125723160 | fmol/mL |
| normalized aldosterone | 1.0 | dimensionless |
| normalized ADH | 1.0 | dimensionless |
| baroreceptor integral | 0.0 | min |
| right-atrial ADH integral | 0.0 | min |
| vascularity | 1.0 | dimensionless |

The eight RAAS concentrations solve the isolated Hallow cascade under the
agreed unit-feedback assumptions. This does not imply the sodium, water,
aldosterone, ADH, and vascularity blocks are all simultaneously at equilibrium.

## 10. Retained pre-calibration CVP-enabled equilibrium

Use this state only for source-reconstruction audits or reproducing the
pre-calibration behavior with `cvp_extended_cardiorenal_parameters()` and
baseline external inputs. New healthy-baseline experiments should use the
calibrated configuration in Section 1.

| State | Equilibrium value |
|---|---:|
| total sodium, mEq | 2078.06772219359 |
| `V_ECF`, L | 14.5538978853956 |
| PRC, pg/mL | 31.2084493147996 |
| AGT, fmol/mL | 471996.215368979 |
| AngI, fmol/mL | 13.6488910514231 |
| AngII, fmol/mL | 7.20193521776991 |
| Ang(1-7), fmol/mL | 23.2983885566133 |
| AngIV, fmol/mL | 2.03474676043154 |
| AT1-bound AngII, fmol/mL | 24.5208630874133 |
| AT2-bound AngII, fmol/mL | 8.10435305431458 |
| normalized aldosterone | 1.39843091679175 |
| normalized ADH | 0.601991304698793 |
| baroreceptor integral, min | 976.073779120935 |
| right-atrial ADH integral, min | -8.9707788694985 |
| vascularity | 1.08021901986482 |

Verified algebraic outputs:

| Output | Value |
|---|---:|
| MAP | 92.4324859581 mmHg |
| `P_ra` | -0.0313977260 mmHg |
| cardiac output | 4.8392085502 L/min |
| blood volume | 4.9395191249 L |
| mean filling pressure | 6.6655587079 mmHg |
| total peripheral resistance | 19.1007444707 mmHg min/L |
| GFR | 69.6730093299 mL/min |
| renal blood flow | 0.4990051603 L/min |
| plasma sodium | 142.7842725404 mEq/L |
| urinary sodium | 0.126 mEq/min |
| urine and water intake | 0.0003 L/min |
| `P_B,CVP` | -0.0046091862 mmHg |

The maximum absolute state derivative is below `5e-13` in native units.

## 11. Corrected CVP-off full equilibrium

Use this state only with `corrected_hallow_cardiorenal_parameters()`.

| State | Equilibrium value |
|---|---:|
| total sodium, mEq | 2078.14079361099 |
| `V_ECF`, L | 14.5544071921238 |
| PRC, pg/mL | 31.2095064306459 |
| AGT, fmol/mL | 471995.285060547 |
| AngI, fmol/mL | 13.6493533768296 |
| AngII, fmol/mL | 7.20217916708525 |
| Ang(1-7), fmol/mL | 23.2991777369884 |
| AngIV, fmol/mL | 2.03481568288976 |
| AT1-bound AngII, fmol/mL | 24.5216936763565 |
| AT2-bound AngII, fmol/mL | 8.1046275709992 |
| normalized aldosterone | 1.39796705524511 |
| normalized ADH | 0.601991304698793 |
| baroreceptor integral, min | 974.990854776849 |
| right-atrial ADH integral, min | -8.96107189235108 |
| vascularity | 1.08012715685293 |

Its MAP, `P_ra`, CO, and GFR are approximately `92.44054 mmHg`, `-0.03136
mmHg`, `4.83939 L/min`, and `69.65606 mL/min`.

## 12. Stored-equilibrium validity rule

A stored equilibrium belongs to one exact combination of:

- equations;
- parameter values;
- enabled mechanisms;
- closure mode; and
- external baseline inputs.

Do not reuse it after changing any of those. Run
`calculate_coupled_equilibrium.py`, verify residuals and balances, and store the
new vector under a new descriptive name.
