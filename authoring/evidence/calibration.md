# Tolerance calibration and wrong routes

Produced with

    O=authoring/evidence/outputs
    python3 authoring/evidence/baselines.py . reference=$O/reference_assessment.json \
        independent=$O/independent_assessment.json \
        stabiliser_variant=$O/stabiliser_variant_assessment.json \
        standard_fit=$O/standard_fit_assessment.json \
        threshold_delay=$O/threshold_delay_assessment.json

Each output file was written by the script of the same name on the shipped data:

| attempt | script | what it does |
|---|---|---|
| reference | `solution/solve.py` | calibration correction, global fit of all runs with stabiliser depletion (zero order, smoothed) inhibiting an autocatalytic step, plus a consecutive step |
| independent | `authoring/evidence/independent_solver.py` | separately written; isothermal conversion curves only, first-order stabiliser depletion with exponential inhibition, hand-built baselines, no second step |
| stabiliser_variant | `authoring/evidence/stabiliser_variant.py` | the reference pipeline with a different stabiliser law (first-order depletion, squared inhibition term) |
| standard_fit | `authoring/evidence/standard_fit.py` | calibration correction and a global two-step autocatalytic fit with no stabiliser |
| threshold_delay | `authoring/evidence/threshold_delay.py` | treats the induction period as a delay equal to the time to 1 % conversion, Arrhenius-extrapolated, then a Kamal-Sourour runaway |

Value and error as a fraction of its tolerance (1.00 is the pass boundary).

| value | truth | reference | independent | stabiliser_variant | standard_fit | threshold_delay |
|---|---|---|---|---|---|---|
| heat J/g | 930.00 | 929.68 (0.01) | 933.24 (0.09) | 929.13 (0.02) | 894.23 (0.96) | 933.24 (0.09) |
| TMRad 110 C, h | 63.78 | 62.45 (0.07) | 51.08 (0.66) | 60.56 (0.17) | 21.46 (2.21) | 35.02 (1.50) |
| TMRad 120 C, h | 20.69 | 20.34 (0.06) | 17.73 (0.48) | 19.94 (0.12) | 9.29 (1.84) | 14.07 (1.07) |
| TMRad 130 C, h | 7.23 | 7.14 (0.04) | 6.55 (0.31) | 7.06 (0.08) | 4.20 (1.40) | 5.92 (0.60) |
| TMRad 140 C, h | 2.72 | 2.70 (0.03) | 2.57 (0.18) | 2.68 (0.05) | 1.97 (0.92) | 2.58 (0.18) |
| TD24 C | 118.64 | 118.48 (0.05) | 117.08 (0.52) | 118.29 (0.12) | 108.70 (3.31) | 114.08 (1.52) |

| attempt | values inside tolerance (of 6) | worst error / tolerance | outcome |
|---|---|---|---|
| reference | 6 | 0.07 | pass |
| independent | 6 | 0.66 | pass |
| stabiliser_variant | 6 | 0.17 | pass |
| standard_fit | 2 | 3.31 | fail |
| threshold_delay | 3 | 1.52 | fail |

Graded values a wrong route clears

- `delta_h_j_per_g` is cleared by every route, including `standard_fit`
  (0.96 of the tolerance). It only needs the dynamic scans integrated after the
  calibration correction. It is graded because the assessment needs it, not to
  separate routes.
- `tmrad_h["140"]` is cleared by `standard_fit` (0.92) and `threshold_delay`
  (0.18). At 140 °C the stabiliser is used up quickly and its kinetics matter
  little.
- `tmrad_h["130"]` is cleared by `threshold_delay` (0.60).
- `tmrad_h["110"]`, `tmrad_h["120"]` and `td24_c` are cleared only by the
  routes that model the stabiliser. They carry the separation.

Why the tolerances are what they are

- The time to maximum rate at 110 °C is extrapolated about 30 K below the
  coolest isothermal run. There, a 5 kJ/mol change in the stabiliser's
  activation energy moves the result by about 16 %. The two correct
  stabiliser laws and the separately written solver span -20 % to -2 % at
  110 °C, so the bound is 30 % and not tighter.
- `td24_c` is bounded at 3.0 K; the correct routes span -1.6 K to -0.2 K, the
  wrong routes are -4.6 K and -9.9 K.
- `delta_h_j_per_g` is bounded at 4 %; the correct routes span 0.0 % to +0.3 %.

Verifier checks run by hand: reference and independent outputs give 4 of 4
tests passing, standard_fit gives 2 of 4, and no file gives 0 of 4.
