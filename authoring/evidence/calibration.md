# Tolerance calibration and wrong routes

Produced with

    O=authoring/evidence/outputs
    python3 authoring/evidence/baselines.py . reference=$O/reference_assessment.json \
        independent=$O/independent_assessment.json \
        stabiliser_variant=$O/stabiliser_variant_assessment.json \
        single_route=$O/single_route_assessment.json \
        standard_fit=$O/standard_fit_assessment.json \
        isothermal_only=$O/isothermal_only_assessment.json \
        onset_shortcut=$O/onset_shortcut_assessment.json

Each output file was written by the script listed for it below, run on the shipped data:

| attempt | script | what it does |
|---|---|---|
| reference | `solution/solve.py` | calibration correction, global fit of all INT-4 runs with a stabiliser used up by two routes (each zero order, smoothed) that inhibits an autocatalytic step, plus a consecutive step |
| independent | `authoring/evidence/independent_solver.py` | separately written (sharing only the reference's two-route start values and TMRad evaluation); decomposition and one consumption route fitted to the isothermal conversion curves, then all parameters with both routes to those curves together with the dynamic-scan onsets (10 % and 30 % of the main peak), first-order consumption with exponential inhibition, hand-built baselines, no second step |
| stabiliser_variant | `authoring/evidence/stabiliser_variant.py` | the reference pipeline with a different stabiliser law (first-order consumption on both routes, squared inhibition term) |
| single_route | `authoring/evidence/single_route.py` | the reference pipeline with one consumption route |
| standard_fit | `authoring/evidence/standard_fit.py` | calibration correction and a global two-step autocatalytic fit with no stabiliser |
| isothermal_only | `authoring/evidence/isothermal_only.py` | kinetics from the isothermal conversion curves alone, one first-order consumption route with exponential inhibition |
| onset_shortcut | `authoring/evidence/onset_shortcut.py` | no kinetic model: an Arrhenius line through the isothermal induction times (signal at 10 % of its peak) read as the TMRad |

Value and error as a fraction of its tolerance (1.00 is the pass boundary).

| value | truth | reference | independent | stabiliser_variant | single_route | standard_fit | isothermal_only | onset_shortcut |
|---|---|---|---|---|---|---|---|---|
| heat J/g | 930.00 | 930.00 (0.00) | 939.03 (0.24) | 930.94 (0.03) | 945.83 (0.43) | 952.35 (0.60) | 939.03 (0.24) | 964.45 (0.93) |
| TMRad 110 C, h | 18.85 | 17.86 (0.17) | 21.70 (0.50) | 17.88 (0.17) | 52.50 (5.95) | 42.72 (4.22) | 52.75 (6.00) | 66.82 (8.48) |
| TMRad 120 C, h | 9.34 | 9.11 (0.08) | 10.33 (0.35) | 9.10 (0.08) | 16.52 (2.57) | 13.80 (1.59) | 16.83 (2.68) | 19.75 (3.72) |
| TMRad 130 C, h | 4.59 | 4.56 (0.02) | 4.78 (0.14) | 4.51 (0.05) | 5.50 (0.66) | 4.72 (0.10) | 5.69 (0.80) | 6.20 (1.17) |
| TMRad 140 C, h | 2.09 | 2.09 (0.00) | 2.08 (0.01) | 2.05 (0.07) | 1.93 (0.25) | 1.70 (0.62) | 2.03 (0.09) | 2.06 (0.04) |
| TD24 C | 106.68 | 105.76 (0.30) | 108.64 (0.65) | 105.74 (0.31) | 116.72 (3.35) | 115.04 (2.79) | 116.84 (3.39) | 118.37 (3.90) |

| attempt | values inside tolerance (of 6) | worst error / tolerance | outcome |
|---|---|---|---|
| reference | 6 | 0.30 | pass |
| independent | 6 | 0.65 | pass |
| stabiliser_variant | 6 | 0.31 | pass |
| single_route | 3 | 5.95 | fail |
| standard_fit | 3 | 4.22 | fail |
| isothermal_only | 3 | 6.00 | fail |
| onset_shortcut | 2 | 8.48 | fail |

## Graded values a wrong route clears

- `delta_h_j_per_g` is cleared by every route; the straight-baseline area of
  `onset_shortcut` comes closest to the bound (0.93). It only needs the dynamic
  scans integrated after the calibration correction. It is graded because the
  assessment needs it, not to separate routes.
- `tmrad_h["140"]` is cleared by every wrong route: `single_route` (0.25),
  `standard_fit` (0.62), `isothermal_only` (0.09) and `onset_shortcut` (0.04).
  140 °C is the temperature of the coolest isothermal run, so its induction
  time is measured rather than extrapolated. It is graded because the
  assessment needs it, and it stops a route that skips the adiabatic
  integration: the time to maximum rate of an isothermal hold at 140 °C is
  36 % longer than the adiabatic one.
- `tmrad_h["130"]` is cleared by `single_route` (0.66), `standard_fit` (0.10)
  and `isothermal_only` (0.80); `onset_shortcut` fails it (1.17).
- `tmrad_h["110"]`, `tmrad_h["120"]` and `td24_c` are cleared by no wrong
  route: the four wrong routes fail them by 1.6 to 8.5 times the tolerance.
  They carry the separation.

## Why the tolerances are what they are

- The time to maximum rate at 110 °C is extrapolated about 30 K below the
  coolest isothermal run. There the stabiliser is used up almost entirely by
  the low-activation-energy route, and a 5 kJ/mol change in that route's
  activation energy moves the result by 8 to 11 %. The three correct routes
  span -5 % to +15 % at 110 °C and -2 % to 0 % at 140 °C, so the bound is
  30 % and not tighter.
- `td24_c` is bounded at 3.0 K; the correct routes span -0.9 K to +2.0 K. The
  wrong routes are +8.4 K (`standard_fit`), +10.0 K (`single_route`), +10.2 K
  (`isothermal_only`) and +11.7 K (`onset_shortcut`).
- `delta_h_j_per_g` is bounded at 4 %; the correct routes span 0.0 % to
  +1.0 %. The as-found heat-flow scale is 0.932, so a heat integrated without
  the indium correction is about 7 % low and fails.

Verifier checks run by hand: the reference, independent and stabiliser_variant
outputs give 4 of 4 tests passing, `single_route`, `standard_fit`,
`isothermal_only` and `onset_shortcut` give 2 of 4, and a missing output file
gives 0 of 4.
