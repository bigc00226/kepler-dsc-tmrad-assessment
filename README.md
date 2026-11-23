# dsc-tmrad-assessment

From the DSC study of a stabilised, autocatalytically decomposing
intermediate, determine the heat of decomposition, the adiabatic time to
maximum rate (TMRad) at four start temperatures and the temperature at which
it is 24 hours (TD24).

## Difficulty

Process safety chemists and engineers do this when a thermal stability
assessment is written for a reactive intermediate. Storage and process
temperatures are set from the TMRad and TD24 that come out of it.

The data are synthetic and disclosed as such. A seeded simulation
(`authoring/provenance/generate.py`) produces the DSC exports of a commercial
grade that carries a stabiliser: the stabiliser holds back an autocatalytic
decomposition while it lasts, and a consecutive step gives a second exotherm.
The stabiliser is used up by two parallel routes, one with a high and one with
a low activation energy, whose rates are about equal near 145 °C. The exports
carry as-found temperature and heat-flow scale errors that the shipped indium
check quantifies, a signal lag from the high-pressure crucibles, run-specific
baselines with a heat-capacity step, a melting endotherm, a settling transient
on the isothermal runs, and noise. The export format, temperature programmes,
crucible weighings and the size of the effects follow real DSC practice.

The difficulty is in the kinetics and the extrapolation:

- The answers are needed up to 33 K below the coolest isothermal run. There
  TMRad is governed by how long the stabiliser lasts (three quarters of TMRad
  at 110 °C), and the stabiliser is used up almost entirely by the
  low-activation-energy route (97 % of it at 110 °C), which is a minor route
  in the 150 to 170 °C runs. The induction time therefore does not follow one
  Arrhenius line: a line through the isothermal induction times puts TMRad at
  110 °C at 67 h against a true 19 h.
- The second route shows in how the onsets of the dynamic scans move with
  heating rate, the slow scans spending hours below the isothermal range, and
  in the curvature of the isothermal induction times. A global fit with one
  consumption route leaves a residual norm 15 times larger than with two, and
  gives TMRad at 110 °C of 52 h and TD24 10 K too high. An autocatalytic model
  without the stabiliser gives 43 h and TD24 8 K too high.
- The adiabatic scenario has to be integrated: the time to maximum rate of an
  isothermal hold at the start temperature is 36 to 77 % longer than the
  adiabatic one.
- The temperature axis and heat flow must be corrected from the indium check
  (onset extrapolated to zero heating rate, peak area against the known heat
  of fusion) before any kinetics.

## Reference solution

`solution/solve.py`, run by `solution/solve.sh`, single-threaded: about an
hour on one core (61 minutes of CPU time measured on a cloud vCPU).

1. From the three indium runs: the temperature offset as the melting onset
   extrapolated to zero heating rate, the heat-flow scale from the peak area,
   and the signal time constant from the decay after the peak. Every INT-4
   run is corrected with them.
2. One kinetic model is fitted to all nine INT-4 runs at once: a stabiliser
   used up by two parallel routes, each at a zero-order rate (smoothed near
   exhaustion), that inhibits an autocatalytic first step by a factor
   1/(1 + Ki s), followed by a consecutive n-th order step, each step with its
   own heat. The model is integrated along each run's measured temperature
   programme (the dynamic scans from 100 °C, above the melting endotherm,
   starting with the full stabiliser, although route b has used about 2 % of it
   by then in the slowest scan), the signal lag is applied to the dynamic
   scans, and each run's baseline (an offset, and on the dynamic scans also a
   slope and a heat-capacity step proportional to conversion) is eliminated
   inside the residual by linear least squares. Start values for the stabiliser
   come from an Arrhenius line through the induction times. The model is fitted
   first with one consumption route; adding the second route, with route a
   started above and route b below the one-route activation energy, then lowers
   the residual norm from 0.068 to 0.0045.
3. The fitted model is integrated adiabatically from material as supplied with
   the stated heat capacity to get TMRad, and TD24 is the root of
   TMRad(T) = 24 h.

## Verification

Ground truth is `tests/truth.json`, obtained by integrating the generator's
own kinetic model: total heat 930.0 J/g, TMRad 18.85, 9.34, 4.59 and 2.09 h
at 110, 120, 130 and 140 °C, TD24 106.68 °C. No fitted value is part of it.

`tests/test_outputs.py` reads only `/app/output/assessment.json` from the
agent's run and compares it with `tests/truth.json`. It requires a regular file
of at most 100 kB with exactly the keys `delta_h_j_per_g`, `tmrad_h` (keys
`"110"`, `"120"`, `"130"`, `"140"`) and `td24_c`, all finite positive numbers;
links, duplicate keys, NaN and infinity are rejected. The heat must be within
4 %, each TMRad within 30 % and TD24 within 3.0 K.

The bounds were set from three correct routes: the reference, a variant of it
with a different stabiliser law (first-order consumption on both routes,
squared inhibition), and `authoring/evidence/independent_solver.py`, written
separately (it shares only the reference's start values for the two-route stage
and its adiabatic TMRad evaluation), which first fits the decomposition and one
consumption route to the isothermal conversion curves, then all parameters with
both routes to those curves together with the onset temperatures of the dynamic
scans, with first-order consumption, exponential inhibition, hand-built
baselines and no second step. TMRad at 110 °C is an extrapolation 30 K below
the coolest isothermal run, where a 5 kJ/mol change in the activation energy of
the low-activation-energy route moves it by 8 to 11 %; the correct routes span
from -5 % to +15 % there and from -2 % to 0 % at 140 °C, hence 30 %. They span
from -0.9 K to +2.0 K in TD24, hence 3.0 K, and from 0.0 % to +1.0 % in the
heat; a heat integrated without the indium heat-flow correction is about 7 %
low, hence 4 %. Error as a fraction of the tolerance, from
`authoring/evidence/calibration.md`:

| route | heat | TMRad 110 | TMRad 120 | TMRad 130 | TMRad 140 | TD24 | outcome |
|---|---|---|---|---|---|---|---|
| reference | 0.00 | 0.17 | 0.08 | 0.02 | 0.00 | 0.30 | pass |
| independent solver | 0.24 | 0.50 | 0.35 | 0.14 | 0.01 | 0.65 | pass |
| stabiliser-law variant | 0.03 | 0.17 | 0.08 | 0.05 | 0.07 | 0.31 | pass |
| one consumption route | 0.43 | 5.95 | 2.57 | 0.66 | 0.25 | 3.35 | fail |
| no stabiliser | 0.60 | 4.22 | 1.59 | 0.10 | 0.62 | 2.79 | fail |
| isothermal runs only | 0.24 | 6.00 | 2.68 | 0.80 | 0.09 | 3.39 | fail |
| line through the induction times | 0.93 | 8.48 | 3.72 | 1.17 | 0.04 | 3.90 | fail |

Graded values a wrong route clears, named one by one:

- The heat of decomposition is cleared by every route; it only needs the
  dynamic scans integrated after the calibration correction.
- TMRad at 140 °C, the temperature of the coolest isothermal run, is cleared by
  every wrong route: there the induction time is measured, not extrapolated.
  It is graded because the assessment needs it, and it stops a route that
  skips the adiabatic integration: the time to maximum rate of an isothermal
  hold is 36 % too long there.
- TMRad at 130 °C is cleared by the one-route fit, the fit without stabiliser
  and the isothermal-only route; the line through the induction times fails
  it.
- TMRad at 110 and 120 °C and TD24 are cleared by no wrong route: all four
  fail them, by 1.6 to 8.5 times the tolerance. They carry the separation.

An empty attempt has no file and fails every test.
