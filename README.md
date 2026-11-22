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
grade that carries a stabiliser: the stabiliser is consumed first and holds
back an autocatalytic decomposition while it lasts, and a consecutive step
gives a second exotherm. The instrument has as-found temperature and heat-flow
scale errors that only the shipped indium check reveals, a signal lag from the
high-pressure crucibles, run-specific baselines with a heat-capacity step, a
melting endotherm, a settling transient on the isothermal runs, and noise. The
export format, temperature programmes, crucible weighings and the size of the
effects follow real DSC practice.

The difficulty is in the kinetics and the extrapolation:

- The answers are needed 10 to 30 K below the coolest isothermal run. In that
  range TMRad is governed by how long the stabiliser lasts, and its depletion
  has a much higher activation energy than the decomposition itself.
- The isothermal runs show a flat induction period followed by an
  autocatalytic peak. An autocatalytic model with an initiation term can
  produce an induction period too, so the stabiliser has to be recognised from
  how the induction time and the peak scale with temperature, and modelled.
  A global two-step autocatalytic fit without it, which reproduces the scans'
  peaks, gives TMRad at 110 °C of 21 h against a true 64 h and TD24 10 K too
  low.
- Taking the induction period as a delay read off a conversion threshold
  mixes the start of the decomposition into it and gives TMRad at 110 °C of
  35 h.
- The temperature axis and heat flow must be corrected from the indium check
  (onset extrapolated to zero heating rate, peak area against the known heat
  of fusion) before any kinetics.

## Reference solution

`solution/solve.py`, run by `solution/solve.sh`, about six minutes on one core.

1. From the three indium runs: the temperature offset as the melting onset
   extrapolated to zero heating rate, the heat-flow scale from the peak area,
   and the signal time constant from the decay after the peak. Every INT-4
   run is corrected with them.
2. One kinetic model is fitted to all nine INT-4 runs at once: a stabiliser
   consumed at a zero-order rate (smoothed near exhaustion) that inhibits an
   autocatalytic first step by a factor 1/(1 + Ki s), followed by a
   consecutive n-th order step, each step with its own heat. The model is
   integrated along each run's measured temperature programme, the signal lag
   is applied to the dynamic scans, and each run's baseline (offset, slope and
   a heat-capacity step proportional to conversion) is eliminated inside the
   residual by linear least squares. Start values for the stabiliser come from
   an Arrhenius line through the induction times.
3. The fitted model is integrated adiabatically from material as supplied with
   the stated heat capacity to get TMRad, and TD24 is the root of
   TMRad(T) = 24 h.

## Verification

Ground truth is `tests/truth.json`, obtained by integrating the generator's
own kinetic model: total heat 930.0 J/g, TMRad 63.78, 20.69, 7.23 and 2.72 h
at 110, 120, 130 and 140 °C, TD24 118.64 °C. No fitted value is part of it.

`tests/test_outputs.py` reads only `/app/output/assessment.json`. It requires a
regular file of at most 100 kB with exactly the keys `delta_h_j_per_g`,
`tmrad_h` (keys `"110"`, `"120"`, `"130"`, `"140"`) and `td24_c`, all finite
positive numbers; links, duplicate keys, NaN and infinity are rejected. The
heat must be within 4 %, each TMRad within 30 % and TD24 within 3.0 K.

The bounds were set from three correct routes: the reference, a variant of it
with a different stabiliser law, and `authoring/evidence/independent_solver.py`,
written separately, which uses only the isothermal conversion curves, a
first-order stabiliser with exponential inhibition, hand-built baselines and no
second step. TMRad at 110 °C is an extrapolation 30 K below the data, where a
5 kJ/mol change in the stabiliser's activation energy moves it by about 16 %;
the correct routes span -20 % to -2 % there, hence 30 %. Error as a fraction of
the tolerance, from `authoring/evidence/calibration.md`:

| route | heat | TMRad 110 | TMRad 120 | TMRad 130 | TMRad 140 | TD24 | outcome |
|---|---|---|---|---|---|---|---|
| reference | 0.01 | 0.07 | 0.06 | 0.04 | 0.03 | 0.05 | pass |
| independent solver | 0.09 | 0.66 | 0.48 | 0.31 | 0.18 | 0.52 | pass |
| stabiliser-law variant | 0.02 | 0.17 | 0.12 | 0.08 | 0.05 | 0.12 | pass |
| two-step fit, no stabiliser | 0.96 | 2.21 | 1.84 | 1.40 | 0.92 | 3.31 | fail |
| delay from a 1 % threshold | 0.09 | 1.50 | 1.07 | 0.60 | 0.18 | 1.52 | fail |

Graded values a wrong route clears, named one by one:

- The heat of decomposition is cleared by every route; it only needs the
  dynamic scans integrated after the calibration correction.
- TMRad at 140 °C is cleared by both wrong routes, and TMRad at 130 °C by the
  threshold delay. Near the measured range the stabiliser matters little.
- TMRad at 110 and 120 °C and TD24 are cleared only by routes that model the
  stabiliser; they carry the separation.

An empty attempt has no file and fails every test.
