# Provenance

All data are synthetic, produced by `generate.py` with the fixed seed in the
file (numpy 2.5.3, scipy 1.18.1). From the task root:

    python3 authoring/provenance/generate.py .

rewrites `environment/data/dsc/`, `environment/data/run_log.csv`,
`tests/truth.json` and `generator_internal.json` (kinetic and instrument
parameters, for audit only).

The generator's kinetics: a stabiliser that is consumed (zero order, smoothed
near exhaustion) and inhibits an autocatalytic first step while it lasts, and
a consecutive n-th order second step. The instrument applies an as-found
temperature offset, a heat-flow scale error, a signal lag, a run-specific
baseline with a heat-capacity step, a melting endotherm and noise.

The ground truth is obtained by integrating the generator's own kinetic model:
the sum of the two step heats, the adiabatic time to maximum rate at the four
start temperatures, and the root of TMRad(T) = 24 h. No fit is involved.
