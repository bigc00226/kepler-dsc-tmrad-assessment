# Provenance

All data are synthetic, produced by `generate.py` with the fixed seed in the
file (numpy 2.5.3, scipy 1.18.1). From the task root:

    python3 authoring/provenance/generate.py .

rewrites `environment/data/dsc/`, `environment/data/run_log.csv`,
`tests/truth.json` and `authoring/provenance/generator_internal.json` (kinetic
and instrument parameters, for audit only).

The generator's kinetics: a stabiliser that is consumed by two parallel routes
(each zero order, smoothed near exhaustion; route a with a high activation
energy, route b with a low one) and inhibits an autocatalytic first step while
it lasts, and a consecutive n-th order second step. The two routes are about
equal near 145 °C; route a dominates above it (89 % of the consumption at
170 °C) and route b below it (62 % at 140 °C, 97 % at 110 °C). The measurement model applies an
as-found temperature offset, a heat-flow scale error, a signal lag, a
run-specific baseline with a heat-capacity step, the sample's melting
endotherm, a settling transient on the isothermal runs and noise. In the
simulation the sample follows each programme in true temperature and the
as-found error sits in the exported temperature channel, so every exported
temperature reads the as-found offset above its `#RANGE` programme.

The ground truth is obtained by integrating the generator's own kinetic model:
the sum of the two step heats, the adiabatic time to maximum rate at the four
start temperatures, and the root of TMRad(T) = 24 h. No fit is involved.
