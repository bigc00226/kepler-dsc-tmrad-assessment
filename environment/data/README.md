# INT-4 DSC study

Notes from the lab on how the runs in `dsc/` were made. `run_log.csv` lists
every run.

## Material

All INT-4 samples are taken from production batch 0427 as supplied. The
commercial grade contains about 0.4 wt% of a phenolic stabiliser; it was not
removed for any run.

## Instrument and crucibles

- Heat-flux DSC, nitrogen purge, empty reference crucible.
- All runs, including the calibration check, used sealed gold-plated
  high-pressure crucibles of the same type. A sealed crucible is a closed
  system.
- Each crucible was weighed with its sample before and after its run. Both
  gross masses are in `run_log.csv`; the balance repeats to about 0.02 mg.

## Calibration check

`IN_check_1Kmin.txt`, `IN_check_5Kmin.txt` and `IN_check_10Kmin.txt` are an
as-found check with indium (melting point 156.60 °C, heat of fusion
28.62 J/g) at 1, 5 and 10 K/min, run immediately before the study. The
instrument was not adjusted after the check, so every export in `dsc/` carries
the as-found temperature and heat-flow scales.

## Export format

Each file is a text export. Header lines start with `#`; `#SAMPLE MASS /mg`
gives the sample mass and `#RANGE` the temperature programme. The data table
follows the line that starts with `##` and has three columns separated by
semicolons: temperature in °C, time in minutes, and the DSC signal in mW per
mg of sample. Exothermic effects are positive.

The signal is the raw instrument signal. It includes the instrument baseline
(an offset and a slow drift that differ from run to run) and, in the dynamic
scans, the heat-capacity contribution of the sample, which changes as the
sample decomposes.

## Dynamic scans

30 °C to 380 °C at 0.5, 1, 2, 4 and 8 K/min. INT-4 melts near 83 °C; the
endotherm there is the melting and is not part of the decomposition.

## Isothermal runs

140, 150, 160 and 170 °C nominal. The sample was heated to the set temperature
at 100 K/min and time zero in the file is the moment the set temperature was
reached. The first two to three minutes of each run are disturbed by the
instrument settling after the heating step.
