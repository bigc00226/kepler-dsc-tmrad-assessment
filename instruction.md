The process safety lab has finished the DSC study of intermediate INT-4, which decomposes exothermically, and the thermal stability assessment for its storage and handling now has to be written. The assessment needs the heat of decomposition and the adiabatic time to maximum rate, worked out from the measurements.

The DSC exports are in `/app/data/dsc/`: five dynamic scans (`INT4_dyn_0p5Kmin.txt`, `INT4_dyn_1Kmin.txt`, `INT4_dyn_2Kmin.txt`, `INT4_dyn_4Kmin.txt`, `INT4_dyn_8Kmin.txt`), four isothermal runs (`INT4_iso_140C.txt`, `INT4_iso_150C.txt`, `INT4_iso_160C.txt`, `INT4_iso_170C.txt`) and three calibration check runs (`IN_check_1Kmin.txt`, `IN_check_5Kmin.txt`, `IN_check_10Kmin.txt`). `/app/data/run_log.csv` lists the runs with their crucible weighings. `/app/data/README.md` describes the instrument, the calibration check, the export format and how the runs were made.

Determine three things. All of them refer to true temperature and true heat, for INT-4 as supplied, kept in a closed system.

`delta_h_j_per_g`: the total heat of decomposition in joules per gram of sample, as a positive number. It covers the whole exothermic decomposition and excludes the melting endotherm.

`tmrad_h`: the adiabatic time to maximum rate, in hours, for start temperatures of 110, 120, 130 and 140 °C. The scenario is fresh material as supplied at the start temperature in a perfectly adiabatic mass with no thermal inertia from the container, a constant specific heat capacity of 1.90 J/(g K), and no loss of material. The time runs from the start until the rate of temperature rise $ dT/dt $ reaches its maximum.

`td24_c`: the start temperature, in °C, at which the adiabatic time to maximum rate in the same scenario is 24 hours.

Write the answers to `/app/output/assessment.json`. It must be a regular file of at most 100 kB holding one JSON object with exactly the keys `delta_h_j_per_g`, `tmrad_h` and `td24_c`. `tmrad_h` must be an object with exactly the keys `"110"`, `"120"`, `"130"` and `"140"`. Every value must be a finite positive JSON number.

Requirements, all of which must hold. `delta_h_j_per_g` must be within 4 % of the true value. Each of the four `tmrad_h` values must be within 30 % of its true value. `td24_c` must be within 3.0 K of the true value. All bounds are inclusive.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
