#!/usr/bin/env python3
"""Wrong route kept as evidence: a model-free shortcut with no kinetic model at all.

This is not a correct method. It reads the induction time of each isothermal run (the signal first
reaching a given fraction of its peak), fits one Arrhenius line through the four times, and takes
the line, extrapolated to each start temperature, as the TMRad; TD24 is where the line crosses
24 h. The heat is the straight-baseline area of the dynamic scans. No stabiliser kinetics and no
adiabatic integration are involved, so the second, low-activation-energy consumption route that
takes over below about 145 C is missed; see authoring/evidence/calibration.md.

Calibration: indium onset from the steepest tangent (lowest over the rates) and the mean peak area.

usage: onset_shortcut.py <data_dir> <out.json> [fraction_percent, default 10]
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq


def load(path):
    rows = [line.split(";") for line in open(path, encoding="utf-8") if line[0] != "#" and line.strip()]
    d = np.array(rows, float)
    return d[:, 0], d[:, 1] * 60.0, d[:, 2]


def indium(data):
    ons, areas = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "IN_check_*.txt"))):
        T, t, y = load(f)
        q = y - np.median(y[T < 154])
        k = int(np.argmin(q))
        best = (0.0, None)
        for j in range(4, k - 4):
            c = np.polyfit(T[j - 4:j + 5], q[j - 4:j + 5], 1)
            if c[0] < best[0]:
                best = (c[0], -c[1] / c[0])
        ons.append(best[1])
        areas.append(-np.trapezoid(q[(T > 156) & (T < 172)], t[(T > 156) & (T < 172)]) / 28.62)
    return float(np.min(ons)) - 156.60, float(np.mean(areas))


def main():
    data, out = sys.argv[1], sys.argv[2]
    frac = (float(sys.argv[3]) if len(sys.argv) > 3 else 10.0) / 100.0
    off, scale = indium(data)
    heats = []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "INT4_dyn_*.txt"))):
        T, t, y = load(f)
        T, y = T - off, y / scale
        w = (T > 110) & (T < 370)
        line = np.polyfit(T[w][[0, -1]], y[w][[0, -1]], 1)
        heats.append(np.trapezoid(y[w] - np.polyval(line, T[w]), t[w]))
    inv_T, ln_t = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "INT4_iso_*.txt"))):
        T, t, y = load(f)
        keep = t > 300.0                                   # skip the settling transient
        T, t, y = T[keep] - off, t[keep], y[keep]
        ys = np.convolve(y, np.ones(15) / 15, mode="same")
        ys = ys - np.median(ys[-200:])
        k = int(np.argmax(ys))
        j = int(np.argmax(ys[:k] > frac * ys[k]))
        inv_T.append(1.0 / (np.median(T) + 273.15))
        ln_t.append(math.log(t[j]))
    slope, icpt = np.polyfit(inv_T, ln_t, 1)

    def hours(T0):
        return math.exp(icpt + slope / (T0 + 273.15)) / 3600.0

    res = {"delta_h_j_per_g": float(np.mean(heats)), "tmrad_h": {str(T0): hours(T0) for T0 in (110, 120, 130, 140)},
           "td24_c": float(brentq(lambda T0: math.log(hours(T0) / 24.0), 60.0, 170.0))}
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
