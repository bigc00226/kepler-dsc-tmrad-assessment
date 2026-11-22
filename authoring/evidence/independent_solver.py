#!/usr/bin/env python3
"""Second, separately written solver used only to calibrate the tolerances.

Differences from solution/solve.py
  * a different stabiliser law: first-order depletion ds/dt = -kS(T) s, and
    exponential inhibition of an autocatalytic step,
    da/dt = kb(T)(z + a^m)(1-a)^n exp(-Ki s)
  * kinetics come from the isothermal runs only, as conversion-time curves
    (integral data), not from a heat-flow fit of all runs; the second exotherm
    is not modelled
  * baselines are built by hand: end-of-run offset for isothermal runs,
    tangential sigmoid for the dynamic scans (used for the total heat only)
  * calibration: steepest-tangent onset of indium, lowest over the rates, and
    the mean peak area; no correction for signal lag
  * the adiabatic prediction uses the first-step heat measured isothermally,
    integrated with scipy LSODA

usage: independent_solver.py <data_dir> <out.json>
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

RG = 8.314462618
CP = 1.90
TREF = 428.15
STARTS = (110, 120, 130, 140)


def load(path):
    rows, rng = [], ""
    for line in open(path, encoding="utf-8"):
        if line.startswith("#RANGE"):
            rng = line
        elif line[0] != "#" and line.strip():
            rows.append([float(v) for v in line.split(";")])
    d = np.array(rows)
    return ("isothermal" in rng), d[:, 0], d[:, 1] * 60.0, d[:, 2]


def total_heat_dynamic(Tc, t, y):
    """Tangential-sigmoid baseline between a pre- and a post-reaction window."""
    pre = (Tc > 100) & (Tc < 116)
    post = (Tc > 352) & (Tc < 378)
    lp = np.polyfit(Tc[pre], y[pre], 1)
    lq = np.polyfit(Tc[post], y[post], 1)
    win = (Tc >= 108) & (Tc <= 365)
    a = np.linspace(0, 1, win.sum())
    for _ in range(8):
        base = (1 - a) * np.polyval(lp, Tc[win]) + a * np.polyval(lq, Tc[win])
        q = y[win] - base
        cum = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t[win]))])
        a = np.clip(cum / cum[-1], 0, 1)
    return cum[-1]


def iso_conversion(t, y):
    off = y[t > 0.97 * t[-1]].mean()
    q = y - off
    early = (t >= 300) & (t <= 900)
    q = np.where(t < 300, np.mean(q[early]), q)
    cum = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
    idx = np.linspace(0, len(t) - 1, 400).astype(int)
    return t[idx], cum[idx] / cum[-1], cum[-1]


def indium_check(data):
    ons, areas = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "IN_check_*.txt"))):
        _, Tc, t, y = load(f)
        q = y - np.median(y[Tc < 154])
        k = int(np.argmin(q))
        Ts = np.polyval(np.polyfit(t, Tc, 1), t)
        best = (0.0, None)
        for j in range(4, k - 4):
            c = np.polyfit(Ts[j - 4:j + 5], q[j - 4:j + 5], 1)
            if c[0] < best[0]:
                best = (c[0], -c[1] / c[0])
        ons.append(best[1])
        w = (Tc > 156) & (Tc < 172)
        areas.append(-np.sum(0.5 * (q[w][1:] + q[w][:-1]) * np.diff(t[w])) / 28.62)
    return float(np.min(ons)) - 156.60, float(np.mean(areas))


def rhs(p, T, y):
    """y = (conversion, stabiliser left). First-order depletion, exponential inhibition."""
    lkb, Eb, lz, m, n, lkS, ES, lKi = p
    kb = math.exp(lkb - Eb / RG * (1 / T - 1 / TREF))
    kS = math.exp(lkS - ES / RG * (1 / T - 1 / TREF))
    a, st = max(y[0], 0.0), max(y[1], 0.0)
    return [kb * (math.exp(lz) + a ** m) * max(1 - a, 0.0) ** n * math.exp(-math.exp(lKi) * st), -kS * st]


def iso_curve(p, T, tt):
    s = solve_ivp(lambda t, y: rhs(p, T, y), (0.0, tt[-1]), [0.0, 1.0], t_eval=tt, method="LSODA",
                  rtol=1e-8, atol=1e-11)
    return s.y[0]


def main():
    data, outp = sys.argv[1], sys.argv[2]
    off, scale = indium_check(data)
    print("offset", round(off, 3), "scale", round(scale, 4))
    heats, iso = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "INT4_*.txt"))):
        is_iso, Tc, t, y = load(f)
        Tc = Tc - off
        y = y / scale
        if is_iso:
            tt, a, h = iso_conversion(t, y)
            iso.append((float(np.median(Tc)) + 273.15, tt, a, h))
        else:
            heats.append(total_heat_dynamic(Tc, t, y))
    dh = float(np.mean(heats))
    h1 = float(np.mean([h for *_, h in iso]))
    # start values: Arrhenius line through the times to 1 % conversion
    Ts = np.array([T for T, *_ in iso])
    d0 = np.array([float(np.interp(0.01, a, tt)) for _, tt, a, _ in iso])
    c0 = np.polyfit(1.0 / Ts, np.log(1.0 / d0), 1)
    ES0, lkS0 = -c0[0] * RG, float(np.polyval(c0, 1.0 / TREF))

    def res(p):
        return np.concatenate([iso_curve(p, T, tt) - a for T, tt, a, _ in iso])

    best = None
    for lKi in (math.log(4.0), math.log(8.0)):
        for kb in (5e-4, 2e-3):
            p0 = [math.log(kb), 90e3, math.log(0.01), 1.0, 1.0, lkS0, ES0, lKi]
            lo = [-15, 3e4, math.log(1e-5), 0.3, 0.3, -25, 5e4, 0.0]
            hi = [0, 2.5e5, math.log(0.3), 3, 3, 0, 3e5, math.log(30.0)]
            r = least_squares(res, np.clip(p0, np.array(lo) + 1e-9, np.array(hi) - 1e-9), bounds=(lo, hi),
                              x_scale="jac")
            if best is None or r.cost < best.cost:
                best = r
    p = best.x
    print("fit", np.round(p, 4), "rms", math.sqrt(2 * best.cost / res(p).size), "H1", round(h1, 1), "dH", round(dh, 1))

    def tmrad(T0c):
        dT = h1 / CP
        ev = lambda t, y: y[0] - 0.999
        ev.terminal = True
        f = lambda t, y: rhs(p, T0c + 273.15 + dT * y[0], y)
        s = solve_ivp(f, (0, 5e8), [0.0, 1.0], method="LSODA", rtol=1e-9, atol=1e-13, events=ev, dense_output=True,
                      max_step=3600)
        te = s.t_events[0][0]
        tt = np.unique(np.clip(np.concatenate([np.linspace(0, te, 4000), te - np.geomspace(1e-3, min(3600, te), 4000)]), 0, te))
        y = s.sol(tt)
        r = [f(0, y[:, i])[0] for i in range(y.shape[1])]
        return tt[int(np.argmax(r))] / 3600.0

    out = {"delta_h_j_per_g": dh, "tmrad_h": {str(T0): float(tmrad(T0)) for T0 in STARTS},
           "td24_c": float(brentq(lambda T0: math.log(tmrad(T0) / 24.0), 60, 170, xtol=1e-3))}
    json.dump(out, open(outp, "w"), indent=1, sort_keys=True)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
