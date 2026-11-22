#!/usr/bin/env python3
"""Second, separately written solver used only to calibrate the tolerances.

Differences from solution/solve.py
  * no stabiliser kinetics: the induction period is treated as a delay. The
    time to 1 % conversion of each isothermal run is fitted with an Arrhenius
    line and extrapolated, and the decomposition after it is modelled on its own
  * kinetics come from the isothermal runs only, as conversion-time curves
    (integral data), not from a heat-flow fit of all runs
  * model after the delay: Kamal-Sourour  da/dt = (ka(T) + kb(T) a^m)(1-a)^n
    with SEPARATE activation energies for initiation and propagation; the
    second exotherm is not modelled at all
  * baselines are built by hand: end-of-run offset for isothermal runs,
    tangential sigmoid for the dynamic scans (used for the total heat only)
  * the adiabatic prediction uses the first-step heat measured isothermally

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
    early = (t >= 300) & (t <= 1200)
    line = np.polyfit(t[early], q[early], 1)
    q = np.where(t < 300, np.polyval(line, t), q)
    cum = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
    n = 300
    idx = np.linspace(0, len(t) - 1, n).astype(int)
    return t[idx], cum[idx] / cum[-1], cum[-1]


A0 = 0.01


def ks_curve(p, T, tt):
    """Conversion from A0 at tt[0] onwards."""
    lka, Ea, lkb, Eb, m, n = p
    ka = math.exp(lka - Ea / RG * (1 / T - 1 / TREF))
    kb = math.exp(lkb - Eb / RG * (1 / T - 1 / TREF))
    s = solve_ivp(lambda t, a: (ka + kb * max(a[0], 0.0) ** m) * max(1 - a[0], 0.0) ** n, (tt[0], tt[-1]), [A0],
                  t_eval=tt, method="LSODA", rtol=1e-8, atol=1e-11)
    return s.y[0]


def indium_check(data):
    """Offset from the steepest-tangent onset (mean over the rates), scale from the peak area."""
    ons, areas = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "IN_check_*.txt"))):
        _, Tc, t, y = load(f)
        q = y - np.median(y[Tc < 154])
        k = int(np.argmin(q))
        Ts = np.polyval(np.polyfit(t, Tc, 1), t)          # smoothed temperature axis
        best = (0.0, None)
        for j in range(4, k - 4):
            c = np.polyfit(Ts[j - 4:j + 5], q[j - 4:j + 5], 1)
            if c[0] < best[0]:
                best = (c[0], -c[1] / c[0])
        ons.append(best[1])
        w = (Tc > 156) & (Tc < 172)
        areas.append(-np.sum(0.5 * (q[w][1:] + q[w][:-1]) * np.diff(t[w])) / 28.62)
    return float(np.min(ons)) - 156.60, float(np.mean(areas))


def main():
    data, outp = sys.argv[1], sys.argv[2]
    off, scale = indium_check(data)
    leaked = set()
    for line in open(os.path.join(data, "run_log.csv")).read().splitlines()[1:]:
        c = line.split(",")
        if float(c[3]) - float(c[4]) > 0.2:
            leaked.add(c[0])
    print("offset", off, "scale", scale, "leaked", sorted(leaked))
    heats, iso = [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "INT4_*.txt"))):
        if os.path.basename(f) in leaked:
            continue
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

    # delay: time to A0 conversion, Arrhenius in 1/T
    Ts, tind, post = [], [], []
    for T, tt, a, _ in iso:
        ti = float(np.interp(A0, a, tt))
        Ts.append(T)
        tind.append(ti)
        keep = tt >= ti
        post.append((T, np.concatenate([[ti], tt[keep]]), np.concatenate([[A0], a[keep]])))
    cd = np.polyfit(1.0 / np.array(Ts), np.log(tind), 1)
    print("delay activation energy %.0f J/mol, delays (min)" % (cd[0] * RG), np.round(np.array(tind) / 60, 1))

    def res(p):
        return np.concatenate([ks_curve(p, T, tt) - a for T, tt, a in post])

    best = None
    for z in (0.003, 0.02):
        for kb in (5e-4, 2e-3):
            p0 = [math.log(z * kb), 90e3, math.log(kb), 90e3, 1.0, 1.0]
            r = least_squares(res, p0, bounds=([-25, 3e4, -15, 3e4, 0.3, 0.3], [0, 2.5e5, 0, 2.5e5, 3, 3]), x_scale="jac")
            if best is None or r.cost < best.cost:
                best = r
    p = best.x
    print("KS fit", p, "rms", math.sqrt(2 * best.cost / res(p).size), "H1", h1, "dH", dh)
    lka, Ea, lkb, Eb, m, n = p

    def tmrad(T0c):
        dT = h1 / CP
        delay = math.exp(np.polyval(cd, 1.0 / (T0c + 273.15)))

        def rate(a, T):
            return (math.exp(lka - Ea / RG * (1 / T - 1 / TREF)) + math.exp(lkb - Eb / RG * (1 / T - 1 / TREF)) * max(a, 0) ** m) * max(1 - a, 0) ** n
        ev = lambda t, a: a[0] - 0.999
        ev.terminal = True
        s = solve_ivp(lambda t, a: [rate(a[0], T0c + 273.15 + dT * a[0])], (0, 5e8), [A0], method="LSODA", rtol=1e-9,
                      atol=1e-13, events=ev, dense_output=True, max_step=3600)
        te = s.t_events[0][0]
        tt = np.unique(np.clip(np.concatenate([np.linspace(0, te, 2000), te - np.geomspace(1e-3, min(3600, te), 4000)]), 0, te))
        a = s.sol(tt)[0]
        r = [rate(x, T0c + 273.15 + dT * x) for x in a]
        return (delay + tt[int(np.argmax(r))]) / 3600.0

    out = {"delta_h_j_per_g": dh, "tmrad_h": {str(T0): float(tmrad(T0)) for T0 in (110, 120, 130, 140)},
           "td24_c": float(brentq(lambda T0: math.log(tmrad(T0) / 24.0), 60, 170, xtol=1e-3))}
    json.dump(out, open(outp, "w"), indent=1, sort_keys=True)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
