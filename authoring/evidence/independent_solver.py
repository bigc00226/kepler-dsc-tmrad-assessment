#!/usr/bin/env python3
"""Second, separately written solver used only to calibrate the tolerances.

Differences from solution/solve.py
  * a different stabiliser law: two first-order consumption routes,
    ds/dt = -(kSa(T) + kSb(T)) s, and exponential inhibition of an
    autocatalytic step, da/dt = kb(T)(z + a^m)(1-a)^n exp(-Ki s)
  * no global heat-flow fit: the isothermal runs enter as conversion-time
    curves (integral data), and the dynamic scans only through the
    temperatures at which their main exotherm first reaches 10 % and 30 % of
    its peak; both consumption routes are fitted to the two together, and the
    second exotherm is not modelled
  * baselines are built by hand: end-of-run offset for the isothermal runs
    (the first 300 s, disturbed by settling and still inside the induction
    period, are set to it), a straight pre-onset line for the dynamic onsets,
    and a tangential sigmoid for the total heat
  * calibration: steepest-tangent onset of indium, lowest over the rates, the
    mean peak area, and the time constant of the decay after the indium peak,
    used only to move each dynamic onset back by heating rate x time constant
  * the adiabatic prediction uses the first-step heat measured isothermally,
    integrated with scipy LSODA

Shared with solution/solve.py: the start values of the two-route stage (the
one-route rate split 0.6 / 0.4, route a started above and route b below the
one-route activation energy) and the adiabatic TMRad evaluation.

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
SIG_A, SIG_T = 0.02, 0.3           # allowed misfit of a conversion point (law differs from the data) and of an onset (K)
FRACS = (0.1, 0.3)                 # onset levels, as fractions of the main exotherm's peak


def load(path):
    rows, rng = [], ""
    for line in open(path, encoding="utf-8"):
        if line.startswith("#RANGE"):
            rng = line
        elif line[0] != "#" and line.strip():
            rows.append([float(v) for v in line.split(";")])
    d = np.array(rows)
    return rng, d[:, 0], d[:, 1] * 60.0, d[:, 2]


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


def dynamic_onset(Tc, y, frac):
    """Temperature where the main exotherm first reaches frac of its peak, above a straight pre-onset line."""
    pre = (Tc > 100) & (Tc < 125)
    q = y - np.polyval(np.polyfit(Tc[pre], y[pre], 1), Tc)
    qs = np.convolve(q, np.ones(9) / 9, mode="same")
    idx = np.where((Tc > 125) & (Tc < 300))[0]
    k = idx[int(np.argmax(qs[idx]))]
    j = idx[int(np.argmax(qs[idx] > frac * qs[k]))]
    return float(Tc[j])


def iso_conversion(t, y):
    off = y[t > 0.97 * t[-1]].mean()
    q = y - off
    q = np.where(t < 300, 0.0, q)      # settling period: still inside the induction period, so at the baseline
    cum = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
    idx = np.linspace(0, len(t) - 1, 200).astype(int)
    return t[idx], cum[idx] / cum[-1], cum[-1]


def indium_check(data):
    ons, areas, taus = [], [], []
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
        after = (np.arange(len(q)) > k) & (q > 0.5 * q[k]) & (q < 0.05 * q[k])
        taus.append(-1.0 / np.polyfit(t[after], np.log(-q[after]), 1)[0])
    return float(np.min(ons)) - 156.60, float(np.mean(areas)), float(np.median(taus))


def unpack(p):
    lkb, Eb, lz, m, n, lkSa, ESa, lKi, lkSb, ESb = p
    return dict(kb=(lkb, Eb), z=math.exp(lz), m=m, n=n, Sa=(lkSa, ESa), Sb=(lkSb, ESb), Ki=math.exp(lKi))


def arr(pair, T):
    return math.exp(pair[0] - pair[1] / RG * (1 / T - 1 / TREF))


def rhs(q, T, y):
    a, st = max(y[0], 0.0), max(y[1], 0.0)
    ra = arr(q["kb"], T) * (q["z"] + a ** q["m"]) * max(1 - a, 0.0) ** q["n"] * math.exp(-q["Ki"] * st)
    return [ra, -(arr(q["Sa"], T) + arr(q["Sb"], T)) * st]


def iso_curve(q, T, tt):
    s = solve_ivp(lambda t, y: rhs(q, T, y), (0.0, tt[-1]), [0.0, 1.0], t_eval=tt, method="LSODA", rtol=1e-8, atol=1e-11)
    return s.y[0]


def model_onsets(q, beta):
    """Onsets (each level in FRACS of the peak rate) of the first step along a scan from 30 C at beta K/s."""
    tend = (380.0 - 30.0) / beta
    tt = np.linspace(0.0, tend, 8001)
    s = solve_ivp(lambda t, y: rhs(q, 303.15 + beta * t, y), (0.0, tend), [0.0, 1.0], t_eval=tt, method="LSODA",
                  rtol=1e-8, atol=1e-11, max_step=60.0)
    T = 303.15 + beta * tt
    a, st = np.clip(s.y[0], 0.0, None), np.clip(s.y[1], 0.0, None)
    kb = np.exp(q["kb"][0] - q["kb"][1] / RG * (1 / T - 1 / TREF))
    r = kb * (q["z"] + a ** q["m"]) * np.clip(1 - a, 0.0, None) ** q["n"] * np.exp(-q["Ki"] * st)
    k = int(np.argmax(r))
    return [float(T[int(np.argmax(r[:k] > f * r[k]))] - 273.15) for f in FRACS]


def main():
    data, outp = sys.argv[1], sys.argv[2]
    off, scale, tau = indium_check(data)
    print("offset", round(off, 3), "scale", round(scale, 4), "tau", round(tau, 1))
    heats, iso, dyn = [], [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "INT4_*.txt"))):
        rng, Tc, t, y = load(f)
        Tc = Tc - off
        y = y / scale
        if "isothermal" in rng:
            tt, a, h = iso_conversion(t, y)
            iso.append((float(np.median(Tc)) + 273.15, tt, a, h))
        else:
            beta = float(np.polyfit(t, Tc, 1)[0])                       # K/s
            heats.append(total_heat_dynamic(Tc, t, y))
            dyn.append((beta, [dynamic_onset(Tc, y, f) - beta * tau for f in FRACS]))   # undo the signal lag
    dh = float(np.mean(heats))
    h1 = float(np.mean([h for *_, h in iso]))
    print("dynamic onsets (C)", [[round(T, 2) for T in Ts] for _, Ts in dyn])

    def res(p, use_dyn):
        q = unpack(p)
        r = [(iso_curve(q, T, tt) - a) / SIG_A for T, tt, a, _ in iso]
        if use_dyn:
            r.append(np.concatenate([(np.array(model_onsets(q, b)) - np.array(Ts)) / SIG_T for b, Ts in dyn]))
        return np.concatenate(r)

    lo = [-15, 3e4, math.log(1e-5), 0.3, 0.3, -25, 5e4, 0.0, -40, 2e4]
    hi = [0, 2.5e5, math.log(0.3), 3, 3, 0, 3e5, math.log(1e4), 0, 3e5]
    # one consumption route from the isothermal curves (route b switched off at its lower bound)
    Ts = np.array([T for T, *_ in iso])
    d0 = np.array([float(np.interp(0.01, a, tt)) for _, tt, a, _ in iso])
    c0 = np.polyfit(1.0 / Ts, np.log(1.0 / d0), 1)
    ES0, lkS0 = -c0[0] * RG, float(np.polyval(c0, 1.0 / TREF))
    best = None
    for lKi in (math.log(4.0), math.log(8.0)):
        for kb in (5e-4, 2e-3):
            p0 = np.clip([math.log(kb), 90e3, math.log(0.01), 1.0, 1.0, lkS0, ES0, lKi, -40, 1e5], np.array(lo) + 1e-9,
                         np.array(hi) - 1e-9)
            free = [0, 1, 2, 3, 4, 5, 6, 7]
            r = least_squares(lambda pf: res(np.r_[pf, p0[8:]], False), p0[free], bounds=(np.array(lo)[free], np.array(hi)[free]),
                              x_scale="jac")
            if best is None or r.cost < best[1]:
                best = (np.r_[r.x, p0[8:]], r.cost)
    p1 = best[0]
    print("one route, isothermal curves:", np.round(p1, 4))
    # two consumption routes from the isothermal curves and the dynamic onsets together
    kS1, ES1 = math.exp(p1[5]), p1[6]
    best = None
    for mb in (0.5, 0.7):
        p0 = p1.copy()
        p0[5], p0[6] = math.log(0.6 * kS1), min(1.2 * ES1, 2.9e5)
        p0[8], p0[9] = math.log(0.4 * kS1), max(mb * ES1, 2.5e4)
        r = least_squares(lambda pf: res(pf, True), p0, bounds=(lo, hi), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    print("two routes:", np.round(p, 4), "rms", math.sqrt(2 * best.cost / best.fun.size), "H1", round(h1, 1), "dH", round(dh, 1))
    q = unpack(p)

    def tmrad(T0c):
        dT = h1 / CP
        ev = lambda t, y: y[0] - 0.999
        ev.terminal = True
        f = lambda t, y: rhs(q, T0c + 273.15 + dT * y[0], y)
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
