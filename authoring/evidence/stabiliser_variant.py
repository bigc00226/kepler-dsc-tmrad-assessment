#!/usr/bin/env python3
"""Variant of the reference with a different stabiliser law (first-order depletion,
squared inhibition term), used to check that the answer does not hinge on the law.

Reference solution: global kinetic model fit of all DSC runs, then adiabatic prediction.

Model: a stabiliser that is consumed first and inhibits the decomposition while
it lasts, an autocatalytic first step with an initiation term, and an n-th order
consecutive step that gives the second exotherm.
Every run is fitted at once. The instrument baseline of each run (offset,
slope and a heat-capacity step proportional to conversion) enters linearly and
is eliminated inside the residual, so no manual baseline construction is needed.
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

DATA = sys.argv[1] if len(sys.argv) > 1 else "/app/data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/app/output/assessment.json"
R = 8.314462618
CP = 1.90
T1REF, T2REF, TSREF = 433.15, 533.15, 423.15
STARTS = [110, 120, 130, 140]
IN_TM, IN_HF = 156.60, 28.62      # indium melting point (C) and heat of fusion (J/g)


def read_run(path):
    mass, rng = None, ""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#SAMPLE MASS"):
                mass = float(line.split(":")[1])
            elif line.startswith("#RANGE"):
                rng = line.split(":", 1)[1]
            elif line and not line.startswith("#"):
                a, b, c = line.split(";")
                rows.append((float(a), float(b) * 60.0, float(c)))
    d = np.array(rows)
    return {"name": os.path.basename(path), "iso": "isothermal" in rng, "T": d[:, 0] + 273.15, "t": d[:, 1],
            "y": d[:, 2], "mass": mass}


def calibration(data):
    """As-found check with indium: temperature offset, heat-flow scale, signal time constant."""
    on, beta, scale, taus = [], [], [], []
    for f in sorted(glob.glob(os.path.join(data, "dsc", "IN_check_*.txt"))):
        r = read_run(f)
        Tc, t, y = r["T"] - 273.15, r["t"], r["y"]
        b = np.polyfit(t, Tc, 1)[0] * 60.0
        base = np.polyfit(Tc[Tc < 155.0], y[Tc < 155.0], 1)
        q = y - np.polyval(base, Tc)
        k = int(np.argmin(q))
        i = np.arange(len(q))
        edge = (i < k) & (q < 0.25 * q[k]) & (q > 0.75 * q[k])
        c = np.polyfit(Tc[edge], q[edge], 1)
        on.append(-c[1] / c[0])
        beta.append(b)
        scale.append(-np.sum(0.5 * (q[1:] + q[:-1]) * np.diff(t)) / IN_HF)
        tail = (i > k) & (q > 0.6 * q[k]) & (q < 0.05 * q[k])
        taus.append(-1.0 / np.polyfit(t[tail], np.log(-q[tail]), 1)[0])
    off = float(np.polyfit(beta, on, 1)[1] - IN_TM)       # onset extrapolated to zero heating rate
    return off, float(np.mean(scale)), float(np.median(taus))


def closed_runs(data):
    """Files whose crucible kept its mass (loss below 1 % of the sample mass)."""
    ok = set()
    with open(os.path.join(data, "run_log.csv")) as fh:
        next(fh)
        for line in fh:
            name, _, m, before, after = line.strip().split(",")
            if float(before) - float(after) < 0.01 * float(m) + 0.05:
                ok.add(name)
    return ok


def prepare(run, npts=700):
    """Cut the analysis window and average the samples into npts bins."""
    T, t, y = run["T"], run["t"], run["y"]
    if run["iso"]:
        keep = t >= 0.0
    else:
        keep = T >= 373.15            # above the melting endotherm
    T, t, y = T[keep], t[keep], y[keep]
    nb = len(t) // npts
    n = nb * npts
    tb = t[:n].reshape(npts, nb).mean(axis=1)
    Tb = T[:n].reshape(npts, nb).mean(axis=1)
    yb = y[:n].reshape(npts, nb).mean(axis=1)
    use = tb >= 300.0 if run["iso"] else np.ones(npts, bool)   # skip the settling transient
    return {"iso": run["iso"], "t": tb, "T": Tb, "y": yb, "use": use, "t0": t[0], "T0": T[0],
            "scale": float(np.max(yb) - np.median(yb))}


def rate_fn(p):
    """p (B,13): step 1, step 2, heats, then stabiliser (lnkS, ES, lnKi)."""
    lk1, E1, n1, m, lz, lk2, E2, n2 = [p[:, i] for i in range(8)]
    lkS, ES, lKi = p[:, 10], p[:, 11], p[:, 12]

    def f(T, x1, x2, st):
        k1 = np.exp(lk1 - E1 / R * (1.0 / T - 1.0 / T1REF))
        k2 = np.exp(lk2 - E2 / R * (1.0 / T - 1.0 / T2REF))
        kS = np.exp(lkS - ES / R * (1.0 / T - 1.0 / TSREF))
        a = np.clip(1.0 - x1, 0.0, 1.0)
        c = np.clip(x1, 0.0, 1.0)
        b = np.clip(x1 - x2, 0.0, 1.0)
        sp = np.clip(st, 0.0, None)
        inh = 1.0 / (1.0 + np.exp(lKi) * sp) ** 2          # variant: squared inhibition term
        return inh * k1 * a ** n1 * (np.exp(lz) + c ** m), k2 * b ** n2, kS * sp           # variant: first-order depletion
    return f


def simulate(run, p):
    """RK4 on the binned grid (with sub-steps). p (B,13) -> x1, r1, r2 arrays (B,N)."""
    f = rate_fn(p)
    B = p.shape[0]
    t, T = run["t"], run["T"]
    tt = np.concatenate([[run["t0"]], t])
    TT = np.concatenate([[run["T0"]], T])
    n = len(t)
    x1 = np.zeros(B)
    x2 = np.zeros(B)
    st = np.ones(B)
    X1 = np.empty((B, n)); R1 = np.empty((B, n)); R2 = np.empty((B, n))
    sub = 4
    for j in range(n):
        h = (tt[j + 1] - tt[j]) / sub
        for q in range(sub):
            Ta = TT[j] + (TT[j + 1] - TT[j]) * q / sub
            Tb = TT[j] + (TT[j + 1] - TT[j]) * (q + 1) / sub
            Tm = 0.5 * (Ta + Tb)
            a1, a2, a3 = f(Ta, x1, x2, st)
            b1, b2, b3 = f(Tm, x1 + 0.5 * h * a1, x2 + 0.5 * h * a2, st - 0.5 * h * a3)
            c1, c2, c3 = f(Tm, x1 + 0.5 * h * b1, x2 + 0.5 * h * b2, st - 0.5 * h * b3)
            d1, d2, d3 = f(Tb, x1 + h * c1, x2 + h * c2, st - h * c3)
            x1 = x1 + h / 6.0 * (a1 + 2 * b1 + 2 * c1 + d1)
            x2 = x2 + h / 6.0 * (a2 + 2 * b2 + 2 * c2 + d2)
            st = np.maximum(st - h / 6.0 * (a3 + 2 * b3 + 2 * c3 + d3), 0.0)
        r1, r2, _ = f(T[j], x1, x2, st)
        X1[:, j], R1[:, j], R2[:, j] = x1, r1, r2
    return X1, R1, R2


TAU = [20.0]


def residuals(runs, P, steps):
    """P (B,10): 8 kinetic parameters, H1, H2. Baselines are projected out per run."""
    B = P.shape[0]
    out = []
    for run in runs:
        X1, R1, R2 = simulate(run, P)
        u = run["use"]
        q = P[:, 8:9] * R1 + (P[:, 9:10] * R2 if steps == 2 else 0.0)
        if not run["iso"]:
            # first-order signal lag of the crucible and sensor
            a = math.exp(-(run["t"][1] - run["t"][0]) / TAU[0])
            ql = np.empty_like(q)
            ql[:, 0] = q[:, 0]
            for j in range(1, q.shape[1]):
                ql[:, j] = a * ql[:, j - 1] + (1 - a) * q[:, j]
            q = ql
        if run["iso"]:
            cols = [np.ones(u.sum())]
        else:
            Tc = run["T"][u] - 500.0
            cols = [np.ones(u.sum()), Tc / 100.0]
        res = np.empty((B, int(u.sum())))
        for b in range(B):
            A = np.array(cols + ([X1[b, u]] if not run["iso"] else [])).T
            d = run["y"][u] - q[b, u]
            c = np.linalg.lstsq(A, d, rcond=None)[0]
            res[b] = (d - A @ c) / run["scale"] / math.sqrt(u.sum())
        # early part of an isothermal run carries the initiation rate: weight it up
        if run["iso"]:
            w = np.where(run["t"][u] < 0.25 * run["t"][u][np.argmax(run["y"][u])], 3.0, 1.0)
            res = res * w[None, :]
        out.append(res)
    return np.concatenate(out, axis=1)


def fit(runs, x0, lo, hi, steps, free):
    x0 = np.array(x0, float)

    def unpack(xf):
        X = np.repeat(x0[None, :], xf.shape[0], axis=0)
        X[:, free] = xf
        return X

    def fun(xf):
        return residuals(runs, unpack(xf[None, :]), steps)[0]

    def jac(xf):
        n = xf.size
        e = 1e-5 * np.maximum(np.abs(xf), 1e-2)
        Xb = np.repeat(xf[None, :], n + 1, axis=0)
        Xb[np.arange(1, n + 1), np.arange(n)] += e
        r = residuals(runs, unpack(Xb), steps)
        return ((r[1:] - r[0]) / e[:, None]).T

    r = least_squares(fun, x0[free], jac=jac, bounds=(np.array(lo)[free], np.array(hi)[free]), x_scale="jac", max_nfev=80)
    x = x0.copy()
    x[free] = r.x
    return x, math.sqrt(2 * r.cost)


def tmrad_h(p, T0c):
    f = rate_fn(p[None, :])
    d1, d2 = p[8] / CP, p[9] / CP

    def rhs(t, y):
        T = T0c + 273.15 + d1 * y[0] + d2 * y[1]
        r1, r2, rs = f(T, np.array([y[0]]), np.array([y[1]]), np.array([y[2]]))
        return [float(r1[0]), float(r2[0]), -float(rs[0])]

    def done(t, y):
        return y[0] - 0.999
    done.terminal = True
    s = solve_ivp(rhs, (0, 5e8), [0.0, 0.0, 1.0], method="LSODA", rtol=1e-9, atol=1e-13, events=done, max_step=3600.0,
                  dense_output=True)
    te = s.t_events[0][0]
    tt = np.unique(np.clip(np.concatenate([np.linspace(0, te, 3000), te - np.geomspace(1e-3, min(3600.0, te), 5000)]), 0, te))
    y = s.sol(tt)
    T = T0c + 273.15 + d1 * y[0] + d2 * y[1]
    r1, r2, _ = f(T, y[0], y[1], y[2])
    return float(tt[int(np.argmax(d1 * r1 + d2 * r2))] / 3600.0)


def main():
    off, scale, TAU[0] = calibration(DATA)
    ok = closed_runs(DATA)
    print("calibration: offset %.2f K, heat scale %.4f, tau %.1f s; closed runs %d" % (off, scale, TAU[0], len(ok)), flush=True)
    files = [f for f in sorted(glob.glob(os.path.join(DATA, "dsc", "INT4_*.txt"))) if os.path.basename(f) in ok]
    runs = []
    for f in files:
        r = read_run(f)
        r["T"] = r["T"] - off          # true temperature
        r["y"] = r["y"] / scale        # true heat flow
        runs.append(prepare(r))
    iso = [r for r in runs if r["iso"]]
    # The isothermal runs show a flat induction period before the autocatalytic peak, and the
    # induction time changes with temperature far more steeply than the peak itself: the
    # as-supplied stabiliser has to be consumed first. Start values for its depletion rate come
    # from an Arrhenius fit of the induction times (signal first reaching 10 % of its peak).
    Ti, ti = [], []
    for r in iso:
        yb = r["y"] - np.median(r["y"][-50:])
        k = int(np.argmax(yb))
        j = int(np.argmax(yb[:k] > 0.1 * yb[k]))
        Ti.append(float(np.median(r["T"])))
        ti.append(float(r["t"][j]))
    sl, ic = np.polyfit(1.0 / np.array(Ti), np.log(1.0 / np.array(ti)), 1)
    ES0 = float(np.clip(-sl * R, 60e3, 250e3))
    lkS0 = float(ic + sl / TSREF)
    print("induction times", np.round(ti, 0), "start ES %.0f" % ES0, flush=True)
    # order: lnk1, E1, n1, m, lnz, lnk2, E2, n2, H1, H2, lnkS, ES, lnKi
    x0 = [math.log(1e-3), 100e3, 1.0, 1.0, math.log(0.03), math.log(3e-3), 130e3, 1.0, 600.0, 300.0, lkS0, ES0, math.log(100.0)]
    lo = [math.log(1e-5), 40e3, 0.3, 0.3, math.log(1e-4), math.log(1e-5), 50e3, 0.3, 100.0, 0.0, math.log(1e-7), 50e3, math.log(1.0)]
    hi = [math.log(1e-1), 200e3, 3.0, 3.0, math.log(0.5), math.log(1e-1), 300e3, 3.0, 2000.0, 2000.0, math.log(1e-1), 300e3, math.log(1e5)]
    # stage A: stabiliser and step 1 on the isothermal runs
    best = None
    for z0 in (0.005, 0.03):
        for k0 in (5e-4, 3e-3):
            xs = list(x0)
            xs[0], xs[4] = math.log(k0), math.log(z0)
            x, c = fit(iso, xs, lo, hi, 1, [0, 1, 2, 3, 4, 8, 10, 11, 12])
            if best is None or c < best[1]:
                best = (x, c)
    x = best[0]
    print("stage A", np.round(x, 4), best[1], flush=True)
    # stage B: step 2 on all runs with step 1 frozen, then everything free
    x, c = fit(runs, x, lo, hi, 2, [5, 6, 7, 9])
    print("stage B", np.round(x, 4), c, flush=True)
    x, c = fit(runs, x, lo, hi, 2, list(range(13)))
    print("stage C", np.round(x, 4), c, flush=True)
    res = {
        "delta_h_j_per_g": float(x[8] + x[9]),
        "tmrad_h": {str(T0): tmrad_h(x, float(T0)) for T0 in STARTS},
        "td24_c": float(brentq(lambda T0: math.log(tmrad_h(x, T0) / 24.0), 60.0, 170.0, xtol=1e-3)),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
