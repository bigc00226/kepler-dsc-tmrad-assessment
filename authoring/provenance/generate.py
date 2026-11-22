#!/usr/bin/env python3
"""Seeded generator for the DSC decomposition study.

Writes
  <out>/environment/data/dsc/*.txt     DSC exports (dynamic and isothermal runs)
  <out>/tests/truth.json               ground truth from the kinetic model itself

Kinetic model (the "plant"):
  stabiliser       consumed        ds/dt  = -kS(T) s/(s+0.02),  s = 1 as supplied
  step 1  A -> B   autocatalytic   dx1/dt = k1(T) (1-x1)^n1 (z + x1^m) / (1 + Ki s)
  step 2  B -> C   n-th order      dx2/dt = k2(T) (x1-x2)^n2
  heat release per gram            q = H1 dx1/dt + H2 dx2/dt
Truth: total heat H1+H2, adiabatic time to maximum rate from fresh material at
four start temperatures and the temperature at which that time is 24 h, all
obtained by integrating the model above (no fitting, no solver output).
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.signal import lfilter

SEED = 7741903
R = 8.314462618
K = dict(E1=96.0e3, k1ref=1.6e-3, T1ref=433.15, n1=1.3, m=0.9, z=0.012, H1=620.0,
         E2=135.0e3, k2ref=3.7e-3, T2ref=538.15, n2=1.1, H2=310.0,
         ES=148.0e3, kSref=4.2e-4, TSref=423.15, Ki=400.0)
CP_MASS = 1.90            # J/(g K), stated in the instruction
TAU_S = 20.0              # signal time constant with the high-pressure crucible, s
T_OFF = 3.6               # as-found temperature scale error, K (indicated - true)
S_HEAT = 0.932            # as-found heat-flow scale (indicated / true)
IN_TM, IN_HF, R_TH = 156.60, 28.62, 45.0   # indium melting point C, heat of fusion J/g, sensor thermal resistance K/W
NOISE_W = 3.0e-6          # W rms per sample
T_MELT, H_FUS = 83.0, 95.0
T_STARTS = [110.0, 120.0, 130.0, 140.0]


def k1(T):
    return K["k1ref"] * np.exp(-K["E1"] / R * (1.0 / T - 1.0 / K["T1ref"]))


def k2(T):
    return K["k2ref"] * np.exp(-K["E2"] / R * (1.0 / T - 1.0 / K["T2ref"]))


def kS(T):
    return K["kSref"] * np.exp(-K["ES"] / R * (1.0 / T - 1.0 / K["TSref"]))


def rates(T, x1, x2, st):
    """st = stabiliser left, as a fraction of the as-supplied content."""
    a = np.clip(1.0 - x1, 0.0, 1.0)
    b = np.clip(x1 - x2, 0.0, 1.0)
    sp = np.clip(st, 0.0, None)
    inh = 1.0 / (1.0 + K["Ki"] * sp)
    r1 = inh * k1(T) * a ** K["n1"] * (K["z"] + np.clip(x1, 0.0, 1.0) ** K["m"])
    r2 = k2(T) * b ** K["n2"]
    rs = kS(T) * sp / (sp + 0.02)
    return r1, r2, rs


def run_program(tgrid, Tfun, leak=None):
    def f(t, y):
        r1, r2, rs = rates(Tfun(t), y[0], y[1], y[2])
        return [r1, r2, -rs]
    s = solve_ivp(f, (tgrid[0], tgrid[-1]), [0.0, 0.0, 1.0], t_eval=tgrid, method="LSODA", rtol=1e-10, atol=1e-13)
    T = Tfun(tgrid)
    r1, r2, _ = rates(T, s.y[0], s.y[1], s.y[2])
    return T, s.y[0], s.y[1], K["H1"] * r1 + K["H2"] * r2


def tmrad(T0c):
    """Adiabatic time to maximum self-heat rate, hours, material as supplied."""
    dT1, dT2 = K["H1"] / CP_MASS, K["H2"] / CP_MASS

    def f(t, y):
        T = T0c + 273.15 + dT1 * y[0] + dT2 * y[1]
        r1, r2, rs = rates(T, y[0], y[1], y[2])
        return [r1, r2, -rs]

    def done(t, y):
        return y[0] - 0.999
    done.terminal = True
    s = solve_ivp(f, (0, 3.0e8), [0.0, 0.0, 1.0], method="LSODA", rtol=1e-10, atol=1e-14, events=done, max_step=3600.0,
                  dense_output=True)
    te = s.t_events[0][0]
    tt = np.concatenate([np.linspace(0, te, 6000), te - np.geomspace(1e-4, min(3600.0, te), 6000)])
    tt = np.unique(np.clip(tt, 0, te))
    y = s.sol(tt)
    T = T0c + 273.15 + dT1 * y[0] + dT2 * y[1]
    r1, r2, _ = rates(T, y[0], y[1], y[2])
    rate = dT1 * r1 + dT2 * r2
    return float(tt[int(np.argmax(rate))] / 3600.0)


def write_export(path, name, mass, rng_line, seg_line, T, t, w_per_g_signal):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        hdr = [
            "#EXPORTTYPE:DATA ALL",
            f"#FILE:{name}.ngb-sd8",
            "#FORMAT:NETZSCH5",
            "#IDENTITY:PS-2291",
            "#DECIMAL:POINT",
            "#SEPARATOR:SEMICOLON",
            "#MTYPE:DSC",
            "#INSTRUMENT:NETZSCH DSC 214",
            "#PROJECT:Intermediate INT-4 thermal stability",
            "#LABORATORY:Process Safety Lab",
            "#SAMPLE:INT-4 batch 0427",
            f"#SAMPLE MASS /mg:{mass:.3f}",
            "#CRUCIBLE:High pressure, gold plated, sealed",
            "#SAMPLE CRUCIBLE MASS /mg:0",
            "#REFERENCE:empty crucible",
            "#PURGE 1 MFC:NITROGEN",
            f"#RANGE:{rng_line}",
            f"#SEGMENT:{seg_line}",
            "#SEG. 1:1/1",
            "##Temp./°C;Time/min;DSC/(mW/mg)",
        ]
        fh.write("\r\n".join(hdr) + "\r\n")
        for a, b, c in zip(T, t / 60.0, w_per_g_signal):
            fh.write(f"{a:.3f};{b:.5f};{c:.6f}\r\n")


def measure(rng, t, T, x1, q, mass_mg, beta_k_s, iso, extra=None):
    """Turn the model heat release (W/g) into the exported signal (mW/mg, exo up)."""
    m = mass_mg * 1e-3                                    # g
    cp_s = 1.60 - 0.35 * x1                               # J/(g K), reactant -> products
    base = -m * beta_k_s * (cp_s - 0.55)                  # W, heat-capacity term
    a = rng.uniform(-25e-6, 25e-6)
    b = rng.uniform(-0.12e-6, 0.12e-6)
    base = base + a + b * (T - 373.15)
    phi = m * q + base
    if extra is not None:
        phi = phi + m * extra
    if not iso:
        sig = 1.0 + 0.15 * beta_k_s * 60.0
        phi = phi - m * H_FUS * beta_k_s * np.exp(-0.5 * ((T - 273.15 - T_MELT - 0.4 * beta_k_s * 60.0) / sig) ** 2) / (sig * math.sqrt(2 * math.pi))
    else:
        phi = phi + rng.choice([-1.0, 1.0]) * rng.uniform(0.25e-3, 0.5e-3) * np.exp(-(t - t[0]) / 25.0)
    dt = t[1] - t[0]
    al = math.exp(-dt / TAU_S)
    phi = lfilter([1 - al], [1, -al], phi, zi=[al * phi[0]])[0]
    phi = S_HEAT * phi + rng.normal(0, NOISE_W, phi.size)
    return phi / m                                        # W/g == mW/mg


def indium_run(rng, beta_min, mass_mg):
    """As-found calibration check: indium melting, same crucible type."""
    b = beta_min / 60.0
    dt = 0.1
    t = np.arange(0, (175.0 - 140.0) / b + dt, dt)
    Ts = 140.0 + b * t                                  # true sample-side temperature, C
    m = mass_mg * 1e-3
    phi = np.zeros_like(t)
    left = m * IN_HF
    for i in range(len(t)):
        if Ts[i] > IN_TM and left > 0:
            p = (Ts[i] - IN_TM) / R_TH
            p = min(p, left / dt)
            phi[i] = -p                                   # endothermic
            left -= p * dt
    phi = phi + rng.uniform(-20e-6, 20e-6) - m * b * 0.24
    al = math.exp(-dt / TAU_S)
    phi = lfilter([1 - al], [1, -al], phi, zi=[al * phi[0]])[0]
    phi = S_HEAT * phi + rng.normal(0, NOISE_W, phi.size)
    step = 5 if beta_min >= 5 else 20
    return (Ts + T_OFF + rng.normal(0, 0.01, t.size))[::step], t[::step], (phi / m)[::step]


def main(out):
    rng = np.random.default_rng(SEED)
    ddir = os.path.join(out, "environment", "data", "dsc")
    os.makedirs(ddir, exist_ok=True)
    os.makedirs(os.path.join(out, "tests"), exist_ok=True)
    info = {}
    log = []

    def log_run(name, crucible, mass, loss_frac):
        tare = float(rng.uniform(1040.0, 1110.0))
        before = tare + mass
        after = before - loss_frac * mass + float(rng.normal(0, 0.012))
        log.append(f"{name}.txt,{crucible},{mass:.3f},{before:.2f},{after:.2f}")

    for q, beta in enumerate([1.0, 5.0, 10.0]):
        mass = float(rng.uniform(5.0, 7.0))
        Tm_, t_, y_ = indium_run(rng, beta, mass)
        name = f"IN_check_{beta:g}Kmin"
        write_export(os.path.join(ddir, name + ".txt"), name, mass, f"140\u00b0C/{beta:.1f}(K/min)/175\u00b0C", "S1/1", Tm_, t_, y_)
        log_run(name, f"HP-{11 + q:02d}", mass, 0.0)
    # dynamic runs
    for q, (beta, step) in enumerate([(0.5, 8), (1.0, 4), (2.0, 2), (4.0, 2), (8.0, 1)]):
        b = beta / 60.0
        dt = 0.5
        t = np.arange(0, (380.0 - 30.0) / b + dt, dt)
        Tfun = lambda tt, b=b: 303.15 + b * tt
        mass = float(rng.uniform(4.0, 9.0))
        leak, extra, loss = None, None, 0.0
        if False:
            tl = (186.0 - 30.0) / b
            leak = (tl, 0.85, 0.70)
            extra = -70.0 * np.exp(-0.5 * ((t - tl - 20.0) / 14.0) ** 2) / (14.0 * math.sqrt(2 * math.pi))
            loss = 0.21
        T, x1, x2, qh = run_program(t, Tfun, leak)
        y = measure(rng, t, T, x1, qh, mass, b, False, extra)
        Tm = T - 273.15 + T_OFF + rng.normal(0, 0.01, T.size)
        name = f"INT4_dyn_{beta:g}Kmin".replace(".", "p")
        write_export(os.path.join(ddir, name + ".txt"), name, mass,
                     f"30\u00b0C/{beta:.1f}(K/min)/380\u00b0C", "S1/1", Tm[::step], t[::step], y[::step])
        info[name] = {"mass_mg": mass, "peak_T": float(T[np.argmax(qh)] - 273.15), "peak_q": float(qh.max()), "leak": leak is not None}
        log_run(name, f"HP-{21 + q:02d}", mass, loss)
    # isothermal runs
    for q, (Tc, hours, step) in enumerate([(140.0, 16.0, 8), (150.0, 9.0, 4), (160.0, 5.0, 4), (170.0, 3.0, 2)]):
        ramp = 100.0 / 60.0
        t_r = (Tc - 100.0) / ramp
        dt = 0.5
        t = np.arange(0, t_r + hours * 3600.0 + dt, dt)
        Tfun = lambda tt, Tc=Tc, t_r=t_r, ramp=ramp: np.minimum(373.15 + ramp * np.asarray(tt), Tc + 273.15)
        mass = float(rng.uniform(10.0, 14.0))
        leak, extra, loss = None, None, 0.0
        if False:
            tl = t_r + 38.0 * 60.0
            leak = (tl, 0.45, 0.80)
            extra = -55.0 * np.exp(-0.5 * ((t - tl - 25.0) / 16.0) ** 2) / (16.0 * math.sqrt(2 * math.pi))
            loss = 0.14
        T, x1, x2, qh = run_program(t, Tfun, leak)
        keep = t >= t_r
        y = measure(rng, t[keep], T[keep], x1[keep], qh[keep], mass, 0.0, True, None if extra is None else extra[keep])
        Tm = T[keep] - 273.15 + T_OFF + rng.normal(0, 0.01, keep.sum())
        name = f"INT4_iso_{Tc:.0f}C"
        write_export(os.path.join(ddir, name + ".txt"), name, mass,
                     f"{Tc:.0f}\u00b0C/{hours:.1f}h isothermal", "S1/1", Tm[::step], (t[keep] - t_r)[::step], y[::step])
        info[name] = {"mass_mg": mass, "leak": leak is not None}
        log_run(name, f"HP-{31 + q:02d}", mass, loss)
    with open(os.path.join(out, "environment", "data", "run_log.csv"), "w") as fh:
        fh.write("file,crucible,sample_mass_mg,gross_mass_before_mg,gross_mass_after_mg\n" + "\n".join(log) + "\n")
    tm = {f"{int(T0)}": tmrad(T0) for T0 in T_STARTS}
    td24 = brentq(lambda T0: math.log(tmrad(T0) / 24.0), 90.0, 150.0, xtol=1e-4)
    truth = {"delta_h_j_per_g": K["H1"] + K["H2"], "tmrad_h": tm, "td24_c": float(td24)}
    with open(os.path.join(out, "tests", "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1, sort_keys=True)
    with open(os.path.join(out, "authoring", "provenance", "generator_internal.json"), "w") as fh:
        json.dump({"kinetics": K, "instrument": {"tau_s": TAU_S, "t_off": T_OFF, "s_heat": S_HEAT}, "runs": info}, fh, indent=1)
    print(json.dumps(truth, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
