#!/usr/bin/env python3
"""Score assessment files with the verifier's own rules, value by value.

usage: baselines.py <task_dir> name=assessment.json [name=assessment.json ...]
"""
import json
import os
import sys

task = os.path.abspath(sys.argv[1])
sys.path.insert(0, os.path.join(task, "tests"))
import test_outputs as V  # noqa: E402

truth = json.load(open(os.path.join(task, "tests", "truth.json")))
att = [a.split("=", 1) for a in sys.argv[2:]]
res = {n: json.load(open(p)) for n, p in att}
rows = [("heat J/g", truth["delta_h_j_per_g"], lambda r: r["delta_h_j_per_g"], lambda v, t: abs(v - t) / (V.REL_DH * t))]
for k in V.STARTS:
    rows.append((f"TMRad {k} C, h", truth["tmrad_h"][k], lambda r, k=k: r["tmrad_h"][k], lambda v, t: abs(v - t) / (V.REL_TMR * t)))
rows.append(("TD24 C", truth["td24_c"], lambda r: r["td24_c"], lambda v, t: abs(v - t) / V.ABS_TD24))
print("Value and error as a fraction of its tolerance (1.00 is the pass boundary).\n")
print("| value | truth | " + " | ".join(n for n, _ in att) + " |")
print("|---|---|" + "---|" * len(att))
for lab, tv, get, err in rows:
    print(f"| {lab} | {tv:.2f} | " + " | ".join(f"{get(res[n]):.2f} ({err(get(res[n]), tv):.2f})" for n, _ in att) + " |")
print("\n| attempt | values inside tolerance (of 6) | worst error / tolerance | outcome |")
print("|---|---|---|---|")
for n, _ in att:
    e = [err(get(res[n]), tv) for _, tv, get, err in rows]
    print(f"| {n} | {sum(x <= 1 for x in e)} | {max(e):.2f} | {'pass' if max(e) <= 1 else 'fail'} |")
