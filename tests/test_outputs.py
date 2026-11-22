"""Verifier for /app/output/assessment.json.

Ground truth (tests/truth.json) comes from integrating the generator's kinetic
model directly: total heat, adiabatic time to maximum rate at four start
temperatures and the temperature at which that time is 24 h. No fitted value
or solver output is part of it.
"""
import json
import math
import os
import stat

RESULT = "/app/output/assessment.json"
TRUTH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "truth.json")
STARTS = ("110", "120", "130", "140")
REL_DH = 0.04
REL_TMR = 0.30
ABS_TD24 = 3.0
MAX_BYTES = 100_000


def _truth():
    with open(TRUTH) as fh:
        return json.load(fh)


def _load():
    st = os.lstat(RESULT)
    assert stat.S_ISREG(st.st_mode), "assessment.json must be a regular file, not a link"
    assert 0 < st.st_size <= MAX_BYTES, "assessment.json is empty or larger than 100 kB"

    def no_dup(pairs):
        keys = [k for k, _ in pairs]
        assert len(keys) == len(set(keys)), "duplicate key in assessment.json"
        return dict(pairs)

    def no_const(name):
        raise AssertionError(f"non-finite constant {name} in assessment.json")

    with open(RESULT) as fh:
        return json.load(fh, object_pairs_hook=no_dup, parse_constant=no_const)


def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def dh_ok(v, true):
    return _num(v) and abs(v - true) <= REL_DH * true


def tmr_ok(v, true):
    return _num(v) and abs(v - true) <= REL_TMR * true


def td24_ok(v, true):
    return _num(v) and abs(v - true) <= ABS_TD24


def test_structure():
    """Exactly the three top-level keys and the four start temperatures, all finite positive numbers."""
    res = _load()
    assert isinstance(res, dict), "top level must be a JSON object"
    assert set(res) == {"delta_h_j_per_g", "tmrad_h", "td24_c"}, f"top-level keys differ: {sorted(res)}"
    assert isinstance(res["tmrad_h"], dict), "tmrad_h must be an object"
    assert set(res["tmrad_h"]) == set(STARTS), f"tmrad_h keys differ: {sorted(res['tmrad_h'])}"
    for v in [res["delta_h_j_per_g"], res["td24_c"]] + [res["tmrad_h"][k] for k in STARTS]:
        assert _num(v), "every value must be a finite number"
        assert v > 0, "every value must be positive"


def test_heat_of_decomposition():
    """delta_h_j_per_g within 4 % of the true total heat."""
    res = _load()
    v = res.get("delta_h_j_per_g") if isinstance(res, dict) else None
    assert dh_ok(v, _truth()["delta_h_j_per_g"]), "delta_h_j_per_g outside tolerance"


def test_time_to_maximum_rate():
    """Each tmrad_h value within 30 % of its true value."""
    res = _load()
    tm = res.get("tmrad_h") if isinstance(res, dict) else None
    tm = tm if isinstance(tm, dict) else {}
    truth = _truth()["tmrad_h"]
    bad = [k for k in STARTS if not tmr_ok(tm.get(k), truth[k])]
    assert not bad, f"tmrad_h outside tolerance at start temperatures: {bad}"


def test_td24():
    """td24_c within 3.0 K of the true value."""
    res = _load()
    v = res.get("td24_c") if isinstance(res, dict) else None
    assert td24_ok(v, _truth()["td24_c"]), "td24_c outside tolerance"
