#!/usr/bin/env python3
"""Recompute the numbers quoted in the article from the per-run tables and compare them with the printed values.

Each line is one statement of the article (abstract, Results, Online Resource 1).  Exit status 0 if all agree.
"""
import csv
import json
import sys

import numpy as np
from scipy import stats

import common as C

D, D5 = C.load_all(), C.load_five()
CONDS, NS = C.CONDS, C.NS
r3 = lambda x: round(x + 0.0, 3)
u = lambda m, c, n, q: float(np.mean(C.vec_all(D, m, c, n, q)))                 # all replicates
up = lambda a, b, c, n, q: tuple(r3(v) for v in C.paired(C.vec_all(D, a, c, n, q), C.vec_all(D, b, c, n, q)))
m5 = lambda c, n, m, q: float(np.mean(C.vec5(D5, c, n, m, q)))                 # replicates 1-5
FAILED = []


def chk(text, ok):
    print(("ok    " if ok else "FAIL  ") + text)
    if not ok:
        FAILED.append(text)


# ---- design ----------------------------------------------------------------------------------------------------
chk("annealed SMC at alpha = 0.99999: 420 runs, about 7,600 CPU-hours",
    sum(1 for k in D if k[0] == "hi") == 420
    and round(sum(v["cpu_s"] for k, v in D.items() if k[0] == "hi") / 3600.0, -2) == 7600)

# ---- accuracy at N = 100 on 50 replicates ----------------------------------------------------------------------
chk("N = 100: mean MCC error, alpha = 0.99999 / BEAST X: 0.046/0.034, 0.021/0.020, 0.062/0.056, 0.015/0.014, 0.055/0.055, 0.064/0.054",
    [(r3(u("hi", c, 100, "nrf_mcc")), r3(u("bx", c, 100, "nrf_mcc"))) for c in CONDS]
    == [(0.046, 0.034), (0.021, 0.020), (0.062, 0.056), (0.015, 0.014), (0.055, 0.055), (0.064, 0.054)])
chk("N = 100: paired differences alpha = 0.99999 minus BEAST X within 0.012; only S1 and R2 exclude zero",
    [up("hi", "bx", c, 100, "nrf_mcc") for c in CONDS]
    == [(0.012, 0.002, 0.021), (0.001, -0.003, 0.005), (0.006, -0.0, 0.012), (0.002, -0.001, 0.005), (0.0, -0.013, 0.013), (0.011, 0.003, 0.018)]
    and [c for c in CONDS if C.paired(C.vec_all(D, "hi", c, 100, "nrf_mcc"), C.vec_all(D, "bx", c, 100, "nrf_mcc"))[1] > 0] == ["S1", "R2"])
_cp = {c: C.paired(C.vec_all(D, "hi", c, 100, "nrf_cons"), C.vec_all(D, "bx", c, 100, "nrf_cons")) for c in CONDS}
chk("N = 100, consensus tree: differences within 0.004; S4 0.003 (0.0003, 0.005) the only interval excluding zero",
    r3(max(abs(v[0]) for v in _cp.values())) == 0.004 and [c for c in CONDS if _cp[c][1] > 0 or _cp[c][2] < 0] == ["S4"]
    and (r3(_cp["S4"][0]), round(_cp["S4"][1], 4), r3(_cp["S4"][2])) == (0.003, 0.0003, 0.005))
chk("N = 100, P(true): alpha = 0.99999 lower in all six conditions (-0.009 to -0.032), no interval excluding zero",
    [up("hi", "bx", c, 100, "p_true") for c in CONDS]
    == [(-0.021, -0.048, 0.006), (-0.019, -0.050, 0.012), (-0.013, -0.027, 0.0), (-0.032, -0.065, 0.001), (-0.022, -0.065, 0.020), (-0.009, -0.025, 0.006)]
    and all(C.paired(C.vec_all(D, "hi", c, 100, "p_true"), C.vec_all(D, "bx", c, 100, "p_true"))[2] > 0 for c in CONDS))

# ---- alpha = 0.99 on 50 replicates ----------------------------------------------------------------------------
chk("alpha = 0.99 minus BEAST X at N = 100: S3 0.558, R2 0.549, S4 0.048, S1 0.047, R1 0.049, S2 0.011",
    [up("lo", "bx", c, 100, "nrf_mcc")[0] for c in ("S3", "R2", "S4", "S1", "R1", "S2")] == [0.558, 0.549, 0.048, 0.047, 0.049, 0.011])
_ratio = [u("lo", c, n, "nrf_mcc") / u("bx", c, n, "nrf_mcc") for c in ("S3", "R2") for n in NS]
chk("alpha = 0.99 on 50 taxa x 1,000 bp: error 4.7 to 11 times that of BEAST X over N (abstract: at least four times)",
    (round(min(_ratio), 1), round(max(_ratio))) == (4.7, 11))
chk("alpha = 0.99 on S3 and R2: P(true) zero in all 50 replicates at every N",
    max(C.vec_all(D, "lo", c, n, "p_true").max() for c in ("S3", "R2") for n in NS) == 0)

# ---- rate variation among lineages ---------------------------------------------------------------------------
_w = lambda m, a, b, n: float(stats.ttest_ind(C.vec_all(D, m, b, n, "nrf_mcc") - C.vec_all(D, "bx", b, n, "nrf_mcc"),
                                              C.vec_all(D, m, a, n, "nrf_mcc") - C.vec_all(D, "bx", a, n, "nrf_mcc"), equal_var=False).pvalue)
chk("rate variation, alpha = 0.99999 at N = 100: p = 0.14 (R1 vs S1) and 0.34 (R2 vs S3)",
    (round(_w("hi", "S1", "R1", 100), 2), round(_w("hi", "S3", "R2", 100), 2)) == (0.14, 0.34))

# ---- cost on replicates 1-5 -----------------------------------------------------------------------------------
_cr = lambda c, n: m5(c, n, "hi", "cpu_s") / m5(c, n, "bx", "cpu_s")
chk("CPU, alpha = 0.99999 / BEAST X: 69-361 at N = 100, 453-1,510 at N = 2,000, 69-1,510 overall",
    (round(min(_cr(c, 100) for c in CONDS)), round(max(_cr(c, 100) for c in CONDS)), round(min(_cr(c, 2000) for c in CONDS)),
     round(max(_cr(c, 2000) for c in CONDS)), round(min(_cr(c, n) for c in CONDS for n in NS)), round(max(_cr(c, n) for c in CONDS for n in NS)))
    == (69, 361, 453, 1510, 69, 1510))
_wr = lambda c, n: m5(c, n, "hi", "wall_s") / m5(c, n, "bx", "wall_s")
chk("wall-clock, alpha = 0.99999 (32 threads) / BEAST X (8 threads): 4 to 91",
    (round(min(_wr(c, n) for c in CONDS for n in NS)), round(max(_wr(c, n) for c in CONDS for n in NS))) == (4, 91))
_lh = [m5(c, n, "hi", "cpu_s") / m5(c, n, "lo", "cpu_s") for c in CONDS for n in NS]
chk("alpha = 0.99 costs 7 to 28 times less CPU than alpha = 0.99999", (round(min(_lh)), round(max(_lh))) == (7, 28))
_st = lambda m: (min(m5(c, n, m, "steps") for c in CONDS for n in NS), max(m5(c, n, m, "steps") for c in CONDS for n in NS))
chk("annealing steps: about 8,500 to 18,400 at alpha = 0.99999, about 320 to 1,400 at alpha = 0.99",
    (round(_st("hi")[0], -2), round(_st("hi")[1], -2)) == (8500, 18400) and (round(_st("lo")[0], -1), round(_st("lo")[1], -2)) == (320, 1400))
chk("S4, N = 2,000: 1,135,473 CPU s against 1,633; 34,526 s wall-clock against 381",
    (round(m5("S4", 2000, "hi", "cpu_s")), round(m5("S4", 2000, "bx", "cpu_s")), round(m5("S4", 2000, "hi", "wall_s")), round(m5("S4", 2000, "bx", "wall_s")))
    == (1135473, 1633, 34526, 381))

# ---- parallel scaling (Fig. 5) --------------------------------------------------------------------------------
_sc = list(csv.DictReader(open(C.table("scaling_runs.tsv")), delimiter="\t"))
_wl = {(x["method"], x["condition"], int(x["threads"])): float(x["wall_s"]) for x in _sc}
_cp_ = {(x["method"], x["condition"], int(x["threads"])): float(x["cpu_s"]) for x in _sc}
_sp = lambda m, t: _wl[(m, "S4", 1)] / _wl[(m, "S4", t)]
chk("annealed SMC on S4: speed-up 1.9, 7.3, 13.9, 23.9, 32.3 at 2, 8, 16, 32, 64 threads",
    [round(_sp("asmc", t), 1) for t in (2, 8, 16, 32, 64)] == [1.9, 7.3, 13.9, 23.9, 32.3])
chk("BEAST X on S4: speed-up 1.70 at 8 threads, 1.30 at 16; CPU 2.1 and 3.4 times the single-thread run",
    (round(_sp("beastx", 8), 2), round(_sp("beastx", 16), 2)) == (1.70, 1.30)
    and (round(_cp_[("beastx", "S4", 8)] / _cp_[("beastx", "S4", 1)], 1), round(_cp_[("beastx", "S4", 16)] / _cp_[("beastx", "S4", 1)], 1)) == (2.1, 3.4))

# ---- chain-length check ---------------------------------------------------------------------------------------
_cc = list(csv.DictReader(open(C.table("conv_check.csv"))))
_f = lambda k: np.array([float(x[k]) for x in _cc])
_dn = _f("nrf_N2000") - _f("nrf_20M")
chk("20,000,000-state chains: 15 alignments, consensus error identical for 12 and one split apart for 3; means 0.028 and 0.026",
    len(_cc) == 15 and int((abs(_dn) < 1e-9).sum()) == 12 and all(abs(abs(x) - 1 / 94.0) < 2e-4 for x in _dn if abs(x) > 1e-9)
    and (r3(_f("nrf_20M").mean()), r3(_f("nrf_N2000").mean())) == (0.028, 0.026))
chk("20,000,000-state chains: P(true) means 0.186 and 0.189", (r3(_f("p_20M").mean()), r3(_f("p_N2000").mean())) == (0.186, 0.189))

# ---- M1809 ----------------------------------------------------------------------------------------------------
_lad = json.load(open(C.table("emp_converge_EMP_M1809.json")))["rungs"]
_e = json.load(open(C.table("emp_converge_EMPA5_M1809.json")))["rungs"]["1000"]
chk("M1809, alpha = 0.99999, N = 1,000, between programs: consensus RF 4 (0.036), MCC RF 6, r = 0.94, |support difference| 0.063, 48 and 46 splits >= 0.95",
    (_e["between"]["cons_rf"], round(_e["between"]["cons_nrf"], 3), _e["between"]["mcc_rf"], round(_e["between"]["pearson"], 2),
     round(_e["between"]["dsupport"], 3), _e["between"]["resolved95_beastx"], _e["between"]["resolved95_asmc"]) == (4, 0.036, 6, 0.94, 0.063, 48, 46))
chk("M1809 between seeds: consensus RF 2.4 for both programs; support difference 0.044 and 0.039; log Z range 23.5",
    (_e["asmc"]["cons_rf_mean"], _e["beastx"]["cons_rf_mean"], round(_e["asmc"]["dsupport_mean"], 3), round(_e["beastx"]["dsupport_mean"], 3),
     round(_e["asmc"]["logZ_spread"], 1)) == (2.4, 2.4, 0.044, 0.039, 23.5))
chk("M1809, alpha = 0.99999: median CPU 333,456 s per run against 290 for BEAST X",
    round(float(np.median(_e["asmc"]["cpu_s"]))) == 333456 and round(float(np.median(_e["beastx"]["cpu_s"]))) == 290)
_rs = sorted(int(k) for k in _lad)
chk("M1809, alpha = 0.99: log Z ranges 860, 864, 719, 587, 189; consensus RF between programs 15-24",
    [round(_lad[str(n)]["asmc"]["logZ_spread"]) for n in _rs] == [860, 864, 719, 587, 189]
    and (min(_lad[str(n)]["between"]["cons_rf"] for n in _rs), max(_lad[str(n)]["between"]["cons_rf"] for n in _rs)) == (15, 24))
chk("M1809, alpha = 0.99: seed-to-seed support difference above the 0.05 criterion at every rung",
    min(_lad[str(n)]["asmc"]["dsupport_mean"] for n in _rs) > 0.05)

print("\n%d failed" % len(FAILED))
sys.exit(1 if FAILED else 0)
