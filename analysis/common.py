"""Loading and summary functions shared by the analysis scripts.

Input tables (written on the server by the scripts in ../pipeline, one row per run):

    pipeline/tables/equal_n_runs.tsv   BEAST X and annealed SMC at alpha = 0.99: 50 replicates of each of the six
                                       conditions at every N (3,000 runs)
    pipeline/tables/alpha_rerun.csv    annealed SMC at alpha = 0.99999: N = 100 on replicates 1-50 and N = 200, 500,
                                       1,000 and 2,000 on replicates 1-5 (420 runs).  Each row also carries, for the
                                       same alignment and N, the BEAST X run and the alpha = 0.99 run, including the
                                       wall-clock runs of replicates 1-5 (annealed SMC at 32 threads, BEAST X at 8
                                       BEAGLE threads)

Methods are keyed "bx" (BEAST X), "hi" (annealed SMC, alpha = 0.99999) and "lo" (annealed SMC, alpha = 0.99).
"""
import csv
import math
import os

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "pipeline", "tables")
OUT_FIG = os.path.join(ROOT, "output", "figures")
OUT_TAB = os.path.join(ROOT, "output", "tables")

CONDS = ["S1", "S2", "S3", "S4", "R1", "R2"]
COND_LABEL = {
    "S1": "S1  20 taxa × 1,000 bp", "S2": "S2  20 taxa × 5,000 bp",
    "S3": "S3  50 taxa × 1,000 bp", "S4": "S4  50 taxa × 5,000 bp",
    "R1": "R1  20 × 1,000, relaxed clock", "R2": "R2  50 × 1,000, relaxed clock",
}
NS = [100, 200, 500, 1000, 2000]
ALL = list(range(1, 51))
FIVE = [1, 2, 3, 4, 5]
METHODS = ["bx", "hi", "lo"]
NAME = {"bx": "BEAST X", "hi": "annealed SMC, α = 0.99999", "lo": "annealed SMC, α = 0.99"}


def table(name):
    return os.path.join(TAB, name)


def mean_sd(x):
    return float(np.mean(x)), float(np.std(x, ddof=1))


def paired(x, y):
    """mean of x - y with its 95% t interval"""
    d = np.asarray(x) - np.asarray(y)
    m, s, k = d.mean(), d.std(ddof=1), len(d)
    h = stats.t.ppf(0.975, k - 1) * s / math.sqrt(k)
    return float(m), float(m - h), float(m + h)


# ---- every replicate each analysis has -------------------------------------------------------------------------
QS_ALL = ["nrf_mcc", "nrf_cons", "p_true", "cpu_s"]


def reps(m, n):
    """the replicates on which method m was run at N = n"""
    return ALL if (m != "hi" or n == 100) else FIVE


def load_all():
    """{(method, condition, N, replicate): {quantity: value}} on all the replicates each analysis has"""
    D = {}
    name = {"beastx": "bx", "asmc": "lo"}
    for r in csv.DictReader(open(table("equal_n_runs.tsv")), delimiter="\t"):
        D[(name[r["method"]], r["condition"], int(r["N"]), int(r["rep"]))] = {q: float(r[q]) for q in QS_ALL}
    for r in csv.DictReader(open(table("alpha_rerun.csv"))):
        D[("hi", r["condition"], int(r["N"]), int(r["rep"]))] = dict(
            {q: float(r[q]) for q in QS_ALL}, wall_s=float(r["wall_s"]), steps=float(r["steps"]))
    for m in METHODS:
        for c in CONDS:
            for n in NS:
                assert all((m, c, n, k) in D for k in reps(m, n)), (m, c, n)
    assert sum(1 for k in D if k[0] == "hi") == 6 * (50 + 4 * 5)
    return D


def vec_all(D, m, c, n, q, rr=None):
    return np.array([D[(m, c, n, k)][q] for k in (rr or reps(m, n))])


# ---- the three analyses on replicates 1-5, with CPU and wall-clock time ---------------------------------------
# column of alpha_rerun.csv per quantity and method ("hi" is the unprefixed column)
COL5 = {
    "nrf_mcc": dict(bx="bx_nrf_mcc", lo="a099_nrf_mcc"),
    "nrf_cons": dict(bx="bx_nrf_cons", lo="a099_nrf_cons"),
    "p_true": dict(bx="bx_p_true", lo="a099_p_true"),
    "cpu_s": dict(bx="bx_cpu_s", lo="a099_cpu32_s"),      # alpha = 0.99: the 32-thread run, as for alpha = 0.99999
    "wall_s": dict(bx="bx_wall_s", lo="a099_wall_s"),
    "steps": dict(lo="a099_steps"),
}
QS5 = ["nrf_mcc", "nrf_cons", "p_true", "cpu_s", "wall_s", "steps"]


def load_five():
    """{(condition, N, replicate): {(method, quantity): value}} on replicates 1-5"""
    out = {}
    for r in csv.DictReader(open(table("alpha_rerun.csv"))):
        if int(r["rep"]) not in FIVE:
            continue
        d = {}
        for q in QS5:
            d[("hi", q)] = float(r[q])
            for m, c in COL5[q].items():
                d[(m, q)] = float(r[c]) if r.get(c) not in ("", None) else float("nan")
        out[(r["condition"], int(r["N"]), int(r["rep"]))] = d
    assert len(out) == 150
    return out


def vec5(D, c, n, m, q):
    return np.array([D[(c, n, k)][(m, q)] for k in FIVE])


def over_n5(D, c, m, q):
    """per-replicate mean over the five N (one value per alignment)"""
    return np.mean([vec5(D, c, n, m, q) for n in NS], axis=0)
