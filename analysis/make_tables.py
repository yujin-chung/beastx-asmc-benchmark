#!/usr/bin/env python3
"""Tables S1-S9 of Online Resource 1, written to output/tables/ as CSV.

    table_S1.csv            N = 100: BEAST X and annealed SMC at alpha = 0.99999.  Accuracy over the 50 replicates;
                            CPU and wall-clock time over replicates 1-5
    table_S2.csv            N = 2,000: BEAST X and annealed SMC at alpha = 0.99.  Accuracy and CPU time over the 50
                            replicates; wall-clock time over replicates 1-5 (pipeline/tables/wall_mean_sd.csv)
    table_S3.csv            annealed SMC at four values of alpha (pipeline/tables/alpha_pilot.csv)
    table_S4a.csv, _S4b     the M1809 alignment: (a) agreement between seeds within each program and cost,
                            (b) agreement between the programs (emp_converge_*.json)
    table_S5a.csv ... _S5d  paired differences in error and support, with 95% t intervals
    table_S6a.csv ... _S6c  error and support across N; alignments with non-zero support for the true topology
    table_S7a.csv ... _S7d  CPU time, wall-clock time, parallel scaling (scaling_runs.tsv), S3 at two settings
    table_S8a.csv, _S8b     strict- against relaxed-clock conditions (Welch tests and intervals)
    table_S9.csv            BEAST X chains of 20,000,000 states (conv_check.csv)

Tables 1 and 2 of the article (the simulation design and the practical guidance) contain no computed values.
Cells of Tables S1 and S2 are "mean (SD)".  Negative values are written with an ASCII minus sign.
"""
import csv
import json
import math
import os

import numpy as np
from scipy import stats

import common as C

COND_NAME = {"S1": "S1 20×1000", "S2": "S2 20×5000", "S3": "S3 50×1000", "S4": "S4 50×5000",
             "R1": "R1 20×1000r", "R2": "R2 50×1000r"}
HEAD = ["Condition", "Method", "nRF MCC", "nRF consensus", "P(true)", "CPU s", "Wall-clock s"]


def small(x):
    return "%.3f (%.3f)" % C.mean_sd(x)


def big(mu, sd):
    return "%s (%s)" % (format(round(mu), ","), format(round(sd), ","))


def write(name, head, rows):
    os.makedirs(C.OUT_TAB, exist_ok=True)
    with open(os.path.join(C.OUT_TAB, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(head)
        w.writerows(rows)
    print("wrote", name, "(%d rows)" % len(rows))


def table_s1(D, D5):
    rows = []
    for c in C.CONDS:
        for m in ("bx", "hi"):
            rows.append([COND_NAME[c], C.NAME[m]] + [small(C.vec_all(D, m, c, 100, q)) for q in ("nrf_mcc", "nrf_cons", "p_true")]
                        + [big(*C.mean_sd(C.vec5(D5, c, 100, m, q))) for q in ("cpu_s", "wall_s")])
    write("table_S1.csv", HEAD, rows)


def table_s2(D):
    wall = {(r["condition"], r["method"]): r for r in csv.DictReader(open(C.table("wall_mean_sd.csv")))
            if r["subset"] == "b" and r["N"] == "2000"}
    rows = []
    for c in C.CONDS:
        for m, key in (("bx", "beastx"), ("lo", "asmc")):
            w = wall[(c, key)]
            rows.append([COND_NAME[c], C.NAME[m]] + [small(C.vec_all(D, m, c, 2000, q)) for q in ("nrf_mcc", "nrf_cons", "p_true")]
                        + [big(*C.mean_sd(C.vec_all(D, m, c, 2000, "cpu_s"))), big(float(w["wall_s_mean"]), float(w["wall_s_sd"]))])
    write("table_S2.csv", HEAD, rows)


def table_s3(D5):
    pil = list(csv.DictReader(open(C.table("alpha_pilot.csv"))))
    mean = lambda rs, q: sum(float(x[q]) for x in rs) / len(rs)
    rows = []
    for cond, n in (("S1", 100), ("S1", 500), ("S3", 100), ("S3", 500), ("S4", 100)):
        sel = lambda a: [x for x in pil if x["condition"] == cond and int(x["N"]) == n and x["alpha"] == a]
        base = sel("0.999")
        reps = sorted(int(x["rep"]) for x in base)
        lo_cpu = sum(D5[(cond, n, r)][("lo", "cpu_s")] for r in reps) / len(reps)      # 32-thread runs, as the others
        rows.append([cond, n, len(reps), "annealed SMC, α = 0.99", "%.3f" % mean(base, "a099_nrf_mcc"),
                     "%.3f" % mean(base, "a099_p_true"), format(round(lo_cpu), ",")])
        for a in ("0.999", "0.9999", "0.99999"):
            rs = sel(a)
            assert sorted(int(x["rep"]) for x in rs) == reps, (cond, n, a)
            rows.append([cond, n, len(reps), "annealed SMC, α = " + a, "%.3f" % mean(rs, "nrf_mcc"), "%.3f" % mean(rs, "p_true"),
                         format(round(mean(rs, "cpu_s")), ",")])
        rows.append([cond, n, len(reps), "BEAST X", "%.3f" % mean(base, "bx_nrf_mcc"), "%.3f" % mean(base, "bx_p_true"),
                     format(round(mean(base, "bx_cpu_s")), ",")])
    write("table_S3.csv", ["Condition", "N", "replicates", "Program and α", "nRF MCC (mean)", "P(true) (mean)", "CPU s (mean)"], rows)


def table_s4():
    lad = json.load(open(C.table("emp_converge_EMP_M1809.json")))["rungs"]               # alpha = 0.99, five rungs
    hi = json.load(open(C.table("emp_converge_EMPA5_M1809.json")))["rungs"]["1000"]       # alpha = 0.99999, N = 1,000

    def within(n, prog, r):
        return [n, r["n_seeds"], prog, "%.1f" % r["cons_rf_mean"], "%.1f" % r["mcc_rf_mean"], "%.3f" % r["dsupport_mean"],
                "" if r.get("logZ_spread") is None else "%.1f" % r["logZ_spread"], format(round(float(np.median(r["cpu_s"]))), ",")]

    def between(label, b):
        return [label, "%.2f" % b["pearson"], b["cons_rf"], b["mcc_rf"], "%.3f" % b["dsupport"],
                "%d / %d" % (b["resolved95_beastx"], b["resolved95_asmc"])]

    a, b = [], []
    for n in sorted(int(k) for k in lad):
        r = lad[str(n)]
        a.append(within(n, C.NAME["bx"], r["beastx"]))
        if n == 1000:
            a.append(within(n, C.NAME["hi"], hi["asmc"]))
            b.append(between("1,000, α = 0.99999", hi["between"]))
        a.append(within(n, C.NAME["lo"], r["asmc"]))
        b.append(between("%s, α = 0.99" % format(n, ","), r["between"]))
    write("table_S4a.csv", ["N", "seeds", "Program", "Consensus RF (mean over seed pairs)", "MCC RF (mean over seed pairs)",
                            "|Δ support| (mean)", "log Z range", "CPU s / run (median)"], a)
    write("table_S4b.csv", ["N and annealed SMC α", "r (split supports)", "Consensus RF", "MCC RF", "|Δ support|",
                            "Splits ≥ 0.95 (BEAST X / annealed SMC)"], b)


def fx(v, k):
    """v with k decimals; a value that rounds to zero is written without a sign"""
    t = "%.*f" % (k, v)
    return t[1:] if t.startswith("-") and float(t) == 0 else t


def ci(t, k):
    return "%s [%s, %s]" % (fx(t[0], k), fx(t[1], k), fx(t[2], k))


def table_s5(D, D5):
    va = lambda m, c, n, q: C.vec_all(D, m, c, n, q)
    write("table_S5a.csv", ["Condition", "Δ MCC nRF [95% CI]", "Δ consensus nRF [95% CI]", "Δ P(true) [95% CI]"],
          [[c] + [ci(C.paired(va("hi", c, 100, q), va("bx", c, 100, q)), 4) for q in ("nrf_mcc", "nrf_cons", "p_true")] for c in C.CONDS])
    rows = []
    for c in C.CONDS:
        for n in (100, 2000):
            p = C.paired(va("lo", c, n, "nrf_mcc"), va("bx", c, n, "nrf_mcc"))
            rows.append([c, n, fx(p[0], 4), "[%s, %s]" % (fx(p[1], 4), fx(p[2], 4))])
    write("table_S5b.csv", ["Condition", "N", "Δ MCC nRF", "95% CI"], rows)
    rows = []
    for c in C.CONDS:
        p = C.paired(va("hi", c, 100, "nrf_mcc"), va("lo", c, 100, "nrf_mcc"))
        rows.append([c, fx(p[0], 4), "[%s, %s]" % (fx(p[1], 4), fx(p[2], 4))])
    write("table_S5c.csv", ["Condition", "Δ MCC nRF", "95% CI"], rows)
    rows = []
    for c in C.CONDS:
        p = C.paired(C.over_n5(D5, c, "hi", "nrf_mcc"), C.over_n5(D5, c, "bx", "nrf_mcc"))
        rows.append([c, fx(p[0], 4), "[%s, %s]" % (fx(p[1], 4), fx(p[2], 4))])
    write("table_S5d.csv", ["Condition", "Δ MCC nRF", "95% CI"], rows)


def table_s6(D, D5):
    m5 = lambda c, n, m, q: float(np.mean(C.vec5(D5, c, n, m, q)))
    mu = lambda m, c, n, q: float(np.mean(C.vec_all(D, m, c, n, q)))
    rg = lambda c, m, q: "%.3f–%.3f" % (min(m5(c, n, m, q) for n in C.NS), max(m5(c, n, m, q) for n in C.NS))
    write("table_S6a.csv", ["Condition", "BEAST X MCC nRF", "SMC, α = 0.99999 MCC nRF", "BEAST X P(true)", "SMC, α = 0.99999 P(true)"],
          [[c, rg(c, "bx", "nrf_mcc"), rg(c, "hi", "nrf_mcc"), rg(c, "bx", "p_true"), rg(c, "hi", "p_true")] for c in C.CONDS])
    write("table_S6b.csv", ["Condition", "N", "BEAST X MCC nRF", "SMC, α = 0.99 MCC nRF", "BEAST X P(true)", "SMC, α = 0.99 P(true)"],
          [[c, n, fx(mu("bx", c, n, "nrf_mcc"), 3), fx(mu("lo", c, n, "nrf_mcc"), 3), fx(mu("bx", c, n, "p_true"), 3), fx(mu("lo", c, n, "p_true"), 3)]
           for c in C.CONDS for n in (100, 200, 2000)])
    write("table_S6c.csv", ["Condition", "BEAST X", "SMC, α = 0.99999", "SMC, α = 0.99"],
          [[c] + ["%d/50" % int((C.vec_all(D, m, c, 100, "p_true") > 0).sum()) for m in ("bx", "hi", "lo")] for c in ("S3", "R2")])


def table_s7(D5):
    m5 = lambda c, n, m, q: float(np.mean(C.vec5(D5, c, n, m, q)))
    com = lambda v: format(int(round(v)), ",")
    head = ["Condition", "N", "BEAST X", "SMC, α = 0.99999", "SMC, α = 0.99"]
    for name, q in (("table_S7a.csv", "cpu_s"), ("table_S7b.csv", "wall_s")):
        write(name, head, [[c, n, com(m5(c, n, "bx", q)), com(m5(c, n, "hi", q)), com(m5(c, n, "lo", q))] for c in C.CONDS for n in C.NS])
    sc = list(csv.DictReader(open(C.table("scaling_runs.tsv")), delimiter="\t"))
    wall = {(x["method"], x["condition"], int(x["threads"])): float(x["wall_s"]) for x in sc}
    cpu = {(x["method"], x["condition"], int(x["threads"])): float(x["cpu_s"]) for x in sc}
    rows = []
    for cond, meth, name in (("S3", "asmc", "SMC, α = 0.99"), ("S4", "asmc", "SMC, α = 0.99"), ("S4", "beastx", "BEAST X")):
        for t in sorted(k[2] for k in wall if k[0] == meth and k[1] == cond):
            sp = wall[(meth, cond, 1)] / wall[(meth, cond, t)]
            rows.append([cond, name, t, "%.2f" % sp, "%.1f" % (100 * sp / t), "%.2f" % (cpu[(meth, cond, t)] / cpu[(meth, cond, 1)])])
    write("table_S7c.csv", ["Condition", "Program", "Threads", "Wall-clock speed-up", "Parallel efficiency (%)", "Relative CPU cost"], rows)
    write("table_S7d.csv", ["α", "N", "Mean MCC nRF", "Mean CPU seconds"],
          [["0.99", 2000, fx(m5("S3", 2000, "lo", "nrf_mcc"), 3), com(m5("S3", 2000, "lo", "cpu_s"))],
           ["0.99999", 100, fx(m5("S3", 100, "hi", "nrf_mcc"), 3), com(m5("S3", 100, "hi", "cpu_s"))]])


def welch(x, y):
    """mean of y minus mean of x for two independent samples, with its 95% Welch interval"""
    d = float(y.mean() - x.mean())
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1))
    h = stats.t.ppf(0.975, df) * math.sqrt(vx + vy)
    return d, d - h, d + h


def table_s8(D):
    va = lambda m, c, n, q: C.vec_all(D, m, c, n, q)
    gap = lambda m, c, n: va(m, c, n, "nrf_mcc") - va("bx", c, n, "nrf_mcc")      # annealed SMC minus BEAST X, per alignment
    rows = []
    for m, a, n in (("hi", "0.99999", 100), ("lo", "0.99", 2000)):
        for taxa, strict, relaxed in ((20, "S1", "R1"), (50, "S3", "R2")):
            p = stats.ttest_ind(gap(m, relaxed, n), gap(m, strict, n), equal_var=False).pvalue
            rows.append([a, n, taxa, fx(float(gap(m, strict, n).mean()), 3), fx(float(gap(m, relaxed, n).mean()), 3), "%.2f" % p])
    write("table_S8a.csv", ["SMC α", "N", "Taxa", "Mean Δ under strict clock", "Mean Δ under relaxed clock", "Welch p-value"], rows)
    rows = []
    for name, m, taxa, strict, relaxed, q, lab, with_ci in (
            ("BEAST X", "bx", 20, "S1", "R1", "nrf_mcc", "MCC", True), ("BEAST X", "bx", 50, "S3", "R2", "nrf_mcc", "MCC", True),
            ("BEAST X", "bx", 20, "S1", "R1", "nrf_cons", "Consensus", True),
            ("SMC, α = 0.99", "lo", 20, "S1", "R1", "nrf_mcc", "MCC", False),         # printed without an interval
            ("SMC, α = 0.99", "lo", 50, "S3", "R2", "nrf_mcc", "MCC", False)):
        x, y = va(m, strict, 2000, q), va(m, relaxed, 2000, q)
        w = welch(x, y)
        rows.append([name, taxa, lab, fx(float(x.mean()), 3), fx(float(y.mean()), 3), ci(w, 3) if with_ci else fx(w[0], 3)])
    write("table_S8b.csv", ["Program", "Taxa", "Tree summary", "Strict-clock mean nRF", "Relaxed-clock mean nRF", "Difference [95% CI]"], rows)


def table_s9():
    cc = sorted(csv.DictReader(open(C.table("conv_check.csv"))), key=lambda x: (x["condition"], int(x["rep"])))
    keys = ("nrf_N2000", "nrf_20M", "p_N2000", "p_20M")
    rows = [[x["condition"], x["rep"]] + ["%.4f" % float(x[k]) for k in keys] for x in cc]
    rows.append(["All", "Mean"] + ["%.4f" % float(np.mean([float(x[k]) for x in cc])) for k in keys])
    write("table_S9.csv", ["Condition", "Replicate", "N = 2,000 consensus nRF", "20M states consensus nRF", "N = 2,000 P(true)", "20M states P(true)"], rows)


if __name__ == "__main__":
    D, D5 = C.load_all(), C.load_five()
    table_s1(D, D5)
    table_s2(D)
    table_s3(D5)
    table_s4()
    table_s5(D, D5)
    table_s6(D, D5)
    table_s7(D5)
    table_s8(D)
    table_s9()
