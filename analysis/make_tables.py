#!/usr/bin/env python3
"""Tables S1-S4 of Online Resource 1, written to output/tables/ as CSV.

    table_S1.csv            N = 100: BEAST X and annealed SMC at alpha = 0.99999.  Accuracy over the 50 replicates;
                            CPU and wall-clock time over replicates 1-5
    table_S2.csv            N = 2,000: BEAST X and annealed SMC at alpha = 0.99.  Accuracy and CPU time over the 50
                            replicates; wall-clock time over replicates 1-5 (pipeline/tables/wall_mean_sd.csv)
    table_S3.csv            annealed SMC at four values of alpha (pipeline/tables/alpha_pilot.csv)
    table_S4a.csv, _S4b     the M1809 alignment: (a) agreement between seeds within each program and cost,
                            (b) agreement between the programs (emp_converge_*.json)

Tables 1 and 2 of the article (the simulation design and the practical guidance) contain no computed values.
Cells are "mean (SD)" unless stated otherwise.
"""
import csv
import json
import os

import numpy as np

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


if __name__ == "__main__":
    D, D5 = C.load_all(), C.load_five()
    table_s1(D, D5)
    table_s2(D)
    table_s3(D5)
    table_s4()
