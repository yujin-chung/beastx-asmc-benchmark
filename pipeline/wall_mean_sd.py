#!/usr/bin/env python3
"""Mean and SD of the wall-clock subsets (replicates 1-5), for the mean/SD
companion file.  score_equal_n.py reports medians and quartiles of the same
runs in tables/wall_summary.csv.

Run in the bench directory on the server:
    python3 wall_mean_sd.py
Writes tables/wall_mean_sd.csv and prints it.  Reads only run.json (wall_s,
cpu_s); nothing is re-run or re-scored.
"""
import csv
import json
import math
import os
import sys

COND = {"S1": (20, 1000), "S2": (20, 5000), "S3": (50, 1000), "S4": (50, 5000),
        "R1": (20, 1000), "R2": (50, 1000)}
NS = [100, 200, 500, 1000, 2000]
LOG_EVERY = 1000


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def wall_of(d):
    try:
        rj = json.load(open(os.path.join(d, "run.json")))
        return float(rj["wall_s"]), float(rj["cpu_s"])
    except Exception:
        return None


def mean_sd(x):
    m = sum(x) / len(x)
    sd = math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1)) if len(x) > 1 else float("nan")
    return m, sd


def main():
    out = []
    for cond in COND:
        for n in NS:
            chain = (500 + n) * LOG_EVERY
            bx = ["runs_beastx/beastx_%s_c%d_i%d_bt8" % (ds_of(cond, r), chain, LOG_EVERY) for r in range(1, 6)]
            arms = [("beastx", "b", "8", bx),
                    ("asmc", "b", "32", ["runs_asmc/asmc_%s_N%d_t32" % (ds_of(cond, r), n) for r in range(1, 6)]),
                    ("beastx", "a", "8", bx),
                    ("asmc", "a", "8", ["runs_asmc/asmc_%s_N%d_t8" % (ds_of(cond, r), n) for r in range(1, 6)])]
            for meth, subset, thr, dirs in arms:
                v = [w for w in (wall_of(d) for d in dirs) if w is not None]
                if not v:
                    continue
                wm, ws = mean_sd([a for a, _ in v])
                cm, cs = mean_sd([b for _, b in v])
                out.append([cond, meth, n, subset, thr, len(v), round(wm, 1), round(ws, 1), round(cm, 1), round(cs, 1)])
    hdr = ["condition", "method", "N", "subset", "threads", "n_reps",
           "wall_s_mean", "wall_s_sd", "cpu_s_mean", "cpu_s_sd"]
    os.makedirs("tables", exist_ok=True)
    with open("tables/wall_mean_sd.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerows(out)
    w = csv.writer(sys.stdout)
    w.writerow(hdr)
    w.writerows(out)
    short = [r for r in out if r[5] < 5]
    sys.stderr.write("wrote tables/wall_mean_sd.csv (%d rows; %d with fewer than 5 replicates)\n" % (len(out), len(short)))


if __name__ == "__main__":
    main()
