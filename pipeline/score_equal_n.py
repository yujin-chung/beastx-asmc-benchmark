#!/usr/bin/env python3
"""
score_equal_n.py — tables for the equal-sample-size comparison on simulated
data, plus the wall-clock subsets.

BEAST X equal-N runs (tags beastx_<ds>_c<chain>_i1000[_bt8]) are re-scored
here with the ladder's fixed convention: drop the trees logged at states
0..500,000 (501 trees), use the next N (= chain/1000 - 500) trees, and compute
    nrf_cons   normalized RF, majority-rule consensus vs true tree
    nrf_mcc    normalized RF, MCC topology vs true tree
    p_true     posterior probability of the true topology
with the shared summariser (metrics.py) and the MCC routine (mcc_rescore.py).
Each run's result is cached in runs_beastx/<tag>/metrics_equaln.json.

aSMC runs with N particles already carry the same three quantities
(metrics_asmc.json from the sweep, metrics_mcc.json from mcc_rescore.py).

Outputs (tables/):
    equal_n_runs.tsv        one row per (dataset, method, N) with all metrics
    equal_n_summary.csv     condition x method x N: medians and quartiles of
                            nrf_cons, nrf_mcc, p_true, cpu_s  (50 replicates)
    wall_summary.csv        condition x method x N x setting: median wall-clock
                            and CPU seconds for subset (b) [aSMC t32, BEAST X bt8]
                            and subset (a) [aSMC t8, BEAST X bt8], replicates 1-5

    cd pipeline && source env.sh
    python3 score_equal_n.py
"""
import csv
import glob
import gzip
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import dendropy
import numpy as np

import metrics as M
from mcc_rescore import mcc_from_equal_weight

COND = {"S1": (20, 1000), "S2": (20, 5000), "S3": (50, 1000), "S4": (50, 5000),
        "R1": (20, 1000), "R2": (50, 1000)}
NS = [100, 200, 500, 1000, 2000]
LOG_EVERY = 1000
BURNIN_TREES = 500000 // LOG_EVERY          # trees at states 0..500,000 -> 501 entries


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def run_json(d):
    try:
        return json.load(open(os.path.join(d, "run.json")))
    except Exception:
        return {}


# ------------------------------------------------------- BEAST X rescoring --
def read_trees(run_dir, tns):
    for name in ("beastx_norm.trees", "beastx.trees", "beastx.trees.gz"):
        p = os.path.join(run_dir, name)
        if os.path.exists(p):
            if p.endswith(".gz"):
                with gzip.open(p, "rt") as fh:
                    return dendropy.TreeList.get(data=fh.read(), schema="nexus",
                                                 taxon_namespace=tns, preserve_underscores=True)
            return dendropy.TreeList.get(path=p, schema="nexus", taxon_namespace=tns,
                                         preserve_underscores=True)
    raise FileNotFoundError("no trees in %s" % run_dir)


def score_beastx_run(args):
    run_dir, ds, n = args
    outj = os.path.join(run_dir, "metrics_equaln.json")
    if os.path.exists(outj):
        try:
            return json.load(open(outj))
        except Exception:
            pass
    try:
        tns = dendropy.TaxonNamespace()
        true_tree = M.load_true(os.path.join("data", ds, "true.nwk"), "newick", tns)
        tl = read_trees(run_dir, tns)
        kept = list(tl[BURNIN_TREES + 1: BURNIN_TREES + 1 + n])
        if len(kept) < n:
            return dict(status="SHORT", n_avail=len(kept))
        summ = M.summarise(kept, true_tree, tns)
        mcc = mcc_from_equal_weight(kept, true_tree, tns)
        rj = run_json(run_dir)
        d = dict(status="OK", nrf_cons=summ["nrf"], nrf_mcc=mcc["nrf_mcc"],
                 p_true=summ["p_true"], n_trees=len(kept),
                 cpu_s=rj.get("cpu_s"), wall_s=rj.get("wall_s"))
        json.dump(d, open(outj, "w"), indent=2)
        return d
    except Exception as e:
        return dict(status="FAIL", err="%s: %s" % (type(e).__name__, e))


# -------------------------------------------------------------- aSMC side --
def asmc_metrics(ds, n):
    """Accuracy of the aSMC run with N particles.  Several run directories can
    exist for one dataset and N (the original sweep run plus the wall-clock
    subset runs at 8 and 32 threads).  Prefer the sweep run; otherwise prefer a
    directory that already carries the MCC score (for replicates 1-5 the sweep
    run itself is tagged _t8), and compute the MCC score on the spot if the
    chosen directory has none, so that every replicate contributes to every
    metric."""
    dirs = sorted(glob.glob("runs_asmc/asmc_%s_N%d_t*/" % (ds, n)))
    dirs = [d for d in dirs if os.path.exists(os.path.join(d, "metrics_asmc.json"))]
    if not dirs:
        return None
    sweep = [d for d in dirs if not (d.rstrip("/").endswith("_t8") or d.rstrip("/").endswith("_t32"))]
    with_mcc = [d for d in dirs if os.path.exists(os.path.join(d, "metrics_mcc.json"))]
    d = (sweep or with_mcc or dirs)[0]
    try:
        ma = json.load(open(os.path.join(d, "metrics_asmc.json")))
    except Exception:
        return None
    mm = {}
    pm = os.path.join(d, "metrics_mcc.json")
    if os.path.exists(pm):
        try:
            mm = json.load(open(pm))
        except Exception:
            mm = {}
    if "nrf_mcc" not in mm:
        try:
            from mcc_rescore import score_asmc
            mm = score_asmc(d.rstrip("/"), os.path.join("data", ds, "true.nwk"))
            json.dump(mm, open(pm, "w"), indent=2)
        except Exception as e:
            print("  !! no MCC score for %s (%s: %s)" % (d, type(e).__name__, e), flush=True)
            mm = {}
    rj = run_json(d)
    return dict(nrf_cons=ma.get("nrf"), nrf_mcc=mm.get("nrf_mcc"), p_true=ma.get("p_true"),
                cpu_s=rj.get("cpu_s"), tag=os.path.basename(d.rstrip("/")))


def wall_of(run_dir):
    rj = run_json(run_dir)
    w, c = rj.get("wall_s"), rj.get("cpu_s")
    try:
        return float(w), float(c)
    except (TypeError, ValueError):
        return None, None


def med_iqr(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], dtype=float)
    if v.size == 0:
        return (float("nan"),) * 3, 0
    return (float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75))), int(v.size)


def main():
    if not os.path.isdir("data"):
        sys.exit("run from the bench directory")
    os.makedirs("tables", exist_ok=True)

    # ---- BEAST X equal-N rescoring (parallel) --------------------------------
    jobs = []
    for cond in COND:
        for rep in range(1, 51):
            ds = ds_of(cond, rep)
            for n in NS:
                chain = (500 + n) * LOG_EVERY
                for suffix in ("", "_bt8"):
                    rd = "runs_beastx/beastx_%s_c%d_i%d%s" % (ds, chain, LOG_EVERY, suffix)
                    if os.path.exists(os.path.join(rd, "run.json")):
                        jobs.append((rd, ds, n))
    print("rescoring %d BEAST X equal-N runs ..." % len(jobs), flush=True)
    bx = {}
    with ProcessPoolExecutor(max_workers=min(32, os.cpu_count() or 8)) as ex:
        for (rd, ds, n), res in zip(jobs, ex.map(score_beastx_run, jobs)):
            bx[rd] = res
            if res.get("status") != "OK":
                print("  %s %s %s" % (res.get("status"), rd, res.get("err", res.get("n_avail", ""))))

    # ---- per-run table --------------------------------------------------------
    rows = []
    for cond in COND:
        for rep in range(1, 51):
            ds = ds_of(cond, rep)
            for n in NS:
                chain = (500 + n) * LOG_EVERY
                rd = "runs_beastx/beastx_%s_c%d_i%d" % (ds, chain, LOG_EVERY)
                r = bx.get(rd)
                if r and r.get("status") == "OK":
                    rows.append(dict(condition=cond, rep=rep, dataset=ds, method="beastx", N=n,
                                     nrf_cons=r["nrf_cons"], nrf_mcc=r["nrf_mcc"],
                                     p_true=r["p_true"], cpu_s=r.get("cpu_s")))
                a = asmc_metrics(ds, n)
                if a:
                    rows.append(dict(condition=cond, rep=rep, dataset=ds, method="asmc", N=n,
                                     nrf_cons=a["nrf_cons"], nrf_mcc=a["nrf_mcc"],
                                     p_true=a["p_true"], cpu_s=a.get("cpu_s")))
    with open("tables/equal_n_runs.tsv", "w") as fh:
        cols = ["condition", "rep", "dataset", "method", "N", "nrf_cons", "nrf_mcc", "p_true", "cpu_s"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r[c]) for c in cols) + "\n")

    # ---- summary: condition x method x N -------------------------------------
    print("\n%-5s %-7s %6s %5s %10s %10s %9s %12s" %
          ("cond", "method", "N", "reps", "nRF cons", "nRF MCC", "P(true)", "CPU s"))
    with open("tables/equal_n_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["condition", "method", "N", "n_reps",
                    "nrf_cons_median", "nrf_cons_q1", "nrf_cons_q3",
                    "nrf_mcc_median", "nrf_mcc_q1", "nrf_mcc_q3",
                    "p_true_median", "p_true_q1", "p_true_q3",
                    "cpu_s_median", "cpu_s_q1", "cpu_s_q3"])
        for cond in COND:
            for meth in ("beastx", "asmc"):
                for n in NS:
                    g = [r for r in rows if r["condition"] == cond and r["method"] == meth and r["N"] == n]
                    if not g:
                        continue
                    (c_m, c_1, c_3), k = med_iqr([r["nrf_cons"] for r in g])
                    (m_m, m_1, m_3), k_mcc = med_iqr([r["nrf_mcc"] for r in g])
                    if k_mcc != k:
                        print("  !! %s %s N=%d: MCC score on %d of %d replicates" % (cond, meth, n, k_mcc, k), flush=True)
                    (p_m, p_1, p_3), _ = med_iqr([r["p_true"] for r in g])
                    (u_m, u_1, u_3), _ = med_iqr([float(r["cpu_s"]) if r["cpu_s"] not in (None, "NA") else None for r in g])
                    w.writerow([cond, meth, n, k, round(c_m, 4), round(c_1, 4), round(c_3, 4),
                                round(m_m, 4), round(m_1, 4), round(m_3, 4),
                                round(p_m, 4), round(p_1, 4), round(p_3, 4),
                                round(u_m, 1), round(u_1, 1), round(u_3, 1)])
                    print("%-5s %-7s %6d %5d %10.4f %10.4f %9.4f %12.1f" %
                          (cond, meth, n, k, c_m, m_m, p_m, u_m))

    # ---- wall-clock subsets --------------------------------------------------
    print("\nwall-clock subsets (replicates 1-5)")
    print("%-5s %-7s %6s %-8s %5s %10s %10s" % ("cond", "method", "N", "setting", "reps", "wall s", "CPU s"))
    with open("tables/wall_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["condition", "method", "N", "subset", "threads", "n_reps",
                    "wall_s_median", "wall_s_q1", "wall_s_q3", "cpu_s_median"])
        for cond in COND:
            for n in NS:
                chain = (500 + n) * LOG_EVERY
                arms = [("beastx", "b", "8", ["runs_beastx/beastx_%s_c%d_i%d_bt8" % (ds_of(cond, r), chain, LOG_EVERY) for r in range(1, 6)]),
                        ("asmc", "b", "32", ["runs_asmc/asmc_%s_N%d_t32" % (ds_of(cond, r), n) for r in range(1, 6)]),
                        ("beastx", "a", "8", ["runs_beastx/beastx_%s_c%d_i%d_bt8" % (ds_of(cond, r), chain, LOG_EVERY) for r in range(1, 6)]),
                        ("asmc", "a", "8", ["runs_asmc/asmc_%s_N%d_t8" % (ds_of(cond, r), n) for r in range(1, 6)])]
                for meth, subset, thr, dirs in arms:
                    ws, cs = [], []
                    for d in dirs:
                        if os.path.exists(os.path.join(d, "run.json")):
                            wv, cv = wall_of(d)
                            if wv is not None:
                                ws.append(wv); cs.append(cv)
                    if not ws:
                        continue
                    (wm, w1, w3), k = med_iqr(ws)
                    (cm, _, _), _ = med_iqr(cs)
                    w.writerow([cond, meth, n, subset, thr, k, round(wm, 1), round(w1, 1), round(w3, 1), round(cm, 1)])
                    print("%-5s %-7s %6d %-8s %5d %10.0f %10.0f" % (cond, meth, n, "(%s) t%s" % (subset, thr), k, wm, cm))

    print("\nwrote tables/equal_n_runs.tsv, tables/equal_n_summary.csv, tables/wall_summary.csv")


if __name__ == "__main__":
    main()
