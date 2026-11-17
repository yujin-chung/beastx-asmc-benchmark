#!/usr/bin/env python3
"""alpha_pilot.py - does annealed SMC do better with the annealing schedule its authors recommend?

All annealed SMC runs so far used -alphaSMCSampler 0.99 (the relative conditional-ESS target that
sets how fine the annealing schedule is).  Wang et al. (2020) define alpha = 1 - 10^(-beta), use
beta = 5 (alpha = 0.99999) by default, 6 for their 50-taxon simulation, and advise "values very
close to one".  A finer schedule means more annealing steps, i.e. more rejuvenation moves per
particle and more CPU.  This pilot re-runs a subset at alpha = 0.999, 0.9999 and 0.99999 and puts
the results next to the existing alpha = 0.99 run and the BEAST X run of the same dataset and N.

    cd pipeline && source env.sh
    python3 alpha_pilot.py --dry-run                       # list the jobs, start nothing
    nohup python3 alpha_pilot.py > alpha_pilot.log 2>&1 &  # run (resumable: finished runs are skipped)
    python3 alpha_pilot.py --score-only --fast             # summary now: consensus nRF, P(true), CPU (no MCC scoring)
    python3 alpha_pilot.py --score-only                    # full summary; MCC scoring runs in parallel and is cached

Default subset (cheapest first): S1 and S3, replicates 1-5, N = 100 and 500; S4, replicates 1-3,
N = 100.  Every run uses 32 threads.  Run directories are runs_asmc/alpha_<ds>_N<N>_a<alpha>/ -
the prefix "alpha_" keeps them out of the globs of score_equal_n.py and the empirical scorers.
Needs run_asmc.sh with the ASMC_ALPHA override (2026-11-13).  Writes tables/alpha_pilot.csv.
"""
import argparse
import csv
import datetime as dt
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial

COND = {"S1": (20, 1000), "S2": (20, 5000), "S3": (50, 1000), "S4": (50, 5000),
        "R1": (20, 1000), "R2": (50, 1000)}
SEED = 452                      # the seed of the equal-N sweep
THREADS = 32
# mean CPU seconds of the alpha = 0.99 sweep runs (N = 100, 500), used only to order and budget jobs
CPU99 = {("S1", 100): 151, ("S1", 500): 618, ("S2", 100): 827, ("S2", 500): 4487,
         ("S3", 100): 848, ("S3", 500): 5486, ("S4", 100): 7984, ("S4", 500): 36435,
         ("R1", 100): 131, ("R1", 500): 595, ("R2", 100): 815, ("R2", 500): 5187}


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def tag_of(ds, n, alpha):
    return "alpha_%s_N%d_a%s" % (ds, n, alpha)


def expected_cpu(cond, n, alpha):
    """Rough guess: the number of annealing steps grows like 1 / sqrt(1 - alpha)."""
    base = CPU99.get((cond, n), CPU99.get((cond, 500), 5000) * n / 500.0)
    return base * (0.01 / (1.0 - float(alpha))) ** 0.5


def now():
    return dt.datetime.now().strftime("%m-%d %H:%M:%S")


def log(msg):
    print("[%s] %s" % (now(), msg), flush=True)


def mem_available_gb():
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) / 1048576.0
    except Exception:
        pass
    return 1e9


def parse_range(s):
    a, _, b = s.partition("-")
    return list(range(int(a), int(b or a) + 1))


def build_jobs(args):
    jobs = []
    for cond in args.conds:
        for rep in args.reps:
            for n in args.ns:
                for alpha in args.alphas:
                    jobs.append(dict(cond=cond, rep=rep, n=n, alpha=alpha))
    if not args.no_s4:
        for rep in [r for r in args.reps if r <= 3]:
            for alpha in args.alphas:
                jobs.append(dict(cond="S4", rep=rep, n=100, alpha=alpha))
    seen, out = set(), []
    for j in jobs:
        j["ds"] = ds_of(j["cond"], j["rep"])
        j["tag"] = tag_of(j["ds"], j["n"], j["alpha"])
        j["cost"] = expected_cpu(j["cond"], j["n"], j["alpha"])
        j["mem"] = 16 if COND[j["cond"]][0] == 20 else (32 if j["cond"] != "S4" else 48)
        if j["tag"] not in seen:
            seen.add(j["tag"])
            out.append(j)
    return sorted(out, key=lambda j: j["cost"])


def finished(j):
    return os.path.exists(os.path.join("runs_asmc", j["tag"], "run.json"))


def run_all(jobs, args):
    if "ASMC_ALPHA" not in open("run_asmc.sh").read():
        sys.exit("run_asmc.sh has no ASMC_ALPHA override - git pull first")
    if not os.environ.get("ASMC_HOME"):
        sys.exit("ASMC_HOME is not set - run 'source env.sh' first")
    os.makedirs("logs", exist_ok=True)
    todo = [j for j in jobs if not finished(j)]
    log("%d jobs, %d already finished, %d to run; up to %d at a time (%d threads each)"
        % (len(jobs), len(jobs) - len(todo), len(todo), args.slots, THREADS))
    log("expected CPU for the remaining jobs (rough guess): %.1f CPU-days"
        % (sum(j["cost"] for j in todo) / 86400.0))
    running, ok, fail = [], 0, 0
    while todo or running:
        for j in list(running):
            rc = j["proc"].poll()
            late = time.time() - j["t0"] > args.timeout_h * 3600
            if rc is None and not late:
                continue
            if rc is None:
                j["proc"].kill()
                j["proc"].wait()
            running.remove(j)
            if finished(j):
                ok += 1
                log("OK      %s  (%.0f min)  [OK %d  FAIL %d  left %d]"
                    % (j["tag"], (time.time() - j["t0"]) / 60.0, ok, fail, len(todo) + len(running)))
            else:
                fail += 1
                log("%s %s  (see logs/%s.log)  [OK %d  FAIL %d  left %d]"
                    % ("TIMEOUT" if late else "FAILED ", j["tag"], j["tag"], ok, fail, len(todo) + len(running)))
        while todo and len(running) < args.slots:
            j = todo[0]
            if running and mem_available_gb() < j["mem"] + 16:
                break
            todo.pop(0)
            env = dict(os.environ)
            env.update(ASMC_ALPHA=j["alpha"], ASMC_XMX="%dg" % j["mem"])
            j["t0"] = time.time()
            j["proc"] = subprocess.Popen(
                ["bash", "run_asmc.sh", j["ds"], str(j["n"]), str(THREADS), j["tag"], str(SEED)],
                stdout=open(os.path.join("logs", j["tag"] + ".log"), "w"), stderr=subprocess.STDOUT, env=env)
            running.append(j)
            log("start   %s  (guess %.1f CPU-h)" % (j["tag"], j["cost"] / 3600.0))
        time.sleep(20)
    log("Done.  OK %d  FAIL %d" % (ok, fail))


# ------------------------------------------------------------------ scoring --
def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def score_run(j, fast=False):
    """Scores of one finished run.  The MCC score is slow (seconds to minutes per run); it is cached in
    metrics_mcc.json, and with fast=True it is only read from the cache, never computed."""
    d = os.path.join("runs_asmc", j["tag"])
    try:
        rj = json.load(open(os.path.join(d, "run.json")))
    except Exception:
        return None
    pm = os.path.join(d, "metrics_mcc.json")
    mm = {}
    if os.path.exists(pm):
        try:
            mm = json.load(open(pm))
        except Exception:
            mm = {}
    if "nrf_mcc" not in mm and not fast:
        try:
            from mcc_rescore import score_asmc
            mm = score_asmc(d, os.path.join("data", j["ds"], "true.nwk"))
            json.dump(mm, open(pm, "w"), indent=2)
        except Exception as e:
            print("  !! no MCC score for %s (%s: %s)" % (d, type(e).__name__, e), flush=True)
            mm = {}
    return dict(nrf_cons=fnum(rj.get("nrf")), nrf_mcc=fnum(mm.get("nrf_mcc")), p_true=fnum(rj.get("p_true")),
                cpu_s=fnum(rj.get("cpu_s")), wall_s=fnum(rj.get("wall_s")), alpha_recorded=rj.get("alpha"))


def reference():
    """alpha = 0.99 annealed SMC and BEAST X runs of the equal-N sweep, keyed (method, condition, rep, N)."""
    ref = {}
    p = os.path.join("tables", "equal_n_runs.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            ref[(r["method"], r["condition"], int(r["rep"]), int(r["N"]))] = {
                k: fnum(r.get(k)) for k in ("nrf_cons", "nrf_mcc", "p_true", "cpu_s")}
    return ref


def mean(x):
    x = [v for v in x if v is not None]
    return sum(x) / len(x) if x else float("nan")


def score(jobs, fast=False, workers=16):
    ref = reference()
    rows = []
    slim = [dict(cond=j["cond"], rep=j["rep"], n=j["n"], alpha=j["alpha"], ds=j["ds"], tag=j["tag"]) for j in jobs if finished(j)]
    print("scoring %d finished runs%s ..." % (len(slim), " (fast: cached MCC scores only)" if fast else
                                             " (MCC scores computed where missing, %d at a time)" % workers), flush=True)
    if fast or len(slim) < 2:
        scored = [score_run(j, fast) for j in slim]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            scored = list(ex.map(partial(score_run, fast=False), slim))
    for j, s in zip(slim, scored):
        if s is None:
            continue
        row = dict(condition=j["cond"], rep=j["rep"], N=j["n"], alpha=j["alpha"], tag=j["tag"])
        row.update(s)
        for meth, pre in (("asmc", "a099_"), ("beastx", "bx_")):
            r = ref.get((meth, j["cond"], j["rep"], j["n"]), {})
            for k in ("nrf_mcc", "p_true", "cpu_s"):
                row[pre + k] = r.get(k)
        rows.append(row)
    if not rows:
        print("no finished runs yet")
        return
    os.makedirs("tables", exist_ok=True)
    with open(os.path.join("tables", "alpha_pilot.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nper condition, N and alpha: means over the replicates finished at that alpha "
          "(reference columns over the same replicates)")
    print("%-4s %5s %-8s %4s | %8s %8s %10s %8s | %8s %8s %10s | %8s %8s %10s"
          % ("cond", "N", "alpha", "reps", "nRF MCC", "P(true)", "CPU s", "wall s",
             "a=.99 nRF", "a=.99 P", "a=.99 CPU", "BX nRF", "BX P", "BX CPU"))
    keys = sorted({(r["condition"], r["N"], float(r["alpha"]), r["alpha"]) for r in rows})
    for cond, n, _, alpha in keys:
        s = [r for r in rows if r["condition"] == cond and r["N"] == n and r["alpha"] == alpha]
        print("%-4s %5d %-8s %4d | %8.3f %8.3f %10.0f %8.0f | %8.3f %8.3f %10.0f | %8.3f %8.3f %10.0f"
              % (cond, n, alpha, len(s), mean([r["nrf_mcc"] for r in s]), mean([r["p_true"] for r in s]),
                 mean([r["cpu_s"] for r in s]), mean([r["wall_s"] for r in s]),
                 mean([r["a099_nrf_mcc"] for r in s]), mean([r["a099_p_true"] for r in s]), mean([r["a099_cpu_s"] for r in s]),
                 mean([r["bx_nrf_mcc"] for r in s]), mean([r["bx_p_true"] for r in s]), mean([r["bx_cpu_s"] for r in s])))
    print("\nper run")
    print("%-4s %3s %5s %-8s | %8s %8s %8s %10s %8s" % ("cond", "rep", "N", "alpha", "nRF MCC", "nRF cons", "P(true)", "CPU s", "wall s"))
    for r in sorted(rows, key=lambda r: (r["condition"], r["N"], r["rep"], float(r["alpha"]))):
        f = lambda v, fmt: (fmt % v) if v is not None else "-"
        print("%-4s %3d %5d %-8s | %8s %8s %8s %10s %8s" % (r["condition"], r["rep"], r["N"], r["alpha"],
              f(r["nrf_mcc"], "%.3f"), f(r["nrf_cons"], "%.3f"), f(r["p_true"], "%.3f"), f(r["cpu_s"], "%.0f"), f(r["wall_s"], "%.0f")))
    # what the sampler leaves behind, so the number of annealing steps can be read off next
    d = os.path.join("runs_asmc", rows[0]["tag"])
    print("\nfiles kept for %s:" % rows[0]["tag"])
    for p in sorted(glob.glob(os.path.join(d, "results", "*")))[:40]:
        print("   %-50s %10d bytes" % (os.path.basename(p), os.path.getsize(p)))
    out = os.path.join(d, "asmc.out")
    if os.path.exists(out):
        lines = open(out, errors="replace").read().splitlines()
        print("asmc.out: %d lines; last 12:" % len(lines))
        for ln in lines[-12:]:
            print("   | " + ln[:160])
    print("\nwrote tables/alpha_pilot.csv (%d runs)" % len(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphas", default="0.999,0.9999,0.99999")
    ap.add_argument("--conds", default="S1,S3")
    ap.add_argument("--reps", default="1-5")
    ap.add_argument("--ns", default="100,500")
    ap.add_argument("--no-s4", action="store_true", help="skip the S4 (50 taxa x 5,000 bp) runs at N = 100, replicates 1-3")
    ap.add_argument("--slots", type=int, default=6, help="concurrent runs, 32 threads each (default 6 = 192 cores)")
    ap.add_argument("--timeout-h", type=float, default=72.0, help="kill a run after this many hours (default 72)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--score-only", action="store_true")
    ap.add_argument("--fast", action="store_true", help="with --score-only: do not compute missing MCC scores")
    ap.add_argument("--score-workers", type=int, default=16, help="parallel MCC scoring processes (default 16)")
    args = ap.parse_args()
    args.alphas = [a.strip() for a in args.alphas.split(",") if a.strip()]
    args.conds = [c.strip() for c in args.conds.split(",") if c.strip()]
    args.reps = parse_range(args.reps)
    args.ns = [int(x) for x in args.ns.split(",")]
    jobs = build_jobs(args)
    if args.dry_run:
        for j in jobs:
            print("%-44s guess %7.1f CPU-h  %s" % (j["tag"], j["cost"] / 3600.0, "(finished)" if finished(j) else ""))
        print("%d jobs; rough total %.1f CPU-days (the guess assumes steps ~ 1/sqrt(1 - alpha); the first runs will show)"
              % (len(jobs), sum(j["cost"] for j in jobs) / 86400.0))
        return
    if not args.score_only:
        run_all(jobs, args)
    score(jobs, fast=args.fast, workers=args.score_workers)


if __name__ == "__main__":
    main()
