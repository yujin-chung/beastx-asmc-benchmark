#!/usr/bin/env python3
"""alpha_rerun.py - annealed SMC at the annealing schedule its authors use by default.

Every annealed SMC run in the 2026-11-13 manuscript used -alphaSMCSampler 0.99.  Wang et al. (2020)
use alpha = 1 - 10^(-5) = 0.99999 by default.  This script re-runs replicates 1-5 of all six simulated
conditions at alpha = 0.99999 for every N of the equal-N design (100, 200, 500, 1,000, 2,000), with
32 threads and the seed of the sweep - the layout of wall-clock subset (b), so each run has an
alpha = 0.99 counterpart (runs_asmc/asmc_<ds>_N<N>_t32) and a BEAST X counterpart.

    cd pipeline && source env.sh
    python3 alpha_rerun.py --dry-run                          # list the jobs, start nothing
    nohup python3 alpha_rerun.py > alpha_rerun.log 2>&1 &     # run (resumable: finished runs are skipped)
    tail -3 alpha_rerun.log
    python3 alpha_rerun.py --score-only --fast                # summary now (no MCC scoring)
    python3 alpha_rerun.py --score-only                       # full summary (MCC scores computed and cached)

150 runs.  Run directories are runs_asmc/alpha_<ds>_N<N>_a0.99999/ - the same names alpha_pilot.py
uses, so the pilot's finished runs at alpha = 0.99999 (same seed, same 32 threads) are reused, and
the prefix "alpha_" keeps the runs out of the globs of score_equal_n.py and the empirical scorers.

Before starting anything the script WAITS until alpha_pilot.py and conv_rerun.py have ended, so it
can be launched while they are still running (--no-wait starts at once; do not use it for timings).
Jobs are started cheapest first, --slots at a time (default 5: the existing wall-clock subsets were
timed with typically 4-5 jobs running).  Writes tables/alpha_rerun.csv.
Needs run_asmc.sh with the ASMC_ALPHA override (2026-11-13).
"""
import argparse
import csv
import datetime as dt
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
# mean CPU seconds of the alpha = 0.99 sweep runs (fifty replicates); used only to order jobs and to
# print rough expectations
CPU99 = {
    "S1": {100: 151, 200: 274, 500: 618, 1000: 1157, 2000: 2637},
    "S2": {100: 827, 200: 1579, 500: 4487, 1000: 8677, 2000: 16743},
    "S3": {100: 848, 200: 2162, 500: 5486, 1000: 11628, 2000: 21571},
    "S4": {100: 7984, 200: 15732, 500: 36435, 1000: 72407, 2000: 139916},
    "R1": {100: 131, 200: 242, 500: 595, 1000: 1178, 2000: 2718},
    "R2": {100: 815, 200: 2160, 500: 5187, 1000: 11344, 2000: 21062},
}
# annealing steps relative to alpha = 0.99, counted in the pilot (lines of essTempDiffAdaptive0.csv):
# measured for S1 and S3; S2 and S4 are guesses (S4 extrapolated from its 0.999 and 0.9999 runs)
STEP_RATIO = {
    "0.999":   {"S1": 2.5, "R1": 2.5, "S2": 2.5, "S3": 2.4, "R2": 2.4, "S4": 2.1},
    "0.9999":  {"S1": 7.6, "R1": 7.6, "S2": 7.0, "S3": 5.9, "R2": 5.9, "S4": 4.8},
    "0.99999": {"S1": 24.5, "R1": 24.5, "S2": 21.0, "S3": 16.9, "R2": 16.9, "S4": 15.0},
}
BUSY_PER_RUN = 28.0             # CPU seconds per wall-clock second of a 32-thread run (measured 25-32)


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def tag_of(ds, n, alpha):
    return "alpha_%s_N%d_a%s" % (ds, n, alpha)


def expected_cpu(cond, n, alpha):
    base = CPU99[cond].get(n, CPU99[cond][2000] * n / 2000.0)
    ratio = STEP_RATIO.get(alpha, {}).get(cond)
    if ratio is None:                       # another alpha: steps grow roughly like 1 / sqrt(1 - alpha)
        ratio = (0.01 / (1.0 - float(alpha))) ** 0.5
    return base * ratio


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


def mem_total_gb():
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
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
                j = dict(cond=cond, rep=rep, n=n, alpha=args.alpha)
                j["ds"] = ds_of(cond, rep)
                j["tag"] = tag_of(j["ds"], n, args.alpha)
                j["cost"] = expected_cpu(cond, n, args.alpha)
                # Java heap limits, as in alpha_pilot.py (so reused pilot runs had the same setting)
                j["mem"] = 16 if COND[cond][0] == 20 else (32 if cond != "S4" else 48)
                jobs.append(j)
    return sorted(jobs, key=lambda j: (j["cost"], j["rep"]))


def finished(j):
    return os.path.exists(os.path.join("runs_asmc", j["tag"], "run.json"))


def processes(text, python_only=False):
    """PIDs of other processes whose command line contains the text (python_only: as a script run by python)."""
    import re
    pattern = re.escape(text)
    if python_only:
        pattern = "python[0-9.]* .*" + pattern
    try:
        out = subprocess.run(["pgrep", "-f", pattern], stdout=subprocess.PIPE, universal_newlines=True).stdout
    except Exception:
        return []
    me = str(os.getpid())
    return [p for p in out.split() if p != me]


def parse_foreign(pgrep_out, own_tags):
    """(tag, heap GB) of every run_asmc.sh job in `pgrep -af` output that is not in own_tags."""
    jobs = {}
    for line in pgrep_out.splitlines():
        a = line.split()
        if "run_asmc.sh" not in a:
            continue
        k = a.index("run_asmc.sh")
        if len(a) < k + 5:
            continue
        ds, tag = a[k + 1], a[k + 4]
        if tag in own_tags:
            continue
        jobs[tag] = 48 if ds.startswith("S4") else (16 if ds[:2] in ("S1", "S2", "R1") or ds.startswith("EMP") else 32)
    return sorted(jobs.items())


def foreign_jobs(own_tags):
    """run_asmc.sh jobs on this machine that this driver did not start - left running by a driver that was
    stopped, or started by another script.  They hold threads and memory, so they are counted against
    --slots and the memory budget; this makes it safe to stop a driver and start a new one at any time."""
    try:
        out = subprocess.run(["pgrep", "-af", "run_asmc[.]sh"], stdout=subprocess.PIPE, universal_newlines=True).stdout
    except Exception:
        return []
    return parse_foreign(out, own_tags)


def drivers_running(names):
    return ["%s (pid %s)" % (n, ",".join(processes(n, python_only=True))) for n in names if processes(n, python_only=True)]


def wait_for_others(names, poll_s=120):
    first = True
    while True:
        busy = drivers_running(names)
        if not busy:
            if not first:
                log("the other drivers have ended - starting")
            return
        if first:
            log("waiting for %s to end before starting (checked every %d s; --no-wait to skip)" % (", ".join(busy), poll_s))
            first = False
        time.sleep(poll_s)


def run_all(jobs, args):
    if "ASMC_ALPHA" not in open("run_asmc.sh").read():
        sys.exit("run_asmc.sh has no ASMC_ALPHA override - git pull first")
    if not os.environ.get("ASMC_HOME"):
        sys.exit("ASMC_HOME is not set - run 'source env.sh' first")
    os.makedirs("logs", exist_ok=True)
    if not args.no_wait:
        wait_for_others(args.wait_for, args.wait_poll_s)
    todo = [j for j in jobs if not finished(j)]
    log("memory: the heap limits of the running jobs may add up to %.0f GB (16 GB per 20-taxon job, 32 GB per "
        "50-taxon job, 48 GB per S4 job)" % args.mem_budget_gb)
    log("%d jobs, %d already finished (reused), %d to run; up to %d at a time, %d threads each; alpha = %s"
        % (len(jobs), len(jobs) - len(todo), len(todo), args.slots, THREADS, args.alpha))
    tot = sum(j["cost"] for j in todo)
    log("rough expectation for the remaining jobs: %.0f CPU-hours, about %.1f days at %d slots"
        % (tot / 3600.0, tot / BUSY_PER_RUN / args.slots / 86400.0, args.slots))
    running, ok, fail, t_start, n_foreign = [], 0, 0, time.time(), -1
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
            left = sum(x["cost"] for x in todo) + sum(x["cost"] for x in running)
            if finished(j):
                ok += 1
                log("OK      %s  (%.0f min)  [OK %d  FAIL %d  left %d; about %.1f h to go]"
                    % (j["tag"], (time.time() - j["t0"]) / 60.0, ok, fail, len(todo) + len(running),
                       left / BUSY_PER_RUN / args.slots / 3600.0))
            else:
                fail += 1
                log("%s %s  (see logs/%s.log)  [OK %d  FAIL %d  left %d]"
                    % ("TIMEOUT" if late else "FAILED ", j["tag"], j["tag"], ok, fail, len(todo) + len(running)))
        foreign = foreign_jobs({x["tag"] for x in running})
        if len(foreign) != n_foreign:
            if foreign or n_foreign > 0:
                log("%d job(s) not started by this driver are running and count towards the slots and the memory "
                    "budget: %s" % (len(foreign), ", ".join(t for t, _ in foreign) or "-"))
            n_foreign = len(foreign)
        i = 0
        while i < len(todo) and len(running) + len(foreign) < args.slots:
            j = todo[i]
            if finished(j):                     # finished meanwhile by another driver
                todo.pop(i)
                continue
            if processes(j["tag"]):             # another driver is running this very job: leave it alone
                i += 1
                continue
            # Memory.  MemAvailable only sees what the running JVMs have claimed SO FAR: five S4 jobs started
            # together each pass that test and then each grow towards their 48 GB heap limit (240 GB on a
            # 188 GB machine).  So the heap limits of the running jobs are reserved in full as well
            # (2026-11-14; a JVM with -Xmx16g was observed at 12.4-13.7 GB resident on M1809, N = 1,000).
            reserved = sum(x["mem"] for x in running) + sum(m for _, m in foreign)
            if (running or foreign) and (reserved + j["mem"] > args.mem_budget_gb or mem_available_gb() < j["mem"] + 16):
                break
            todo.pop(i)
            env = dict(os.environ)
            env.update(ASMC_ALPHA=j["alpha"], ASMC_XMX="%dg" % j["mem"])
            j["t0"] = time.time()
            j["proc"] = subprocess.Popen(
                ["bash", "run_asmc.sh", j["ds"], str(j["n"]), str(THREADS), j["tag"], str(SEED)],
                stdout=open(os.path.join("logs", j["tag"] + ".log"), "w"), stderr=subprocess.STDOUT, env=env)
            running.append(j)
            log("start   %s  (guess %.1f CPU-h, %.1f h)" % (j["tag"], j["cost"] / 3600.0, j["cost"] / BUSY_PER_RUN / 3600.0))
        time.sleep(args.poll_s)
    log("Done.  OK %d  FAIL %d  (%.1f h)" % (ok, fail, (time.time() - t_start) / 3600.0))


# ------------------------------------------------------------------ scoring --
def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_json(path):
    try:
        return json.load(open(path))
    except Exception:
        return {}


def n_steps(run_dir):
    """Annealing steps = lines of results/essTempDiffAdaptive0.csv (load-independent cost measure)."""
    try:
        with open(os.path.join(run_dir, "results", "essTempDiffAdaptive0.csv"), errors="replace") as fh:
            return sum(1 for _ in fh)
    except Exception:
        return None


def score_run(j, fast=False):
    """Scores of one finished run.  The MCC score is slow; it is cached in metrics_mcc.json, and with
    fast=True it is only read from the cache, never computed."""
    d = os.path.join("runs_asmc", j["tag"])
    rj = read_json(os.path.join(d, "run.json"))
    if not rj:
        return None
    pm = os.path.join(d, "metrics_mcc.json")
    mm = read_json(pm) if os.path.exists(pm) else {}
    if "nrf_mcc" not in mm and not fast:
        try:
            from mcc_rescore import score_asmc
            mm = score_asmc(d, os.path.join("data", j["ds"], "true.nwk"))
            json.dump(mm, open(pm, "w"), indent=2)
        except Exception as e:
            print("  !! no MCC score for %s (%s: %s)" % (d, type(e).__name__, e), flush=True)
            mm = {}
    return dict(nrf_cons=fnum(rj.get("nrf")), nrf_mcc=fnum(mm.get("nrf_mcc")), p_true=fnum(rj.get("p_true")),
                cpu_s=fnum(rj.get("cpu_s")), wall_s=fnum(rj.get("wall_s")), steps=n_steps(d),
                alpha_recorded=rj.get("alpha"))


def reference():
    """alpha = 0.99 annealed SMC and BEAST X runs of the equal-N sweep, keyed (method, condition, rep, N)."""
    ref = {}
    p = os.path.join("tables", "equal_n_runs.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            ref[(r["method"], r["condition"], int(r["rep"]), int(r["N"]))] = {
                k: fnum(r.get(k)) for k in ("nrf_cons", "nrf_mcc", "p_true", "cpu_s")}
    return ref


def wall_reference(j):
    """Wall-clock subset (b) counterparts of one job: annealed SMC at alpha = 0.99 with 32 threads and
    BEAST X with 8 BEAGLE threads (run.json of those runs; None where a run is missing)."""
    a = read_json(os.path.join("runs_asmc", "asmc_%s_N%d_t32" % (j["ds"], j["n"]), "run.json"))
    b = read_json(os.path.join("runs_beastx", "beastx_%s_c%d_i1000_bt8" % (j["ds"], (500 + j["n"]) * 1000), "run.json"))
    return dict(a099_wall_s=fnum(a.get("wall_s")), a099_cpu32_s=fnum(a.get("cpu_s")),
                a099_steps=n_steps(os.path.join("runs_asmc", "asmc_%s_N%d_t32" % (j["ds"], j["n"]))),
                bx_wall_s=fnum(b.get("wall_s")))


def mean(x):
    x = [v for v in x if v is not None]
    return sum(x) / len(x) if x else float("nan")


def score(jobs, fast=False, workers=16):
    ref = reference()
    slim = [dict(cond=j["cond"], rep=j["rep"], n=j["n"], alpha=j["alpha"], ds=j["ds"], tag=j["tag"])
            for j in jobs if finished(j)]
    print("scoring %d of %d runs%s ..." % (len(slim), len(jobs), " (fast: cached MCC scores only)" if fast else
                                           " (MCC scores computed where missing, %d at a time)" % workers), flush=True)
    if fast or len(slim) < 2:
        scored = [score_run(j, fast) for j in slim]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            scored = list(ex.map(partial(score_run, fast=False), slim))
    rows = []
    for j, s in zip(slim, scored):
        if s is None:
            continue
        row = dict(condition=j["cond"], rep=j["rep"], N=j["n"], alpha=j["alpha"], tag=j["tag"])
        row.update(s)
        for meth, pre in (("asmc", "a099_"), ("beastx", "bx_")):
            r = ref.get((meth, j["cond"], j["rep"], j["n"]), {})
            for k in ("nrf_mcc", "nrf_cons", "p_true", "cpu_s"):
                row[pre + k] = r.get(k)
        row.update(wall_reference(j))
        rows.append(row)
    if not rows:
        print("no finished runs yet")
        return
    os.makedirs("tables", exist_ok=True)
    with open(os.path.join("tables", "alpha_rerun.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    order = list(COND)
    keys = sorted({(r["condition"], r["N"]) for r in rows}, key=lambda k: (order.index(k[0]), k[1]))
    alpha = rows[0]["alpha"]

    def sub(cond, n):
        return [r for r in rows if r["condition"] == cond and r["N"] == n]

    print("\nACCURACY - means over the replicates finished at alpha = %s; reference columns are the same "
          "replicates at alpha = 0.99 and BEAST X" % alpha)
    print("%-4s %5s %4s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s"
          % ("cond", "N", "reps", "nRF MCC", "nRF cons", "P(true)", ".99 MCC", ".99 cons", ".99 P", "BX MCC", "BX cons", "BX P"))
    for cond, n in keys:
        s = sub(cond, n)
        print("%-4s %5d %4d | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f"
              % (cond, n, len(s), mean([r["nrf_mcc"] for r in s]), mean([r["nrf_cons"] for r in s]), mean([r["p_true"] for r in s]),
                 mean([r["a099_nrf_mcc"] for r in s]), mean([r["a099_nrf_cons"] for r in s]), mean([r["a099_p_true"] for r in s]),
                 mean([r["bx_nrf_mcc"] for r in s]), mean([r["bx_nrf_cons"] for r in s]), mean([r["bx_p_true"] for r in s])))
    print("\nCOST - same replicates.  steps = annealing steps; CPU and wall-clock in seconds (32 threads; the "
          ".99 columns are the _t32 runs, BEAST X wall-clock is the 8-thread run, BEAST X CPU the sweep run)")
    print("%-4s %5s %4s | %7s %10s %8s | %7s %10s %8s | %8s %8s"
          % ("cond", "N", "reps", "steps", "CPU s", "wall s", ".99 stp", ".99 CPU", ".99 wall", "BX CPU", "BX wall"))
    for cond, n in keys:
        s = sub(cond, n)
        print("%-4s %5d %4d | %7.0f %10.0f %8.0f | %7.0f %10.0f %8.0f | %8.0f %8.0f"
              % (cond, n, len(s), mean([r["steps"] for r in s]), mean([r["cpu_s"] for r in s]), mean([r["wall_s"] for r in s]),
                 mean([r["a099_steps"] for r in s]), mean([r["a099_cpu32_s"] for r in s]), mean([r["a099_wall_s"] for r in s]),
                 mean([r["bx_cpu_s"] for r in s]), mean([r["bx_wall_s"] for r in s])))
    print("\nper run")
    print("%-4s %3s %5s | %8s %8s %8s %7s %10s %8s" % ("cond", "rep", "N", "nRF MCC", "nRF cons", "P(true)", "steps", "CPU s", "wall s"))
    f = lambda v, fmt: (fmt % v) if v is not None else "-"
    for r in sorted(rows, key=lambda r: (order.index(r["condition"]), r["N"], r["rep"])):
        print("%-4s %3d %5d | %8s %8s %8s %7s %10s %8s" % (r["condition"], r["rep"], r["N"],
              f(r["nrf_mcc"], "%.3f"), f(r["nrf_cons"], "%.3f"), f(r["p_true"], "%.3f"), f(r["steps"], "%d"),
              f(r["cpu_s"], "%.0f"), f(r["wall_s"], "%.0f")))
    print("\nwrote tables/alpha_rerun.csv (%d of %d runs)" % (len(rows), len(jobs)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.99999")
    ap.add_argument("--conds", default="S1,S2,S3,S4,R1,R2")
    ap.add_argument("--reps", default="1-5")
    ap.add_argument("--ns", default="100,200,500,1000,2000")
    ap.add_argument("--slots", type=int, default=5, help="concurrent runs, 32 threads each (default 5)")
    ap.add_argument("--timeout-h", type=float, default=120.0, help="kill a run after this many hours (default 120)")
    ap.add_argument("--mem-budget-gb", type=float, default=None,
                    help="largest sum of Java heap limits (16/32/48 GB per job) allowed to run at once "
                         "(default: total RAM minus 24 GB, i.e. 164 on the R840: 5 x 32 GB or 3 x 48 GB)")
    ap.add_argument("--wait-for", default="alpha_pilot.py,conv_rerun.py",
                    help="driver scripts that must have ended before anything is started")
    ap.add_argument("--no-wait", action="store_true", help="start at once even if those drivers are running")
    ap.add_argument("--poll-s", type=int, default=20)
    ap.add_argument("--wait-poll-s", type=int, default=120)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--score-only", action="store_true")
    ap.add_argument("--fast", action="store_true", help="with --score-only: do not compute missing MCC scores")
    ap.add_argument("--score-workers", type=int, default=16, help="parallel MCC scoring processes (default 16)")
    args = ap.parse_args()
    args.alpha = args.alpha.strip()
    args.conds = [c.strip() for c in args.conds.split(",") if c.strip()]
    for c in args.conds:
        if c not in COND:
            sys.exit("unknown condition %r" % c)
    args.reps = parse_range(args.reps)
    args.ns = [int(x) for x in args.ns.split(",")]
    args.wait_for = [x.strip() for x in args.wait_for.split(",") if x.strip()]
    if args.mem_budget_gb is None:
        args.mem_budget_gb = mem_total_gb() - 24.0
    jobs = build_jobs(args)
    if args.dry_run:
        for j in jobs:
            print("%-44s guess %8.1f CPU-h %6.1f h  %s" % (j["tag"], j["cost"] / 3600.0, j["cost"] / BUSY_PER_RUN / 3600.0,
                                                          "(finished - reused)" if finished(j) else ""))
        todo = [j for j in jobs if not finished(j)]
        tot = sum(j["cost"] for j in todo)
        print("%d jobs, %d finished, %d to run; rough expectation %.0f CPU-hours = about %.1f days at %d slots"
              % (len(jobs), len(jobs) - len(todo), len(todo), tot / 3600.0, tot / BUSY_PER_RUN / args.slots / 86400.0, args.slots))
        busy = drivers_running(args.wait_for)
        if busy:
            print("still running: %s - a real run waits for them to end" % ", ".join(busy))
        return
    if not args.score_only:
        run_all(jobs, args)
    score(jobs, fast=args.fast, workers=args.score_workers)


if __name__ == "__main__":
    main()
