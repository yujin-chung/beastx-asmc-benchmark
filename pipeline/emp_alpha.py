#!/usr/bin/env python3
"""emp_alpha.py - the empirical M1809 alignment with annealed SMC at the settings of its authors.

The convergence ladder in the 2026-11-13 manuscript ran annealed SMC with -alphaSMCSampler 0.99.
Wang et al. (2020) analysed M1809 with K = 1,000 particles and beta = 5, i.e. alpha = 1 - 10^-5 =
0.99999 (and the JC69 model; we keep our model, HKY with kappa fixed at 2 and equal base frequencies,
because the BEAST X runs it is compared with use it).  This script repeats rung N = 1,000 of the
ladder at alpha = 0.99999 with the ladder's five seeds, 32 threads each.

    cd pipeline && source env.sh
    python3 emp_alpha.py --dry-run
    nohup python3 emp_alpha.py > emp_alpha.log 2>&1 &
    tail -3 emp_alpha.log
    python3 emp_alpha.py --status          # annealing steps, CPU, wall-clock, log Z of the finished seeds
    python3 emp_alpha.py --score           # seed-to-seed consistency and comparison with BEAST X

Run directories are runs_asmc/asmc_EMPA5_M1809_N<N>_s<seed>/.  The alias EMPA5_M1809 ("alpha with
beta = 5") is only a name for the run directories - the data are data/EMP_M1809 - and is chosen so
that (i) the existing scorer, which globs runs_*/*EMP_M1809*, never mixes these runs with the
alpha = 0.99 ladder, and (ii) the same scorer can score them unchanged when it is called with
--dataset EMPA5_M1809.  For (ii) --score first makes symbolic links
runs_beastx/beastx_EMPA5_M1809_c<chain>_i1000_s<seed> -> the existing BEAST X runs of the same rung
(BEAST X does not depend on alpha and is not re-run), then calls emp_rung_score.py, which writes
tables/emp_converge_EMPA5_M1809.csv/.json.

Before starting, the script waits until alpha_pilot.py and conv_rerun.py have ended, and refuses to
start while simulated-data jobs of alpha_rerun.py are running (--no-wait skips both checks).  To
queue the simulated-data rerun behind this one, start this script first and then:
    nohup python3 alpha_rerun.py --wait-for emp_alpha.py > alpha_rerun.log 2>&1 &
Higher rungs: --ns 1000,5000 (cost grows steeply; decide after seeing N = 1,000).
"""
import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import time

DS = "EMP_M1809"                # data/EMP_M1809
ALIAS = "EMPA5_M1809"           # name used in the run-directory tags (see the docstring)
ASMC_SEEDS = [452, 1453, 2454, 3455, 4456]      # the seeds of emp_converge.py
BEASTX_SEEDS = [42, 1042, 2042, 3042, 4042]
LOG_EVERY = 1000
# mean CPU seconds per run of the alpha = 0.99 ladder; used only for rough expectations
CPU99 = {1000: 31604, 5000: 96246, 10000: 161995, 20000: 351902, 50000: 912944}
HEAP_GB = {1000: 16, 5000: 24, 10000: 40, 20000: 72, 50000: 160}
STEP_RATIO = 17.0               # annealing steps at alpha = 0.99999 relative to 0.99 (pilot: 12-25)
BUSY_PER_RUN = 28.0             # CPU seconds per wall-clock second of a 32-thread run


def now():
    return dt.datetime.now().strftime("%m-%d %H:%M:%S")


def log(msg):
    print("[%s] %s" % (now(), msg), flush=True)


def tag_of(n, seed):
    return "asmc_%s_N%d_s%d" % (ALIAS, n, seed)


def run_dir(n, seed):
    return os.path.join("runs_asmc", tag_of(n, seed))


def finished(n, seed):
    return os.path.exists(os.path.join(run_dir(n, seed), "run.json"))


def mem_available_gb():
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) / 1048576.0
    except Exception:
        pass
    return 1e9


def processes(text, python_only=False, raw_pattern=None):
    pattern = raw_pattern or re.escape(text)
    if python_only:
        pattern = "python[0-9.]* .*" + pattern
    try:
        out = subprocess.run(["pgrep", "-f", pattern], stdout=subprocess.PIPE, universal_newlines=True).stdout
    except Exception:
        return []
    me = str(os.getpid())
    return [p for p in out.split() if p != me]


def drivers_running(names):
    return ["%s (pid %s)" % (n, ",".join(processes(n, True))) for n in names if processes(n, True)]


def wait_for_others(names, poll_s):
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


def expected_cpu(n):
    base = CPU99.get(n, CPU99[1000] * n / 1000.0)
    return base * STEP_RATIO


def run_all(args):
    if "ASMC_ALPHA" not in open("run_asmc.sh").read():
        sys.exit("run_asmc.sh has no ASMC_ALPHA override - git pull first")
    if not os.environ.get("ASMC_HOME"):
        sys.exit("ASMC_HOME is not set - run 'source env.sh' first")
    if not os.path.exists(os.path.join("data", DS, "aln.fasta")):
        sys.exit("data/%s/aln.fasta not found - run from the bench directory" % DS)
    os.makedirs("logs", exist_ok=True)
    if not args.no_wait:
        wait_for_others(args.wait_for, args.wait_poll_s)
        # alpha_rerun.py may be queued behind this script (--wait-for emp_alpha.py), so it is not waited
        # for; but if its simulated-data jobs are already running, starting now would spoil both timings
        sim = processes("", raw_pattern=r"run_asmc\.sh .* alpha_")
        if sim:
            sys.exit("%d simulated-data jobs of alpha_rerun.py are running.  Either stop it first\n"
                     "    pkill -f alpha_rerun.py; pkill -f 'run_asmc.sh .* alpha_'\n"
                     "(finished runs are kept; restart it later with --wait-for emp_alpha.py), or start this\n"
                     "script with --no-wait to run alongside it (timings of both will be inflated)." % len(sim))
    jobs = [dict(n=n, seed=s, tag=tag_of(n, s), cost=expected_cpu(n), mem=HEAP_GB.get(n, 160))
            for n in args.ns for s in ASMC_SEEDS[:args.seeds]]
    todo = [j for j in jobs if not finished(j["n"], j["seed"])]
    log("%d jobs, %d already finished, %d to run; up to %d at a time, %d threads each; alpha = %s"
        % (len(jobs), len(jobs) - len(todo), len(todo), args.slots, args.threads, args.alpha))
    tot = sum(j["cost"] for j in todo)
    log("rough expectation: %.0f CPU-hours; about %.1f h per N = %d run"
        % (tot / 3600.0, expected_cpu(args.ns[0]) / BUSY_PER_RUN / 3600.0, args.ns[0]))
    running, ok, fail, t_start = [], 0, 0, time.time()
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
            if finished(j["n"], j["seed"]):
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
            env.update(ASMC_ALPHA=args.alpha, ASMC_XMX="%dg" % j["mem"])
            j["t0"] = time.time()
            j["proc"] = subprocess.Popen(
                ["bash", "run_asmc.sh", DS, str(j["n"]), str(args.threads), j["tag"], str(j["seed"])],
                stdout=open(os.path.join("logs", j["tag"] + ".log"), "w"), stderr=subprocess.STDOUT, env=env)
            running.append(j)
            log("start   %s  (guess %.0f CPU-h, %.1f h)" % (j["tag"], j["cost"] / 3600.0, j["cost"] / BUSY_PER_RUN / 3600.0))
        time.sleep(args.poll_s)
    log("Done.  OK %d  FAIL %d  (%.1f h)" % (ok, fail, (time.time() - t_start) / 3600.0))


# ------------------------------------------------------------------ reporting --
def read_json(path):
    try:
        return json.load(open(path))
    except Exception:
        return {}


def n_lines(path):
    try:
        with open(path, errors="replace") as fh:
            return sum(1 for _ in fh)
    except Exception:
        return None


def raw(path, limit=200):
    try:
        return " | ".join(x.strip() for x in open(path, errors="replace").read().strip().splitlines())[:limit]
    except Exception:
        return "-"


def status(args):
    """What the finished runs say without any tree scoring: annealing steps (Wang et al. report 17,639
    for M1809), CPU, wall-clock and the log marginal-likelihood file of each seed, next to the
    alpha = 0.99 run of the same N and seed."""
    for n in args.ns:
        print("\nN = %d, alpha = %s   (reference: the alpha = 0.99 ladder run of the same seed)" % (n, args.alpha))
        print("%-6s %8s %12s %10s | %8s %12s | %s" % ("seed", "steps", "CPU s", "wall s", ".99 stp", ".99 CPU s", "logZout.csv (new run)"))
        for s in ASMC_SEEDS[:args.seeds]:
            d, d0 = run_dir(n, s), os.path.join("runs_asmc", "asmc_%s_N%d_s%d" % (DS, n, s))
            rj, r0 = read_json(os.path.join(d, "run.json")), read_json(os.path.join(d0, "run.json"))
            f = lambda v, fmt="%.0f": (fmt % float(v)) if v is not None else "-"
            state = "" if rj else ("  (running)" if processes(tag_of(n, s)) else "  (not started)")
            print("%-6d %8s %12s %10s | %8s %12s | %s%s"
                  % (s, f(n_lines(os.path.join(d, "results", "essTempDiffAdaptive0.csv")), "%d"), f(rj.get("cpu_s")), f(rj.get("wall_s")),
                     f(n_lines(os.path.join(d0, "results", "essTempDiffAdaptive0.csv")), "%d"), f(r0.get("cpu_s")),
                     raw(os.path.join(d, "results", "logZout.csv")), state))


def link_beastx(args):
    made = 0
    for n in args.ns:
        chain = (500 + n) * 1000
        for s in BEASTX_SEEDS[:args.seeds]:
            src = "beastx_%s_c%d_i%d_s%d" % (DS, chain, LOG_EVERY, s)
            dst = os.path.join("runs_beastx", "beastx_%s_c%d_i%d_s%d" % (ALIAS, chain, LOG_EVERY, s))
            if not os.path.isdir(os.path.join("runs_beastx", src)):
                print("  !! missing BEAST X run runs_beastx/%s" % src)
                continue
            if not os.path.lexists(dst):
                os.symlink(src, dst)            # relative link inside runs_beastx/
                made += 1
    print("BEAST X links for the alias %s: %d new" % (ALIAS, made))


def score(args):
    done = [(n, s) for n in args.ns for s in ASMC_SEEDS[:args.seeds] if finished(n, s)]
    if len(done) < 2:
        print("fewer than two finished runs - nothing to compare yet")
        return
    link_beastx(args)
    cmd = [sys.executable, "emp_rung_score.py", "--dataset", ALIAS]
    print("running: %s\n" % " ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        print("\nemp_rung_score.py exited with status %d.  If it asked for data/%s, make the link\n"
              "    ln -s %s data/%s\nand run 'python3 emp_alpha.py --score' again." % (rc, ALIAS, DS, ALIAS))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.99999")
    ap.add_argument("--ns", default="1000", help="particle counts (default 1000, the K of Wang et al. 2020)")
    ap.add_argument("--seeds", type=int, default=5, help="number of seeds (default 5, the ladder's)")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--slots", type=int, default=5, help="concurrent runs (default 5)")
    ap.add_argument("--timeout-h", type=float, default=240.0)
    ap.add_argument("--wait-for", default="alpha_pilot.py,conv_rerun.py")
    ap.add_argument("--no-wait", action="store_true")
    ap.add_argument("--poll-s", type=int, default=20)
    ap.add_argument("--wait-poll-s", type=int, default=120)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--score", action="store_true")
    args = ap.parse_args()
    args.alpha = args.alpha.strip()
    args.ns = [int(x) for x in args.ns.split(",")]
    args.wait_for = [x.strip() for x in args.wait_for.split(",") if x.strip()]
    if args.dry_run:
        for n in args.ns:
            for s in ASMC_SEEDS[:args.seeds]:
                print("%-34s guess %6.0f CPU-h %5.1f h  %s" % (tag_of(n, s), expected_cpu(n) / 3600.0,
                      expected_cpu(n) / BUSY_PER_RUN / 3600.0, "(finished)" if finished(n, s) else ""))
        busy = drivers_running(args.wait_for)
        if busy:
            print("still running: %s - a real run waits for them to end" % ", ".join(busy))
        return
    if args.status:
        status(args)
        return
    if args.score:
        score(args)
        return
    run_all(args)
    status(args)
    print("\nnext: python3 emp_alpha.py --score")


if __name__ == "__main__":
    main()
