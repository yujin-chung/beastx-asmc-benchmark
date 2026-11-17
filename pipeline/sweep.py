#!/usr/bin/env python3
"""sweep.py — the budget-sweep harness.

    cd pipeline && source env.sh
    python3 sweep.py --dry-run                       # plan + time estimate
    python3 sweep.py --until "2026-08-29 22:00"      # run it

WHAT A SWEEP IS, AND WHY IT REPLACED MATCHING
---------------------------------------------
The original design matched the two methods on core-seconds: pick a budget,
tune each method until it costs that much, compare accuracy at that one point.
That cannot be done here.  aSMC's cost is dominated by R, the number of
annealing steps, which the adaptive temperature schedule CHOOSES at run time
from the ESS trajectory - it is an output, not a setting.  Matching would mean
iterating runs until the cost happened to land on target, and would still yield
a single point.

So we sweep instead.  Each method is run at several of its own natural
settings, the ACTUAL cost of every run is measured, and accuracy is plotted
against measured cost.  The comparison is then between two curves rather than
two points, no calibration loop is needed, and the reader can read off whichever
budget they care about.

    BEAST X   chain length in {250k, 500k, 1M, 2M}
    aSMC      particles N in {100, 200, 500, 1000, 2000}

The upper end of the BEAST X ladder is 2M and not higher because the
chain-length check settled it: at 20M states the summaries do not move
(|dnrf| = 0.0000, |dp_true| <= 0.0053 across six datasets).  The interesting
question is where each method BREAKS DOWN as the budget shrinks, not what it
does with more.

COST IS MEASURED, NEVER INFERRED
--------------------------------
Every run records User+System CPU seconds from /usr/bin/time -v.  Wall x threads
would overstate BEAST X, which does not keep its threads busy, and so would
manufacture the paper's conclusion instead of testing it.

PHASES, IN PRIORITY ORDER
-------------------------
  scal    thread ladders, ONE JOB AT A TIME so the timings are clean.  This is
          Figure 4 and nothing else may contend for memory bandwidth while it
          runs, so it goes first, when the time budget is certain.
  beastx  the BEAST X arm of the sweep, heavily parallel.
  asmc    the aSMC arm of the sweep, 8 threads per job.
  extra   spare time: the ladders extended upward.

Each phase checks the clock before starting each job and skips rather than
starts anything that cannot finish, so whatever exists at the deadline is the
most useful subset.  Everything is resumable: a job whose run.json already
exists is skipped, so re-running after a crash continues instead of restarting.
"""

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time

# condition -> (taxa, sites, mode, asmc_cpu_s_at_N200, beastx_cpu_s_at_2M)
# The two cost numbers are measured for S1-S4 and copied by size for R1/R2.
# They are used ONLY to decide whether a job fits before the deadline; the
# recorded cost is always the measured one.
COND = {
    "S1": (20, 1000, "clock",     165,  250),
    "S2": (20, 5000, "clock",     885, 1100),
    "S3": (50, 1000, "clock",     960,  700),
    "S4": (50, 5000, "clock",    5300, 2500),
    "R1": (20, 1000, "nonclock",  165,  250),
    "R2": (50, 1000, "nonclock",  960,  700),
}
SWEEP_CONDS = ["S1", "S2", "S3", "S4", "R1", "R2"]

# Empirical data.  M1809 (= DS7 of Lakner et al. 2008) is the primary: 59 taxa
# puts it in the same regime as S4, which is exactly where the simulated result
# says annealed SMC's posterior stops concentrating, so it is the case that can
# corroborate or contradict that finding on real sequences.  It is also one of
# the two datasets Wang, Wang & Bouchard-Cote (2020) analysed, and we feed aSMC
# the .msf from their own experiment repository, so our numbers are checkable
# against theirs.  Cost estimates scaled from S3 by taxa x sites.
COND["EMP_M1809"] = (59, 1824, "empirical", 2130, 420)
COND["EMP_M336"] = (27, 1949, "empirical", 700, 200)
# Repeated seeds are the point of the empirical design: with no true tree, the
# only thing measurable without knowing the answer is whether a method agrees
# with ITSELF across seeds.
EMP_SEEDS_ASMC = [452, 1453, 2454, 3455, 4456]
EMP_SEEDS_BEASTX = [42, 1042, 2042, 3042, 4042]
EMP_DATASETS = ["EMP_M1809"]     # M336 available; add it here to include it
REPS = 10
ASMC_N = [100, 200, 500, 1000, 2000]
BEASTX_CHAIN = [250000, 500000, 1000000, 2000000]
ASMC_THREADS = 8
BEASTX_CORES = 4          # accounting estimate, not a flag we pass

TABLE = "runs.tsv"
STATUS = "sweep_status.txt"
COLUMNS = ["phase", "tag", "method", "dataset", "condition", "taxa", "sites",
           "mode", "budget", "threads", "status", "wall_s", "cpu_s",
           "nrf", "p_true", "n_units", "finished_at"]


# --------------------------------------------------------------------- job --
class Job:
    def __init__(self, phase, tag, method, dataset, cond, budget, threads,
                 cmd, result, cores, est_s, mem_gb=4, env=None, done_glob=None):
        self.done_glob = done_glob
        self.phase, self.tag, self.method = phase, tag, method
        self.dataset, self.cond, self.budget = dataset, cond, budget
        self.threads, self.cmd, self.result = threads, cmd, result
        self.cores, self.est_s = cores, est_s
        self.mem_gb = mem_gb
        self.env = env or {}
        self.proc = None
        self.started = None
        self.logpath = "logs/%s.log" % tag

    def done(self):
        if os.path.exists(self.result):
            return True
        # Resume must survive a change of thread count.  The tag encodes
        # threads, so a job already run at 8 threads would otherwise look
        # missing once the scheduler starts choosing threads adaptively, and
        # the whole ten-replicate sweep would be recomputed for nothing.
        # Throughput settings are not part of the experimental condition.
        if self.done_glob:
            import glob as _g
            return bool(_g.glob(self.done_glob))
        return False


def asmc_mem_gb(cond, n):
    """Heap a job needs, in GB.

    Memory, not cores, is what actually caps concurrency here.  A 50-taxon
    5000-site run at N = 2000 wants about 16 GB; running enough of those to
    fill 224 cores at 4 threads each would ask for far more than the machine's
    188 GB, so the scheduler has to budget both resources or it will thrash or
    be OOM-killed.  Scaled from the one measured point by particles x taxa x
    sites, floored at 2 GB and capped at 24.
    """
    taxa, sites = COND[cond][0], COND[cond][1]
    scale = (n * taxa * sites) / (2000.0 * 50 * 5000)
    return int(max(2, min(24, round(2 + 14 * scale))))


def auto_threads(mem_gb, cores_per_gb=216 / 150.0):
    """Threads for a job, chosen so cores and memory run out together.

    The machine has 224 cores and 188 GB, so about 1.4 usable cores per usable
    GB.  A 16 GB aSMC job at 8 threads would let memory cap concurrency at nine
    jobs and leave 144 cores idle; the same job at 24 threads fills them.  A
    3 GB job wants the opposite treatment.  Clamped to [4, 32]: below 4 the
    tail job runs for most of a day, and above 32 the measured parallel
    overhead starts to bite (CPU x1.22 at 32 threads, x1.49 at 64).
    """
    return int(max(4, min(32, round(mem_gb * cores_per_gb))))


def asmc_job(phase, cond, rep, n, threads=None, tag=None, ds=None, seed=452):
    ds = ds or "%s_t%d_L%d_rep%d" % (cond, COND[cond][0], COND[cond][1], rep)
    if threads is None:
        threads = (auto_threads(asmc_mem_gb(cond, n)) if ASMC_THREADS == 0
                   else ASMC_THREADS)
    tag = tag or "asmc_%s_N%d_t%d" % (ds, n, threads)
    est = COND[cond][3] * (n / 200.0) / (threads * 0.9)
    mem = asmc_mem_gb(cond, n)
    return Job(phase, tag, "asmc", ds, cond, n, threads,
               ["bash", "run_asmc.sh", ds, str(n), str(threads), tag, str(seed)],
               "runs_asmc/%s/run.json" % tag, threads, est,
               mem_gb=mem, env={"ASMC_XMX": "%dg" % mem},
               done_glob=("runs_asmc/asmc_%s_N%d_t*/run.json" % (ds, n)
                          if phase not in ("scal", "emp", "smoke") else None))


def beastx_job(phase, cond, rep, chain, bthreads="auto", tag=None, ds=None,
               seed=42):
    ds = ds or "%s_t%d_L%d_rep%d" % (cond, COND[cond][0], COND[cond][1], rep)
    tag = tag or "beastx_%s_c%d" % (ds, chain)
    if bthreads != "auto":
        tag = "beastx_%s_c%d_bt%s" % (ds, chain, bthreads)
    est = COND[cond][4] * (chain / 2000000.0) / 1.5
    cores = BEASTX_CORES if bthreads == "auto" else max(1, int(bthreads))
    return Job(phase, tag, "beastx", ds, cond, chain, bthreads,
               ["bash", "run_beastx.sh", ds, str(chain), str(bthreads), tag,
                str(seed)],
               "runs_beastx/%s/run.json" % tag, cores, est, mem_gb=4)


# ------------------------------------------------------------------ phases --
def build_phases(reps):
    """reps: int, or a dict condition -> replicate count.

    Per-condition replicate counts are deliberate, not opportunistic.  The
    primary comparison - cost against accuracy - is settled at ten replicates:
    on S4 the gap between the two methods is nearly three orders of magnitude
    in CPU seconds, which no plausible amount of between-replicate variation
    touches.  The MARGINAL question is the relaxed-clock contrast, where the
    sign of the effect currently disagrees between the 20- and 50-taxon pairs,
    and those four conditions (S1/R1, S3/R2) are also the cheap ones.  So the
    rule, fixed before the runs and reported in Methods, is: fifty replicates
    where the answer is still open and affordable, ten where it is not in
    doubt.  Choosing the counts after seeing which way the result fell would be
    a different and much worse thing.
    """
    if isinstance(reps, int):
        reps = {c: reps for c in SWEEP_CONDS}
    nrep = lambda c: reps.get(c, 0)
    P = {}

    # --- scal: thread ladders, run one at a time -----------------------------
    scal = []
    for cond, ladder in (("S3", [1, 2, 4, 8, 16, 32]),
                         ("S4", [1, 2, 4, 8, 16, 32, 64])):
        for th in ladder:
            scal.append(asmc_job("scal", cond, 1, 200, threads=th,
                                 tag="scal_asmc_%s_t%d" % (cond, th)))
    for bt in ["1", "2", "4", "8", "16"]:
        scal.append(beastx_job("scal", "S4", 1, 2000000, bthreads=bt))
    P["scal"] = scal

    # --- beastx: the full sweep ---------------------------------------------
    bx = []
    for cond in SWEEP_CONDS:
        for rep in range(1, nrep(cond) + 1):
            for ch in BEASTX_CHAIN:
                bx.append(beastx_job("beastx", cond, rep, ch))
    P["beastx"] = bx

    # --- asmc: the full sweep ------------------------------------------------
    am = []
    for cond in SWEEP_CONDS:
        for rep in range(1, nrep(cond) + 1):
            for n in ASMC_N:
                am.append(asmc_job("asmc", cond, rep, n))
    P["asmc"] = am

    # --- extra: spare time, ladders extended upward --------------------------
    ex = []
    for cond in ["S1", "S2", "S3"]:
        for rep in range(1, min(6, nrep(cond) + 1)):
            ex.append(asmc_job("extra", cond, rep, 5000))
    for cond in SWEEP_CONDS:
        for rep in range(1, min(6, nrep(cond) + 1)):
            for ch in [4000000, 8000000]:
                ex.append(beastx_job("extra", cond, rep, ch))
    P["extra"] = ex

    # --- emp: empirical data, no true tree ----------------------------------
    # Two arms.  A budget ladder at one seed, so the cost axis is comparable
    # with the simulated figures; and repeated seeds at one budget, which is
    # the only reliability measurement available without a true tree.
    em = []
    for ds in EMP_DATASETS:
        if not os.path.isdir(os.path.join("data", ds)):
            continue
        for ch in BEASTX_CHAIN + [4000000, 8000000]:
            em.append(beastx_job("emp", ds, 1, ch, ds=ds,
                                 tag="beastx_%s_c%d_s%d" % (ds, ch, EMP_SEEDS_BEASTX[0]),
                                 seed=EMP_SEEDS_BEASTX[0]))
        for sd in EMP_SEEDS_BEASTX[1:]:
            em.append(beastx_job("emp", ds, 1, 2000000, ds=ds,
                                 tag="beastx_%s_c2000000_s%d" % (ds, sd), seed=sd))
        for n in ASMC_N:
            em.append(asmc_job("emp", ds, 1, n, ds=ds,
                               tag="asmc_%s_N%d_s%d" % (ds, n, EMP_SEEDS_ASMC[0]),
                               seed=EMP_SEEDS_ASMC[0]))
        for sd in EMP_SEEDS_ASMC[1:]:
            em.append(asmc_job("emp", ds, 1, 1000, ds=ds,
                               tag="asmc_%s_N1000_s%d" % (ds, sd), seed=sd))
    P["emp"] = em

    # --- sweep: beastx + asmc + extra in ONE pool ---------------------------
    # Running them as separate phases leaves the machine nearly idle for the
    # whole BEAST X phase, which is 1% of the work.  Mixed, and scheduled
    # longest-first under both a core and a memory budget, the small jobs fill
    # the gaps the big ones leave and the tail shortens.
    P["sweep"] = bx + am + ex
    # everything, empirical included - the empirical runs are under 1% of the
    # total work, so they cost no wall time when mixed into the same pool
    P["all"] = bx + am + ex + em

    return P


# --------------------------------------------------------------- recording --
def ensure_table():
    os.makedirs("logs", exist_ok=True)
    if not os.path.exists(TABLE) or os.path.getsize(TABLE) == 0:
        with open(TABLE, "w") as fh:
            fh.write("\t".join(COLUMNS) + "\n")


def record(job, status, wall, extra=None):
    taxa, sites, mode = COND[job.cond][0], COND[job.cond][1], COND[job.cond][2]
    d = extra or {}
    row = [job.phase, job.tag, job.method, job.dataset, job.cond, taxa, sites,
           mode, job.budget, job.threads, status,
           d.get("wall_s", wall), d.get("cpu_s", ""),
           d.get("nrf", ""), d.get("p_true", ""),
           d.get("_n_particles", d.get("n_trees_logged", "")),
           dt.datetime.now().isoformat(timespec="seconds")]
    with open(TABLE, "a") as fh:
        fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def refresh_status(phase, left, running, queued, done, failed):
    with open(STATUS, "w") as fh:
        fh.write("sweep.py   updated %s\n" % dt.datetime.now().isoformat(timespec="seconds"))
        fh.write("phase      %s\n" % phase)
        fh.write("time left  %s\n" % hhmm(left))
        fh.write("running %d   queued %d   done %d   failed %d\n\n"
                 % (running, queued, done, failed))
        try:
            out = subprocess.run(["column", "-t", "-s", "\t", TABLE],
                                 capture_output=True, text=True, timeout=20)
            fh.write(out.stdout[-200000:])
        except Exception:
            fh.write(open(TABLE).read()[-200000:])


def hhmm(s):
    s = max(0, int(s))
    return "%dh%02dm" % (s // 3600, (s % 3600) // 60)


# --------------------------------------------------------------- scheduler --
def run_phase(name, jobs, deadline, max_cores, max_mem_gb=150, exclusive=False,
              reserve=600, poll=5):
    todo = [j for j in jobs if not j.done()]
    skipped = len(jobs) - len(todo)
    # Longest first, so short jobs fill the gaps at the end instead of one long
    # job straggling alone on an otherwise idle machine.  With BEAST X and aSMC
    # in the same pool this matters more than it looks: BEAST X is about 1% of
    # the total work, so as its own phase it leaves the machine nearly idle,
    # whereas mixed in it costs no wall time at all.
    todo.sort(key=lambda j: -j.est_s)

    print("\n" + "=" * 70)
    print("PHASE %s   %d jobs (%d already done)   est %s of core time"
          % (name, len(todo), skipped,
             hhmm(sum(j.est_s * j.cores for j in todo) / max(1, max_cores))))
    print("=" * 70, flush=True)

    running, done, failed, cores_used, mem_used = [], 0, 0, 0, 0
    limit_cores = 1 if exclusive else max_cores

    while todo or running:
        left = deadline - time.time()

        # launch: walk the queue and take the first job that fits BOTH budgets,
        # rather than head-of-line blocking on one big job
        while todo:
            pick = None
            for i, j in enumerate(todo):
                if exclusive:
                    pick = 0 if not running else None
                    break
                if (cores_used + j.cores <= limit_cores
                        and mem_used + j.mem_gb <= max_mem_gb):
                    pick = i
                    break
            if pick is None:
                break
            j = todo[pick]
            if left < j.est_s + reserve:
                todo.pop(pick)
                record(j, "SKIPPED_TIME", "")
                continue
            todo.pop(pick)
            os.makedirs("logs", exist_ok=True)
            env = dict(os.environ); env.update(j.env)
            j.proc = subprocess.Popen(j.cmd, stdout=open(j.logpath, "w"),
                                      stderr=subprocess.STDOUT, env=env)
            j.started = time.time()
            running.append(j)
            cores_used += j.cores
            mem_used += j.mem_gb
            print("  [%s] start %-44s cores=%-3d mem=%-3dG est=%s left=%s"
                  % (dt.datetime.now().strftime("%H:%M:%S"), j.tag, j.cores,
                     j.mem_gb, hhmm(j.est_s), hhmm(left)), flush=True)
            if exclusive:
                break

        if not running:
            break

        time.sleep(poll)

        # reap
        for j in list(running):
            if j.proc.poll() is None:
                continue
            running.remove(j)
            cores_used -= j.cores
            mem_used -= j.mem_gb
            wall = int(time.time() - j.started)
            if j.proc.returncode == 0 and os.path.exists(j.result):
                try:
                    extra = json.load(open(j.result))
                except Exception:
                    extra = {}
                record(j, "OK", wall, extra)
                done += 1
                print("  [%s] OK    %-46s wall=%-7s cpu=%s nrf=%s p_true=%s"
                      % (dt.datetime.now().strftime("%H:%M:%S"), j.tag,
                         "%ds" % wall, extra.get("cpu_s", "?"),
                         extra.get("nrf", "?"), extra.get("p_true", "?")),
                      flush=True)
            else:
                record(j, "FAILED", wall)
                failed += 1
                print("  [%s] FAIL  %-46s rc=%s  see %s"
                      % (dt.datetime.now().strftime("%H:%M:%S"), j.tag,
                         j.proc.returncode, j.logpath), flush=True)
        refresh_status(name, deadline - time.time(), len(running), len(todo),
                       done, failed)

    print("PHASE %s finished: %d ok, %d failed, %d not started"
          % (name, done, failed, len(todo)), flush=True)
    return done, failed


# -------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default=None,
                    help='deadline, e.g. "2026-08-29 22:00"')
    ap.add_argument("--phases", default="scal,beastx,asmc,extra")
    ap.add_argument("--reps", type=int, default=REPS,
                    help="replicates for every condition (default 10)")
    ap.add_argument("--reps-map", default=None,
                    help='per-condition override, e.g. '
                         '"S1=50,R1=50,S3=50,R2=50,S2=10,S4=10"')
    ap.add_argument("--max-cores", type=int, default=216,
                    help="of the machine's 224")
    ap.add_argument("--max-mem-gb", type=int, default=150,
                    help="total JVM heap budget, of the machine's 188 GB")
    ap.add_argument("--asmc-threads", type=int, default=0,
                    help="threads per aSMC job; 0 (default) picks per job so "
                         "that cores and memory run out together")
    ap.add_argument("--drop", default="",
                    help='exclude budgets, e.g. "asmc:S4:1000,asmc:S4:2000" '
                         "to skip the two most expensive points")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="one tiny job of each kind, then stop")
    a = ap.parse_args()

    if not os.path.isdir("data"):
        sys.exit("!! no data/ - run gen_data.sh first")
    if not os.environ.get("ASMC_HOME"):
        sys.exit("!! ASMC_HOME unset - run 'source env.sh' first")

    if a.until:
        deadline = dt.datetime.fromisoformat(a.until).timestamp()
    else:
        deadline = time.time() + 24 * 3600

    ensure_table()
    global ASMC_THREADS
    ASMC_THREADS = a.asmc_threads

    drop = set()
    for item in a.drop.split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) != 3:
            sys.exit('!! --drop entries look like "asmc:S4:2000"')
        drop.add((parts[0], parts[1], float(parts[2])))

    if a.reps_map:
        rmap = {c: a.reps for c in SWEEP_CONDS}
        for item in a.reps_map.split(","):
            k, _, v = item.partition("=")
            k = k.strip()
            if k not in SWEEP_CONDS:
                sys.exit("!! unknown condition in --reps-map: %s" % k)
            rmap[k] = int(v)
        print("replicates: " + "  ".join("%s=%d" % (c, rmap[c]) for c in SWEEP_CONDS))
        phases = build_phases(rmap)
    else:
        phases = build_phases(a.reps)

    if drop:
        for p in list(phases):
            if p == "scal":
                continue
            before = len(phases[p])
            phases[p] = [j for j in phases[p]
                         if (j.method, j.cond, float(j.budget)) not in drop]
            if before != len(phases[p]):
                print("dropped %d jobs from phase %s" % (before - len(phases[p]), p))

    if a.smoke:
        jobs = [beastx_job("smoke", "S1", 1, 250000, tag="smoke_beastx"),
                asmc_job("smoke", "S1", 1, 100, tag="smoke_asmc")]
        for j in jobs:
            if os.path.exists(j.result):
                os.remove(j.result)
        ok, bad = run_phase("smoke", jobs, time.time() + 3600, a.max_cores)
        sys.exit(0 if bad == 0 else 1)

    order = [p.strip() for p in a.phases.split(",") if p.strip()]

    if a.dry_run:
        print("deadline : %s" % dt.datetime.fromtimestamp(deadline))
        print("budget   : %s\n" % hhmm(deadline - time.time()))
        grand = 0.0
        for p in order:
            jobs = [j for j in phases[p] if not j.done()]
            if p == "scal":
                cost = sum(j.est_s for j in jobs)          # serial
            else:
                core_bound = sum(j.est_s * j.cores for j in jobs) / a.max_cores
                mem_bound = sum(j.est_s * j.mem_gb for j in jobs) / a.max_mem_gb
                cost = max(core_bound, mem_bound)
            grand += cost
            print("%-8s %4d jobs   est wall %-8s  (cores %s | memory %s)"
                  % (p, len(jobs), hhmm(cost),
                     hhmm(sum(j.est_s * j.cores for j in jobs) / a.max_cores),
                     hhmm(sum(j.est_s * j.mem_gb for j in jobs) / a.max_mem_gb)))
        print("\ntotal estimated wall: %s   budget: %s"
              % (hhmm(grand), hhmm(deadline - time.time())))
        print("\nnote: estimates are rough and ignore contention; the schedule"
              "\nis deadline-driven, so overruns cost jobs at the tail, not"
              "\ncorrectness of what completed.")
        return

    print("sweep.py start %s" % dt.datetime.now().isoformat(timespec="seconds"))
    print("deadline       %s  (budget %s)"
          % (dt.datetime.fromtimestamp(deadline), hhmm(deadline - time.time())))
    print("max cores      %d of 224      max heap %d GB of 188"
          % (a.max_cores, a.max_mem_gb))
    print("asmc threads   %s" % ("auto (matched to each job's memory)"
                                 if ASMC_THREADS == 0 else str(ASMC_THREADS)))

    for p in order:
        if time.time() > deadline - 600:
            print("out of time - skipping phase %s" % p)
            continue
        run_phase(p, phases[p], deadline, a.max_cores,
                  max_mem_gb=a.max_mem_gb, exclusive=(p == "scal"))

    print("\nDONE %s" % dt.datetime.now().isoformat(timespec="seconds"))
    print("table:  %s" % TABLE)
    print("status: %s" % STATUS)


if __name__ == "__main__":
    main()
