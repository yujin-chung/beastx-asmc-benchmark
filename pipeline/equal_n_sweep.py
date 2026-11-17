#!/usr/bin/env python3
"""
equal_n_sweep.py — simulated-data runs for the EQUAL-SAMPLE-SIZE comparison
and the wall-clock subsets.  Safe to run alongside emp_converge.py: instead of
a fixed budget it looks at the machine's actual free memory and load before
starting each job.

Jobs (in priority order; anything whose run.json exists is skipped):

  1. BEAST X equal-N sweep        6 conditions x 50 replicates x N in
     {100, 200, 500, 1000, 2000}; chain = (500 + N) x 1,000 states, logged
     every 1,000 states (BEASTX_LOGEVERY=1000), so after the fixed 500,000-
     state burn-in each chain holds exactly N trees - the same number as the
     aSMC run with N particles it is compared with.  Default BEAGLE threading.
       tag  beastx_<ds>_c<chain>_i1000

  2. Wall-clock subset (b), main text: each method at its best thread count,
     replicates 1-5 of every condition, same N grid.
       aSMC    32 threads         tag  asmc_<ds>_N<N>_t32
       BEAST X  8 BEAGLE threads  tag  beastx_<ds>_c<chain>_i1000_bt8
     (the BEAST X runs double as the 8-thread arm of subset (a))

  3. Wall-clock subset (a), supplementary: both methods at 8 threads.
       aSMC     8 threads         tag  asmc_<ds>_N<N>_t8   (reused if present)

    cd pipeline && source env.sh
    nohup python3 equal_n_sweep.py > equal_n_sweep.log 2>&1 &
    tail -n 30 equal_n_sweep.log
"""
import datetime as dt
import os
import subprocess
import time

COND = {           # condition -> (taxa, sites)
    "S1": (20, 1000), "S2": (20, 5000), "S3": (50, 1000), "S4": (50, 5000),
    "R1": (20, 1000), "R2": (50, 1000),
}
NS = [100, 200, 500, 1000, 2000]
REPS_ALL = range(1, 51)
REPS_WALL = range(1, 6)
LOG_EVERY = 1000
BEASTX_SEED, ASMC_SEED = 42, 452

MAX_LOAD = 210           # 1-min load average the machine may reach (224 cores)
MEM_HEADROOM_GB = 12     # keep at least this much MemAvailable after launching
MAX_JOBS = 60
POLL_S = 45


def asmc_mem_gb(cond, n):
    taxa, sites = COND[cond]
    scale = (n * taxa * sites) / (2000.0 * 50 * 5000)
    return int(max(2, min(24, round(2 + 14 * scale))))


class Job:
    def __init__(self, kind, tag, cmd, cores, mem_gb, env=None, phase=""):
        self.kind, self.tag, self.cmd, self.cores, self.mem_gb = kind, tag, cmd, cores, mem_gb
        self.env, self.phase = env or {}, phase
        self.result = os.path.join("runs_beastx" if kind == "beastx" else "runs_asmc",
                                   tag, "run.json")
        self.proc = self.started = None

    def done(self):
        return os.path.exists(self.result)


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def build_jobs():
    jobs = []
    # 1. BEAST X equal-N sweep
    for cond in COND:
        for rep in REPS_ALL:
            ds = ds_of(cond, rep)
            for n in NS:
                chain = (500 + n) * LOG_EVERY
                tag = "beastx_%s_c%d_i%d" % (ds, chain, LOG_EVERY)
                jobs.append(Job("beastx", tag,
                                ["bash", "run_beastx.sh", ds, str(chain), "auto", tag, str(BEASTX_SEED)],
                                cores=4, mem_gb=4, env={"BEASTX_LOGEVERY": str(LOG_EVERY)},
                                phase="equalN"))
    # 2. wall subset (b): aSMC 32 threads + BEAST X 8 BEAGLE threads
    for cond in COND:
        for rep in REPS_WALL:
            ds = ds_of(cond, rep)
            for n in NS:
                chain = (500 + n) * LOG_EVERY
                tag = "beastx_%s_c%d_i%d_bt8" % (ds, chain, LOG_EVERY)
                jobs.append(Job("beastx", tag,
                                ["bash", "run_beastx.sh", ds, str(chain), "8", tag, str(BEASTX_SEED)],
                                cores=8, mem_gb=4, env={"BEASTX_LOGEVERY": str(LOG_EVERY)},
                                phase="wall_b"))
                mem = asmc_mem_gb(cond, n)
                tag = "asmc_%s_N%d_t32" % (ds, n)
                jobs.append(Job("asmc", tag,
                                ["bash", "run_asmc.sh", ds, str(n), "32", tag, str(ASMC_SEED)],
                                cores=32, mem_gb=mem + 2, env={"ASMC_XMX": "%dg" % mem},
                                phase="wall_b"))
    # 3. wall subset (a): aSMC 8 threads
    for cond in COND:
        for rep in REPS_WALL:
            ds = ds_of(cond, rep)
            for n in NS:
                mem = asmc_mem_gb(cond, n)
                tag = "asmc_%s_N%d_t8" % (ds, n)
                jobs.append(Job("asmc", tag,
                                ["bash", "run_asmc.sh", ds, str(n), "8", tag, str(ASMC_SEED)],
                                cores=8, mem_gb=mem + 2, env={"ASMC_XMX": "%dg" % mem},
                                phase="wall_a"))
    return jobs


def mem_available_gb():
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 / 1024
    return 0.0


def now():
    return dt.datetime.now().strftime("%m-%d %H:%M:%S")


def log(msg):
    print("[%s] %s" % (now(), msg), flush=True)


def main():
    if not os.path.isdir("data"):
        raise SystemExit("run from the bench directory")
    os.makedirs("logs", exist_ok=True)
    jobs = build_jobs()
    todo = [j for j in jobs if not j.done()]
    by_phase = {}
    for j in todo:
        by_phase[j.phase] = by_phase.get(j.phase, 0) + 1
    log("%d jobs, %d already done, %d to run %s" % (len(jobs), len(jobs) - len(todo),
                                                    len(todo), by_phase))
    running = []
    reserved_mem = 0.0      # memory of jobs started in the last poll, not yet visible
    while todo or running:
        started_now = 0
        for j in list(todo):
            if len(running) >= MAX_JOBS:
                break
            load = os.getloadavg()[0] + sum(x.cores for x in running if time.time() - x.started < 120)
            if load + j.cores > MAX_LOAD:
                break
            if mem_available_gb() - reserved_mem - j.mem_gb < MEM_HEADROOM_GB:
                break
            env = dict(os.environ); env.update(j.env)
            j.proc = subprocess.Popen(j.cmd, stdout=open(os.path.join("logs", j.tag + ".log"), "w"),
                                      stderr=subprocess.STDOUT, env=env)
            j.started = time.time()
            running.append(j); todo.remove(j)
            reserved_mem += j.mem_gb; started_now += 1
            log("start %-48s %-7s cores=%-2d mem=%2dGB  (running %d, todo %d, free %.0f GB)"
                % (j.tag, j.phase, j.cores, j.mem_gb, len(running), len(todo), mem_available_gb()))
            if started_now >= 8:
                break                      # let the machine settle before more
        time.sleep(POLL_S)
        reserved_mem = 0.0
        for j in list(running):
            if j.proc.poll() is None:
                continue
            running.remove(j)
            ok = j.proc.returncode == 0 and j.done()
            log("%s %-48s wall=%ds" % ("OK   " if ok else "FAIL ", j.tag, int(time.time() - j.started)))
    log("equal-N sweep and wall-clock subsets finished.  Next: python3 score_equal_n.py")


if __name__ == "__main__":
    main()
