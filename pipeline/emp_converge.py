#!/usr/bin/env python3
"""
emp_converge.py — convergence ladder on the empirical M1809 alignment, with
the two methods matched on the NUMBER OF TREES per run at every rung.

    rung N   aSMC particles   BEAST X chain = (500 + N) x 1,000 states      seeds
     1,000       1,000 (exists)  1,500,000                                  5
     5,000       5,000           5,500,000                                  5
    10,000      10,000          10,500,000                                  5
    20,000      20,000          20,500,000                                  5
    50,000      50,000          50,500,000                                  3
Logged every 1,000 states with the first 500,000 states discarded as
burn-in, every chain leaves exactly N trees, the same number as the aSMC
particles it is compared with.

BEAST X's sampling interval is pinned to 1,000 states for every chain
(BEASTX_LOGEVERY=1000, honoured by run_beastx.sh) and the burn-in is a fixed
500,000 states, so thinning and burn-in are identical across rungs and only
the chain length grows.  BEAST X tags
carry "_i1000" so they are not confused with the earlier 8M run that used
run_beastx.sh's default interval.

Rungs run in order.  After each aSMC rung finishes, emp_rung_score.py
measures seed-to-seed consistency; once the mean absolute difference in
split support across seeds falls below STOP_THRESHOLD (0.05) the ladder
stops and higher rungs (both methods) are skipped.

Resumable: a job whose runs_*/<tag>/run.json exists is skipped.

    cd pipeline && source env.sh
    nohup python3 emp_converge.py > emp_converge.log 2>&1 &
    tail -n 40 emp_converge.log          # (no -f: tail -f never exits by itself)
"""
import datetime as dt
import json
import os
import subprocess
import sys
import time

DS = "EMP_M1809"
CORE_BUDGET = 200
MEM_BUDGET_GB = 180
STOP_THRESHOLD = 0.05
POLL_S = 60
LOG_EVERY = 1000

ASMC_SEEDS = [452, 1453, 2454, 3455, 4456]
BEASTX_SEEDS = [42, 1042, 2042, 3042, 4042]

# (N, chain, n_seeds, asmc threads, asmc heap GB, asmc scheduler GB)
RUNGS = [
    (1000,   1500000, 5,  7,   8,  10),
    (5000,   5500000, 5, 24,  20,  22),
    (10000, 10500000, 5, 32,  36,  38),
    (20000, 20500000, 5, 32,  70,  72),
    (50000, 50500000, 3, 48, 160, 165),
]


class Job:
    def __init__(self, kind, tag, cmd, cores, mem_gb, rung, env=None):
        self.kind, self.tag, self.cmd = kind, tag, cmd
        self.cores, self.mem_gb, self.rung, self.env = cores, mem_gb, rung, env
        self.result = os.path.join("runs_beastx" if kind == "beastx" else "runs_asmc",
                                   tag, "run.json")
        self.proc = None
        self.started = None

    def done(self):
        return os.path.exists(self.result)


def build_jobs():
    jobs = []
    for n, chain, nseeds, thr, xmx, mem in RUNGS:
        for s in BEASTX_SEEDS[:nseeds]:
            tag = "beastx_%s_c%d_i%d_s%d" % (DS, chain, LOG_EVERY, s)
            jobs.append(Job("beastx", tag,
                            ["bash", "run_beastx.sh", DS, str(chain), "auto", tag, str(s)],
                            cores=8, mem_gb=6, rung=n,
                            env={"BEASTX_LOGEVERY": str(LOG_EVERY)}))
        for s in ASMC_SEEDS[:nseeds]:
            tag = "asmc_%s_N%d_s%d" % (DS, n, s)
            jobs.append(Job("asmc", tag,
                            ["bash", "run_asmc.sh", DS, str(n), str(thr), tag, str(s)],
                            cores=thr, mem_gb=mem, rung=n,
                            env={"ASMC_XMX": "%dg" % xmx}))
    return jobs


def now():
    return dt.datetime.now().strftime("%m-%d %H:%M:%S")


def log(msg):
    print("[%s] %s" % (now(), msg), flush=True)


def score_group(method, budget):
    try:
        out = subprocess.run([sys.executable, "emp_rung_score.py", "--dataset", DS,
                              "--group", method, str(budget)],
                             capture_output=True, text=True, timeout=14400)
        js = [l for l in out.stdout.strip().splitlines() if l.startswith("{")]
        if not js:
            log("score %s %d: no JSON\n%s" % (method, budget, out.stderr[-2000:]))
            return None
        return json.loads(js[-1])
    except Exception as e:
        log("score %s %d failed: %s" % (method, budget, e))
        return None


def main():
    if not os.path.isdir("data/%s" % DS):
        sys.exit("data/%s not found - run from the bench directory" % DS)
    if "BEASTX_LOGEVERY" not in open("run_beastx.sh").read():
        sys.exit("run_beastx.sh does not honour BEASTX_LOGEVERY yet - git pull first")
    os.makedirs("logs", exist_ok=True)

    jobs = build_jobs()
    todo = [j for j in jobs if not j.done()]
    log("%d jobs, %d already done, %d to run" % (len(jobs), len(jobs) - len(todo), len(todo)))

    running, cores, mem = [], 0, 0
    stop_at = None            # rung N at which the criterion was met
    scored = set()

    while todo or running:
        for j in list(todo):
            if stop_at is not None and j.rung > stop_at:
                log("skip %s (ladder stopped at N=%d)" % (j.tag, stop_at))
                todo.remove(j); continue
            # aSMC rungs are sequential (the stop test needs the smaller rung first);
            # BEAST X chains are cheap and may run whenever resources allow
            if j.kind == "asmc" and any(x.kind == "asmc" and x.rung < j.rung
                                        for x in running + todo):
                continue
            if cores + j.cores > CORE_BUDGET or mem + j.mem_gb > MEM_BUDGET_GB:
                continue
            env = dict(os.environ); env.update(j.env or {})
            j.proc = subprocess.Popen(j.cmd, stdout=open(os.path.join("logs", j.tag + ".log"), "w"),
                                      stderr=subprocess.STDOUT, env=env)
            j.started = time.time()
            running.append(j); todo.remove(j)
            cores += j.cores; mem += j.mem_gb
            log("start %-44s cores=%-3d mem=%3dGB  (in use %d cores, %d GB)"
                % (j.tag, j.cores, j.mem_gb, cores, mem))

        time.sleep(POLL_S)

        for j in list(running):
            if j.proc.poll() is None:
                continue
            running.remove(j); cores -= j.cores; mem -= j.mem_gb
            ok = j.proc.returncode == 0 and j.done()
            log("%s %-44s wall=%ds" % ("OK   " if ok else "FAIL ", j.tag,
                                       int(time.time() - j.started)))

        # score each aSMC rung once all of its runs are finished
        for n, chain, nseeds, *_ in RUNGS:
            if n in scored:
                continue
            tags = ["asmc_%s_N%d_s%d" % (DS, n, s) for s in ASMC_SEEDS[:nseeds]]
            if any(x.tag in tags for x in running + todo):
                continue
            finished = [t for t in tags if os.path.exists(os.path.join("runs_asmc", t, "run.json"))]
            if len(finished) < 2:
                scored.add(n); continue
            scored.add(n)
            d = score_group("asmc", n)
            if d and d.get("n_seeds", 0) >= 2:
                log("rung N=%d  aSMC seeds=%d  |dsupport|=%.4f  cons RF %.1f  MCC RF %.1f  "
                    "logZ spread %s"
                    % (n, d["n_seeds"], d["dsupport_mean"], d["cons_rf_mean"],
                       d["mcc_rf_mean"], d.get("logZ_spread")))
                if d["dsupport_mean"] < STOP_THRESHOLD and stop_at is None:
                    stop_at = n
                    log("STOP: aSMC reaches |dsupport| < %.2f at N=%d; higher rungs skipped"
                        % (STOP_THRESHOLD, n))

    log("ladder finished.  Now run:  python3 emp_rung_score.py --dataset %s" % DS)


if __name__ == "__main__":
    main()
