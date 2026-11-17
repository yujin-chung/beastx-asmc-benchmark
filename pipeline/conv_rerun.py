#!/usr/bin/env python3
"""conv_rerun.py - BEAST X chain-length check on the CURRENT alignments.

The 20,000,000-state runs in runs_long/ were made before the simulated data were regenerated, so they
cannot validate the current design (conv_check.py showed them disagreeing with every current run).
This script repeats the check on the current data: for each chosen dataset it runs one chain of
20,000,000 states, logged every 10,000, with the first 25% (5,000,000 states) discarded - an
independent chain (seed 43; the sweep used 42), ten times the length of the longest equal-N chain -
and compares the majority-rule consensus nRF and P(true) with the equal-N run of the same dataset at
N = 2,000 (2,500,000 states, 500,000-state burn-in) and at N = 100.

    cd pipeline && source env.sh
    python3 conv_rerun.py --dry-run
    nohup python3 conv_rerun.py > conv_rerun.log 2>&1 &
    python3 conv_rerun.py --score-only

Default: S3, S4, R2 (the 50-taxon conditions), replicates 1-5 -> 15 runs, roughly 20 CPU-hours in
all.  Run directories are runs_beastx/long20M_<ds>/ (not matched by the globs of score_equal_n.py).
Writes tables/conv_check.csv (replacing the table made from the stale runs_long/).
"""
import argparse
import csv
import datetime as dt
import glob
import json
import os
import re
import subprocess
import sys
import time

COND = {"S1": (20, 1000), "S2": (20, 5000), "S3": (50, 1000), "S4": (50, 5000),
        "R1": (20, 1000), "R2": (50, 1000)}
CHAIN, LOGEVERY, SEED = 20000000, 10000, 43


def ds_of(cond, rep):
    taxa, sites = COND[cond]
    return "%s_t%d_L%d_rep%d" % (cond, taxa, sites, rep)


def tag_of(ds):
    return "long20M_%s" % ds


def log(msg):
    print("[%s] %s" % (dt.datetime.now().strftime("%m-%d %H:%M:%S"), msg), flush=True)


def finished(ds):
    return os.path.exists(os.path.join("runs_beastx", tag_of(ds), "run.json"))


def first_sequence(ds):
    txt = open(os.path.join("data", ds, "aln.fasta")).read()
    return "".join(txt.split(">")[1].split("\n", 1)[1].split()).upper()


def xml_has_alignment(xml_path, ds):
    xml = re.sub(r"\s", "", open(xml_path, errors="replace").read()).upper()
    return first_sequence(ds)[:200] in xml


def alignment_control(dsets):
    """The test that called runs_long/ stale needs a positive control: the XML of a CURRENT run
    must contain the current alignment."""
    ok = tot = 0
    for ds in dsets:
        for x in sorted(glob.glob(os.path.join("runs_beastx", "beastx_%s_c2500000_i1000*" % ds, "beastx.xml")))[:1]:
            tot += 1
            ok += xml_has_alignment(x, ds)
    stale = [ds for ds in dsets if os.path.exists(os.path.join("runs_long", ds, "beastx.xml"))
             and not xml_has_alignment(os.path.join("runs_long", ds, "beastx.xml"), ds)]
    n_long = sum(os.path.exists(os.path.join("runs_long", ds, "beastx.xml")) for ds in dsets)
    print("alignment control: %d of %d current equal-N runs contain the current alignment "
          "(must be all of them for the test to mean anything)" % (ok, tot))
    print("runs_long/: %d of %d old 20M runs were made on a different alignment" % (len(stale), n_long))


def run_all(dsets, args):
    if not os.environ.get("BEAST_HOME"):
        sys.exit("BEAST_HOME is not set - run 'source env.sh' first")
    os.makedirs("logs", exist_ok=True)
    todo = [ds for ds in dsets if not finished(ds)]
    log("%d datasets, %d already finished, %d to run; up to %d at a time" % (len(dsets), len(dsets) - len(todo), len(todo), args.slots))
    running, ok, fail = [], 0, 0
    while todo or running:
        for ds, proc, t0 in list(running):
            if proc.poll() is None:
                continue
            running.remove((ds, proc, t0))
            if finished(ds):
                ok += 1
                log("OK     %s  (%.0f min)  [OK %d  FAIL %d  left %d]" % (tag_of(ds), (time.time() - t0) / 60.0, ok, fail, len(todo) + len(running)))
            else:
                fail += 1
                log("FAILED %s  (see logs/%s.log)" % (tag_of(ds), tag_of(ds)))
        while todo and len(running) < args.slots:
            ds = todo.pop(0)
            env = dict(os.environ)
            env.update(BEASTX_LOGEVERY=str(LOGEVERY))
            proc = subprocess.Popen(["bash", "run_beastx.sh", ds, str(CHAIN), "auto", tag_of(ds), str(SEED)],
                                    stdout=open(os.path.join("logs", tag_of(ds) + ".log"), "w"),
                                    stderr=subprocess.STDOUT, env=env)
            running.append((ds, proc, time.time()))
            log("start  %s" % tag_of(ds))
        time.sleep(20)
    log("Done.  OK %d  FAIL %d" % (ok, fail))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mean(x):
    return sum(x) / len(x) if x else float("nan")


def score(dsets):
    eq = {}
    for r in csv.DictReader(open(os.path.join("tables", "equal_n_runs.tsv")), delimiter="\t"):
        if r["method"] == "beastx":
            eq[(r["condition"], int(r["rep"]), int(r["N"]))] = (fnum(r["nrf_cons"]), fnum(r["p_true"]))
    rows = []
    for ds in dsets:
        try:
            rj = json.load(open(os.path.join("runs_beastx", tag_of(ds), "run.json")))
        except Exception:
            continue
        m = re.match(r"([A-Z]\d)_t\d+_L\d+_rep(\d+)$", ds)
        cond, rep = m.group(1), int(m.group(2))
        e2k, e100 = eq.get((cond, rep, 2000), (None, None)), eq.get((cond, rep, 100), (None, None))
        rows.append(dict(dataset=ds, condition=cond, rep=rep, nrf_20M=fnum(rj.get("nrf")), p_20M=fnum(rj.get("p_true")),
                         n_trees_logged=rj.get("n_trees_logged"), cpu_s=fnum(rj.get("cpu_s")),
                         nrf_N2000=e2k[0], p_N2000=e2k[1], nrf_N100=e100[0], p_N100=e100[1]))
    if not rows:
        print("no finished 20M runs yet")
        return
    os.makedirs("tables", exist_ok=True)
    with open(os.path.join("tables", "conv_check.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n%-22s %8s %8s | %9s %9s | %8s %8s" % ("dataset", "nRF 20M", "P 20M", "nRF N=2k", "P N=2k", "nRF N=100", "P N=100"))
    f = lambda v: ("%.4f" % v) if v is not None else "-"
    for r in rows:
        print("%-22s %8s %8s | %9s %9s | %8s %8s" % (r["dataset"], f(r["nrf_20M"]), f(r["p_20M"]), f(r["nrf_N2000"]), f(r["p_N2000"]), f(r["nrf_N100"]), f(r["p_N100"])))
    for label, kn, kp in (("20M-state chain vs equal-N run at N = 2,000", "nrf_N2000", "p_N2000"),
                          ("20M-state chain vs equal-N run at N = 100", "nrf_N100", "p_N100")):
        print("\n%s (consensus nRF; P(true)); change = equal-N minus 20M" % label)
        for cond in sorted({r["condition"] for r in rows}) + ["all"]:
            s = [r for r in rows if (cond == "all" or r["condition"] == cond) and r[kn] is not None and r["nrf_20M"] is not None]
            if not s:
                continue
            dn = [r[kn] - r["nrf_20M"] for r in s]
            dp = [r[kp] - r["p_20M"] for r in s]
            print("  %-4s n=%2d  nRF mean %.4f vs %.4f (mean change %+.4f, largest |change| %.4f)   "
                  "P(true) mean %.4f vs %.4f (mean change %+.4f, largest |change| %.4f)"
                  % (cond, len(s), mean([r[kn] for r in s]), mean([r["nrf_20M"] for r in s]), mean(dn), max(abs(x) for x in dn),
                     mean([r[kp] for r in s]), mean([r["p_20M"] for r in s]), mean(dp), max(abs(x) for x in dp)))
    print("\nwrote tables/conv_check.csv (%d datasets)" % len(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conds", default="S3,S4,R2")
    ap.add_argument("--reps", default="1-5")
    ap.add_argument("--slots", type=int, default=8, help="concurrent BEAST X runs (default 8)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--score-only", action="store_true")
    args = ap.parse_args()
    a, _, b = args.reps.partition("-")
    reps = list(range(int(a), int(b or a) + 1))
    dsets = [ds_of(c.strip(), r) for c in args.conds.split(",") if c.strip() for r in reps]
    if os.path.isdir("data"):
        try:
            alignment_control(dsets)
        except Exception as e:
            print("alignment control skipped (%s: %s)" % (type(e).__name__, e))
    if args.dry_run:
        for ds in dsets:
            print("%-28s chain %d  logEvery %d  seed %d  %s" % (tag_of(ds), CHAIN, LOGEVERY, SEED, "(finished)" if finished(ds) else ""))
        return
    if not args.score_only:
        run_all(dsets, args)
    score(dsets)


if __name__ == "__main__":
    main()
