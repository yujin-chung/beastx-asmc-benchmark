#!/usr/bin/env python3
"""
mcc_rescore.py — add a maximum clade credibility (MCC) accuracy metric to
every accuracy-bearing run of the sweep, using the SAME shared machinery as
the published consensus scores.

What it does, per run listed in runs.tsv (phases beastx / asmc / extra):
    1. reload the run's posterior sample
         BEAST X : runs_beastx/<tag>/beastx.trees   (nexus, 25% burn-in)
         aSMC    : runs_asmc/<tag>/particles.tsv    (weights + root repair,
                                                     exactly as parse_asmc.py)
    2. reduce to an equal-weight tree list (metrics.to_equal_weight — the one
       method-specific step, identical to the published pipeline)
    3. compute per-split posterior frequencies, pick the SAMPLED topology that
       maximizes the sum of log split frequencies (the MCC topology; ties go
       to the more frequent topology)
    4. write runs_*/<tag>/metrics_mcc.json with
         nrf_mcc        normalized RF(MCC topology, true topology)
         p_mcc          posterior mass on the MCC topology
         mcc_equals_consensus_nrf   convenience echo of the stored consensus nrf

Resumable: a run whose metrics_mcc.json already exists is skipped.
Output table: mcc.tsv (one row per run) and mcc_summary.csv
(condition x method x budget medians, same grouping as sweep_summary.csv).

Run on the analysis server:
    cd pipeline && source env.sh
    nohup python3 mcc_rescore.py > mcc_rescore.log 2>&1 &
    tail -f mcc_rescore.log
"""

import csv
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import dendropy

import metrics as M
import parse_asmc as PA

BURNIN = 0.25
N_RESAMPLE = 10000
RESAMPLE_SEED = 0
PHASES = ("beastx", "asmc", "extra")
TABLE = "mcc.tsv"
SUMMARY = "mcc_summary.csv"
MAX_WORKERS = 48


# ---------------------------------------------------------------- core MCC --
def mcc_from_equal_weight(trees, true_tree, tns):
    """Equal-weight tree list -> (nrf_mcc, p_mcc, n_trees, n_distinct)."""
    n_tax = len(tns)
    M._unroot(true_tree)
    true_key = M._topology_key(true_tree)

    split_counts = {}
    topo_counts = {}
    n = 0
    keys = []
    for t in trees:
        M._unroot(t)
        k = M._topology_key(t)
        keys.append(k)
        topo_counts[k] = topo_counts.get(k, 0) + 1
        n += 1
    for k, c in topo_counts.items():
        for s in k:
            split_counts[s] = split_counts.get(s, 0) + c

    if n == 0:
        raise ValueError("empty tree list")

    best_key, best_score, best_count = None, -float("inf"), -1
    for k, c in topo_counts.items():
        score = sum(math.log(split_counts[s] / n) for s in k)
        if (score > best_score
                or (score == best_score and c > best_count)):
            best_key, best_score, best_count = k, score, c

    rf = len(best_key ^ true_key)
    denom = 2 * (n_tax - 3) if n_tax > 3 else 1
    return dict(
        nrf_mcc=rf / denom,
        p_mcc=topo_counts[best_key] / n,
        n_trees=n,
        n_distinct_topologies=len(topo_counts),
    )


# ------------------------------------------------------------- per-run I/O --
def _load_beastx_trees(run_dir, tns):
    """Load the posterior tree sample, transparently handling the gzipped
    form: run_beastx.sh gzips beastx.trees for every simulated run below the
    largest budget (chain < 2,000,000) once the summary is computed."""
    norm = os.path.join(run_dir, "beastx_norm.trees")
    plain = os.path.join(run_dir, "beastx.trees")
    gz = plain + ".gz"
    if os.path.exists(norm):
        return M.load_treelist(norm, "nexus", tns, burnin_frac=BURNIN)
    if os.path.exists(plain):
        return M.load_treelist(plain, "nexus", tns, burnin_frac=BURNIN)
    if os.path.exists(gz):
        import gzip as _gz
        with _gz.open(gz, "rt") as fh:
            data = fh.read()
        tl = dendropy.TreeList.get(data=data, schema="nexus",
                                   taxon_namespace=tns)
        start = int(math.floor(len(tl) * BURNIN))
        out = dendropy.TreeList(taxon_namespace=tns)
        for t in tl[start:]:
            out.append(t)
        return out
    raise FileNotFoundError("no beastx(.norm).trees or .trees.gz in %s"
                            % run_dir)


def score_beastx(run_dir, true_path):
    tns = dendropy.TaxonNamespace()
    true_tree = M.load_true(true_path, "newick", tns)
    tl = _load_beastx_trees(run_dir, tns)
    equal = M.to_equal_weight(tl, weights=None)
    return mcc_from_equal_weight(equal, true_tree, tns)


def score_asmc(run_dir, true_path):
    particles = os.path.join(run_dir, "particles.tsv")
    wraw, wlog, raw_trees, _, _ = PA.load(particles, "newick_nonclock")
    if not raw_trees:
        raise ValueError("no particles in %s" % particles)
    w, _uniform = PA.normalise_weights(wraw, wlog)
    fixed = [PA.repair_root_leaf(nwk)[0] for nwk in raw_trees]

    tns = dendropy.TaxonNamespace()
    true_tree = M.load_true(true_path, "newick", tns)
    trees = [dendropy.Tree.get(data=nwk, schema="newick",
                               taxon_namespace=tns) for nwk in fixed]
    if len(tns) != len({x.label for x in tns}):
        raise ValueError("duplicate taxon labels after repair")
    equal = M.to_equal_weight(trees, weights=w,
                              n_resample=N_RESAMPLE, seed=RESAMPLE_SEED)
    return mcc_from_equal_weight(equal, true_tree, tns)


def score_one(row):
    tag, method = row["tag"], row["method"]
    run_dir = os.path.join("runs_beastx" if method == "beastx"
                           else "runs_asmc", tag)
    outj = os.path.join(run_dir, "metrics_mcc.json")
    base = dict(tag=tag, method=method, dataset=row["dataset"],
                condition=row["condition"], phase=row["phase"],
                budget=row["budget"], threads=row["threads"],
                nrf_consensus=row["nrf"])

    if os.path.exists(outj):
        try:
            d = json.load(open(outj))
            base.update(status="CACHED", **{k: d[k] for k in
                        ("nrf_mcc", "p_mcc", "n_trees",
                         "n_distinct_topologies") if k in d})
            return base
        except Exception:
            pass  # fall through and recompute

    true_path = os.path.join("data", row["dataset"], "true.nwk")
    if not os.path.isdir(run_dir):
        base.update(status="NO_RUN_DIR")
        return base
    if not os.path.exists(true_path):
        base.update(status="NO_TRUTH")
        return base
    try:
        d = score_beastx(run_dir, true_path) if method == "beastx" \
            else score_asmc(run_dir, true_path)
    except Exception as e:
        base.update(status="FAIL", err="%s: %s" % (type(e).__name__, e))
        return base

    json.dump(d, open(outj, "w"), indent=2)
    base.update(status="OK", **d)
    return base


# ------------------------------------------------------------------ driver --
def load_runs():
    rows = []
    with open("runs.tsv") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["status"] != "OK":
                continue
            if r["phase"] not in PHASES:
                continue
            rows.append(r)
    return rows


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return float("nan")
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def summarize(rows):
    groups = {}
    for r in rows:
        if r["status"] not in ("OK", "CACHED"):
            continue
        key = (r["condition"], r["method"], int(r["budget"]))
        groups.setdefault(key, []).append(r)
    with open(SUMMARY, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["condition", "method", "budget", "n_reps",
                    "nrf_mcc_median", "nrf_consensus_median",
                    "p_mcc_median"])
        for key in sorted(groups):
            g = groups[key]
            w.writerow([key[0], key[1], key[2], len(g),
                        round(median([float(x["nrf_mcc"]) for x in g]), 4),
                        round(median([float(x["nrf_consensus"]) for x in g]), 4),
                        round(median([float(x["p_mcc"]) for x in g]), 4)])
    print("summary written: %s" % SUMMARY)


def main():
    if not os.path.exists("runs.tsv"):
        sys.exit("runs.tsv not found - run from the bench directory")
    rows = load_runs()
    print("MCC re-scoring %d runs (phases %s) ..."
          % (len(rows), "/".join(PHASES)))

    results = []
    with ProcessPoolExecutor(max_workers=min(MAX_WORKERS,
                                             os.cpu_count() or 8)) as ex:
        for i, r in enumerate(ex.map(score_one, rows), 1):
            results.append(r)
            if r["status"] == "FAIL":
                print("  FAIL  %-50s %s" % (r["tag"], r.get("err", "")),
                      flush=True)
            if i % 100 == 0:
                print("  ... %d / %d done" % (i, len(rows)), flush=True)

    with open(TABLE, "w") as fh:
        cols = ["tag", "method", "dataset", "condition", "phase", "budget",
                "threads", "status", "nrf_mcc", "p_mcc", "nrf_consensus",
                "n_trees", "n_distinct_topologies"]
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print("table written: %s  (%d rows)" % (TABLE, len(results)))

    bad = [r for r in results if r["status"] not in ("OK", "CACHED")]
    print("status: %d ok/cached, %d problems"
          % (len(results) - len(bad), len(bad)))
    for r in bad[:20]:
        print("   %-10s %s" % (r["status"], r["tag"]))

    summarize(results)

    # quick on-screen comparison at the representative budgets
    print()
    print("MCC vs consensus nRF, representative budgets "
          "(BEAST X 2,000,000; aSMC N = 2,000):")
    print("%-5s %-8s %12s %12s" % ("cond", "method", "nrf_mcc", "nrf_cons"))
    for key in [(c, m, b) for c in ("S1", "S2", "S3", "S4", "R1", "R2")
                for m, b in (("beastx", 2000000), ("asmc", 2000))]:
        g = [r for r in results if r["status"] in ("OK", "CACHED")
             and (r["condition"], r["method"], int(r["budget"])) == key]
        if g:
            print("%-5s %-8s %12.4f %12.4f"
                  % (key[0], key[1],
                     median([float(x["nrf_mcc"]) for x in g]),
                     median([float(x["nrf_consensus"]) for x in g])))


if __name__ == "__main__":
    main()
