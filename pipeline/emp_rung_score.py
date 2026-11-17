#!/usr/bin/env python3
"""
emp_rung_score.py — scoring for the empirical convergence ladder, with the
two methods compared on the SAME NUMBER OF TREES.

Design (one rung = one aSMC particle count N):
    aSMC     : each run contributes its N particles (weighted), i.e. exactly
               the trees its consensus / MCC are computed from.
    BEAST X  : each chain is logged every 1,000 states and the first 500,000
               states are discarded as burn-in whatever the chain length, so
               a chain of (500 + N) x 1,000 states leaves EXACTLY N trees:

                    N        chain length    trees after burn-in
                    1,000     1,500,000        1,000
                    5,000     5,500,000        5,000
                   10,000    10,500,000       10,000
                   20,000    20,500,000       20,000
                   50,000    50,500,000       50,000
               Sampling interval and burn-in never change between rungs; only
               the chain length is set from N.  (500,000 states is the burn-in
               the existing 2,000,000-state chains used, and the
               20,000,000-state check showed no summary statistic moving
               beyond that point.)

Per rung the script reports, for each method, seed-to-seed consistency on
unrooted topologies (consensus RF, MCC RF, mean |delta split support|, and
for aSMC the spread of the log marginal likelihood estimates), then pools
the seeds of each method (equal run weights, no thinning) and compares the
two pooled posteriors, which now contain the same number of trees.

    python3 emp_rung_score.py --dataset EMP_M1809            # all rungs
    python3 emp_rung_score.py --dataset EMP_M1809 --group asmc 5000
        -> one JSON line with that group's seed consistency (used by emp_converge.py)

Writes tables/emp_converge_<dataset>.json and .csv
"""
import argparse
import csv
import glob
import itertools
import json
import math
import os
import re
import sys

import dendropy
import numpy as np

import gzip

import metrics as M
import parse_asmc as PA

LOG_EVERY = 1000
BURNIN_STATES = 500000
BURNIN_TREES = BURNIN_STATES // LOG_EVERY          # 500 logged trees discarded
RUNG_CHAIN = {n: (500 + n) * LOG_EVERY for n in (1000, 5000, 10000, 20000, 50000)}


# ------------------------------------------------------------------ loading --
def _read_beastx_all(run_dir, tns):
    """All logged trees of a chain, in order (plain or gzipped file)."""
    for name in ("beastx_norm.trees", "beastx.trees", "beastx.trees.gz"):
        p = os.path.join(run_dir, name)
        if os.path.exists(p):
            if p.endswith(".gz"):
                with gzip.open(p, "rt") as fh:
                    data = fh.read()
                return dendropy.TreeList.get(data=data, schema="nexus",
                                             taxon_namespace=tns,
                                             preserve_underscores=True)
            return dendropy.TreeList.get(path=p, schema="nexus", taxon_namespace=tns,
                                         preserve_underscores=True)
    raise FileNotFoundError("no beastx trees in %s" % run_dir)


def keys_beastx(run_dir, tns):
    tl = _read_beastx_all(run_dir, tns)
    keys = []
    # logged states are 0, 1000, ..., chain: drop every tree with state <= 500,000
    # (the first BURNIN_TREES + 1 entries), which leaves chain/1000 - 500 = N trees
    for t in tl[BURNIN_TREES + 1:]:
        M._unroot(t)
        keys.append(M._topology_key(t))
    return keys                                     # equal weights, in chain order


def keys_asmc(run_dir, tns):
    wraw, wlog, raw, _, _ = PA.load(os.path.join(run_dir, "particles.tsv"),
                                    "newick_nonclock")
    w, _ = PA.normalise_weights(wraw, wlog)
    out = []
    for nwk, wi in zip(raw, w):
        if wi <= 0:
            continue
        fixed, _ = PA.repair_root_leaf(nwk)
        t = dendropy.Tree.get(data=fixed, schema="newick", taxon_namespace=tns,
                              preserve_underscores=True)
        M._unroot(t)
        out.append((M._topology_key(t), wi))
    s = sum(wi for _, wi in out)
    return [(k, wi / s) for k, wi in out]


# -------------------------------------------------------------- summaries --
def summarise(sample):
    """sample: list of (key, weight), weights summing to 1."""
    sup, topo = {}, {}
    for k, w in sample:
        topo[k] = topo.get(k, 0.0) + w
        for s in k:
            sup[s] = sup.get(s, 0.0) + w
    cons = frozenset(s for s, v in sup.items() if v > 0.5)
    best, best_sc, best_w = None, -float("inf"), -1.0
    for k, w in topo.items():
        sc = sum(math.log(sup[s]) for s in k)
        if sc > best_sc or (sc == best_sc and w > best_w):
            best, best_sc, best_w = k, sc, w
    return dict(support=sup, consensus=cons, mcc=best, p_mcc=topo[best],
                n_distinct=len(topo))


def pair_stats(a, b):
    keys = set(a["support"]) | set(b["support"])
    x = np.array([a["support"].get(k, 0.0) for k in keys])
    y = np.array([b["support"].get(k, 0.0) for k in keys])
    r = float(np.corrcoef(x, y)[0, 1]) if len(keys) > 2 and x.std() > 0 and y.std() > 0 else float("nan")
    return dict(cons_rf=len(a["consensus"] ^ b["consensus"]),
                mcc_rf=len(a["mcc"] ^ b["mcc"]),
                dsupport=float(np.mean(np.abs(x - y))),
                pearson=r)


def resolution(summary, thr=0.95):
    return sum(1 for v in summary["support"].values() if v >= thr)


def pool(samples):
    n = len(samples)
    return [(k, w / n) for smp in samples for k, w in smp]


def equal_weight(keys):
    n = len(keys)
    return [(k, 1.0 / n) for k in keys]


# --------------------------------------------------------------- metadata --
def budget_of(tag):
    m = re.search(r"_(?:c|N)(\d+)_", tag + "_")
    return int(m.group(1)) if m else -1


def seed_of(tag):
    m = re.search(r"_s(\d+)$", tag)
    return int(m.group(1)) if m else -1


def run_json(run_dir):
    j = os.path.join(run_dir, "run.json")
    try:
        return json.load(open(j))
    except Exception:
        return {}


def cpu_of(run_dir):
    v = run_json(run_dir).get("cpu_s")
    if v not in (None, "NA"):
        try:
            return float(v)
        except Exception:
            pass
    t = os.path.join(run_dir, "time.txt")
    if os.path.exists(t):
        u = s = None
        for line in open(t):
            if "User time" in line: u = float(line.split(":")[-1])
            if "System time" in line: s = float(line.split(":")[-1])
        if u is not None:
            return u + (s or 0.0)
    return float("nan")


LOGZ_RE = re.compile(r"Estimate of log\(Z\):\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
LOGZ_FALLBACK_RE = re.compile(r"(log\s*\(?Z\)?\w*|logNorm\w*|marginal\s*likelihood\w*|logEvidence)"
                              r"[^-\d]{0,20}(-?\d+\.\d+)", re.I)


def logZ_of(tag):
    """The sampler prints 'Estimate of log(Z): <value>' to stdout at the end of a
    run (kept in runs_asmc/<tag>/asmc.out) and, in newer runs, also leaves
    results/logZout.csv.  Returns (value or None, where-found)."""
    run_dir = os.path.join("runs_asmc", tag)
    out = os.path.join(run_dir, "asmc.out")
    if os.path.exists(out):
        vals = []
        with open(out, errors="replace") as fh:
            for line in fh:
                m = LOGZ_RE.search(line)
                if m:
                    vals.append(float(m.group(1)))
        if vals:
            return vals[-1], "asmc.out"
    csvp = os.path.join(run_dir, "results", "logZout.csv")
    if os.path.exists(csvp):
        nums = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", open(csvp, errors="replace").read())
        if nums:
            return float(nums[-1]), "results/logZout.csv"
    if os.path.exists(out):                       # last resort: any Z-like number
        cands = []
        with open(out, errors="replace") as fh:
            for line in fh:
                for m in LOGZ_FALLBACK_RE.finditer(line):
                    cands.append(float(m.group(2)))
        if cands:
            return cands[-1], "asmc.out (fallback pattern)"
    return None, None


def annealing_steps(tag):
    """Number of annealing steps R from the ESS/temperature trajectory, if kept."""
    for name in ("essTempDiffAdaptive0.csv",):
        p = os.path.join("runs_asmc", tag, "results", name)
        if os.path.exists(p):
            with open(p) as fh:
                return max(0, sum(1 for _ in fh) - 1)
    return None


# --------------------------------------------------------------- loading --
CACHE = "topokeys_v1.pkl"


def canonical_labels(ds):
    """Sorted taxon labels, read from the first aSMC particle of any run of the
    dataset.  Every worker builds its TaxonNamespace in this order, so split
    bitmasks are identical across runs and across processes (cacheable)."""
    for p in sorted(glob.glob("runs_asmc/*%s*/particles.tsv" % ds)):
        wraw, wlog, raw, _, _ = PA.load(p, "newick_nonclock")
        if raw:
            fixed, _ = PA.repair_root_leaf(raw[0])
            t = dendropy.Tree.get(data=fixed, schema="newick", preserve_underscores=True)
            return sorted(x.label for x in t.taxon_namespace)
    raise SystemExit("no aSMC particles found for %s" % ds)


def _load_one(args):
    """Worker: returns (tag, record) with topology keys, using a per-run cache."""
    import pickle
    tag, method, run_dir, labels = args
    cache = os.path.join(run_dir, CACHE)
    rec = dict(method=method, budget=budget_of(tag), seed=seed_of(tag), cpu_s=cpu_of(run_dir))
    try:
        if os.path.exists(cache):
            c = pickle.load(open(cache, "rb"))
            if c.get("labels") == labels:
                rec.update(c["data"])
                return tag, rec
        tns = dendropy.TaxonNamespace(labels)
        if method == "beastx":
            data = dict(keys=keys_beastx(run_dir, tns))
        else:
            data = dict(sample=keys_asmc(run_dir, tns))
        if len(tns) != len(labels):
            extra = sorted(set(x.label for x in tns) - set(labels))
            return tag, dict(error="labels not in the canonical set: %s" % extra[:5])
        pickle.dump(dict(labels=labels, data=data), open(cache, "wb"))
        rec.update(data)
        return tag, rec
    except Exception as e:
        return tag, dict(error="%s: %s" % (type(e).__name__, e))


def load_runs(ds, tns, want=None, workers=4):
    """want: None (all), or ('asmc', N) / ('beastx', chain).  `tns` is filled
    with the canonical labels so that the caller's bitmask arithmetic matches."""
    from concurrent.futures import ProcessPoolExecutor
    labels = canonical_labels(ds)
    for lab in labels:
        tns.require_taxon(label=lab)
    jobs = []
    for p in sorted(glob.glob("runs_beastx/*%s*/" % ds)):
        tag = os.path.basename(p.rstrip("/"))
        rj = run_json(p)
        if not rj or int(rj.get("log_every", 0) or 0) != LOG_EVERY:
            continue
        if want and (want[0] != "beastx" or budget_of(tag) != want[1]):
            continue
        jobs.append((tag, "beastx", p, labels))
    for p in sorted(glob.glob("runs_asmc/*%s*/" % ds)):
        tag = os.path.basename(p.rstrip("/"))
        if want and (want[0] != "asmc" or budget_of(tag) != want[1]):
            continue
        if not os.path.exists(os.path.join(p, "particles.tsv")):
            continue
        jobs.append((tag, "asmc", p, labels))
    runs = {}
    print("  loading %d runs (%d parallel workers; cached runs are instant) ..."
          % (len(jobs), workers), file=sys.stderr, flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for tag, rec in ex.map(_load_one, jobs):
            if "error" in rec:
                print("  skip %s (%s)" % (tag, rec["error"]), file=sys.stderr, flush=True)
                continue
            if rec["method"] == "asmc":
                lz, where = logZ_of(tag)
                rec.update(logZ=lz, logZ_where=where, R=annealing_steps(tag))
            else:
                rec.update(logZ=None)
            runs[tag] = rec
            print("  loaded %s" % tag, file=sys.stderr, flush=True)
    return runs


def sample_for(run, n_trees):
    """The trees a run contributes at a rung of size n_trees."""
    if run["method"] == "asmc":
        return run["sample"]                          # its N weighted particles
    if len(run["keys"]) < n_trees:
        raise ValueError("only %d retained trees, need %d" % (len(run["keys"]), n_trees))
    return equal_weight(run["keys"][:n_trees])        # first N trees after burn-in


def group_stats(tags, runs, n_tax, n_trees):
    sums = {t: summarise(sample_for(runs[t], n_trees)) for t in tags}
    pairs = [pair_stats(sums[a], sums[b]) for a, b in itertools.combinations(tags, 2)]
    denom = 2 * (n_tax - 3)
    lz = [runs[t]["logZ"] for t in tags if runs[t].get("logZ") is not None]
    d = dict(n_seeds=len(tags), n_trees_per_run=n_trees, tags=tags,
             cpu_s=[runs[t]["cpu_s"] for t in tags],
             n_distinct=[sums[t]["n_distinct"] for t in tags],
             p_mcc=[sums[t]["p_mcc"] for t in tags])
    if pairs:
        for k in ("cons_rf", "mcc_rf", "dsupport"):
            v = [p[k] for p in pairs]
            d[k + "_mean"] = float(np.mean(v)); d[k + "_max"] = float(np.max(v))
        d["cons_nrf_mean"] = d["cons_rf_mean"] / denom
        d["mcc_nrf_mean"] = d["mcc_rf_mean"] / denom
    d["logZ"] = lz
    d["annealing_steps"] = [runs[t].get("R") for t in tags if runs[t].get("R") is not None]
    d["logZ_spread"] = (max(lz) - min(lz)) if len(lz) >= 2 else None
    d["logZ_sd"] = float(np.std(lz, ddof=1)) if len(lz) >= 2 else None
    return d, sums


# ------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="EMP_M1809")
    ap.add_argument("--group", nargs=2, metavar=("METHOD", "BUDGET"))
    args = ap.parse_args()
    ds = args.dataset
    tns = dendropy.TaxonNamespace()

    if args.group:
        meth, bud = args.group[0], int(args.group[1])
        runs = load_runs(ds, tns, want=(meth, bud))
        tags = sorted(runs)
        if len(tags) < 2:
            print(json.dumps(dict(n_seeds=len(tags)))); return
        n_trees = bud if meth == "asmc" else next(
            (n for n, c in RUNG_CHAIN.items() if c == bud), len(runs[tags[0]]["keys"]))
        d, _ = group_stats(tags, runs, len(tns), n_trees)
        d.pop("tags", None)
        print(json.dumps(d))
        return

    runs = load_runs(ds, tns)
    if len(tns) < 4:
        sys.exit("no runs loaded")
    n_tax = len(tns); denom = 2 * (n_tax - 3)
    # every run must carry the same taxon set, or between-method RF is meaningless
    full = tns.all_taxa_bitmask()
    for t, r in runs.items():
        keys = r["keys"] if r["method"] == "beastx" else [k for k, _ in r["sample"]]
        seen = 0
        for k in keys[:50]:
            for s_ in k:
                seen |= s_
        if seen & ~full:
            sys.exit("taxon bitmask overflow in %s" % t)
    labels = sorted(x.label for x in tns)
    print("taxa (%d): %s%s" % (n_tax, ", ".join(labels[:6]), " ..." if n_tax > 6 else ""))
    if n_tax != 59 and ds == "EMP_M1809":
        sys.exit("!! expected 59 taxa for M1809 but the runs define %d distinct labels - "
                 "taxon labels differ between programs; refusing to compare" % n_tax)
    print("%s: %d taxa, %d runs loaded (BEAST X chains logged every %d states, "
          "first %d states discarded as burn-in)" % (ds, n_tax, len(runs), LOG_EVERY,
                                                     BURNIN_STATES))

    groups = {}
    for t, r in runs.items():
        groups.setdefault((r["method"], r["budget"]), []).append(t)

    out = dict(dataset=ds, n_taxa=n_tax, rungs={})
    rows = []
    print("\n%-7s %-9s %-8s %5s %12s %12s %10s %11s %10s" %
          ("rung N", "method", "setting", "seeds", "cons RF", "MCC RF",
           "|dsupp|", "logZ spread", "CPU s/run"))
    for N in sorted(RUNG_CHAIN):
        chain = RUNG_CHAIN[N]
        a_tags = sorted(groups.get(("asmc", N), []))
        b_tags = sorted(groups.get(("beastx", chain), []))
        if not a_tags and not b_tags:
            continue
        rung = dict(N=N, chain=chain)
        sums = {}
        for meth, tags, setting in (("beastx", b_tags, "%dM" % (chain // 1000000)),
                                    ("asmc", a_tags, "N=%d" % N)):
            if not tags:
                continue
            try:
                d, s = group_stats(tags, runs, n_tax, N)
            except ValueError as e:
                print("%-7d %-9s %-8s   %s" % (N, meth, setting, e)); continue
            sums[meth] = s
            rung[meth] = d
            if d["n_seeds"] >= 2:
                print("%-7d %-9s %-8s %5d %5.1f / %-4d %5.1f / %-4d %10.4f %11s %10.0f" %
                      (N, meth, setting, d["n_seeds"], d["cons_rf_mean"], d["cons_rf_max"],
                       d["mcc_rf_mean"], d["mcc_rf_max"], d["dsupport_mean"],
                       ("%.2f" % d["logZ_spread"]) if d["logZ_spread"] is not None else "-",
                       float(np.nanmedian(d["cpu_s"]))))
                rows.append([N, meth, setting, d["n_seeds"], N,
                             round(d["cons_rf_mean"], 2), d["cons_rf_max"],
                             round(d["mcc_rf_mean"], 2), d["mcc_rf_max"],
                             round(d["dsupport_mean"], 4),
                             "" if d["logZ_spread"] is None else round(d["logZ_spread"], 3),
                             round(float(np.nanmedian(d["cpu_s"])), 1),
                             int(np.median(d["n_distinct"]))])
            else:
                print("%-7d %-9s %-8s %5d   (single seed)" % (N, meth, setting, d["n_seeds"]))
        # pooled between-method comparison at equal tree counts
        if "beastx" in sums and "asmc" in sums and len(b_tags) >= 1 and len(a_tags) >= 1:
            k = min(len(b_tags), len(a_tags))        # same number of runs on both sides
            pb = pool([sample_for(runs[t], N) for t in b_tags[:k]])
            pa = pool([sample_for(runs[t], N) for t in a_tags[:k]])
            sb, sa = summarise(pb), summarise(pa)
            st = pair_stats(sb, sa)
            rung["between"] = dict(runs_each=k, trees_each=k * N,
                                   pearson=st["pearson"], cons_rf=st["cons_rf"],
                                   mcc_rf=st["mcc_rf"], cons_nrf=st["cons_rf"] / denom,
                                   mcc_nrf=st["mcc_rf"] / denom, dsupport=st["dsupport"],
                                   resolved95_beastx=resolution(sb),
                                   resolved95_asmc=resolution(sa),
                                   pooled_distinct=[sb["n_distinct"], sa["n_distinct"]])
            print("        between methods, %d runs x %d trees each side: r=%.3f  "
                  "cons RF %d (%.3f)  MCC RF %d (%.3f)  |dsupp| %.3f  "
                  "resolved>=0.95 %d vs %d"
                  % (k, N, st["pearson"], st["cons_rf"], st["cons_rf"] / denom,
                     st["mcc_rf"], st["mcc_rf"] / denom, st["dsupport"],
                     resolution(sb), resolution(sa)))
        out["rungs"][str(N)] = rung

    found = [(t, r["logZ_where"]) for t, r in runs.items() if r.get("logZ") is not None]
    print()
    if found:
        print("logZ read from e.g. %s -> %s" % found[0])
    else:
        print("logZ: no marginal-likelihood estimate found in asmc.out / asmc_work "
              "(adjust LOGZ_RE once the output format is known)")

    os.makedirs("tables", exist_ok=True)
    jp = os.path.join("tables", "emp_converge_%s.json" % ds)
    json.dump(out, open(jp, "w"), indent=2)
    cp = os.path.join("tables", "emp_converge_%s.csv" % ds)
    with open(cp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rung_N", "method", "setting", "n_seeds", "trees_per_run",
                    "cons_rf_mean", "cons_rf_max", "mcc_rf_mean", "mcc_rf_max",
                    "dsupport_mean", "logZ_spread", "cpu_s_median", "n_distinct_median"])
        w.writerows(rows)
    print("wrote %s and %s" % (jp, cp))


if __name__ == "__main__":
    main()
