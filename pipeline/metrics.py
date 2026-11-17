#!/usr/bin/env python3
"""
metrics.py — SHARED SUMMARISER for the BEAST X vs annealed-SMC benchmark.

The whole point: reduce EITHER method's output to a set of weighted tree
samples, then run ONE identical summarisation path.  This removes
"summarisation bias" (differences that come from how each package summarises
its own posterior rather than from the algorithms themselves).

Pipeline per method:
    posterior trees (+ optional weights)
        -> unroot (topology comparison is on UNROOTED bipartitions)
        -> [resample to equal weight if weights are unequal]
        -> SplitDistribution (identical code for both methods)
        -> metrics vs the true tree

Topology metrics (all unrooted):
    rf                 Robinson-Foulds( majority-rule consensus , true )
    nrf                normalized RF = rf / (2*(n-3))
    split_precision/recall/f1   consensus splits vs true splits
    p_true             posterior probability mass on the true topology
    cred95_size        # distinct topologies to reach 95% posterior mass
    mean_rf_to_truth   posterior-weighted E[RF(tree, true)]

Usage as a library:
    from metrics import summarise_treefile
    m = summarise_treefile("run.trees", "true.nwk", schema="nexus",
                            burnin_frac=0.25)

CLI:
    python3 metrics.py --trees run.trees --true true.nwk --schema nexus \
        --burnin 0.25 --out metrics.json
"""
import argparse, json, math, random
import dendropy
from dendropy.calculate import treecompare


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_treelist(path, schema, tns, burnin_frac=0.0, thin=1):
    """Load posterior trees into a TreeList sharing taxon namespace `tns`."""
    tl = dendropy.TreeList.get(path=path, schema=schema, taxon_namespace=tns)
    n = len(tl)
    start = int(math.floor(n * burnin_frac))
    kept = tl[start::thin]
    out = dendropy.TreeList(taxon_namespace=tns)
    for t in kept:
        out.append(t)
    return out


def load_true(path, schema, tns):
    return dendropy.Tree.get(path=path, schema=schema, taxon_namespace=tns)


# ---------------------------------------------------------------------------
# core: reduce to equal-weight, then summarise
# ---------------------------------------------------------------------------
def _unroot(tree):
    tree.is_rooted = False
    tree.collapse_basal_bifurcation()
    tree.encode_bipartitions(suppress_unifurcations=True)
    return tree


def to_equal_weight(trees, weights=None, n_resample=10000, seed=0):
    """
    Return a list of trees with EQUAL weight.
    - weights None/uniform  -> return trees unchanged (BEAST X case).
    - weights given/unequal -> multinomial resample to n_resample (aSMC case).
    This preprocessing is the ONLY method-specific step; everything after is
    identical for both methods.
    """
    if weights is None:
        return list(trees)
    w = list(weights)
    s = sum(w)
    if s <= 0:
        raise ValueError("non-positive total weight")
    w = [x / s for x in w]
    if max(w) - min(w) < 1e-12:          # already uniform
        return list(trees)
    rng = random.Random(seed)
    cum, acc = [], 0.0
    for x in w:
        acc += x
        cum.append(acc)
    out = []
    for _ in range(n_resample):
        r = rng.random()
        lo, hi = 0, len(cum) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if r <= cum[mid]:
                hi = mid
            else:
                lo = mid + 1
        out.append(trees[lo])
    return out


def _topology_key(tree):
    """Canonical unrooted-topology key from the set of nontrivial bipartitions."""
    return frozenset(b.split_bitmask for b in tree.bipartition_encoding
                     if not b.is_trivial())


def summarise(treelist, true_tree, tns):
    """Compute all metrics from an EQUAL-WEIGHT list of trees vs the true tree."""
    n_tax = len(tns)
    # prep true tree
    _unroot(true_tree)
    true_key = _topology_key(true_tree)
    true_splits = set(true_key)

    # split distribution (weighted counting via equal-weight multiset)
    sd = dendropy.SplitDistribution(taxon_namespace=tns)
    topo_counts = {}
    rf_sum = 0.0
    ntrees = 0
    for t in treelist:
        _unroot(t)
        sd.count_splits_on_tree(t, is_bipartitions_updated=True)
        k = _topology_key(t)
        topo_counts[k] = topo_counts.get(k, 0) + 1
        rf_sum += treecompare.symmetric_difference(true_tree, t)
        ntrees += 1

    # majority-rule consensus (>50%)
    consensus = sd.consensus_tree(min_freq=0.5)
    consensus.taxon_namespace = tns
    _unroot(consensus)
    cons_splits = set(_topology_key(consensus))

    # RF(consensus, true), unrooted, and normalization
    rf = treecompare.symmetric_difference(true_tree, consensus)
    denom = 2 * (n_tax - 3) if n_tax > 3 else 1
    nrf = rf / denom if denom else 0.0

    # split precision/recall/F1 (consensus splits vs true splits)
    tp = len(cons_splits & true_splits)
    prec = tp / len(cons_splits) if cons_splits else 0.0
    rec = tp / len(true_splits) if true_splits else 0.0
    f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0

    # P(true topology) and 95% credible set size
    p_true = topo_counts.get(true_key, 0) / ntrees if ntrees else 0.0
    probs = sorted((c / ntrees for c in topo_counts.values()), reverse=True)
    cum, cred95 = 0.0, 0
    for p in probs:
        cum += p
        cred95 += 1
        if cum >= 0.95:
            break

    mean_rf = rf_sum / ntrees if ntrees else float("nan")

    return dict(
        n_taxa=n_tax, n_trees=ntrees,
        rf=rf, nrf=nrf,
        split_precision=prec, split_recall=rec, split_f1=f1,
        p_true=p_true, cred95_size=cred95,
        n_distinct_topologies=len(topo_counts),
        mean_rf_to_truth=mean_rf,
    )


# ---------------------------------------------------------------------------
# convenience wrappers
# ---------------------------------------------------------------------------
def summarise_treefile(trees_path, true_path, schema="nexus",
                       burnin_frac=0.25, thin=1, true_schema="newick"):
    tns = dendropy.TaxonNamespace()
    true_tree = load_true(true_path, true_schema, tns)
    tl = load_treelist(trees_path, schema, tns, burnin_frac=burnin_frac, thin=thin)
    equal = to_equal_weight(tl, weights=None)
    return summarise(equal, true_tree, tns)


def summarise_weighted(newick_weight_pairs, true_path, true_schema="newick",
                       n_resample=10000, seed=0):
    """For aSMC: input is a list of (newick_string, weight)."""
    tns = dendropy.TaxonNamespace()
    true_tree = load_true(true_path, true_schema, tns)
    trees, weights = [], []
    for nwk, w in newick_weight_pairs:
        trees.append(dendropy.Tree.get(data=nwk, schema="newick",
                                        taxon_namespace=tns))
        weights.append(float(w))
    equal = to_equal_weight(trees, weights=weights, n_resample=n_resample, seed=seed)
    return summarise(equal, true_tree, tns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--true", required=True)
    ap.add_argument("--schema", default="nexus")
    ap.add_argument("--true-schema", default="newick")
    ap.add_argument("--burnin", type=float, default=0.25)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    m = summarise_treefile(args.trees, args.true, schema=args.schema,
                           burnin_frac=args.burnin, thin=args.thin,
                           true_schema=args.true_schema)
    txt = json.dumps(m, indent=2)
    print(txt)
    if args.out:
        open(args.out, "w").write(txt)


if __name__ == "__main__":
    main()
