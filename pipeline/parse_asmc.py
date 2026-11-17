#!/usr/bin/env python3
"""
parse_asmc.py — connect the patched annealedSMC particle dump to the SHARED
SUMMARISER, so both methods are scored by identical code.

Input   results/particles_<tree>.tsv   (written by the patch in patch_asmc.py)
            weight <tab> newick_unrooted <tab> newick_nonclock
Output  metrics JSON, same schema as metrics.py's BEAST X output.

TWO THINGS THIS FILE EXISTS TO HANDLE
--------------------------------------
1. pty.UnrootedTree.toNewick() roots the tree at an ARBITRARY LEAF and writes
   that leaf as a ROOT LABEL:

        ((T7:...,(...)internal_1:0.0309)T4;

   Standard parsers (DendroPy included) read this as a SEVEN-taxon tree with a
   root labelled "T4" - and raise no error.  Scored naively against the true
   eight-taxon tree, every RF would come out at its maximum and nothing would
   look wrong.  We rewrite it to the equivalent proper form

        (T7:..., internal_4:..., T4:0.0309);

   which preserves both the topology and every branch length.

2. metrics.py is the validated shared summariser and is deliberately NOT
   modified: normalisation happens here, on the input strings.

Usage
    python3 parse_asmc.py --particles runs_asmc/<id>/particles.tsv \
        --true data/<id>/true.nwk --out runs_asmc/<id>/metrics_asmc.json
"""
import argparse, json, math, os, sys

import dendropy
import metrics as M  # the shared summariser


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load(path, column):
    """Accept a 'weight' column (already normalised, written by the tdp.process
    hook) or a legacy 'logWeight' column."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if column not in header:
            sys.exit("column %r not in header %r" % (column, header))
        ci = header.index(column)
        if "weight" in header:
            wi, wlog = header.index("weight"), False
        elif "logWeight" in header:
            wi, wlog = header.index("logWeight"), True
        else:
            wi, wlog = None, False
        other = [c for c in header if c.startswith("newick_") and c != column]
        oi = header.index(other[0]) if other else None

        wraw, trees, others = [], [], []
        for ln in fh:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) <= ci:
                continue
            trees.append(parts[ci])
            wraw.append(float(parts[wi]) if wi is not None else (0.0 if wlog else 1.0))
            if oi is not None and len(parts) > oi:
                others.append(parts[oi])
    return wraw, wlog, trees, others, (other[0] if other else None)


def normalise_weights(wraw, is_log):
    n = len(wraw)
    if is_log:
        finite = [v for v in wraw if math.isfinite(v)]
        if not finite:
            return [1.0 / n] * n, True
        mx = max(finite)
        raw = [math.exp(v - mx) if math.isfinite(v) else 0.0 for v in wraw]
    else:
        raw = [v if math.isfinite(v) and v > 0 else 0.0 for v in wraw]
    s = sum(raw)
    if s <= 0:
        return [1.0 / n] * n, True
    w = [x / s for x in raw]
    return w, (max(w) - min(w)) < 1e-12


# ---------------------------------------------------------------------------
# the root-label repair
# ---------------------------------------------------------------------------
def repair_root_leaf(nwk):
    """(X:d)LABEL;  ->  (..children of X.., LABEL:d);

    Returns (newick_string, repaired_label_or_None)."""
    t = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True)
    seed = t.seed_node
    kids = seed.child_nodes()
    label = seed.label
    repaired = None

    if label and len(kids) == 1:
        child = kids[0]
        d = child.edge.length
        seed.remove_child(child)
        taxon = t.taxon_namespace.require_taxon(label=label)
        leaf = dendropy.Node(taxon=taxon)
        leaf.edge.length = d
        child.add_child(leaf)
        child.edge.length = None
        child.label = None
        t.seed_node = child
        repaired = label
    elif label and len(kids) > 1:
        # a labelled root with several children: the label is a real internal
        # label, nothing to repair
        seed.label = None

    out = t.as_string(schema="newick", suppress_rooting=True,
                      unquoted_underscores=True,
                      suppress_internal_node_labels=True).strip()
    return out, repaired


def taxa_of(nwk):
    t = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True)
    return {x.label for x in t.taxon_namespace}


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--particles", required=True)
    ap.add_argument("--true", dest="true_path", required=True)
    ap.add_argument("--column", default="newick_nonclock",
                    choices=["newick_nonclock", "newick_unrooted"])
    ap.add_argument("--n-resample", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    wraw, wlog, trees, others, other_name = load(args.particles, args.column)
    if not trees:
        sys.exit("no particles found in %s" % args.particles)
    print("particles: %d   column: %s   weights: %s"
          % (len(trees), args.column, "log" if wlog else "normalised"))

    # --- do the two tree columns agree? --------------------------------------
    if others:
        same = sum(1 for a, b in zip(trees, others) if a == b)
        print("column agreement (%s vs %s): %d/%d identical"
              % (args.column, other_name, same, len(trees)))

    # --- weights -------------------------------------------------------------
    w, uniform = normalise_weights(wraw, wlog)
    print("weights: %s" % ("UNIFORM" if uniform else
                           "non-uniform - applied (resampleAndPack collapsed "
                           "duplicates and summed their weights)"))
    print("  min=%.4g  max=%.4g  sum=%.6f" % (min(w), max(w), sum(w)))

    # --- repair the root-label leaf -----------------------------------------
    fixed, repaired_labels = [], set()
    for nwk in trees:
        f, rep = repair_root_leaf(nwk)
        fixed.append(f)
        if rep:
            repaired_labels.add(rep)
    if repaired_labels:
        print("root-label leaf restored on %d/%d particles: %s"
              % (sum(1 for a, b in zip(trees, fixed) if a != b), len(trees),
                 sorted(repaired_labels)))
    print("distinct particle topologies (as strings): %d / %d"
          % (len(set(fixed)), len(fixed)))

    # --- taxon labels --------------------------------------------------------
    t_true = taxa_of(open(args.true_path).read())
    t_asmc = taxa_of(fixed[0])
    if t_true == t_asmc:
        print("taxon labels: MATCH (%d taxa)" % len(t_true))
    else:
        print("taxon labels: *** MISMATCH ***")
        print("  true only: %s" % sorted(t_true - t_asmc)[:8])
        print("  aSMC only: %s" % sorted(t_asmc - t_true)[:8])
        sys.exit("refusing to score with mismatched labels "
                 "(every RF would be maximal and no error would be raised)")

    # every particle must carry the full taxon set, not just the first
    bad = [i for i, f in enumerate(fixed) if taxa_of(f) != t_true]
    if bad:
        sys.exit("%d particle(s) have a different taxon set, first at index %d"
                 % (len(bad), bad[0]))

    # --- score with the SHARED summariser ------------------------------------
    m = M.summarise_weighted(list(zip(fixed, w)), args.true_path,
                             n_resample=args.n_resample, seed=args.seed)
    m["_source"] = os.path.basename(args.particles)
    m["_column"] = args.column
    m["_n_particles"] = len(fixed)
    m["_weights_uniform"] = bool(uniform)
    m["_root_label_repaired"] = sorted(repaired_labels)

    txt = json.dumps(m, indent=2)
    print()
    print(txt)
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        open(args.out, "w").write(txt)
        print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
