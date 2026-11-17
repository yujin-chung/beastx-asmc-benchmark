#!/usr/bin/env python3
"""make_nj_ref.py — a neighbour-joining tree from an alignment, in the exact
newick dialect aSMC's parser accepts.

WHY THIS EXISTS
    SMCSamplerExperiments takes -refTree and hands it to TreeDistancesProcessor,
    which reports each particle's distance to it.  On simulated data that is the
    true tree.  On EMPIRICAL data there is no true tree, but the option still has
    to be given something, so we supply a neighbour-joining tree built from the
    same alignment.

    It is a placeholder for aSMC's own reporting only.  It never enters the
    sampler: the transition kernel, the proposal set, the temperature schedule
    and the weights are all functions of the alignment and the model, not of
    -refTree.  Nothing we report about empirical data is scored against it -
    empirical scoring is method-versus-method in score_empirical.py.

    (Said plainly because a reviewer will ask, and the honest answer is that
    handing a sampler a starting-quality tree could bias it if it were used for
    anything.  It is not.)

Usage:
    python3 make_nj_ref.py --fasta data/EMP_M1809/aln.fasta \
                           --out runs_asmc/<tag>/ref_nj.nwk
"""
import argparse
import re
import sys

LABEL_OK = re.compile(r"^[A-Za-z0-9_]+$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="identity",
                    help="distance model for Bio.Phylo (identity is fine here)")
    a = ap.parse_args()

    try:
        from Bio import AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        from Bio import Phylo
    except ImportError:
        sys.exit("!! biopython missing:  python3 -m pip install --user biopython")
    import io
    import dendropy

    aln = AlignIO.read(a.fasta, "fasta")
    for rec in aln:
        if not LABEL_OK.match(rec.id):
            sys.exit("!! label %r is not plain alphanumeric; aSMC's parser will "
                     "reject the tree built from it" % rec.id)

    dm = DistanceCalculator(a.model).get_distance(aln)
    tree = DistanceTreeConstructor().nj(dm)
    # NJ puts "Inner<N>" labels on internal nodes and can emit negative branch
    # lengths; both break downstream parsers, so clear and clamp them.
    for cl in tree.find_clades():
        if not cl.is_terminal():
            cl.name = None
        if cl.branch_length is not None and cl.branch_length < 0:
            cl.branch_length = 0.0

    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")

    t = dendropy.Tree.get(data=buf.getvalue(), schema="newick",
                          preserve_underscores=True)
    t.seed_node.edge.length = None
    for nd in t.preorder_node_iter():
        if nd.taxon is None:
            nd.label = None
    s = t.as_string(schema="newick", suppress_rooting=True,
                    unquoted_underscores=True,
                    suppress_internal_node_labels=True,
                    suppress_internal_taxon_labels=True,
                    real_value_format_specifier=".12f").strip()
    s = s.replace("'", "").replace('"', "")
    if not s.endswith(";"):
        s += ";"
    for ch in "'\"[]":
        if ch in s:
            sys.exit("!! character %r survived normalisation" % ch)

    with open(a.out, "w") as fh:
        fh.write(s + "\n")
    print("NJ reference tree (placeholder, not used for scoring): %s" % a.out)
    print("  %d taxa" % len(t.taxon_namespace))


if __name__ == "__main__":
    main()
