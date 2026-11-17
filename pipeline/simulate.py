#!/usr/bin/env python3
"""
simulate.py — generate a true tree and a fixed-kappa HKY nucleotide alignment.

Design (BEAST X vs annealed SMC benchmark, version 2):
  - Substitution model: HKY with FIXED kappa (default 2.0), equal base freqs.
    Rationale: the annealedSMC reference implementation is locked to a fixed
    kappa / equal-frequency / non-clock space, so we match that exactly.
  - Gamma rate heterogeneity: OFF by default (sensitivity only).
  - Tree:
      clock    -> pure-birth (Yule) ultrametric tree, scaled so root-to-tip
                  height = --height subs/site.  BEAST X (strict clock) is then
                  well-specified.
      nonclock -> same topology but branch lengths perturbed by lognormal rate
                  multipliers with MEAN ONE (non-ultrametric).
  - kappa <-> seq-gen ts/tv:  with EQUAL base frequencies, seq-gen's -t
    (expected ts/tv ratio) equals kappa/2.  We default to equal freqs so
    ts/tv = kappa/2.  (Documented; override with --tstv.)

Outputs (per dataset), written to <outdir>/<dataset_id>/ :
  true.nwk         true tree (newick, branch lengths in subs/site)
  aln.fasta        alignment (FASTA)
  aln.nex          alignment (NEXUS, for BEAST X / downstream)
  aln.phylip       alignment (relaxed PHYLIP)
  meta.json        all parameters + seed for reproducibility

Usage:
  python3 simulate.py --taxa 20 --length 1000 --replicate 1 \
      --condition S1 --seed 12345 --outdir data


CHANGES IN VERSION 2 - four defects found while debugging the first run
=======================================================================

(1) LOGNORMAL RATE MULTIPLIERS DID NOT HAVE MEAN ONE.
    v1:  mult = exp(gauss(0.0, sd))          E[mult] = exp(sd^2/2)
    At the default sd = 0.3 that is 1.046, so every 'nonclock' tree was
    silently inflated by ~4.6% in total tree length, and the inflation grew
    with sd.  The rate-variation sensitivity condition is supposed to change
    only the SHAPE of rate variation, holding the expected amount of evolution
    fixed; otherwise BEAST X's strict-clock misspecification is confounded with
    a longer tree.  Fixed by centring the normal at -sd^2/2, the standard
    mean-one normalisation (the same one BEAST's UCLN model uses).

(2) THE NEWICK CARRIED A ROOT STEM LENGTH.
    DendroPy writes the seed node's edge length, producing "...):0.0;".
    aSMC's ma.newick.NewickParser expects ';' straight after the final ')' and
    threw a ParseException.  run_asmc.sh had to strip it at run time.  Fixed at
    the source: the seed edge length is cleared before writing.

(3) TWO-DIGIT LABELS WERE QUOTED, AND THE QUOTES BECAME PART OF THE NAME.
    DendroPy emitted 'T10', 'T20', ... (every label that is a multiple of ten).
    Because run_seqgen fed that same newick to seq-gen, the quotes ended up
    INSIDE the taxon names in the FASTA (">'T10'"), then inside BEAST X's XML,
    and finally in the .trees file as "'T10'".  The scoring step then saw
    'T10' and T10 as two different taxa and every 20- and 50-taxon dataset
    failed to score - 12 of 13 runs in the first overnight batch.  Fixed by
    routing every newick through newick_string(), which strips quotes and then
    ASSERTS that all labels match ^[A-Za-z0-9_]+$, so this cannot recur
    silently.

(4) THE TREE WAS NOT A FAIR SAMPLE OF YULE TREES CONDITIONED ON N TIPS.
    birth_death_tree(num_extant_tips=n) stops the process at the instant the
    n-th lineage appears, so the final speciation always sits exactly at the
    present and the last-born pair always has two zero-length terminal
    branches.  Every v1 tree contains such a cherry (e.g. "(T5:0.0,T18:0.0)"),
    which is not a property of the Yule process - it is an artefact of the
    stopping rule.  Fixed with DendroPy's generalised sampling approach
    (gsa_ntax), which simulates past n tips and samples a tree from the correct
    conditional distribution (Hartmann, Wong & Stadler 2010, Syst Biol 59:465).
    Pass --gsa-multiple 0 to reproduce the v1 behaviour.
"""
import argparse, json, os, re, subprocess, sys, random, math
import dendropy

LABEL_OK = re.compile(r"^[A-Za-z0-9_]+$")


def newick_string(tree, drop_root_stem=True):
    """One-line newick that BOTH seq-gen and aSMC's parser accept.

    DendroPy's defaults are wrong for us in three ways, each of which cost a
    failed run: it writes a root stem length, it quotes some labels, and it can
    emit branch lengths in scientific notation.  Everything is normalised here,
    in one place, so no downstream script has to repair it.
    """
    if drop_root_stem:
        tree.seed_node.edge.length = None
    s = tree.as_string(schema="newick",
                       suppress_rooting=True,
                       unquoted_underscores=True,
                       suppress_internal_node_labels=True,
                       suppress_internal_taxon_labels=True,
                       real_value_format_specifier=".12f").strip()
    s = s.replace("'", "").replace('"', "")
    s = " ".join(s.split())
    if not s.endswith(";"):
        s += ";"

    # If a label ever needs quoting for a real reason, stripping the quotes
    # would corrupt it.  Refuse rather than write a file that scores wrongly.
    bad = [t.label for t in tree.taxon_namespace
           if not t.label or not LABEL_OK.match(t.label)]
    if bad:
        raise RuntimeError("taxon labels are not plain alphanumerics: %s" % bad)
    for ch in "'\"[]":
        if ch in s:
            raise RuntimeError("character %r survived normalisation" % ch)
    return s


def build_tree(n_taxa, mode, height, seed, rate_sd=0.3, gsa_multiple=4):
    """Return a dendropy Tree with branch lengths in subs/site."""
    rng = random.Random(seed)
    taxa = dendropy.TaxonNamespace(["T%d" % (i + 1) for i in range(n_taxa)])

    kw = dict(birth_rate=1.0, death_rate=0.0,
              num_extant_tips=n_taxa, taxon_namespace=taxa,
              rng=rng, repeat_until_success=True)
    # Generalised sampling approach: simulate past n tips and sample back, so
    # the tree is a fair draw conditioned on n tips instead of one frozen at
    # the instant of the n-th speciation (which always yields a zero-length
    # cherry).  gsa_ntax must exceed num_extant_tips.
    if gsa_multiple and gsa_multiple > 1:
        kw["gsa_ntax"] = max(n_taxa + 5, int(n_taxa * gsa_multiple))
    tree = dendropy.simulate.treesim.birth_death_tree(**kw)
    tree.is_rooted = True

    # Scale so root-to-tip (tree height) == height subs/site.
    current = tree.max_distance_from_root()
    if current <= 0:
        raise RuntimeError("degenerate tree height")
    scale = height / current
    for edge in tree.preorder_edge_iter():
        if edge.length is not None:
            edge.length *= scale

    if mode == "nonclock":
        # Lognormal rate multipliers with MEAN ONE.  Centring the normal at
        # -sd^2/2 makes E[exp(X)] = 1, so only the shape of rate variation
        # changes and the expected amount of evolution is held fixed.  v1
        # centred at 0 and inflated every tree by exp(sd^2/2).
        mu = -0.5 * rate_sd * rate_sd
        for edge in tree.preorder_edge_iter():
            if edge.length is not None and edge.length > 0:
                edge.length *= math.exp(rng.gauss(mu, rate_sd))
    return tree


def tree_diagnostics(tree):
    lens = [e.length for e in tree.preorder_edge_iter()
            if e.length is not None and e.tail_node is not None]
    zero = sum(1 for x in lens if x <= 0.0)
    return dict(n_edges=len(lens), total_length=sum(lens),
                min_edge=min(lens) if lens else None,
                zero_length_edges=zero)


def write_alignment(nex_path, fasta_path, phy_path, taxon_labels, seqs):
    with open(fasta_path, "w") as fh:
        for name, s in zip(taxon_labels, seqs):
            fh.write(">%s\n%s\n" % (name, s))
    with open(phy_path, "w") as fh:
        fh.write(" %d %d\n" % (len(seqs), len(seqs[0])))
        for name, s in zip(taxon_labels, seqs):
            fh.write("%s  %s\n" % (name, s))
    with open(nex_path, "w") as fh:
        fh.write("#NEXUS\nbegin data;\n")
        fh.write("  dimensions ntax=%d nchar=%d;\n" % (len(seqs), len(seqs[0])))
        fh.write("  format datatype=dna missing=? gap=-;\n  matrix\n")
        for name, s in zip(taxon_labels, seqs):
            fh.write("    %s  %s\n" % (name, s))
        fh.write("  ;\nend;\n")


def run_seqgen(newick, length, tstv, freqs, gamma_shape, seed):
    """Call seq-gen on the (already normalised) newick string."""
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as tf:
        tf.write(newick + "\n")
        tree_path = tf.name
    cmd = ["seq-gen", "-mHKY", "-t%g" % tstv, "-l%d" % length, "-z%d" % seed,
           "-f", "%g" % freqs[0], "%g" % freqs[1], "%g" % freqs[2], "%g" % freqs[3]]
    if gamma_shape is not None:
        cmd += ["-a%g" % gamma_shape]
    cmd.append(tree_path)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(tree_path)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError("seq-gen failed")
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    ntax = int(lines[0].split()[0])
    names, seqs = [], []
    for ln in lines[1:1 + ntax]:
        parts = ln.split()
        names.append(parts[0])
        seqs.append(parts[1])
    return names, seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxa", type=int, required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--condition", default="S")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--mode", choices=["clock", "nonclock"], default="clock")
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--tstv", type=float, default=None,
                    help="seq-gen ts/tv; default kappa/2 (equal freqs)")
    ap.add_argument("--freqs", default="0.25,0.25,0.25,0.25")
    ap.add_argument("--gamma-shape", type=float, default=None)
    ap.add_argument("--height", type=float, default=0.10,
                    help="root-to-tip distance in subs/site (clock)")
    ap.add_argument("--rate-sd", type=float, default=0.3,
                    help="lognormal sd for mean-one nonclock rate multipliers")
    ap.add_argument("--gsa-multiple", type=float, default=4,
                    help="gsa_ntax = taxa * this (0 or 1 disables GSA)")
    ap.add_argument("--outdir", default="data")
    args = ap.parse_args()

    freqs = tuple(float(x) for x in args.freqs.split(","))
    tstv = args.tstv if args.tstv is not None else args.kappa / 2.0

    dataset_id = "%s_t%d_L%d_rep%d" % (args.condition, args.taxa, args.length,
                                       args.replicate)
    dpath = os.path.join(args.outdir, dataset_id)
    os.makedirs(dpath, exist_ok=True)

    tree = build_tree(args.taxa, args.mode, args.height, args.seed,
                      rate_sd=args.rate_sd, gsa_multiple=args.gsa_multiple)
    diag = tree_diagnostics(tree)

    # ONE normalised newick, used for BOTH the stored true tree and seq-gen's
    # input.  In v1 these were produced separately, which is how quoted labels
    # got into the alignment but not into the scoring namespace.
    nwk = newick_string(tree)
    with open(os.path.join(dpath, "true.nwk"), "w") as fh:
        fh.write(nwk + "\n")

    names, seqs = run_seqgen(nwk, args.length, tstv, freqs,
                             args.gamma_shape, args.seed)

    bad = [n for n in names if not LABEL_OK.match(n)]
    if bad:
        raise RuntimeError("seq-gen returned labels we cannot use: %s" % bad[:5])
    if set(names) != {t.label for t in tree.taxon_namespace}:
        raise RuntimeError("alignment and true tree disagree on the taxon set")

    write_alignment(os.path.join(dpath, "aln.nex"),
                    os.path.join(dpath, "aln.fasta"),
                    os.path.join(dpath, "aln.phylip"),
                    names, seqs)

    meta = dict(dataset_id=dataset_id, taxa=args.taxa, length=args.length,
                replicate=args.replicate, condition=args.condition,
                seed=args.seed, mode=args.mode, kappa=args.kappa, tstv=tstv,
                freqs=freqs, gamma_shape=args.gamma_shape, height=args.height,
                rate_sd=args.rate_sd, gsa_multiple=args.gsa_multiple,
                model="HKY", seqgen_version_note="fixed-kappa, equal-freq",
                simulate_version=2,
                lognormal_mean_one=True, root_stem_written=False,
                tree_diagnostics=diag)
    with open(os.path.join(dpath, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print("wrote %s  (ntax=%d nchar=%d mode=%s tstv=%.3f "
          "tree_len=%.4f zero_edges=%d)" %
          (dpath, len(names), len(seqs[0]), args.mode, tstv,
           diag["total_length"], diag["zero_length_edges"]))


if __name__ == "__main__":
    main()
