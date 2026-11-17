#!/usr/bin/env python3
"""Validate metrics.py on controlled posteriors with known answers."""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import dendropy
import metrics as M


def make_true_and_alt(n=8, seed=1):
    """Return (true_tree, alt_tree) on the same taxa, different unrooted topo."""
    import random
    tns = dendropy.TaxonNamespace(["T%d" % (i + 1) for i in range(n)])
    rng = random.Random(seed)
    true = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.0, num_extant_tips=n,
        taxon_namespace=tns, rng=rng, repeat_until_success=True)
    # alt: a clearly different random topology on same taxa
    rng2 = random.Random(seed + 999)
    alt = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.0, num_extant_tips=n,
        taxon_namespace=tns, rng=rng2, repeat_until_success=True)
    for t in (true, alt):
        for e in t.preorder_edge_iter():
            if e.length is not None:
                e.length = 0.05
    return tns, true, alt


def write_posterior(path, tns, trees):
    tl = dendropy.TreeList(taxon_namespace=tns)
    for t in trees:
        tl.append(t)
    tl.write(path=path, schema="nexus", suppress_rooting=True)


def approx(a, b, tol=1e-6):
    return abs(a - b) <= tol


def run():
    tns, true, alt = make_true_and_alt(n=8, seed=1)
    d = tempfile.mkdtemp()
    true_path = os.path.join(d, "true.nwk")
    true.write(path=true_path, schema="newick", suppress_rooting=True)

    # verify true and alt are actually different unrooted topologies
    from dendropy.calculate import treecompare
    tt = dendropy.Tree.get(path=true_path, schema="newick", taxon_namespace=tns)
    tt.is_rooted = False; tt.encode_bipartitions(suppress_unifurcations=True)
    a2 = alt.clone(depth=1); a2.is_rooted = False; a2.encode_bipartitions(suppress_unifurcations=True)
    d_ta = treecompare.symmetric_difference(tt, a2)
    print("sanity: RF(true, alt) =", d_ta, "(must be > 0)")
    assert d_ta > 0

    results = {}

    # Case A: posterior = 100% true  -> rf=0, f1=1, p_true=1, cred95=1
    p = os.path.join(d, "A.trees")
    write_posterior(p, tns, [true.clone(depth=1) for _ in range(100)])
    mA = M.summarise_treefile(p, true_path, schema="nexus", burnin_frac=0.0)
    results["A_all_true"] = mA
    assert mA["rf"] == 0 and approx(mA["split_f1"], 1.0) and approx(mA["p_true"], 1.0)
    assert mA["cred95_size"] == 1

    # Case B: 90 true + 10 alt -> p_true=0.9, consensus==true (majority), rf=0
    p = os.path.join(d, "B.trees")
    trees = [true.clone(depth=1) for _ in range(90)] + [alt.clone(depth=1) for _ in range(10)]
    write_posterior(p, tns, trees)
    mB = M.summarise_treefile(p, true_path, schema="nexus", burnin_frac=0.0)
    results["B_90true_10alt"] = mB
    assert approx(mB["p_true"], 0.90, tol=1e-6), mB["p_true"]
    assert mB["rf"] == 0, mB["rf"]                     # 90% true -> consensus true
    assert mB["n_distinct_topologies"] == 2
    assert mB["cred95_size"] in (1, 2)

    # Case C: 40 true + 60 alt -> consensus should be alt-driven; p_true=0.4
    p = os.path.join(d, "C.trees")
    trees = [true.clone(depth=1) for _ in range(40)] + [alt.clone(depth=1) for _ in range(60)]
    write_posterior(p, tns, trees)
    mC = M.summarise_treefile(p, true_path, schema="nexus", burnin_frac=0.0)
    results["C_40true_60alt"] = mC
    assert approx(mC["p_true"], 0.40, tol=1e-6), mC["p_true"]
    assert mC["rf"] > 0, "consensus should differ from true when true is minority"

    # Case D: weighted (aSMC-like) -> true weight 0.9, alt weight 0.1
    tnsW = dendropy.TaxonNamespace(t.label for t in tns)
    pairs = [
        (true.as_string(schema="newick", suppress_rooting=True).strip(), 0.9),
        (alt.as_string(schema="newick", suppress_rooting=True).strip(), 0.1),
    ]
    mD = M.summarise_weighted(pairs, true_path, n_resample=20000, seed=0)
    results["D_weighted_0.9true"] = mD
    assert 0.87 <= mD["p_true"] <= 0.93, mD["p_true"]   # resample noise band
    assert mD["rf"] == 0

    print("\n=== metric results ===")
    for k, v in results.items():
        print(f"[{k}] rf={v['rf']} nrf={v['nrf']:.3f} f1={v['split_f1']:.3f} "
              f"p_true={v['p_true']:.3f} cred95={v['cred95_size']} "
              f"nTopo={v['n_distinct_topologies']} meanRF={v['mean_rf_to_truth']:.2f}")
    print("\nALL ASSERTIONS PASSED ✅")


if __name__ == "__main__":
    run()
