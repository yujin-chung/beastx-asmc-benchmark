#!/usr/bin/env python3
"""
patch_asmc.py — make the annealedSMC reference implementation write out its final
weighted particle population.

WHY
    The implementation only ever writes its OWN consensus tree
    (`consensus_annealing.newick`) plus summary CSVs.  Scoring that against
    BEAST X summarised by our code would reinstate exactly the summarisation
    bias the shared summariser exists to remove.

WHY *THIS* HOOK
    A first attempt dumped `pc.getSamples()` after `pc.sample(ppk, tdp)` and got
    zero particles.  SMCSampler.sample() ends with:

        if ((kernel.isLastIter()) && schedule != null) {
            if (resampleLastRound) {
                Pair<List<S>, double[]> resampled =
                        resampleAndPack(samples, normalizedWeights, rand);
                samples          = resampled.getFirst();
                normalizedWeights = resampled.getSecond();
            }
            for (int n = 0; n < normalizedWeights.length; n++)
                if (samples.get(n) != null) {
                    tdp.process(samples.get(n), normalizedWeights[n]);
                    samples.set(n, null);          // <-- freed straight away
                }
        }

    So the only moment the final particles and their normalised weights coexist
    is inside that tdp.process() call - and it runs on the last iteration only.
    resampleAndPack() collapses duplicates and sums their weights, so what comes
    through is DISTINCT particles with proper weights, not an equal-weight
    multiset.

WHAT IT DOES
    ev.poi.processors.TreeDistancesProcessor is a non-final public class whose
    `process(Object, double)` is non-final, so we substitute an anonymous
    subclass that calls super and additionally appends to

        results/particles_<treeName>.tsv      weight  newick_unrooted  newick_nonclock

    Nothing else about the run changes.  Two edits, both inside the chosen enum
    branch (ANNEALING by default), since both anchors occur in several branches.

USAGE
    python3 patch_asmc.py --asmc-home $ASMC_HOME --list
    python3 patch_asmc.py --asmc-home $ASMC_HOME --revert   # if already patched
    python3 patch_asmc.py --asmc-home $ASMC_HOME

THEN
    cd $ASMC_HOME && ./gradlew installDist
"""
import argparse, os, re, shutil, sys

ANCHOR_TDP = "TreeDistancesProcessor tdp = new TreeDistancesProcessor();"
ANCHOR_RUN = "pc.sample(ppk, tdp);"
MARKER = "BENCH PATCH"
BRANCH_RE = re.compile(r'^[ \t]+([A-Z][A-Z_0-9]*)[ \t]*\{[ \t]*$')

BLOCK_TDP = r'''				// ==== BENCH PATCH: capture the final weighted particle population ====
				final String benchTag = treeName.replaceAll("[^A-Za-z0-9_.-]", "_");
				final java.io.PrintWriter benchOut = IOUtils.openOutEasy(
						new File(instance.output, "particles_" + benchTag + ".tsv"));
				benchOut.println("weight\tnewick_unrooted\tnewick_nonclock");
				final int[] benchCount = new int[]{0};
				TreeDistancesProcessor tdp = new TreeDistancesProcessor() {
					@Override
					public void process(Object benchState, double benchWeight) {
						super.process(benchState, benchWeight);
						try {
							if (benchState instanceof UnrootedTreeState) {
								UnrootedTreeState benchSt = (UnrootedTreeState) benchState;
								benchOut.println(benchWeight + "\t"
										+ benchSt.getUnrootedTree().toNewick() + "\t"
										+ benchSt.getNonClockTree().toNewick());
								benchCount[0]++;
							}
						} catch (Exception benchEx) {
							LogInfo.logsForce("BENCH: dump error: " + benchEx);
						}
					}
				};
				// ==== END BENCH PATCH ===='''

BLOCK_CLOSE = r'''				// ==== BENCH PATCH: close particle dump ====
				benchOut.close();
				LogInfo.logsForce("BENCH: dumped " + benchCount[0]
						+ " particles -> particles_" + benchTag + ".tsv");
				// ==== END BENCH PATCH ===='''


def branch_of(lines, i):
    for j in range(i - 1, -1, -1):
        m = BRANCH_RE.match(lines[j])
        if m:
            return m.group(1)
    return "?"


def occurrences(lines, anchor):
    return [(i, branch_of(lines, i)) for i, ln in enumerate(lines) if anchor in ln]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asmc-home", required=True)
    ap.add_argument("--branch", default="ANNEALING")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--revert", action="store_true")
    args = ap.parse_args()

    target = os.path.join(args.asmc_home, "src", "main", "java", "smcsampler",
                          "SMCSamplerExperiments.java")
    backup = target + ".bak"
    if not os.path.exists(target):
        sys.exit("not found: %s" % target)

    if args.revert:
        if not os.path.exists(backup):
            sys.exit("no backup at %s" % backup)
        shutil.copy2(backup, target)
        print("reverted from %s" % backup)
        print("rebuild:  cd %s && ./gradlew installDist" % args.asmc_home)
        return

    src = open(target, encoding="utf-8", errors="replace").read()
    lines = src.split("\n")

    occ_t = occurrences(lines, ANCHOR_TDP)
    occ_r = occurrences(lines, ANCHOR_RUN)
    print("TreeDistancesProcessor construction sites:")
    for i, b in occ_t:
        print("  line %-5d branch=%s" % (i + 1, b))
    print("pc.sample(ppk, tdp) sites:")
    for i, b in occ_r:
        print("  line %-5d branch=%s" % (i + 1, b))
    if args.list:
        return

    if MARKER in src:
        sys.exit("\nalready patched - run with --revert first, then re-apply")

    sel_t = [i for i, b in occ_t if b == args.branch]
    sel_r = [i for i, b in occ_r if b == args.branch]
    if len(sel_t) != 1 or len(sel_r) != 1:
        sys.exit("\nbranch %r: found %d tdp site(s) and %d sample site(s); "
                 "need exactly 1 of each" % (args.branch, len(sel_t), len(sel_r)))
    i_t, i_r = sel_t[0], sel_r[0]
    if i_r <= i_t:
        sys.exit("\nunexpected order: pc.sample at %d precedes tdp at %d"
                 % (i_r + 1, i_t + 1))

    if not os.path.exists(backup):
        shutil.copy2(target, backup)
        print("\nbackup -> %s" % backup)

    # edit from the bottom up so the earlier index stays valid
    lines.insert(i_r + 1, BLOCK_CLOSE)
    lines[i_t] = BLOCK_TDP

    open(target, "w", encoding="utf-8").write("\n".join(lines))
    print("patched branch %s: processor at line %d, close after line %d"
          % (args.branch, i_t + 1, i_r + 1))
    print()
    print("next:")
    print("  cd %s && ./gradlew installDist" % args.asmc_home)


if __name__ == "__main__":
    main()
