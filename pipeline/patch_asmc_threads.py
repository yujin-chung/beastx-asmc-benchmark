#!/usr/bin/env python3
"""patch_asmc_threads.py - make AnnealingKernel.nextProposal thread-safe.

THE DEFECT
----------
SMCSampler.propagateAndComputeWeights shares ONE kernel instance across every
Parallelizer worker thread:

    Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(nThreads);
    parallelizer.process(CollUtils.ints(N), new Parallelizer.Processor<Integer>() {
        public void process(Integer x, int _i, int _n, boolean log) {
            Random rand = new Random(seeds[x]);
            final Pair<S, Double> current = kernel.next(rand, samples.get(x));
            ...

and AnnealingKernel fills its proposal list lazily, without synchronisation:

    private ProposalDistribution nextProposal(Random rand) {
        if (proposalDistributions.isEmpty())
            proposalDistributions.addAll(ProposalDistribution.Util.proposalList(...));
        return proposalDistributions.get(rand.nextInt(proposalDistributions.size()));
    }

At t = 1 every worker sees an empty LinkedList at the same moment and they all
call addAll on it concurrently.  LinkedList is not thread-safe, so its size and
links go inconsistent and the run dies with IndexOutOfBoundsException inside
LinkedList.get (line 481) or LinkedList.addAll (line 391).

Observed in our benchmark: 6 of 27 aSMC runs died this way, always within the
first ~15 seconds, and always at higher thread counts (16, 32) or higher
particle counts (N >= 500) - i.e. exactly when more workers touch the empty
list simultaneously.  The list is filled only once and is read-only afterwards,
so any run that survives its first iteration is unaffected.

THE FIX
-------
Double-checked locking on a volatile flag, so the one-time initialisation
happens exactly once.

WHY IT CANNOT CHANGE THE ANSWER
-------------------------------
  * The list contents are deterministic - proposalList() is a function of
    proposalOptions and the initial tree, both fixed before sampling starts.
  * Each particle keeps drawing with its own Random(seeds[x]); the sampler
    already pre-draws per-particle seeds precisely so that results are
    "#-of-threads-invariant" (its own comment).
  * After initialisation the list is never modified, so the added volatile
    read is the entire remaining cost.
Neutrality is verified empirically by verify_thread_patch.sh, which reproduces
an unpatched run bit-for-bit and shows thread-count invariance.

USAGE
-----
    python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME"          # apply
    python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME" --list   # inspect
    python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME" --revert # undo
Then rebuild:  bash rebuild_asmc.sh
"""

import argparse
import os
import re
import shutil
import sys

REL = "src/main/java/smcsampler/AnnealingKernel.java"
MARK = "BENCH THREAD PATCH"
SIG = "private ProposalDistribution nextProposal(Random rand)"
BACKUP_SUFFIX = ".bench_threads_orig"

NEW_METHOD = """\
\t// ---- BENCH THREAD PATCH -------------------------------------------------
\t// One kernel instance is shared by every Parallelizer worker, and the
\t// original code filled this list lazily without synchronisation, so all
\t// threads raced to fill it at t=1 and corrupted the LinkedList
\t// (IndexOutOfBoundsException in LinkedList.get / addAll).  Double-checked
\t// locking makes the one-time initialisation safe.  The sampled distribution
\t// is unchanged: the list contents are deterministic, each particle still
\t// draws with its own Random(seeds[x]), and the list is read-only afterwards.
\tprivate volatile boolean benchProposalsReady = false;

\tprivate ProposalDistribution nextProposal(Random rand) {
\t\tif (!benchProposalsReady) {
\t\t\tsynchronized (this) {
\t\t\t\tif (!benchProposalsReady) {
\t\t\t\t\tif (proposalDistributions.isEmpty())
\t\t\t\t\t\tproposalDistributions.addAll(ProposalDistribution.Util
\t\t\t\t\t\t\t\t.proposalList(proposalOptions,
\t\t\t\t\t\t\t\t\t\tinitial.getNonClockTree(), rand));
\t\t\t\t\tbenchProposalsReady = true;
\t\t\t\t}
\t\t\t}
\t\t}
\t\treturn proposalDistributions.get(rand.nextInt(proposalDistributions
\t\t\t\t.size()));
\t}
\t// ---- END BENCH THREAD PATCH ---------------------------------------------
"""


def find_method(src):
    """Return (start, end) character offsets of the whole nextProposal method,
    including its leading indentation, by brace counting."""
    i = src.find(SIG)
    if i < 0:
        return None
    # walk back to the start of the line (keep the indentation)
    start = src.rfind("\n", 0, i) + 1
    ob = src.find("{", i)
    if ob < 0:
        return None
    depth = 0
    j = ob
    while j < len(src):
        c = src[j]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    else:
        return None
    end = j + 1
    # swallow the trailing newline so the replacement lines up
    if end < len(src) and src[end] == "\n":
        end += 1
    return start, end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asmc-home", default=os.environ.get("ASMC_HOME", ""))
    ap.add_argument("--list", action="store_true",
                    help="show the current state of the method, change nothing")
    ap.add_argument("--revert", action="store_true",
                    help="restore the file from the backup made when patching")
    a = ap.parse_args()

    if not a.asmc_home:
        sys.exit("!! set ASMC_HOME or pass --asmc-home")
    path = os.path.join(a.asmc_home, REL)
    if not os.path.isfile(path):
        sys.exit("!! not found: %s" % path)
    backup = path + BACKUP_SUFFIX

    src = open(path).read()
    patched = MARK in src

    if a.list:
        print("file    : %s" % path)
        print("patched : %s" % ("YES" if patched else "no"))
        print("backup  : %s" % ("present" if os.path.exists(backup) else "none"))
        span = find_method(src)
        if span is None:
            print("!! nextProposal() not located - the source may have changed")
        else:
            s, e = span
            line = src.count("\n", 0, s) + 1
            print("nextProposal() at line %d, %d chars" % (line, e - s))
            print("-" * 70)
            print(src[s:e].rstrip())
            print("-" * 70)
        # show how the kernel is shared, for the record
        smc = os.path.join(a.asmc_home, "src/main/java/smcsampler/SMCSampler.java")
        if os.path.isfile(smc):
            t = open(smc).read()
            m = re.search(r"new Parallelizer<Integer>\([^)]*\)", t)
            if m:
                print("SMCSampler shares the kernel via: %s" % m.group(0))
        return

    if a.revert:
        if not os.path.exists(backup):
            sys.exit("!! no backup at %s - nothing to revert" % backup)
        shutil.copyfile(backup, path)
        print("reverted %s from %s" % (path, backup))
        print("now rebuild:  bash rebuild_asmc.sh")
        return

    if patched:
        print("already patched - nothing to do (use --revert to undo)")
        return

    span = find_method(src)
    if span is None:
        sys.exit("!! could not locate nextProposal() - refusing to guess")
    s, e = span
    old = src[s:e]

    # sanity: the body must be the lazy-init form we expect
    if "proposalDistributions.isEmpty()" not in old or "addAll" not in old:
        sys.exit("!! nextProposal() does not look like the lazy-init version:\n%s" % old)

    if not os.path.exists(backup):
        shutil.copyfile(path, backup)
        print("backup written: %s" % backup)

    new = src[:s] + NEW_METHOD + src[e:]
    open(path, "w").write(new)

    print("patched: %s" % path)
    print("  replaced %d chars at line %d" % (e - s, src.count("\n", 0, s) + 1))
    print("  markers now present: %d" % new.count(MARK))
    print()
    print("--- removed ---")
    print(old.rstrip())
    print("--- inserted ---")
    print(NEW_METHOD.rstrip())
    print()
    print("next:  bash rebuild_asmc.sh   then   bash verify_thread_patch.sh")


if __name__ == "__main__":
    main()
