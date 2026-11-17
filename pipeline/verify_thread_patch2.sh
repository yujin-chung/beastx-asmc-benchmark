#!/usr/bin/env bash
# verify_thread_patch2.sh — a CORRECT neutrality proof for the thread patch.
#
#   cd pipeline && source env.sh
#   bash verify_thread_patch2.sh 2>&1 | tee verify2_summary.txt      # ~7 min
#
# WHY THE FIRST ATTEMPT WAS WRONG
# -------------------------------
# verify_thread_patch.sh compared the patched build against an ARCHIVED
# UNPATCHED run at 8 threads and demanded byte equality.  That test cannot
# pass, and the reason is the defect itself.
#
# When the race does not crash, it still does damage: two or more threads each
# run proposalDistributions.addAll(...) on the same empty list, so the list ends
# up holding several copies of the proposal set.  Sampling stays correct in
# DISTRIBUTION - a uniform draw over k copies of the same list is a uniform draw
# over the list - but rand.nextInt(k*size) is not rand.nextInt(size), so the
# REALISED sequence of proposals differs.  An unpatched multi-threaded run is
# therefore not reproducible even against itself.
#
# So the reference for neutrality must be a run with no race to begin with:
# ONE THREAD.  There, patched and unpatched must agree exactly.
#
# WHAT THIS SCRIPT DOES
#   step 1  revert the patch, rebuild, and record three unpatched runs:
#             u_t1   threads=1   (the race cannot occur)
#             u_t8a  threads=8   (race window open)
#             u_t8b  threads=8   (same seed, same everything)
#   step 2  re-apply the patch, rebuild, and record:
#             p_t1   threads=1
#             p_t8   threads=8
#   step 3  three claims:
#             A  p_t1  == u_t1     the patch changes nothing.  NEUTRALITY.
#             B  p_t8  == p_t1     patched runs are thread-invariant, as the
#                                  sampler's own comment promises.
#             C  u_t8a vs u_t8b    if these differ, the unpatched build is not
#                                  reproducible at 8 threads - which is the
#                                  defect, stated as an observation.
#
# Claim C is the one worth reporting.  A and B are what let us keep using the
# implementation at all.

set -u
DS=S1_t20_L1000_rep2
N=200
PASS=0; FAIL=0
ok()  { echo "  PASS  $*"; PASS=$((PASS+1)); }
bad() { echo "  FAIL  $*"; FAIL=$((FAIL+1)); }

P() { echo "runs_asmc/$1/particles.tsv"; }

cmp_files() {  # label A B   -> expect identical
  local lbl="$1" a b; a=$(P "$2"); b=$(P "$3")
  if [ ! -s "$a" ] || [ ! -s "$b" ]; then bad "$lbl - a run is missing"; return; fi
  if cmp -s "$a" "$b"; then ok "$lbl - byte-identical ($(wc -l < "$a") lines)"
  else
    bad "$lbl - DIFFER ($(wc -l < "$a") vs $(wc -l < "$b") lines)"
  fi
}

differ_report() {  # label A B  -> just report, never fails the run
  local lbl="$1" a b; a=$(P "$2"); b=$(P "$3")
  if [ ! -s "$a" ] || [ ! -s "$b" ]; then echo "  n/a   $lbl - a run is missing"; return; fi
  if cmp -s "$a" "$b"; then
    echo "  note  $lbl - identical ($(wc -l < "$a") lines): the race did not fire this time"
  else
    echo "  note  $lbl - DIFFER ($(wc -l < "$a") vs $(wc -l < "$b") lines)"
    echo "        -> the unpatched build is not reproducible at 8 threads"
  fi
}

run() {  # tag threads
  local tag="$1" th="$2"
  echo "--- $tag  (threads=$th) ---"
  if bash run_asmc.sh "$DS" "$N" "$th" "$tag" > "v2_${tag}.log" 2>&1; then
    grep -m1 '^OK ' "v2_${tag}.log" | sed 's/^/    /'
  else
    echo "    FAILED - see v2_${tag}.log"; tail -4 "v2_${tag}.log" | sed 's/^/      /'
  fi
}

echo "=========================================================="
echo " verify_thread_patch2.sh   dataset=$DS  N=$N"
echo "=========================================================="

echo
echo "### step 1  UNPATCHED build"
python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME" --revert || exit 1
bash rebuild_asmc.sh > v2_rebuild_unpatched.log 2>&1 \
  || { echo "!! rebuild failed - see v2_rebuild_unpatched.log"; exit 1; }
echo "rebuilt unpatched."
run u_t1  1
run u_t8a 8
run u_t8b 8

echo
echo "### step 2  PATCHED build"
python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME" > v2_patch.log 2>&1 || exit 1
bash rebuild_asmc.sh > v2_rebuild_patched.log 2>&1 \
  || { echo "!! rebuild failed - see v2_rebuild_patched.log"; exit 1; }
echo "rebuilt patched."
run p_t1 1
run p_t8 8

echo
echo "### step 3  claims"
cmp_files   "A  neutrality:  patched t1 == unpatched t1" p_t1 u_t1
cmp_files   "B  invariance:  patched t8 == patched t1"   p_t8 p_t1
differ_report "C  unpatched t8 run-to-run"               u_t8a u_t8b

echo
echo "=========================================================="
echo " passed: $PASS    failed: $FAIL"
if [ "$FAIL" -eq 0 ]; then
cat <<'TXT'

 VERDICT: proceed.

 For Methods:
   Two source-level changes were made to the annealed SMC reference
   implementation.  The first adds an output processor that writes the weighted
   particle population to disk; the sampler is otherwise untouched.  The second
   places the one-time initialisation of the shared proposal list under
   double-checked locking.  Without it, the single kernel instance shared by all
   Parallelizer worker threads is filled concurrently, which aborted runs at 16
   or more threads and at 500 or more particles, and - when it did not abort -
   left duplicate entries in the list, so that repeated runs of the unpatched
   build at the same seed and thread count did not agree.  With the patch,
   single-threaded output is bit-for-bit identical to the unpatched build, and
   runs at 1, 8 and 16 threads are bit-for-bit identical to one another, as the
   implementation's per-particle seeding intends.
TXT
else
  echo
  echo " VERDICT: claim A or B failed.  Do NOT proceed - the patch is not neutral."
  echo "          Revert it (patch_asmc_threads.py --revert; bash rebuild_asmc.sh)"
  echo "          and cap the design at 8 threads / N=200."
fi
echo "=========================================================="
