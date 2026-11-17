#!/usr/bin/env bash
# verify_thread_patch.sh — prove that the thread-safety patch changes nothing
# about the sampled posterior, and that it fixes what it was meant to fix.
#
#   cd pipeline && source env.sh
#   bash verify_thread_patch.sh            # ~4 minutes (stages 1-3)
#   bash verify_thread_patch.sh --full     # + stage 4, the S4 N=500 case (~35 min)
#
# WHY THIS EXISTS
#   The Methods section claims the sampler's transition kernel is unmodified in
#   substance.  A reviewer is entitled to disbelieve that.  So we do not argue
#   it, we show it: the patched build must reproduce an UNPATCHED run
#   bit-for-bit on the same dataset, seed and thread count.
#
# WHAT IS COMPARED
#   stage 1  patched, threads=8  vs  the archived UNPATCHED threads=8 run
#            -> must be byte-identical.  This is the neutrality proof.
#   stage 2  patched, threads=1  vs  patched, threads=8
#            -> must be byte-identical.  aSMC pre-draws per-particle seeds
#               ("so that result for a given randomization is
#                #-of-threads-invariant"), so this is the implementation
#               honouring its own contract once the race is removed.
#   stage 3  patched, threads=16
#            -> unpatched this ALWAYS died in ~10 s.  Must now finish, and must
#               match the threads=1 result.
#   stage 4  patched, S4 N=500, threads=8  (--full only)
#            -> unpatched this died in 17 s.  Must now finish.  This is the run
#               the particle-count axis of the paper depends on.

set -u
FULL=0
[ "${1:-}" = "--full" ] && FULL=1

REF_TAG=asmc_S1_t20_L1000_rep2_N200          # archived unpatched run, threads=8
REF="runs_asmc_arch/$REF_TAG/particles.tsv"
DS=S1_t20_L1000_rep2
N=200

PASS=0; FAIL=0
ok()   { echo "  PASS  $*"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL  $*"; FAIL=$((FAIL+1)); }

echo "=========================================================="
echo " verify_thread_patch.sh"
echo "=========================================================="
echo "patch state:"
python3 patch_asmc_threads.py --asmc-home "$ASMC_HOME" --list | head -4
echo

cmp_files() {   # $1 label  $2 fileA  $3 fileB
  local lbl="$1" a="$2" b="$3"
  if [ ! -s "$a" ]; then bad "$lbl - missing $a"; return; fi
  if [ ! -s "$b" ]; then bad "$lbl - missing $b"; return; fi
  local na nb
  na=$(wc -l < "$a"); nb=$(wc -l < "$b")
  if cmp -s "$a" "$b"; then
    ok "$lbl - byte-identical ($na lines)"
  else
    bad "$lbl - DIFFER ($na vs $nb lines)"
    echo "        first difference:"
    diff "$a" "$b" | head -6 | sed 's/^/        /'
    # a weaker check, in case only ordering changed
    if cmp -s <(sort "$a") <(sort "$b"); then
      echo "        (note: same multiset of lines, only the order differs)"
    fi
  fi
}

run() {  # $1 tag  $2 threads  [$3 dataset] [$4 N]
  local tag="$1" th="$2" ds="${3:-$DS}" n="${4:-$N}"
  echo "--- running $tag  (dataset=$ds N=$n threads=$th) ---"
  if bash run_asmc.sh "$ds" "$n" "$th" "$tag" > "verify_${tag}.log" 2>&1; then
    echo "    finished: $(grep -m1 '^OK ' "verify_${tag}.log")"
    return 0
  else
    echo "    FAILED (see verify_${tag}.log)"
    tail -5 "verify_${tag}.log" | sed 's/^/        /'
    return 1
  fi
}

# ------------------------------------------------------------------ stage 1 --
echo "### stage 1  neutrality: patched t8 must reproduce the unpatched t8 run"
if [ ! -s "$REF" ]; then
  echo "  SKIP - no archived unpatched run at $REF"
  echo "         (that archive is the only unpatched reference we have;"
  echo "          without it neutrality cannot be shown on this machine)"
else
  run verify_t8 8 && cmp_files "stage 1  patched t8 vs UNPATCHED t8" \
      "runs_asmc/verify_t8/particles.tsv" "$REF"
fi
echo

# ------------------------------------------------------------------ stage 2 --
echo "### stage 2  thread invariance: t1 must equal t8"
run verify_t1 1 && cmp_files "stage 2  t1 vs t8" \
    "runs_asmc/verify_t1/particles.tsv" "runs_asmc/verify_t8/particles.tsv"
echo

# ------------------------------------------------------------------ stage 3 --
echo "### stage 3  the fix works: t16 used to die in ~10 s"
if run verify_t16 16; then
  ok "stage 3  threads=16 completed"
  cmp_files "stage 3  t16 vs t1" \
      "runs_asmc/verify_t16/particles.tsv" "runs_asmc/verify_t1/particles.tsv"
else
  bad "stage 3  threads=16 still fails - the patch did not take effect"
  echo "        did you rebuild?   bash rebuild_asmc.sh"
fi
echo

# ------------------------------------------------------------------ stage 4 --
if [ "$FULL" -eq 1 ]; then
  echo "### stage 4  the particle-count axis: S4 N=500 used to die in 17 s"
  if run verify_S4_N500 8 S4_t50_L5000_rep1 500; then
    ok "stage 4  S4 N=500 completed"
    python3 -c "
import json;d=json.load(open('runs_asmc/verify_S4_N500/run.json'))
print('        particles=%s  nrf=%s  p_true=%s  wall=%ss  cpu=%ss'%(
  d.get('_n_particles'),d.get('nrf'),d.get('p_true'),d.get('wall_s'),d.get('cpu_s')))"
  else
    bad "stage 4  S4 N=500 still fails"
  fi
  echo
else
  echo "### stage 4 skipped (pass --full to include the ~35 min S4 N=500 case)"
  echo
fi

# -------------------------------------------------------------------- verdict --
echo "=========================================================="
echo " passed: $PASS    failed: $FAIL"
if [ "$FAIL" -eq 0 ]; then
  echo
  echo " VERDICT: the patch is neutral and effective."
  echo "   Methods can state: the annealed SMC transition kernel, proposal set,"
  echo "   temperature schedule and resampling rule are those of the reference"
  echo "   implementation; the only source changes are (i) an output processor"
  echo "   that writes the weighted particle population to disk and (ii)"
  echo "   double-checked locking on the one-time initialisation of the shared"
  echo "   proposal list, without which runs at 16 or more threads, or with 500"
  echo "   or more particles, abort during their first iteration.  Both were"
  echo "   verified to reproduce unpatched runs bit-for-bit."
else
  echo
  echo " VERDICT: do NOT proceed to the sweep.  Resolve the failures above first."
fi
echo "=========================================================="
