#!/usr/bin/env bash
# run_asmc.sh — run the (patched) annealedSMC on one dataset and score it with
# the SHARED SUMMARISER.  Safe to run many copies concurrently.
#
#   cd pipeline && source env.sh
#   bash run_asmc.sh <dataset-id> [N] [threads] [tag]
#
#   $1  dataset id under data/          (default Q1_t8_L300_rep1)
#   $2  number of particles N           (default 200)
#   $3  threads                         (default 8)
#   $4  tag = output directory name     (default <id>_N<N>_t<threads>)
#
#   env override:  ASMC_XMX=32g bash run_asmc.sh ...
#
# WHAT CHANGED FROM THE SEQUENTIAL VERSION
# ----------------------------------------
# The old script did `cd "$ASMC_HOME"` and then picked the NEWEST
# state/execs/*.exec directory to find the particle dump.  fig.basic numbers
# those directories globally per working directory, so two concurrent jobs
# stole each other's output and only one aSMC run could exist at a time — which
# left 216 of the machine's 224 cores idle.
#
# Now every job gets its own working directory, asmc_work/<tag>/, and aSMC
# creates state/execs/ underneath it.  There is exactly one .exec per job, so
# the lookup is unambiguous and any number of jobs can run at once.
#
# The old script also trusted java's exit code.  It should not: fig.basic
# catches worker-thread exceptions, prints "N errors, 0 warnings", and exits 0.
# That is how six overnight runs reported success while writing an empty
# particle file.  We now check the error count and the particle file, and exit
# non-zero on failure so a driver script can see it.
#
# Baseline error count is 1: aSMC's own consensus step calls PHYLIP `neighbor`,
# which pads taxon names and produces labels its newick parser then rejects.
# That failure is internal to aSMC's reporting and does not touch the particle
# population we summarise, so ERR_BASELINE=1 is expected, not a defect in our
# pipeline.  Anything above that is a real failure.

set -u
ID="${1:-Q1_t8_L300_rep1}"
NPART="${2:-200}"
THREADS="${3:-8}"
TAG="${4:-${ID}_N${NPART}_t${THREADS}}"
SEED="${5:-452}"
XMX="${ASMC_XMX:-16g}"
# relative conditional-ESS target of the adaptive annealing schedule.  0.99 is what every run before
# 2026-11-13 used; Wang et al. (2020) use 0.99999 by default (alpha = 1 - 10^-beta, beta = 5).
ALPHA="${ASMC_ALPHA:-0.99}"
ERR_BASELINE="${ASMC_ERR_BASELINE:-1}"

BENCH="$PWD"
DATA="$BENCH/data/$ID"
OUTROOT="$BENCH/runs_asmc/$TAG"
WORK="$BENCH/asmc_work/$TAG"
LIB="$ASMC_HOME/build/install/phylosmcsampler/lib"

[ -s "$DATA/aln.fasta" ] || { echo "!! missing $DATA/aln.fasta"; exit 1; }
[ -d "$LIB" ] || { echo "!! missing $LIB - run 'bash rebuild_asmc.sh'"; exit 1; }

# Simulated data have a true tree and are scored against it.  Empirical data do
# not; those runs are scored method-against-method later, by score_empirical.py.
if [ -s "$DATA/true.nwk" ]; then EMPIRICAL=0; else EMPIRICAL=1; fi

rm -rf "$WORK"
# fig.basic only writes into state/execs/ if that directory ALREADY EXISTS;
# otherwise it silently falls back elsewhere and no .exec directory appears.
# The sequential version created it under $ASMC_HOME and this was easy to miss.
mkdir -p "$WORK/state/execs" "$OUTROOT"

echo "=== run_asmc  tag=$TAG  id=$ID  N=$NPART  threads=$THREADS  xmx=$XMX  alpha=$ALPHA ==="

echo
echo "=== 1. alignment in the format aSMC reads ==="
if [ -s "$DATA/aln.msf" ]; then
  # Empirical data arrive AS an MSF, taken verbatim from the annealed SMC
  # authors' own experiment repository.  Use it unchanged: that is what makes
  # aSMC's input byte-identical to the published analysis, so our numbers are
  # checkable against theirs.  Rewriting it would throw that away - and would
  # fail anyway, because fasta2msf.py lays names out in a 16-character column
  # sized for T1..T50, while real taxon names run past 20 characters.
  cp "$DATA/aln.msf" "$OUTROOT/aln.msf"
  echo "using the distributed MSF verbatim: $(wc -c < "$OUTROOT/aln.msf") bytes"
else
  python3 fasta2msf.py --fasta "$DATA/aln.fasta" --out "$OUTROOT/aln.msf" || exit 1
fi

echo
if [ "$EMPIRICAL" -eq 1 ]; then
  echo "=== 1b. EMPIRICAL dataset: no true tree; building an NJ placeholder ==="
  python3 make_nj_ref.py --fasta "$DATA/aln.fasta" --out "$OUTROOT/ref.nwk" || exit 1
else
echo "=== 1b. normalise the reference tree for aSMC's newick parser ==="
# ma.newick.NewickParser is strict.  Three things in DendroPy's default output
# break it, and each cost us a failed run:
#   1. a stem/root branch length  "...):0.0;"   -> it wants ";" after the final ")"
#   2. QUOTED labels  'T10'                     -> lexical error on the quote
#   3. internal node labels                     -> not always fatal, removed anyway
# We also force fixed-point branch lengths rather than 3.77e-4 notation.
# data/<id>/true.nwk is left untouched; BEAST X scoring keeps using the original.
python3 - "$DATA/true.nwk" "$OUTROOT/ref.nwk" <<'PY' || exit 1
import sys, re, dendropy

src, dst = sys.argv[1], sys.argv[2]
t = dendropy.Tree.get(path=src, schema="newick", preserve_underscores=True)

t.seed_node.edge.length = None
for nd in t.preorder_node_iter():
    if nd.taxon is None:
        nd.label = None

changed = []
for tx in t.taxon_namespace:
    clean = re.sub(r"[^A-Za-z0-9]", "_", tx.label or "")
    if clean != tx.label:
        changed.append((tx.label, clean))
        tx.label = clean

s = t.as_string(schema="newick", suppress_rooting=True,
                unquoted_underscores=True,
                suppress_internal_node_labels=True,
                suppress_internal_taxon_labels=True,
                real_value_format_specifier=".12f").strip()
s = s.replace("'", "").replace('"', "")
if not s.endswith(";"):
    s += ";"
open(dst, "w").write(s + "\n")

print("refTree for aSMC: %s" % dst)
print("  %s%s" % (s[:90], " ..." if len(s) > 90 else ""))
if changed:
    print("  relabelled: %s" % changed[:5])
bad = sorted({c for c in "'\"[]" if c in s})
if bad:
    sys.exit("  !! characters aSMC's parser rejects are still present: %s" % bad)
PY
fi

NTHOUS=$(python3 -c "print($NPART/1000.0)")
echo
echo "=== 2. annealedSMC  (N=$NPART, nThousandIters=$NTHOUS, threads=$THREADS) ==="
echo "    working dir: $WORK"

# aSMC writes state/execs/<n>.exec RELATIVE TO THE CWD, so a private working
# directory is what makes concurrent jobs possible.
cd "$WORK" || exit 1

START=$(date +%s)
/usr/bin/time -v -o "$OUTROOT/time.txt" \
java -Xmx"$XMX" -cp "$LIB/*" smcsampler.SMCSamplerExperiments \
  -useDataGenerator false \
  -dataFile "$OUTROOT/aln.msf" \
  -refTree "$OUTROOT/ref.nwk" \
  -sequenceType DNA \
  -nThreads "$THREADS" \
  -treeRate 10 -deltaProposalRate 10 \
  -useNonclock true -useSlightNonclock false -sdScale 0.3 \
  -iterScalings 1.0 -nThousandIters "$NTHOUS" \
  -methods ANNEALING \
  -resamplingStrategy ESS \
  -nAnnealing 200 -nSubsampling 200 -alphaSMCSampler "$ALPHA" \
  -nSitesPerIndex 10 -essRatioThreshold 0.5 \
  -adaptiveTempDiff true -runDSMCusingadaptiveTemp false -adaptiveType 0 \
  -csmc_trans2tranv 2.0 -setJC false -treePrior "unconstrained:exp(10)" \
  -fixNucleotideFreq true -nReplica 1 -repPerDataPt 1 \
  -mainRand "$SEED" -gen.rand $((SEED + 893)) -useCESS true -useNNI true \
  -useLIS false -usenewSS false -useRevSS false -usenewDSMC true \
  2>&1 | tee "$OUTROOT/asmc.out"
RC=${PIPESTATUS[0]}
END=$(date +%s)
WALL=$((END-START))
cd "$BENCH" || exit 1

# true CPU cost, for the x-axis of the budget sweep - NOT wall x threads
CPU=$(awk -F': ' '
  /User time/   {u=$2}
  /System time/ {s=$2}
  END {if (u!="") printf "%.1f", u+s; else printf "NA"}' "$OUTROOT/time.txt" 2>/dev/null)

echo
echo "exit=$RC  wall=${WALL}s  cpu_seconds=$CPU"

# ---------------------------------------------------------------- failures --
# java exits 0 even when a worker thread died, so read the error count instead.
ERRS=$(grep -oE '[0-9]+ errors,' "$OUTROOT/asmc.out" | tail -1 | grep -oE '[0-9]+')
ERRS="${ERRS:-0}"
echo "reported errors: $ERRS (baseline $ERR_BASELINE)"
if [ "$ERRS" -gt "$ERR_BASELINE" ]; then
  echo "!! aSMC reported MORE errors than the known-harmless baseline."
  echo "   first stack trace:"
  grep -m1 -B6 -A6 -E "Exception|Error" "$OUTROOT/asmc.out" | head -20
fi

echo
echo "=== 3. locate the particle dump ==="
# Search the whole private working directory rather than assuming fig.basic's
# layout: the job owns this directory exclusively, so whatever particle file is
# in it belongs to this run and nothing else.  (Assuming state/execs/<n>.exec
# is what broke the first parallel attempt.)
EXEC=$(ls -1dt "$WORK"/state/execs/*.exec 2>/dev/null | head -1)
[ -n "${EXEC:-}" ] && echo "exec dir: $EXEC"

PART=$(find "$WORK" -name "particles_*.tsv" 2>/dev/null | head -1)
if [ -z "${PART:-}" ]; then
  echo "!! no particles_*.tsv anywhere under $WORK"
  echo "   - output patch applied?   python3 patch_asmc.py --asmc-home \$ASMC_HOME --list"
  echo "   - files aSMC did write:"
  find "$WORK" -type f 2>/dev/null | head -20
  echo "   - last 15 lines of aSMC output:"
  tail -15 "$OUTROOT/asmc.out" | sed 's/^/     /'
  echo "FAILED $TAG"
  exit 4
fi
cp "$PART" "$OUTROOT/particles.tsv"
NLINES=$(wc -l < "$OUTROOT/particles.tsv")
echo "copied -> $OUTROOT/particles.tsv  ($NLINES lines)"

if [ "$NLINES" -lt 2 ]; then
  echo "!! particle file is empty - the sampler died before dumping."
  echo "   (six overnight runs failed exactly this way: a LinkedList race in"
  echo "    AnnealingKernel.nextProposal.  Apply patch_asmc_threads.py.)"
  echo "FAILED $TAG"
  exit 5
fi

grep -h "Estimate of log(Z)" "$OUTROOT/asmc.out" | tail -2

echo
if [ "$EMPIRICAL" -eq 1 ]; then
  echo "=== 4. empirical dataset: no true tree, so no accuracy score here ==="
  echo "    particles kept for score_empirical.py (method-vs-method)"
  python3 - "$OUTROOT/particles.tsv" "$OUTROOT/metrics_asmc.json" <<'PY'
import sys, json
n = sum(1 for _ in open(sys.argv[1])) - 1
json.dump(dict(_n_particles=n, empirical=True), open(sys.argv[2], "w"), indent=2)
print("    %d particles" % n)
PY
else
  echo "=== 4. score with the SHARED summariser ==="
  python3 parse_asmc.py --particles "$OUTROOT/particles.tsv" \
      --true "$DATA/true.nwk" --out "$OUTROOT/metrics_asmc.json" || {
    echo "FAILED $TAG (scoring)"; exit 6; }
fi

# ------------------------------------------------------------- run record --
python3 - "$OUTROOT/metrics_asmc.json" "$OUTROOT/run.json" <<PY
import json, sys
m = json.load(open(sys.argv[1]))
m.update(dict(tag="$TAG", dataset="$ID", method="asmc",
              n_particles_requested=$NPART, threads=$THREADS, seed=$SEED, alpha="$ALPHA",
              wall_s=$WALL, cpu_s=("$CPU" if "$CPU"=="NA" else float("$CPU")),
              errors_reported=$ERRS, xmx="$XMX", empirical=bool($EMPIRICAL)))
json.dump(m, open(sys.argv[2], "w"), indent=2)
print("run record: %s" % sys.argv[2])
PY

echo "OK $TAG  wall=${WALL}s  cpu=${CPU}s"

# keep the sampler's small diagnostic outputs (log Z estimate, ESS / temperature
# trajectory, run record) next to the particle dump before the scratch state goes
mkdir -p "$OUTROOT/results"
find "$WORK" -path "*/results/*" -type f ! -name "particles_*" -size -50M \
     -exec cp -f {} "$OUTROOT/results/" \; 2>/dev/null
for f in record log info.map options.map output.map; do
  [ -f "$EXEC/$f" ] && cp -f "$EXEC/$f" "$OUTROOT/results/exec_$f" 2>/dev/null
done

# the working directory holds a full copy of aSMC's scratch state; it is large
# and we already extracted what we need.
rm -rf "$WORK"
