#!/usr/bin/env bash
# run_beastx.sh — one BEAST X run at one chain length, scored with the SHARED
# summariser, with its true CPU cost recorded.  Safe to run many at once.
#
#   cd pipeline && source env.sh
#   bash run_beastx.sh <dataset-id> [chain-length] [beagle-threads] [tag]
#
#   $1  dataset id under data/            (required)
#   $2  chain length                      (default 2000000)
#   $3  beagle threads, or "auto"         (default auto)
#   $4  tag = output directory name       (default <id>_c<chain>)
#
# ON FLAGS
#   -threads N is NOT used: the scaling probe showed that naming it makes BEAST
#   X about 45% SLOWER than leaving the default alone.  Thread control for the
#   scalability figure goes through -beagle_threads, which is the knob that
#   actually moves (1.22x at 5000 bp, 2.91x at 20000 bp, saturating by 8).
#   -beagle_instances is inert here: this BEAGLE build exposes a single CPU
#   resource, so multiple instances have nothing to spread across.
#
# ON COST
#   Core-seconds are taken from /usr/bin/time -v (User + System), never from
#   wall x threads.  BEAST X does not keep every thread busy, so wall x threads
#   would overstate its cost and flatter the SMC side of the comparison - the
#   exact bias this paper exists to avoid.

set -u
ID="${1:?dataset id required}"
CHAIN="${2:-2000000}"
BTHREADS="${3:-auto}"
TAG="${4:-${ID}_c${CHAIN}}"
SEED="${5:-42}"

BENCH="$PWD"
DATA="$BENCH/data/$ID"
OUT="$BENCH/runs_beastx/$TAG"

[ -s "$DATA/aln.fasta" ] || { echo "!! missing $DATA/aln.fasta"; exit 1; }
command -v beast >/dev/null || { echo "!! beast not on PATH - source env.sh"; exit 1; }

# Simulated data are scored against the true tree; empirical data are scored
# method-against-method afterwards, by score_empirical.py.
if [ -s "$DATA/true.nwk" ]; then EMPIRICAL=0; else EMPIRICAL=1; fi

rm -rf "$OUT"; mkdir -p "$OUT"

# ~2000 posterior samples whatever the chain length, so the summariser sees a
# comparable number of trees at every budget and the sweep compares sampling
# effort rather than output size.
LOGEVERY="${BEASTX_LOGEVERY:-$(python3 -c "print(max(100, int($CHAIN/2000)))")}"   # env override pins the interval
# Tree prior: "yule" (default, every main-text run) or "expbrlen" (prior-
# sensitivity arm: the annealed SMC program's iid Exp(10) branch-length prior).
TREEPRIOR="${BEASTX_TREEPRIOR:-yule}"

echo "=== run_beastx  tag=$TAG  id=$ID  chain=$CHAIN  logEvery=$LOGEVERY  beagle_threads=$BTHREADS  treePrior=$TREEPRIOR ==="

python3 make_beastx_xml.py --fasta "$DATA/aln.fasta" \
    --out "$OUT/beastx.xml" --file-stem beastx \
    --chain-length "$CHAIN" --log-every "$LOGEVERY" --kappa 2.0 \
    --tree-prior "$TREEPRIOR" || exit 1

BFLAGS=(-beagle_SSE -overwrite -seed "$SEED")
if [ "$BTHREADS" != "auto" ]; then
  if [ "$BTHREADS" = "1" ]; then
    BFLAGS+=(-beagle_threading_off -beagle_instances 1)
  else
    BFLAGS+=(-beagle_threads "$BTHREADS")
  fi
fi

START=$(date +%s)
( cd "$OUT" && /usr/bin/time -v -o time.txt \
    beast "${BFLAGS[@]}" beastx.xml > beastx.out 2>&1 )
RC=$?
END=$(date +%s); WALL=$((END-START))

CPU=$(awk -F': ' '/User time/{u=$2} /System time/{s=$2}
      END{if(u!="") printf "%.1f", u+s; else printf "NA"}' "$OUT/time.txt" 2>/dev/null)

echo "exit=$RC  wall=${WALL}s  cpu_seconds=$CPU"

if [ ! -s "$OUT/beastx.trees" ]; then
  echo "!! no beastx.trees produced"
  tail -15 "$OUT/beastx.out" | sed 's/^/     /'
  echo "FAILED $TAG"; exit 3
fi
NTREES=$(grep -c "^tree " "$OUT/beastx.trees" 2>/dev/null || echo 0)
echo "trees logged: $NTREES"

if [ "$EMPIRICAL" -eq 1 ]; then
  echo "=== empirical dataset: no true tree, so no accuracy score here ==="
  echo "    trees kept for score_empirical.py (method-vs-method)"
  echo "{\"empirical\": true}" > "$OUT/metrics_beastx.json"
else
  echo "=== score with the SHARED summariser ==="
  python3 metrics.py --trees "$OUT/beastx.trees" --true "$DATA/true.nwk" \
      --schema nexus --burnin 0.25 --out "$OUT/metrics_beastx.json" || {
    echo "!! scoring failed.  If this mentions taxon labels, the data were built"
    echo "   by simulate.py v1 - regenerate with gen_data.sh."
    echo "FAILED $TAG"; exit 4; }
fi

python3 - "$OUT/metrics_beastx.json" "$OUT/run.json" <<PY
import json, sys
m = json.load(open(sys.argv[1]))
m.update(dict(tag="$TAG", dataset="$ID", method="beastx",
              chain_length=$CHAIN, log_every=$LOGEVERY, seed=$SEED,
              tree_prior="$TREEPRIOR",
              beagle_threads="$BTHREADS", n_trees_logged=$NTREES,
              wall_s=$WALL, cpu_s=("$CPU" if "$CPU"=="NA" else float("$CPU")),
              empirical=bool($EMPIRICAL)))
json.dump(m, open(sys.argv[2], "w"), indent=2)
PY

# The .trees file is the bulky part and we have the summary, so keep it only
# where it is still needed: the largest simulated budget, and every empirical
# run (score_empirical.py reads them).
if [ "$EMPIRICAL" -eq 0 ] && [ "$CHAIN" -lt 2000000 ]; then
  gzip -f "$OUT/beastx.trees"
fi

echo "OK $TAG  wall=${WALL}s  cpu=${CPU}s"
