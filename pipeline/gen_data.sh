#!/usr/bin/env bash
# gen_data.sh — regenerate every simulated dataset with simulate.py v2.
#
#   cd pipeline && source env.sh && bash gen_data.sh
#
# The v1 data are archived, not deleted: data_v1_<date>/.  Every result we have
# so far was computed from them, so they must stay readable if a number is ever
# questioned.
#
# DESIGN
#   4 clock conditions      BEAST X's strict clock is WELL specified
#     S1  20 taxa x 1000 bp
#     S2  20 taxa x 5000 bp
#     S3  50 taxa x 1000 bp
#     S4  50 taxa x 5000 bp
#   2 nonclock conditions   BEAST X's strict clock is MISspecified, aSMC is not
#     R1  20 taxa x 1000 bp   lognormal rate multipliers, mean one, sd 0.3
#     R2  50 taxa x 1000 bp
#   x 10 replicates = 60 datasets.
#
# WHY TEN REPLICATES
#   Across the three replicates we had, p_true for BEAST X on S1 ran 0.437,
#   0.995, 0.999, and on S2 it ran 0.169, 0.009, 0.999.  Between-replicate
#   variation is larger than any difference we hope to attribute to the methods,
#   so three replicates cannot support a claim.  Ten is the minimum that lets us
#   report a median and an interquartile range rather than three points.

set -u
BASE_SEED=20260828
REP_START="${REP_START:-1}"
REPS="${REPS:-10}"
OUT=data

# REP_START>1 means "add replicates to what is already there".  Seeds are a
# pure function of (condition index, replicate number), so replicate 7 built
# today is byte-identical to replicate 7 built yesterday - adding replicates
# never disturbs the ones the sweep has already used.

CONDS=(
  "S1 20 1000 clock"
  "S2 20 5000 clock"
  "S3 50 1000 clock"
  "S4 50 5000 clock"
  "R1 20 1000 nonclock"
  "R2 50 1000 nonclock"
)

if [ -d "$OUT" ] && [ -z "${KEEP_OLD:-}" ] && [ "$REP_START" -eq 1 ]; then
  ARCH="data_v1_$(date +%Y%m%d_%H%M)"
  echo "archiving existing $OUT -> $ARCH"
  mv "$OUT" "$ARCH"
fi
mkdir -p "$OUT"

echo "generating $(( ${#CONDS[@]} * (REPS - REP_START + 1) )) datasets"
echo "  replicates $REP_START..$REPS   simulate.py v2"
echo

idx=0
fail=0
for spec in "${CONDS[@]}"; do
  set -- $spec
  COND=$1; TAXA=$2; LEN=$3; MODE=$4
  idx=$((idx+1))
  for rep in $(seq "$REP_START" "$REPS"); do
    SEED=$(( BASE_SEED + idx*1000 + rep ))
    if ! python3 simulate.py --taxa "$TAXA" --length "$LEN" --replicate "$rep" \
            --condition "$COND" --seed "$SEED" --mode "$MODE" \
            --outdir "$OUT"; then
      echo "!! FAILED: $COND rep$rep"
      fail=$((fail+1))
    fi
  done
done

echo
echo "=== verification ==="
# The whole point of v2 is that these two checks can never fail silently again.
bad_quote=$(grep -l "'" "$OUT"/*/true.nwk 2>/dev/null | wc -l)
bad_fasta=$(grep -h "^>" "$OUT"/*/aln.fasta 2>/dev/null | grep -c "['\"]" || true)
bad_stem=$(grep -lE '\):[0-9.eE+-]+;[[:space:]]*$' "$OUT"/*/true.nwk 2>/dev/null | wc -l)
zero_edges=$(python3 - "$OUT" <<'PY'
import glob, json, os, sys
tot = 0
for m in glob.glob(os.path.join(sys.argv[1], "*", "meta.json")):
    d = json.load(open(m))
    tot += d.get("tree_diagnostics", {}).get("zero_length_edges", 0)
print(tot)
PY
)

printf "  datasets written        : %d\n" "$(ls -1d "$OUT"/*/ 2>/dev/null | wc -l)"
printf "  simulate.py failures    : %d   (expect 0)\n" "$fail"
printf "  true.nwk with quotes    : %d   (expect 0 - v1 bug 3)\n" "$bad_quote"
printf "  FASTA labels with quotes: %d   (expect 0 - v1 bug 3)\n" "$bad_fasta"
printf "  true.nwk with root stem : %d   (expect 0 - v1 bug 2)\n" "$bad_stem"
printf "  zero-length edges total : %s   (expect 0 - v1 bug 4)\n" "$zero_edges"

echo
if [ "$fail" -eq 0 ] && [ "$bad_quote" -eq 0 ] && [ "$bad_fasta" -eq 0 ] \
   && [ "$bad_stem" -eq 0 ] && [ "$zero_edges" = "0" ]; then
  echo "ALL CHECKS PASSED - data are ready for the sweep."
else
  echo "!! CHECKS FAILED - do not start the sweep."
  exit 1
fi
