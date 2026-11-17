# bench/env.sh — source this FIRST in EVERY server session:
#     cd ~/bench && source env.sh
#
# Purpose: guarantee BEAST X loads the system BEAGLE 4.0.1 from /usr/local/lib
# and NOT a conda copy.  A silently different BEAGLE makes every timing in the
# scalability figure meaningless, so this runs before anything else.
#
# NOTE: we deliberately do NOT strip conda from PATH — python3/dendropy/biopython
# may legitimately live in conda.  Only LD_LIBRARY_PATH is forced.

# ---- EDIT THESE TWO LINES ----------------------------------------------------
export BEAST_HOME="${BEAST_HOME:-/opt/beastx-10.5.0}"
export ASMC_HOME="${ASMC_HOME:-$HOME/annealedSMC}"
# ------------------------------------------------------------------------------

export PATH="$BEAST_HOME/bin:$PATH"

# EXACTLY this — do not append anything (conda paths here break BEAGLE)
export LD_LIBRARY_PATH="/usr/local/lib"

# annealedSMC refuses to start unless these exist
mkdir -p "$PWD/state/execs"

echo "=============================== env check ==============================="
printf '%-12s %s\n' "BEAST_HOME" "$BEAST_HOME"
printf '%-12s %s\n' "ASMC_HOME"  "$ASMC_HOME"
printf '%-12s %s\n' "LD_LIB"     "$LD_LIBRARY_PATH"
printf '%-12s %s\n' "beast"      "$(command -v beast || echo '*** NOT FOUND ***')"
printf '%-12s %s\n' "java"       "$(command -v java  || echo '*** NOT FOUND ***')"
printf '%-12s %s\n' "seq-gen"    "$(command -v seq-gen || echo '*** NOT FOUND ***')"
printf '%-12s %s\n' "python3"    "$(command -v python3 || echo '*** NOT FOUND ***')"
printf '%-12s %s\n' "BEAGLE"     "$(ls /usr/local/lib/libhmsbeagle*.so* 2>/dev/null | head -3 | tr '\n' ' ')"
printf '%-12s %s\n' "cores"      "$(nproc)"

# java from conda is a common silent problem — warn but do not act
case "$(command -v java)" in
  *conda*|*miniconda*|*anaconda*)
    echo "!! WARNING: java resolves inside conda ($(command -v java))."
    echo "!! BEAST may then load a conda BEAGLE. Prefer the system JDK." ;;
esac

python3 - <<'PY'
for mod, attr in (("dendropy", "__version__"), ("Bio", "__version__")):
    try:
        m = __import__(mod)
        print("%-12s %s" % (mod, getattr(m, attr, "?")))
    except Exception as e:
        print("%-12s *** MISSING *** (%s)" % (mod, e))
PY
echo "========================================================================="
