#!/usr/bin/env bash
# rebuild_asmc.sh — rebuild the patched annealedSMC.
#
#   cd pipeline && source env.sh && bash rebuild_asmc.sh
#
# Gradle's installDist refuses to write into an existing install directory
# ("neither empty nor does it contain an installation for 'phylosmcsampler'"),
# so we clear it first.  build/libs/phylosmcsampler.jar is untouched and the
# lib/ tree is regenerated from scratch.

set -u
: "${ASMC_HOME:?run 'source env.sh' first}"

echo "clearing $ASMC_HOME/build/install/phylosmcsampler"
rm -rf "$ASMC_HOME/build/install/phylosmcsampler"

cd "$ASMC_HOME" || exit 1
./gradlew installDist || { echo "!! build failed"; exit 1; }

N=$(ls "$ASMC_HOME"/build/install/phylosmcsampler/lib/*.jar 2>/dev/null | wc -l)
M=$(grep -c "BENCH PATCH" "$ASMC_HOME/src/main/java/smcsampler/SMCSamplerExperiments.java")
echo
echo "jars installed : $N        (expect 113)"
echo "patch markers  : $M        (expect 4 = BEGIN/END x 2 sites, 0 = unpatched)"
[ "$M" -eq 4 ] || echo "!! patch not present - run: python3 patch_asmc.py --asmc-home \$ASMC_HOME"
