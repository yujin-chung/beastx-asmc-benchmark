#!/usr/bin/env bash
# get_empirical.sh — fetch the empirical benchmark alignments and put them into
# our data/<id>/ layout with IDENTICAL taxon labels for both methods.
#
#   cd pipeline && source env.sh && bash get_empirical.sh
#
# WHICH DATA AND WHY
#   M336 (= DS1)   27 taxa x 1949 bp   18S rRNA, tetrapods        Hedges et al. 1990
#   M1809 (= DS7)  59 taxa x 1824 bp   mtDNA COII+cytb, primates  Yoder & Yang 2004
#   Both are from Table 1 of Lakner et al. (2008) and are the two real datasets
#   analysed by Wang, Wang & Bouchard-Cote (2020).  Using the .msf files from
#   their own experiment repository means aSMC receives byte-identical input to
#   the published analysis, so our numbers can be checked against theirs - the
#   empirical analysis doubles as a validation of the patched build + summariser.
#
# THE LABEL TRAP
#   The MSF and FASTA distributions of these alignments do NOT use the same taxon
#   labels (Microcebus__rufus_1_138 vs Microcebus__rufus_1.138).  If aSMC reads
#   the MSF and BEAST X reads someone else's FASTA, the two posteriors are over
#   different taxon sets and every comparison is silently wrong.  So we derive
#   the FASTA FROM the MSF here: one source, one label set, by construction.

set -u
BASE="https://raw.githubusercontent.com/shijiaw/annealingExperiment/master"
mkdir -p data

fetch_one() {   # $1 = folder, $2 = name, $3 = our dataset id
  local folder="$1" name="$2" id="$3"
  local dir="data/$id"
  mkdir -p "$dir"
  echo "=== $id  ($name.msf from $folder/) ==="
  if [ ! -s "$dir/aln.msf" ]; then
    curl -fsSL "$BASE/$folder/$name.msf" -o "$dir/aln.msf" \
      || { echo "!! download failed - no network, or the repo moved."; return 1; }
  fi
  echo "  msf: $(wc -c < "$dir/aln.msf") bytes"

  python3 - "$dir" "$id" "$name" <<'PY' || return 1
import json, re, sys
d, ds_id, src = sys.argv[1], sys.argv[2], sys.argv[3]

raw = open("%s/aln.msf" % d, errors="replace").read().splitlines()

# --- header: taxon order (some Name: lines lack the space after the colon) ---
names = []
for ln in raw:
    m = re.match(r"\s*Name:\s*(\S+)", ln)
    if m:
        names.append(m.group(1))
    if ln.strip() == "//":
        break
if not names:
    sys.exit("!! no Name: records found in the MSF header")

# --- body: accumulate sequence blocks, tolerating GCG coordinate lines --------
seqs = {n: [] for n in names}
seen_sep = False
for ln in raw:
    if not seen_sep:
        seen_sep = ln.strip() == "//"
        continue
    parts = ln.split()
    if len(parts) < 2 or parts[0] not in seqs:
        continue                      # blank line or coordinate ruler
    seqs[parts[0]].append("".join(parts[1:]))

joined = {n: "".join(v).replace(".", "-").upper() for n, v in seqs.items()}
lengths = {len(v) for v in joined.values()}
if len(lengths) != 1:
    sys.exit("!! ragged alignment: lengths %s" % sorted(lengths))
n_sites = lengths.pop()

# --- FASTA derived from the SAME source, so labels match by construction -----
with open("%s/aln.fasta" % d, "w") as fh:
    for n in names:
        fh.write(">%s\n%s\n" % (n, joined[n]))

json.dump(dict(dataset_id=ds_id, source_msf="%s.msf" % src,
               source_repo="shijiaw/annealingExperiment",
               taxa=len(names), sites=n_sites, empirical=True,
               note="labels taken verbatim from the MSF; FASTA derived from it"),
          open("%s/meta.json" % d, "w"), indent=2)

print("  taxa=%d  sites=%d" % (len(names), n_sites))
print("  first labels: %s" % ", ".join(names[:3]))
gaps = sum(v.count("-") for v in joined.values())
print("  gap fraction: %.3f" % (gaps / float(len(names) * n_sites)))
PY
  echo
}

fetch_one Table5 M336  EMP_M336
fetch_one Table6 M1809 EMP_M1809

echo "=== summary ==="
for id in EMP_M336 EMP_M1809; do
  [ -s "data/$id/meta.json" ] && python3 -c "
import json;d=json.load(open('data/$id/meta.json'))
print('  %-12s %3d taxa x %5d bp' % (d['dataset_id'], d['taxa'], d['sites']))"
done
echo
echo "NOTE: there is no true.nwk for empirical data.  Scoring is method-vs-method"
echo "      (consensus distance, split-support agreement, run-to-run consistency),"
echo "      not accuracy against a known tree."
