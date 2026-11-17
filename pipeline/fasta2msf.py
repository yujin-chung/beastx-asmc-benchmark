#!/usr/bin/env python3
"""
fasta2msf.py — convert a FASTA alignment to the GCG/PileUp ".msf" flavour that
the annealedSMC reference implementation reads.

WHY: SMCSamplerExperiments globs the data directory for "*.msf" only
     (`IO.ls(dataDir, "msf")`).  A .phylip or .fasta file is copied in and then
     silently ignored - the run completes having analysed nothing, with no error.

The layout below was reverse-engineered byte-for-byte (via `cat -A`) from a file
written by the implementation's own data generator, so it should be accepted
verbatim.

VERIFY IT before trusting it - round-trip the generator's own output:

    python3 fasta2msf.py \
        --fasta $ASMC_HOME/state/execs/1.exec/output/sim--2142322983.fasta \
        --out /tmp/roundtrip.msf
    diff /tmp/roundtrip.msf \
         $ASMC_HOME/state/execs/1.exec/output/sim--2142322983.msf \
      && echo "FORMAT MATCH"

Usage:
    python3 fasta2msf.py --fasta data/<id>/aln.fasta --out data/<id>/aln.msf
"""
import argparse, os, sys


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
                names.append(line[1:].split()[0])
            elif line.strip():
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if len(names) != len(seqs):
        sys.exit("malformed FASTA: %d names, %d sequences" % (len(names), len(seqs)))
    return names, seqs


def to_msf(names, seqs, name_width=16, wrap=0):
    """wrap=0 reproduces the generator's own output (one long line per taxon)."""
    n = len(seqs[0])
    for i, s in enumerate(seqs):
        if len(s) != n:
            sys.exit("sequence %s has length %d, expected %d" % (names[i], len(s), n))

    out = []
    out.append("PileUp")
    out.append("")
    out.append("")
    out.append("")
    out.append("   MSF:  %d  Type: P    Check:  666   .. " % n)
    out.append("")
    for nm in names:
        out.append(" Name: %s oo  Len:   %d  Check:  666  Weight:  1.00" % (nm, n))
    out.append("")
    out.append("//")
    # the generator emits three blank lines after the "//" separator
    out.extend(["", "", ""])
    if wrap and wrap > 0:
        for start in range(0, n, wrap):
            for nm, s in zip(names, seqs):
                out.append(nm.ljust(name_width) + s[start:start + wrap])
            out.append("")
    else:
        for nm, s in zip(names, seqs):
            out.append(nm.ljust(name_width) + s)
    # ... and three trailing blank lines at end of file
    out.extend(["", "", ""])
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name-width", type=int, default=16)
    ap.add_argument("--wrap", type=int, default=0,
                    help="0 (default) = one line per taxon, matching the "
                         "reference generator. Set e.g. 60 only if the reader "
                         "chokes on very long lines.")
    args = ap.parse_args()

    names, seqs = read_fasta(args.fasta)
    if len(set(names)) != len(names):
        sys.exit("duplicate taxon names in FASTA")
    too_long = [x for x in names if len(x) >= args.name_width]
    if too_long:
        sys.exit("taxon name(s) >= name width %d, would collide with the "
                 "sequence column: %s" % (args.name_width, too_long[:5]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(to_msf(names, seqs, args.name_width, args.wrap))
    print("wrote %s  (ntax=%d, nchar=%d, names=%s%s)"
          % (args.out, len(names), len(seqs[0]), ", ".join(names[:4]),
             " ..." if len(names) > 4 else ""))


if __name__ == "__main__":
    main()
