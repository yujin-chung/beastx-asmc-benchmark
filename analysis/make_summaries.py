#!/usr/bin/env python3
"""Means, standard deviations and paired differences behind the numbers quoted in the text.

Writes to output/tables/:
    summary_all_replicates.csv   mean and SD per condition, N and method on all the replicates each analysis has
                                 (BEAST X and annealed SMC at alpha = 0.99: 50; annealed SMC at alpha = 0.99999:
                                 50 at N = 100 and replicates 1-5 at N > 100)
    paired_n100.csv              paired differences between methods at N = 100 over the 50 replicates, with 95%
                                 t intervals
    summary_replicates_1-5.csv   mean and SD of accuracy, CPU time, wall-clock time and annealing steps for the three
                                 analyses on replicates 1-5
    paired_replicates_1-5.csv    paired differences on replicates 1-5, per N and for the per-alignment mean over N
"""
import csv
import math
import os

import common as C


def fmt(v):
    return "" if (v is None or (isinstance(v, float) and math.isnan(v))) else ("%.6g" % v if isinstance(v, float) else v)


def main():
    os.makedirs(C.OUT_TAB, exist_ok=True)
    D, D5 = C.load_all(), C.load_five()
    with open(os.path.join(C.OUT_TAB, "summary_all_replicates.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["condition", "N", "method", "replicates"] + [q + s for q in C.QS_ALL for s in ("_mean", "_sd")])
        for c in C.CONDS:
            for n in C.NS:
                for m in C.METHODS:
                    w.writerow([c, n, C.NAME[m], len(C.reps(m, n))] +
                               ["%.6g" % v for q in C.QS_ALL for v in C.mean_sd(C.vec_all(D, m, c, n, q))])
    with open(os.path.join(C.OUT_TAB, "paired_n100.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["condition", "contrast", "quantity", "replicates", "diff", "lo95", "hi95"])
        for c in C.CONDS:
            for q in ("nrf_mcc", "nrf_cons", "p_true"):
                for a, b in (("hi", "bx"), ("lo", "bx"), ("hi", "lo")):
                    w.writerow([c, "%s minus %s" % (C.NAME[a], C.NAME[b]), q, 50] +
                               ["%.4f" % v for v in C.paired(C.vec_all(D, a, c, 100, q), C.vec_all(D, b, c, 100, q))])
    keys = ["condition", "N", "method"] + [q + s for q in C.QS5 for s in ("_mean", "_sd")]
    with open(os.path.join(C.OUT_TAB, "summary_replicates_1-5.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(keys)
        for c in C.CONDS:
            for n in C.NS:
                for m in C.METHODS:
                    row = [c, n, C.NAME[m]]
                    for q in C.QS5:
                        if (m, q) in D5[(c, n, 1)]:
                            row += [fmt(v) for v in C.mean_sd(C.vec5(D5, c, n, m, q))]
                        else:
                            row += ["", ""]
                    w.writerow(row)
    with open(os.path.join(C.OUT_TAB, "paired_replicates_1-5.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["condition", "N", "contrast", "quantity", "diff", "lo95", "hi95"])
        for c in C.CONDS:
            for q in ("nrf_mcc", "nrf_cons", "p_true"):
                for a, b in (("hi", "bx"), ("lo", "bx"), ("hi", "lo")):
                    lab = "%s minus %s" % (C.NAME[a], C.NAME[b])
                    for n in C.NS:
                        w.writerow([c, n, lab, q] + ["%.4f" % v for v in C.paired(C.vec5(D5, c, n, a, q), C.vec5(D5, c, n, b, q))])
                    w.writerow([c, "all", lab, q] + ["%.4f" % v for v in C.paired(C.over_n5(D5, c, a, q), C.over_n5(D5, c, b, q))])
    print("wrote four summary tables to", os.path.relpath(C.OUT_TAB, C.ROOT))
    print("\nN = 100, 50 replicates: annealed SMC at alpha = 0.99999 minus BEAST X, normalized RF of the MCC tree (95% CI)")
    for c in C.CONDS:
        print("  %s  %+.3f (%+.4f, %+.4f)" % ((c,) + C.paired(C.vec_all(D, "hi", c, 100, "nrf_mcc"), C.vec_all(D, "bx", c, 100, "nrf_mcc"))))


if __name__ == "__main__":
    main()
