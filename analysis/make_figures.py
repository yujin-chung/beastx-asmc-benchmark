#!/usr/bin/env python3
"""Figures 1-6 of the article and Figs. S1-S3 of Online Resource 1.

    Figs. 1, 2, S1   error and posterior support against N: BEAST X and annealed SMC at alpha = 0.99 on 50 replicates
                     at every N; annealed SMC at alpha = 0.99999 on 50 replicates at N = 100 (filled square, with its
                     mean drawn across the panel as a reference line)
    Figs. 3, 4       CPU time, wall-clock time and error against cost for the three analyses on replicates 1-5
    Figs. S2, S3     error and posterior support against N for the three analyses on replicates 1-5
    Fig. 5           thread-scaling ladders (pipeline/tables/scaling_runs.tsv)
    Fig. 6           the M1809 alignment (pipeline/tables/emp_converge_EMP_M1809.json: every rung at alpha = 0.99;
                     emp_converge_EMPA5_M1809.json: N = 1,000 at alpha = 0.99999)

Points are means and bars one standard deviation.  Each figure is drawn at its printed size and written to
output/figures/ as PDF and as PNG (600 dpi).  No colour.
"""
import csv
import json
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

import common as C

INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d9d9d9"
plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300,
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "legend.frameon": False, "savefig.bbox": "standard", "text.hinting": "no_hinting",
})
# The figures of the article are lettered in Inter; without that typeface matplotlib's default is used and only
# the lettering differs.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = ["Inter", "DejaVu Sans"]

BX = dict(color="#000000", marker="o", ls="-", mfc="#000000", mec="#000000")
HI = dict(color="#4d4d4d", marker="s", mec="#4d4d4d")
LO = dict(color="#9a9a9a", marker="^", ls=":", mfc="#ffffff", mec="#9a9a9a")
KW = dict(lw=1.4, ms=4.0, mew=0.9, capsize=1.8, elinewidth=0.7)
H_BX = Line2D([], [], lw=1.4, ms=4.0, **BX)
H_LO = Line2D([], [], lw=1.4, ms=4.0, mew=0.9, **LO)
H_HI_REF = Line2D([], [], color=HI["color"], marker="s", ls=(0, (5, 3)), mfc=HI["color"], mec=HI["mec"], lw=0.9, ms=4.0)
H_HI = Line2D([], [], color=HI["color"], marker="s", ls="--", mfc="#ffffff", mec=HI["mec"], lw=1.4, ms=4.0, mew=0.9)
LEGEND3 = [C.NAME["bx"], C.NAME["hi"], C.NAME["lo"]]

D = C.load_all()
D5 = C.load_five()


def new_fig(nrows, ncols, size_cm, **kw):
    return plt.subplots(nrows, ncols, figsize=(size_cm[0] / 2.54, size_cm[1] / 2.54), **kw)


def finish(fig, stem, handles=None, labels=None, ncol=3, **tl):
    """Legend in a strip at the top inside the canvas; the file has exactly the size of the canvas."""
    top = 1.0
    if handles:
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=ncol, fontsize=6.8,
                   columnspacing=1.0, handlelength=2.6, handletextpad=0.5, borderaxespad=0.25)
        top = 1.0 - 0.25 / fig.get_figheight()
    kw = dict(pad=0.3, h_pad=0.7, w_pad=0.9)
    kw.update(tl)
    fig.tight_layout(rect=[0, 0, 1, top], **kw)
    os.makedirs(C.OUT_FIG, exist_ok=True)
    fig.savefig(os.path.join(C.OUT_FIG, stem + ".pdf"))
    fig.savefig(os.path.join(C.OUT_FIG, stem + ".png"), dpi=600)
    print("wrote %s  %.1f x %.1f cm" % (stem, fig.get_figwidth() * 2.54, fig.get_figheight() * 2.54))
    plt.close(fig)


def panel_axes(ax, title, logx_ns=True):
    if logx_ns:
        ax.set_xscale("log")
        ax.set_xticks(C.NS)
        ax.set_xticklabels(["100", "200", "500", "1k", "2k"])
        ax.minorticks_off()
    ax.set_title(title, loc="left", color=INK, fontsize=7.8, pad=3)
    ax.grid(True, which="major", axis="both", alpha=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---- Figs. 1, 2, S1 ------------------------------------------------------------------------------------------
def series(c, m, q, ns):
    mu = np.array([C.mean_sd(C.vec_all(D, m, c, n, q))[0] for n in ns])
    sd = np.array([C.mean_sd(C.vec_all(D, m, c, n, q))[1] for n in ns])
    return mu, np.minimum(sd, mu), np.minimum(sd, 1 - mu)      # bars truncated at 0 and 1


def facet_vs_n(q, ylabel, stem, size_cm):
    dodge = {"bx": 0.94, "hi": 1.0, "lo": 1.06}
    fig, axes = new_fig(3, 2, size_cm)
    for pi, (ax, c) in enumerate(zip(axes.ravel(), C.CONDS)):
        top = 0.0
        mu, lo, hi = series(c, "lo", q, C.NS)
        ax.errorbar([n * dodge["lo"] for n in C.NS], mu, yerr=[lo, hi], zorder=2, **LO, **KW)
        top = max(top, float((mu + hi).max()))
        mu, lo, hi = series(c, "bx", q, C.NS)
        ax.errorbar([n * dodge["bx"] for n in C.NS], mu, yerr=[lo, hi], zorder=3, **BX, **KW)
        top = max(top, float((mu + hi).max()))
        mu, lo, hi = series(c, "hi", q, [100])                    # 50 replicates, N = 100 only
        ax.axhline(float(mu[0]), color=HI["color"], ls=(0, (5, 3)), lw=0.9, zorder=1)
        ax.errorbar([100 * dodge["hi"]], mu, yerr=[lo, hi], ls="none", mfc=HI["color"], zorder=4, **HI, **KW)
        top = max(top, float((mu + hi).max()))
        panel_axes(ax, "(%s) %s" % ("abcdef"[pi], C.COND_LABEL[c]))
        ax.set_ylim(-0.035 * top, top * 1.06)
        ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top * 1.06])
    for ax in axes[2]:
        ax.set_xlabel("number of posterior trees, N (log)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    finish(fig, stem, [H_BX, H_HI_REF, H_LO], LEGEND3)


# ---- Figs. 3, 4, S2, S3 (replicates 1-5; three series) ---------------------------------------------------------
S5 = {"bx": dict(BX), "hi": dict(HI, ls="--", mfc="#ffffff"), "lo": dict(LO)}
DODGE5 = {"bx": 0.92, "hi": 1.0, "lo": 1.08}


def ms5(c, n, m, q):
    return C.mean_sd(C.vec5(D5, c, n, m, q))


def fig_cost(stem, size_cm):
    fig, axes = new_fig(6, 2, size_cm, sharex=True)
    for ri, c in enumerate(C.CONDS):
        for ci, (q, lab) in enumerate((("cpu_s", "CPU"), ("wall_s", "wall-clock"))):
            ax = axes[ri, ci]
            for m in ("lo", "hi", "bx"):
                mu = np.array([ms5(c, n, m, q)[0] for n in C.NS])
                sd = np.array([ms5(c, n, m, q)[1] for n in C.NS])
                ax.errorbar([n * DODGE5[m] for n in C.NS], mu, yerr=[np.minimum(sd, mu * 0.999), sd],
                            zorder=3 if m == "bx" else 2, **S5[m], **KW)
            panel_axes(ax, "(%s) %s, %s" % ("abcdefghijkl"[2 * ri + ci],
                                            C.COND_LABEL[c].replace(", relaxed clock", ", relaxed"), lab))
            ax.title.set_fontsize(7.2)
            ax.set_yscale("log")
            ax.set_ylim((10, 3e6) if q == "cpu_s" else (3, 1e5))
            ax.set_yticks([10, 1e3, 1e5])
            ax.minorticks_off()
            if ci == 0:
                ax.set_ylabel("seconds (log)", fontsize=7)
    for ax in axes[5]:
        ax.set_xlabel("number of posterior trees, N (log)")
    finish(fig, stem, [H_BX, H_HI, H_LO], LEGEND3, h_pad=0.45)


def facet5(q, ylabel, stem, size_cm):
    fig, axes = new_fig(3, 2, size_cm)
    for pi, (ax, c) in enumerate(zip(axes.ravel(), C.CONDS)):
        top = 0.0
        for m in ("lo", "hi", "bx"):
            mu = np.array([ms5(c, n, m, q)[0] for n in C.NS])
            sd = np.array([ms5(c, n, m, q)[1] for n in C.NS])
            lo, hi = np.minimum(sd, mu), np.minimum(sd, 1 - mu)      # bars truncated at 0 and 1
            ax.errorbar([n * DODGE5[m] for n in C.NS], mu, yerr=[lo, hi], zorder=3 if m == "bx" else 2, **S5[m], **KW)
            top = max(top, float((mu + hi).max()))
        panel_axes(ax, "(%s) %s" % ("abcdef"[pi], C.COND_LABEL[c]))
        ax.set_ylim(-0.035 * top, top * 1.06)
        ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top * 1.06])
    for ax in axes[2]:
        ax.set_xlabel("number of posterior trees, N (log)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    finish(fig, stem, [H_BX, H_HI, H_LO], LEGEND3)


def fig_acc_vs_cost(stem, size_cm):
    fig, axes = new_fig(3, 2, size_cm)
    for pi, (ax, c) in enumerate(zip(axes.ravel(), C.CONDS)):
        top = 0.0
        for m in ("lo", "hi", "bx"):
            x = np.array([ms5(c, n, m, "cpu_s")[0] for n in C.NS])
            mu = np.array([ms5(c, n, m, "nrf_mcc")[0] for n in C.NS])
            sd = np.array([ms5(c, n, m, "nrf_mcc")[1] for n in C.NS])
            ax.errorbar(x, mu, yerr=[np.minimum(sd, mu), sd], zorder=3 if m == "bx" else 2, **S5[m], **KW)
            top = max(top, float((mu + sd).max()))
        panel_axes(ax, "(%s) %s" % ("abcdef"[pi], C.COND_LABEL[c]), logx_ns=False)
        ax.set_xscale("log")
        ax.set_xlim(10, 3e6)
        ax.set_xticks([10, 1e2, 1e3, 1e4, 1e5, 1e6])
        ax.minorticks_off()
        ax.set_ylim(-0.035 * top, top * 1.06)
        ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top * 1.06])
    for ax in axes[2]:
        ax.set_xlabel("measured cost (CPU seconds, log)")
    for ax in axes[:, 0]:
        ax.set_ylabel("normalized RF, MCC tree")
    finish(fig, stem, [H_BX, H_HI, H_LO], LEGEND3)


# ---- Fig. 5 (scaling ladders) ------------------------------------------------------------------------------
def fig_scaling(stem, size_cm):
    rows = [r for r in csv.DictReader(open(C.table("scaling_runs.tsv")), delimiter="\t") if r["status"] == "OK"]
    lad = defaultdict(dict)
    for r in rows:
        lad[(r["method"], r["condition"])][int(r["threads"])] = r
    fig, (a1, a2) = new_fig(2, 1, size_cm)
    ser = [("asmc", "S4", "annealed SMC  (S4)", "--", dict(color="#6b6b6b", marker="s", mfc="#ffffff", mec="#6b6b6b")),
           ("asmc", "S3", "annealed SMC  (S3)", ":", dict(color="#6b6b6b", marker="s", mfc="#ffffff", mec="#6b6b6b")),
           ("beastx", "S4", "BEAST X  (S4)", "-", dict(color="#000000", marker="o", mfc="#000000", mec="#000000"))]
    allth = sorted({t for k in lad for t in lad[k]})
    a1.plot(allth, allth, color=MUTED, ls=(0, (2, 2)), lw=0.9, label="linear speed-up", zorder=1)
    for meth, cond, lbl, ls, st in ser:
        d = lad[(meth, cond)]
        th = sorted(d)
        sp = [float(d[th[0]]["wall_s"]) / float(d[t]["wall_s"]) for t in th]
        cpu = [float(d[t]["cpu_s"]) / float(d[th[0]]["cpu_s"]) for t in th]
        kw = dict(ls=ls, lw=1.6, ms=4.5, mew=0.9, **st)
        a1.plot(th, sp, label=lbl, zorder=3, **kw)
        a2.plot(th, cpu, label=lbl, **kw)
        a1.annotate("%.1f×" % sp[-1], (th[-1], sp[-1]), textcoords="offset points", xytext=(5, -2), color=INK, fontsize=7)
    for ax in (a1, a2):
        ax.set_xscale("log", base=2)
        ax.set_xticks(allth)
        ax.set_xticklabels([str(t) for t in allth])
        ax.set_xlabel("threads")
        ax.grid(True, alpha=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    a1.set_yscale("log", base=2)
    a1.set_ylabel("speed-up over one thread")
    a1.set_title("(a) strong scaling", loc="left", color=INK)
    a2.set_ylabel("CPU seconds relative to one thread")
    a2.set_title("(b) cost of that scaling", loc="left", color=INK)
    a2.axhline(1.0, color=MUTED, ls=(0, (2, 2)), lw=0.9)
    a1.legend(loc="upper left")
    finish(fig, stem, h_pad=1.2)


# ---- Fig. 6 (M1809) ----------------------------------------------------------------------------------------
def fig_ladder(stem, size_cm):
    lad = json.load(open(C.table("emp_converge_EMP_M1809.json")))["rungs"]              # alpha = 0.99, five rungs
    hi = json.load(open(C.table("emp_converge_EMPA5_M1809.json")))["rungs"]["1000"]["asmc"]      # alpha = 0.99999
    N = sorted(int(k) for k in lad)
    bx = [lad[str(n)]["beastx"] for n in N]
    lo = [lad[str(n)]["asmc"] for n in N]
    med = lambda r: float(np.median(r["cpu_s"]))
    kw = dict(lw=1.4, ms=4.2, mew=0.9)
    pt = dict(color=HI["color"], marker="s", ls="none", mfc=HI["color"], mec=HI["color"], ms=6.0, zorder=5)
    fig, axes = new_fig(2, 2, size_cm)
    ax = axes[0, 0]
    ax.plot(N, [r["dsupport_mean"] for r in bx], **BX, **kw)
    ax.plot(N, [r["dsupport_mean"] for r in lo], **LO, **kw)
    ax.plot([1000], [hi["dsupport_mean"]], **pt)
    ax.axhline(0.05, color=MUTED, ls=(0, (2, 2)), lw=0.9)
    ax.text(N[1] * 1.05, 0.062, "0.05 criterion", fontsize=6.5, color=MUTED)
    ax.set_ylim(-0.01, 0.33)
    ax.set_ylabel("mean |Δ split support|\nbetween seeds", fontsize=7.5)
    ax.set_title("(a) seed-to-seed support difference", loc="left", color=INK, fontsize=7.8)
    ax = axes[0, 1]
    ax.plot(N, [r["cons_rf_mean"] for r in bx], **BX, **kw)
    ax.plot(N, [r["cons_rf_mean"] for r in lo], **LO, **kw)
    ax.plot([1000], [hi["cons_rf_mean"]], **pt)
    ax.set_ylim(-2, 62)
    ax.set_ylabel("RF between consensus trees\nof seeds (max 112)", fontsize=7.5)
    ax.set_title("(b) consensus-tree disagreement", loc="left", color=INK, fontsize=7.8)
    ax = axes[1, 0]
    ax.plot(N, [r["logZ_spread"] for r in lo], **LO, **kw)
    ax.plot([1000], [hi["logZ_spread"]], **pt)
    ax.set_ylim(-20, 950)
    ax.set_ylabel("range of log Z estimates\nacross seeds", fontsize=7.5)
    ax.set_title("(c) marginal-likelihood spread", loc="left", color=INK, fontsize=7.8)
    ax = axes[1, 1]
    ax.plot(N, [med(r) for r in bx], **BX, **kw)
    ax.plot(N, [med(r) for r in lo], **LO, **kw)
    ax.plot([1000], [med(hi)], **pt)
    ax.set_yscale("log")
    ax.set_ylim(100, 3e6)
    ax.set_ylabel("CPU seconds per run (log)", fontsize=7.5)
    ax.set_title("(d) cost per run", loc="left", color=INK, fontsize=7.8)
    for ax in axes.ravel():
        ax.set_xscale("log")
        ax.set_xticks(N)
        ax.set_xticklabels(["1k", "5k", "10k", "20k", "50k"])
        ax.minorticks_off()
        ax.grid(True, which="major", axis="both", alpha=0.6)
        ax.set_axisbelow(True)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("trees per run, N (log)")
    finish(fig, stem, [H_BX, H_LO, Line2D([], [], **pt)],
           [C.NAME["bx"], C.NAME["lo"], "annealed SMC, α = 0.99999, N = 1,000"], ncol=3)


if __name__ == "__main__":
    facet_vs_n("nrf_mcc", "normalized RF, MCC tree", "fig1_accuracy_vs_n", (13.0, 13.98))
    facet_vs_n("p_true", "posterior probability\nof true topology", "fig2_ptrue_vs_n", (13.0, 15.70))
    fig_cost("fig3_cost_vs_n", (13.0, 14.84))
    fig_acc_vs_cost("fig4_accuracy_vs_cost", (13.0, 14.84))
    fig_scaling("fig5_scaling", (13.0, 14.84))
    fig_ladder("fig6_empirical_ladder", (13.0, 10.54))
    facet_vs_n("nrf_cons", "normalized RF, majority-rule consensus", "figS1_accuracy_vs_n_consensus", (16.0, 17.5))
    facet5("nrf_mcc", "normalized RF, MCC tree", "figS2_accuracy_vs_n_reps1-5", (14.7, 16.6))
    facet5("p_true", "posterior probability\nof true topology", "figS3_ptrue_vs_n_reps1-5", (16.0, 18.0))
