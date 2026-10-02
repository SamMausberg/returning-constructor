#!/usr/bin/env python3
"""Render the certified data as separate column-width vector plots.

All plotted values/enclosures are computed by certified_figures.py. Lines join
integer samples for readability. See data/*.json for exact rational endpoints.
"""

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source/figures"
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.15,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }
)
BLACK = "0.05"
GRAY = "0.42"
PALE = "0.83"


def canvas(w=3.35, h=2.25):
    f, a = plt.subplots(figsize=(w, h))
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.tick_params(direction="out", length=3, width=0.6)
    return f, a


def save(f, name):
    f.savefig(OUT / f"{name}.pdf")
    f.savefig(OUT / f"{name}.png", dpi=220)
    plt.close(f)


def ratios():
    d = pd.read_csv(ROOT / "data/ratio_reservoir.csv")
    e = pd.read_csv(ROOT / "data/ratio_uses.csv")
    f, a = canvas()
    f.subplots_adjust(left=0.18, right=0.955, bottom=0.22, top=0.91)
    a.fill_between(
        d.M, d.forward_bound_lo, d.forward_bound_hi, color=PALE, label=r"$p\to\tau$: enclosure"
    )
    a.plot(d.M, d.forward_bound_lo, color=GRAY, lw=0.5)
    a.plot(d.M, d.forward_bound_hi, color=GRAY, lw=0.5)
    exact = d.dropna(subset=["forward_ratio"])
    a.plot(
        exact.M, exact.forward_ratio, "o", ms=3, color=BLACK, label=r"$p\to\tau$: certified values"
    )
    a.plot(
        d.M,
        d.reverse_ratio,
        "--",
        color=BLACK,
        marker="s",
        ms=2.5,
        markevery=4,
        label=r"$\tau\to p$: certified values",
    )
    a.set_yscale("log")
    a.set_ylim(1e-16, 3)
    a.set_xlim(0.5, 32.5)
    a.set_yticks([1, 1e-4, 1e-8, 1e-12, 1e-16])
    a.set_xticks([1, 8, 16, 24, 32])
    a.set_xlabel(r"Reservoir size $M$")
    a.set_ylabel(r"Current-error ratio $R_{M,2}$")
    a.text(0.02, 1.05, r"(a) Two uses, $\eta=\pi/4$", transform=a.transAxes)
    a.legend(
        loc="lower left",
        frameon=False,
        handlelength=1.7,
        borderpad=0,
        labelspacing=0.25,
        fontsize=7.4,
    )
    save(f, "ratio_reservoir")
    f, a = canvas()
    f.subplots_adjust(left=0.18, right=0.955, bottom=0.22, top=0.91)
    a.axhline(2, color=GRAY, lw=0.7, ls=":")
    a.plot(e.n, e.forward_ratio, color=BLACK, label=r"$p\to\tau$", marker="o", ms=2.3, markevery=6)
    a.plot(
        e.n,
        e.reverse_ratio,
        "--",
        color=GRAY,
        label=r"$\tau\to p$",
        marker="s",
        ms=2.3,
        markevery=6,
    )
    a.set_xlim(1, 48)
    a.set_ylim(0, 2.18)
    a.set_xticks([1, 12, 24, 36, 48])
    a.set_yticks([0, 0.5, 1, 1.5, 2])
    a.set_xlabel(r"Number of uses $n$")
    a.set_ylabel(r"Current-error ratio $R_{2,n}$")
    a.text(0.02, 1.05, r"(b) Two cells, $\eta=\pi/4$", transform=a.transAxes)
    a.text(28, 1.79, r"Common limit $2^{M-1}=2$", fontsize=8)
    a.legend(loc="lower right", frameon=False)
    save(f, "ratio_uses")
    # A single compact axes lets a one-page Comment show both slices.
    f, a = canvas(3.35, 1.25)
    f.subplots_adjust(left=0.15, right=0.965, bottom=0.29, top=0.86)
    a.fill_between(d.M, d.forward_bound_lo, d.forward_bound_hi, color=PALE)
    a.plot(d.M, d.forward_bound_lo, color=GRAY, lw=0.45)
    a.plot(d.M, d.forward_bound_hi, color=GRAY, lw=0.45)
    a.plot(exact.M, exact.forward_ratio, "o", color=BLACK, ms=2)
    a.plot(d.M, d.reverse_ratio, "--", color=BLACK, lw=0.9)
    e = e[e.n <= 32]
    a.plot(e.n, e.forward_ratio, color=BLACK, lw=0.9)
    a.plot(e.n, e.reverse_ratio, "--", color=GRAY, lw=0.9)
    a.set_yscale("log")
    a.set_ylim(2e-17, 20)
    a.set_yticks([1, 1e-8, 1e-16])
    a.set_xticks([1, 8, 16, 24, 32])
    a.set_xlim(1, 32)
    a.tick_params(labelsize=7.5, pad=1.5, length=2)
    a.set_ylabel(r"$R$", labelpad=2, fontsize=8)
    a.set_xlabel(r"$k$", labelpad=0, fontsize=8)
    a.text(0.36, 1.045, r"fixed reservoir $M=2$, uses $n=k$", transform=a.transAxes, fontsize=7.7)
    a.text(1.8, 3e-15, r"fixed uses $n=2$, reservoir $M=k$", fontsize=7.7)
    save(f, "comment_ratio")


def singlet():
    d = pd.read_csv(ROOT / "data/singlet.csv")
    f, a = canvas(3.35, 2.55)
    f.subplots_adjust(left=0.18, right=0.955, bottom=0.19, top=0.97)
    a.plot(d.M, d.joint_error, color=BLACK, label="Joint error", marker="o", ms=2.5, markevery=3)
    a.plot(d.M, d.singlet_lower, ":", color=GRAY, label="Singlet lower bound", lw=1.0)
    a.plot(
        d.M,
        d.marginal_error,
        "--",
        color=BLACK,
        label="Largest marginal error",
        marker="s",
        ms=2.2,
        markevery=3,
    )
    a.plot(d.M, d.prefix_return_upper, "-.", color=GRAY, label="Prefix-return upper bound")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_ylim(2e-15, 1.2)
    a.set_xlim(1, 4e7)
    a.set_xticks([1, 1e2, 1e4, 1e6, 1e8])
    a.set_xlim(1, 4e7)
    a.set_yticks([1, 1e-4, 1e-8, 1e-12])
    a.set_xlabel(r"Reservoir size $M=k^5$")
    a.set_ylabel("Trace distance")
    a.text(5e4, 0.39, r"$1/4$", fontsize=8)
    a.legend(loc="lower left", frameon=False, fontsize=7.7, handlelength=2.3, labelspacing=0.3)
    save(f, "singlet_data")


def costs():
    d = pd.read_csv(ROOT / "data/cost_horizon.csv")
    f, a = canvas()
    f.subplots_adjust(left=0.18, right=0.955, bottom=0.22, top=0.91)
    a.fill_between(
        d.n, d.reverse_lower, d.reverse_upper, color=PALE, label="Purification: lower / upper"
    )
    a.plot(d.n, d.reverse_lower, ":", color=GRAY, lw=0.9)
    a.plot(d.n, d.reverse_upper, "--", color=BLACK, lw=1)
    a.plot(d.n, d.forward, color=BLACK, label="Randomization: exact optimum")
    a.set_yscale("log")
    a.set_xlim(1, 128)
    a.set_ylim(0.9, 1800)
    a.set_xticks([1, 32, 64, 96, 128])
    a.set_yticks([1, 10, 100, 1000])
    a.set_xlabel(r"Horizon $n$")
    a.set_ylabel("Memory qubits $m$")
    a.text(0.02, 1.05, r"(a) Return tolerance $\delta=0.1$", transform=a.transAxes)
    a.legend(loc="lower right", frameon=False, fontsize=7.5)
    save(f, "cost_horizon")
    d = pd.read_csv(ROOT / "data/cost_tolerance.csv").sort_values("delta")
    f, a = canvas()
    f.subplots_adjust(left=0.18, right=0.955, bottom=0.22, top=0.91)
    a.fill_between(d.delta, d.reverse_lower, d.reverse_upper, color=PALE)
    a.plot(d.delta, d.reverse_lower, ":", color=GRAY, lw=0.9)
    a.plot(d.delta, d.reverse_upper, "--", color=BLACK, lw=1)
    a.plot(d.delta, d.forward, color=BLACK)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(0.01, 0.25)
    a.set_ylim(10, 4000)
    a.set_xticks([0.01, 0.025, 0.05, 0.1, 0.25])
    a.set_xticklabels(["0.01", "0.025", "0.05", "0.1", "0.25"])
    a.set_yticks([10, 100, 1000])
    a.set_xlabel(r"Return tolerance $\delta$")
    a.set_ylabel("Memory qubits $m$")
    a.text(0.02, 1.05, r"(b) Horizon $n=32$", transform=a.transAxes)
    a.text(0.03, 1700, r"Purification $\sim n/\delta$", fontsize=8)
    a.text(0.019, 23, r"Randomization $=\lceil n/2\rceil$", fontsize=8)
    save(f, "cost_tolerance")


def main():
    OUT.mkdir(exist_ok=True)
    ratios()
    singlet()
    costs()
    print("Created six vector plots and their inspection PNGs from certified data.")


if __name__ == "__main__":
    main()
