#!/usr/bin/env python3
"""Summarise HTSeq-count output: % of reads assigned to genes vs. each
unassigned category (__no_feature, __ambiguous, __too_low_aQual, ...).

Usage:
    python scripts/07_plot_count_summary.py results/06_counts/*.counts.tsv \
        --out summary.png --table summary.tsv
"""
import argparse
import os
from pathlib import Path

# Force a non-interactive backend before importing matplotlib: notebooks such as
# Colab export MPLBACKEND=module://matplotlib_inline..., which breaks plain Python.
os.environ["MPLBACKEND"] = "Agg"
import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

CATEGORIES = ["assigned", "no_feature", "ambiguous", "too_low_aQual",
              "not_aligned", "alignment_not_unique"]
COLORS = ["#2E75B6", "#C0392B", "#E67E22", "#27AE60", "#8E44AD", "#7F8C8D"]


def summarise(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=None, names=["feature", "count"])
    special = df["feature"].str.startswith("__")
    counts = {"assigned": df.loc[~special, "count"].sum()}
    for feat, n in df.loc[special].itertuples(index=False):
        counts[feat.lstrip("_")] = n
    s = pd.Series(counts).reindex(CATEGORIES, fill_value=0)
    return s / s.sum() * 100


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("counts", nargs="+", type=Path, help="HTSeq-count .tsv files")
    ap.add_argument("--out", type=Path, default=Path("htseq_assignment_summary.png"))
    ap.add_argument("--table", type=Path, default=None)
    args = ap.parse_args()

    pct = pd.DataFrame({p.name.split(".")[0]: summarise(p) for p in args.counts}).T
    pct = pct.loc[:, (pct > 0).any()]          # drop categories absent everywhere
    if args.table:
        pct.round(2).to_csv(args.table, sep="\t")
    print(pct.round(1).to_string())

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(pct) + 4), 5))
    colors = [COLORS[CATEGORIES.index(c)] for c in pct.columns]
    if len(pct) == 1:                          # one sample: simple bar chart
        bars = ax.bar(pct.columns, pct.iloc[0], color=colors)
        ax.bar_label(bars, fmt="%.1f%%", fontsize=9)
        ax.set_title(f"HTSeq-count read assignment — {pct.index[0]}")
    else:                                      # several samples: stacked bars
        pct.plot(kind="bar", stacked=True, color=colors, ax=ax, width=0.8)
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        ax.set_title("HTSeq-count read assignment per sample")
    ax.set_ylabel("Reads (%)")
    ax.tick_params(axis="x", rotation=30)
    for lbl in ax.get_xticklabels():
        lbl.set_horizontalalignment("right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
