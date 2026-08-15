from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# --------------------------------------------------------------------------- #
# Style                                                                       #
# --------------------------------------------------------------------------- #

COL_W  = 3.5
DCOL_W = 7.16

# High-contrast palette (ColorBrewer Set1 + complements): each method sits
# on its own hue family so the active methods are unambiguous even at small
# sizes and on greyscale-ish printouts.
PALETTE = {
    "mink_unrooted_greedy":       "#984EA3",  # purple
    "mink_unrooted_ilp":          "#A65628",  # brown
    "mink_unrooted_annealing":    "#F781BF",  # pink
    "mink_rooted_greedy":         "#377EB8",  # blue
    "mink_rooted_ILP":            "#4DAF4A",  # green
    "mink_rooted_annealing":      "#E41A1C",  # red
    "mink_rooted_annealing_best": "#FF7F00",  # orange
}

METHOD_ORDER = [
    # "mink_unrooted_greedy",
    # "mink_unrooted_ilp",
    # "mink_unrooted_annealing",
    "mink_rooted_greedy",
    "mink_rooted_ILP",
    "mink_rooted_annealing",
    # "mink_rooted_annealing_best",
]

METHOD_LABEL = {
    "mink_unrooted_greedy":       r"U-Greedy",
    "mink_unrooted_ilp":          r"U-ILP",
    "mink_unrooted_annealing":    r"U-SA",
    "mink_rooted_greedy":         r"R-Greedy",
    "mink_rooted_ILP":            r"R-ILP",
    "mink_rooted_annealing":      r"R-SA",
    "mink_rooted_annealing_best": r"R-SA$^\star$",
}

MARKERS = {
    "mink_unrooted_greedy":       "o",
    "mink_unrooted_ilp":          "D",
    "mink_unrooted_annealing":    "s",
    "mink_rooted_greedy":         "v",
    "mink_rooted_ILP":            "P",
    "mink_rooted_annealing":      "^",
    "mink_rooted_annealing_best": "X",
}

ILP_METHODS = {"mink_unrooted_ilp", "mink_rooted_ILP"}

# Anomalous nPath / nPathBound ratio above which we treat the ILP row as
# a failed solve (typically a timeout that emitted a bogus incumbent).
ANOMALY_RATIO = 2.0


def set_rtss_style() -> None:
    mpl.rcParams.update({
        "font.family":       "serif",
        "font.serif":        ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset":  "stix",
        "font.size":         10,
        "axes.labelsize":    10,
        "axes.titlesize":    10,
        "legend.fontsize":   9,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "axes.linewidth":    0.6,
        "grid.linewidth":    0.4,
        "lines.linewidth":   1.1,
        "lines.markersize":  4.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size":  2.5,
        "ytick.major.size":  2.5,
        "legend.frameon":    False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.0,
        "pdf.fonttype":      42,
        "ps.fonttype":       42,
        "savefig.dpi":       300,
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.02,
    })


# --------------------------------------------------------------------------- #
# Data loading                                                                #
# --------------------------------------------------------------------------- #

def load(csv_path: Path) -> pd.DataFrame:
    """Load the benchmark CSV by column POSITION (not name), then aggregate.

    Reads the first 12 fields of each row, ignoring any trailing `path k`
    columns. Drops:
      - rows where the ILP solver reported IsValid=false
      - ILP rows where nPath > ANOMALY_RATIO * nPathBound (timeout/bug)
    """
    col_names = ["Method", "Dataset", "Budget", "w_max", "w_acc",
                 "nTiles", "CoveredTiles", "nPath", "nPathBound",
                 "TilingTime", "PlanningTime", "IsValid"]

    df = pd.read_csv(
        csv_path,
        header=None,
        skiprows=1,
        names=col_names,
        usecols=range(len(col_names)),
        engine="c",
        on_bad_lines="skip",
    )

    for c in ("Budget", "nTiles", "CoveredTiles", "nPath", "nPathBound"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("w_max", "w_acc", "TilingTime", "PlanningTime"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["Budget", "nPath", "nPathBound", "PlanningTime"])
    for c in ("Budget", "nTiles", "CoveredTiles", "nPath", "nPathBound"):
        df[c] = df[c].astype(int)

    iv = df["IsValid"]
    if iv.dtype == bool:
        is_valid = iv
    else:
        is_valid = iv.astype(str).str.strip().str.lower().eq("true")

    is_ilp = df["Method"].isin(ILP_METHODS)

    # 1. Drop ILP rows with IsValid=false
    n0 = len(df)
    df = df[~(is_ilp & ~is_valid)].copy()
    if (drop_iv := n0 - len(df)) > 0:
        print(f"[load] dropped {drop_iv} ILP rows with IsValid=false")

    # 2. Drop anomalous ILP rows (nPath >> nPathBound)
    is_ilp = df["Method"].isin(ILP_METHODS)
    anomaly = is_ilp & (df["nPathBound"] > 0) & \
              (df["nPath"] > ANOMALY_RATIO * df["nPathBound"])
    if anomaly.any():
        bad = df[anomaly][["Method", "Dataset", "Budget", "nPath", "nPathBound"]]
        print(f"[load] dropped {anomaly.sum()} anomalous ILP rows "
              f"(nPath > {ANOMALY_RATIO}x bound):")
        for _, r in bad.iterrows():
            print(f"       {r['Method']:>26}  {r['Dataset']:>30}  "
                  f"Budget={r['Budget']:<3}  nPath={r['nPath']:<3}  "
                  f"bound={r['nPathBound']}")
        df = df[~anomaly]

    if df.empty:
        raise RuntimeError(f"No rows left after filtering in {csv_path}")

    agg = (df.groupby(["Method", "Dataset", "Budget"], as_index=False)
             .agg(nPath=("nPath", "min"),
                  nPathBound=("nPathBound", "max"),
                  PlanningTime=("PlanningTime", "mean"),
                  TilingTime=("TilingTime", "mean"),
                  nTiles=("nTiles", "first")))
    return agg


_TILE_RE = re.compile(r"_(\d+)\.txt$|_(\d+)$")


def _tile_count(name: str) -> int:
    """Pull the trailing tile-count int out of a dataset filename."""
    m = _TILE_RE.search(name)
    if not m:
        return 0
    return int(m.group(1) or m.group(2))


def short_dataset(name: str) -> str:
    """Compact axis label, e.g. 'GW191105_143521_59.txt' -> '191105(59)'."""
    s = name.replace("GW", "").replace(".txt", "")
    parts = s.split("_")
    if len(parts) >= 3:
        return f"{parts[0]}({parts[-1]})"
    return s[:12]


def _sorted_datasets(df: pd.DataFrame) -> list[str]:
    return sorted(df["Dataset"].unique(), key=_tile_count)


def plot_npath_by_dataset(df: pd.DataFrame, budget: int, out: Path) -> plt.Figure:
    """Grouped bars of nPath per (dataset, method) at one Budget.

    ILP bars are split: solid base of height nPathBound + hatched cap of
    height (nPath - nPathBound). Total bar height = nPath.
    """
    sub = df[df["Budget"] == budget].copy()
    if sub.empty:
        raise ValueError(f"No rows at Budget={budget}")

    datasets = _sorted_datasets(sub)
    n_methods = len(METHOD_ORDER)

    fig, ax = plt.subplots(figsize=(DCOL_W, 2.0))
    x = np.arange(len(datasets))
    total_w = 0.84
    bar_w = total_w / n_methods

    for i, m in enumerate(METHOD_ORDER):
        offsets = x - total_w / 2 + (i + 0.5) * bar_w
        npaths, bounds = [], []
        for d in datasets:
            row = sub[(sub["Dataset"] == d) & (sub["Method"] == m)]
            if row.empty:
                npaths.append(np.nan); bounds.append(np.nan)
            else:
                npaths.append(float(row["nPath"].iloc[0]))
                bounds.append(float(row["nPathBound"].iloc[0]))

        if m in ILP_METHODS:
            ax.bar(offsets, bounds, bar_w,
                   color=PALETTE[m], edgecolor="black", linewidth=0.3,
                   label=METHOD_LABEL[m])
            gap = [max(0.0, n - b) if not (np.isnan(n) or np.isnan(b)) else 0.0
                   for n, b in zip(npaths, bounds)]
            for xo, lb, g in zip(offsets, bounds, gap):
                if g > 0:
                    ax.bar(xo, g, bar_w, bottom=lb,
                           facecolor="white", edgecolor=PALETTE[m],
                           hatch="////", linewidth=0.5)
        else:
            ax.bar(offsets, npaths, bar_w,
                   color=PALETTE[m], edgecolor="black", linewidth=0.3,
                   label=METHOD_LABEL[m])

    ax.set_xticks(x)
    ax.set_xticklabels([short_dataset(d) for d in datasets],
                       rotation=25, ha="right")
    ax.set_ylabel(r"$n_{\mathrm{path}}$")
    ax.set_xlabel("GW event (tile count)")
    # ax.set_title(fr"Budget = {budget}")
    ax.grid(axis="y", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    # Zoom y-axis: don't waste space below the smallest bar
    all_vals = sub["nPath"].astype(float).tolist()
    if all_vals:
        # ymin = max(0.0, min(all_vals)-1)
        ymin = 0.0
        ymax = max(all_vals)
        ax.set_ylim(ymin, ymax)

    # Legend at the bottom, inside the axes box
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor="white", edgecolor=PALETTE["mink_rooted_ILP"],
                         hatch="////", linewidth=0.5))
    labels.append(r"ILP gap to bound")
    ax.legend(handles, labels,
              loc="lower center",
              ncol=len(handles),
              fontsize=8.5,
              columnspacing=1.2,
              handlelength=1.6,
              handletextpad=0.5,
              frameon=True,
              framealpha=0.9,
              edgecolor="0.7")

    fig.savefig(out)
    return fig



def plot_runtime_by_dataset(df: pd.DataFrame, budget: int,
                            out: Path) -> plt.Figure:
    """Line chart of PlanningTime per (dataset, method) at one Budget.

    Datasets are ordered left-to-right by tile count, so the line slope
    is a visual proxy for how each method scales with problem size. Y is
    log because the dynamic range across methods spans 6-7 orders of
    magnitude (sub-millisecond greedy through 30-minute ILP cap).
    """
    sub = df[df["Budget"] == budget].copy()
    if sub.empty:
        raise ValueError(f"No rows at Budget={budget}")

    datasets = _sorted_datasets(sub)
    x = np.arange(len(datasets))

    fig, ax = plt.subplots(figsize=(DCOL_W, 2.0))

    for m in METHOD_ORDER:
        times = []
        for d in datasets:
            row = sub[(sub["Dataset"] == d) & (sub["Method"] == m)]
            times.append(float(row["PlanningTime"].iloc[0]) if not row.empty
                         else np.nan)
        # Skip methods that have no data at this budget
        if not np.any(np.isfinite(times)):
            continue
        ax.plot(x, times, marker=MARKERS[m], color=PALETTE[m],
                label=METHOD_LABEL[m], linewidth=1.2, markersize=4.5)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([short_dataset(d) for d in datasets],
                       rotation=25, ha="right")
    ax.set_ylabel("Planning time (s)")
    ax.set_xlabel("GW event (tile count)")
    # ax.set_title(fr"Budget = {budget}")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    ax.axhline(7200, color="grey", linestyle="--", linewidth=0.5)
    ax.text(x[-1], 7200, " 120 min ILP cap",
            color="grey", fontsize=8, va="bottom", ha="right")
    # ax.axhline(3600, color="grey", linestyle="--", linewidth=0.5)
    # ax.text(x[-1], 3600, " 60 min ILP cap",
    #         color="grey", fontsize=8, va="bottom", ha="right")

    # Legend at the bottom, inside the axes box
    ax.legend(loc="lower center", ncol=len(METHOD_ORDER), fontsize=8.5,
              columnspacing=1.2, handlelength=1.6, handletextpad=0.5,
              frameon=True, framealpha=0.9, edgecolor="0.7")

    fig.savefig(out)
    return fig


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    csv_path   = Path("mink_large.csv")
    csv_path   = Path("mink_small.csv")
    outdir     = Path("figs")
    budget_bar = 100      # which budget to feature in the per-dataset bar chart

    set_rtss_style()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load(csv_path)

    figs = [
        plot_npath_by_dataset             (df, budget_bar, outdir / f"mink_small_npath_budget_{budget_bar}.png"),
        plot_runtime_by_dataset           (df, budget_bar, outdir / f"mink_small_runtime_budget_{budget_bar}.png"),
        # plot_npath_by_dataset_all_budgets (df,             outdir / "mink_npath_by_dataset_all_budgets.png"),
        # plot_npath_by_dataset_stacked     (df,             outdir / "mink_npath_by_dataset_stacked.png"),
        # plot_runtime_by_dataset_per_budget(df,             outdir / "mink_runtime_by_dataset_per_budget.png"),
        # plot_npath_heatmap                (df,             outdir / "npath_heatmap.png"),
        # plot_runtime_vs_budget            (df,             outdir / "runtime_vs_budget.png"),
        # plot_gap_to_lb                    (df,             outdir / "gap_to_lb.png"),
    ]
    print(f"Wrote {len(figs)} figures to {outdir}/")

    plt.show()