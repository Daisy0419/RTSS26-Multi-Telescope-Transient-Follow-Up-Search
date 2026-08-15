from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

COL_W = 3.5
DCOL_W = 7.16

PALETTE = {
    "greedy+gcp":       "#984EA3",  # purple
    "greedy+annealing": "#377EB8",  # blue
    "greedy+ilp":       "#4DAF4A",  # green
    "top_annealing":    "#FF7F00",  # orange
    "top_PathCover":    "#A65628",  # brown
    "top_ILP":          "#E41A1C",  # red
    "top_greedy":         "#984EA3",  # purple
}

METHOD_ORDER = [
    # "greedy+gcp",
    "greedy+annealing",
    "greedy+ilp",
    # "top_PathCover",
    "top_annealing",
    "top_ILP",
    "top_greedy",
]

METHOD_LABEL = {
    "greedy+gcp":       r"Greedy+GCP",
    "greedy+annealing": r"Greedy+SA",
    "greedy+ilp":       r"Greedy+ILP",
    "top_PathCover":    r"TOP-PathCover",
    "top_annealing":    r"TOP-SA",
    "top_ILP":          r"TOP-ILP",
    "top_greedy":       r"TilePy",
}

MARKERS = {
    "greedy+gcp":       "o",
    "greedy+annealing": "s",
    "greedy+ilp":       "D",
    "top_annealing":    "^",
    "top_PathCover":    "v",
    "top_ILP":          "P",
    "top_greedy":       ".",
}

ILP_METHODS = {"greedy+ilp", "top_ILP"}


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
        "lines.linewidth":   1.0,
        "lines.markersize":  3.5,
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
    """Load the benchmark CSV, keeping both Budget and nPath.

    CSV columns include:
      Method, Dataset, Budget, w_max, w_acc, nPath, mapTiles,
      SumProb, SumProbBound, TilingTime, PlanningTime, IsValid, path 0, ...
    """
    keep = ["Method", "Dataset", "Budget", "nPath", "mapTiles",
            "SumProb", "SumProbBound", "TilingTime", "PlanningTime", "IsValid"]

    df = pd.read_csv(csv_path, usecols=keep, engine="python",
                     on_bad_lines="skip")

    # Robust boolean parsing for IsValid
    if df["IsValid"].dtype == bool:
        is_valid = df["IsValid"]
    else:
        is_valid = df["IsValid"].astype(str).str.strip().str.lower().eq("true")
    df = df[is_valid].copy()

    for c in ["Budget", "nPath", "mapTiles"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["SumProb", "SumProbBound", "TilingTime", "PlanningTime"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["Budget", "nPath", "mapTiles",
                           "SumProb", "SumProbBound", "PlanningTime"])

    df["Budget"] = df["Budget"].astype(int)
    df["nPath"] = df["nPath"].astype(int)
    df["mapTiles"] = df["mapTiles"].astype(int)

    agg = (df.groupby(["Method", "Dataset", "Budget", "nPath"], as_index=False)
             .agg(SumProb=("SumProb", "max"),
                  SumProbBound=("SumProbBound", "max"),
                  PlanningTime=("PlanningTime", "mean"),
                  TilingTime=("TilingTime", "mean"),
                  mapTiles=("mapTiles", "first")))
    
    # mean_sumprob = (
    #     agg.groupby("Method", as_index=False)["SumProb"]
    #     .mean()
    #     .sort_values("SumProb", ascending=False)
    # )

    # print(mean_sumprob.to_string(index=False))    

    # agg["Method"] = agg["Method"].replace(METHOD_LABEL)

    # mean_sumprob = (
    #     agg.groupby("Method", as_index=False)["SumProb"]
    #     .mean()
    #     .sort_values("SumProb", ascending=False)
    # )

    # print(mean_sumprob.to_string(index=False))    

    return agg


def _line_styles() -> dict[str, str]:
    METHOD_LINESTYLES = ["-", "--", ":", "-."]
    return {m: METHOD_LINESTYLES[i % len(METHOD_LINESTYLES)]
            for i, m in enumerate(METHOD_ORDER)}


def _legend_with_ilp_bound(ax: plt.Axes, ncol: int | None = None) -> None:
    handles, labels = ax.get_legend_handles_labels()

    handles.append(Patch(facecolor="white",
                         edgecolor="black",
                         hatch="////",
                         linewidth=0.5))
    labels.append("ILP bound")

    by_label = dict(zip(labels, handles))

    ax.legend(by_label.values(),
              by_label.keys(),
              loc="upper left",
              fontsize=8,
              ncol=ncol if ncol is not None else len(by_label),
              columnspacing=0.8,
              handlelength=1.6,
              frameon=True,
              framealpha=0.9,
              edgecolor="0.7")


def _legend_methods_only(ax: plt.Axes, ncol: int | None = None) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    by_label = dict(zip(labels, handles))

    ax.legend(by_label.values(),
              by_label.keys(),
              loc="upper left",
              fontsize=8,
              ncol=ncol if ncol is not None else len(by_label),
              columnspacing=0.8,
              handlelength=1.8,
              frameon=True,
              framealpha=0.9,
              edgecolor="0.7")


def _sorted_datasets(df: pd.DataFrame) -> list[str]:
    """Datasets sorted ascending by mapTiles (problem size proxy).

    Same dataset appears across multiple nPath rows; we take the first
    mapTiles value per dataset since the column should be constant
    within a dataset.
    """
    tile_lookup = (df.groupby("Dataset")["mapTiles"].first()
                     .sort_values()
                     .index.tolist())
    return tile_lookup


def short_dataset(name: str, tiles: int | None = None) -> str:
    """'GW191219_163120_117.txt' → '191219(117)' when tiles given."""
    s = name.replace("GW", "").replace(".txt", "")
    head = s.split("_")[0]
    return f"{head}({tiles})" if tiles is not None else head

def plot_detection_prob_fixed(df: pd.DataFrame,
                              n_path: int,
                              budget: int,
                              out: Path,
                              shared_y: bool = True) -> plt.Figure:
    """Grouped bar chart of detection probability at fixed (nPath, Budget).

    One bar per (dataset, method). ILP-bound caps are overlaid as hatched
    extensions for ILP methods. Legend sits below the axes box.
    """
    sub = df[(df["nPath"] == n_path) & (df["Budget"] == budget)].copy()
    if sub.empty:
        raise ValueError(f"No rows at nPath={n_path}, Budget={budget}")

    datasets = _sorted_datasets(sub)
    tiles_by_ds = dict(sub.groupby("Dataset")["mapTiles"].first())
    n_m = len(METHOD_ORDER)

    fig, ax = plt.subplots(figsize=(DCOL_W, 2.0))

    x = np.arange(len(datasets))
    total_w = 0.84
    bar_w = total_w / n_m

    for i, m in enumerate(METHOD_ORDER):
        offsets = x - total_w / 2 + (i + 0.5) * bar_w
        color = PALETTE[m]

        vals, bounds = [], []
        for d in datasets:
            row = sub[(sub["Dataset"] == d) & (sub["Method"] == m)]
            if row.empty:
                vals.append(np.nan)
                bounds.append(np.nan)
            else:
                vals.append(float(row["SumProb"].iloc[0]))
                bounds.append(float(row["SumProbBound"].iloc[0]))

        if not np.any(np.isfinite(vals)):
            continue

        ax.bar(offsets, vals, bar_w,
               color=color, edgecolor="black", linewidth=0.3,
               label=METHOD_LABEL[m])

        if m in ILP_METHODS:
            for xo, lo, hi in zip(offsets, vals, bounds):
                if (not np.isnan(lo) and not np.isnan(hi)
                        and hi > lo + 1e-6):
                    ax.bar(xo, hi - lo, bar_w, bottom=lo,
                           facecolor="white", edgecolor=color,
                           hatch="////", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([short_dataset(d, tiles_by_ds.get(d))
                        for d in datasets],
                       rotation=25, ha="right")
    ax.set_ylabel("Detection probability")
    ax.set_xlabel("GW event (tile count)")
    # ax.set_title(fr"$n_{{\mathrm{{path}}}}={n_path}$, Budget = {budget}")

    if shared_y:
        ymax = max(sub["SumProbBound"].max(), sub["SumProb"].max())
        ymax = min(1.02, ymax * 1.05) if ymax <= 1.0 else ymax * 1.05
        ax.set_ylim(0, ymax)

    # LO, HI_PAD = 0.5, 1.05
    # tight_y = False

    # if tight_y:
    #     data_max = max(sub["SumProbBound"].max(), sub["SumProb"].max())
    #     ax.set_ylim(LO, min(HI_PAD, data_max * 1.03))
    # else:
    #     ax.set_ylim(LO, HI_PAD)

    #     ax.axhline(1.0, color="grey", linestyle="--", linewidth=0.5)
    #     ax.grid(axis="y", linestyle=":", alpha=0.6)
    #     ax.set_axisbelow(True)

    # Legend at the bottom, inside the axes box
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor="white", edgecolor=PALETTE["top_ILP"],
                         hatch="////", linewidth=0.5))
    labels.append("ILP bound")
    ax.legend(handles, labels,
              loc="lower center",
              ncol=len(handles),
              frameon=True,
              framealpha=0.9,
              edgecolor="0.7",
              columnspacing=1.2,
              handlelength=1.6,
              handletextpad=0.5)

    fig.savefig(out)
    return fig


def plot_runtime_fixed(df: pd.DataFrame,
                       n_path: int,
                       budget: int,
                       out: Path,
                       timeout: float = 3600.0) -> plt.Figure:
    """Single-panel line chart of PlanningTime at fixed (nPath, Budget).

    One line per method, log y, datasets sorted by mapTiles. Legend sits
    below the axes box.
    """
    sub = df[(df["nPath"] == n_path) & (df["Budget"] == budget)].copy()
    if sub.empty:
        raise ValueError(f"No rows at nPath={n_path}, Budget={budget}")

    datasets = _sorted_datasets(sub)
    tiles_by_ds = dict(sub.groupby("Dataset")["mapTiles"].first())
    linestyle_for_method = _line_styles()

    fig, ax = plt.subplots(figsize=(DCOL_W, 2.0))
    x = np.arange(len(datasets))

    for m in METHOD_ORDER:
        times = []
        for d in datasets:
            row = sub[(sub["Dataset"] == d) & (sub["Method"] == m)]
            times.append(float(row["PlanningTime"].iloc[0])
                         if not row.empty else np.nan)
        if not np.any(np.isfinite(times)):
            continue
        ax.plot(x, times,
                marker=MARKERS[m],
                color=PALETTE[m],
                linestyle=linestyle_for_method[m],
                linewidth=1.2,
                markersize=4.5,
                label=METHOD_LABEL[m])

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([short_dataset(d, tiles_by_ds.get(d))
                        for d in datasets],
                       rotation=25, ha="right")
    ax.set_ylabel("Planning time (s)")
    ax.set_xlabel("GW event (tile count)")
    # ax.set_title(fr"$n_{{\mathrm{{path}}}}={n_path}$, Budget = {budget}")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    ax.axhline(timeout, color="grey", linestyle="--", linewidth=0.5)
    ax.text(x[-1], timeout, " 1 h timeout", color="grey",
            fontsize=8, va="bottom", ha="right")

    # Legend at the bottom, inside the axes box
    # ax.legend(ncol=len(METHOD_ORDER),
    #           frameon=True,
    #           framealpha=0.9,)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 0.15),  # increase 0.08 to move it higher
        ncol=len(METHOD_ORDER),
        frameon=True,
        framealpha=0.9,
    )
    # ax.legend(loc="lower center",
    #           ncol=len(METHOD_ORDER),
    #           frameon=True,
    #           framealpha=0.9,
    #           edgecolor="0.7",
    #           columnspacing=1.2,
    #           handlelength=2.0,
    #           handletextpad=0.5)

    fig.savefig(out)
    return fig


if __name__ == "__main__":
    csv_path = Path("maxp_large.csv")
    csv_path = Path("maxp_small.csv")
    # csv_path = Path("maxp_6.9x6.9_tiling.csv")
    outdir = Path("figs")

    fixed_npath = 4
    fixed_budget = 100

    set_rtss_style()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load(csv_path)

    figs = [
        plot_detection_prob_fixed(
            df,
            fixed_npath,
            fixed_budget,
            outdir / f"maxp_small_detprob_npath_{fixed_npath}_budget_{fixed_budget}.png",
            shared_y=True,
        ),

        plot_runtime_fixed(
            df,
            fixed_npath,
            fixed_budget,
            outdir / f"maxp_small_runtime_npath_{fixed_npath}_budget_{fixed_budget}.png",
        ),

    ]

    print(f"Wrote {len(figs)} figures to {outdir}/")
    plt.show()