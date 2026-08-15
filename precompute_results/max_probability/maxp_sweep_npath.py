from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# Style                                                                       #
# --------------------------------------------------------------------------- #

COL_W = 3.5
DCOL_W = 7.16

PALETTE = {
    "greedy+gcp":       "#984EA3",  # purple
    "greedy+annealing": "#377EB8",  # blue
    "greedy+ilp":       "#4DAF4A",  # green
    "top_annealing":    "#FF7F00",  # orange
    "top_PathCover":    "#A65628",  # brown
    "top_ILP":          "#E41A1C",  # red
}

METHOD_ORDER = [
    # "greedy+gcp",
    "greedy+annealing",
    "greedy+ilp",
    # "top_PathCover",
    "top_annealing",
    "top_ILP",
]

METHOD_LABEL = {
    "greedy+gcp":       r"Greedy+GCP",
    "greedy+annealing": r"Greedy+SA",
    "greedy+ilp":       r"Greedy+ILP",
    "top_PathCover":    r"TOP-PathCover",
    "top_annealing":    r"TOP-SA",
    "top_ILP":          r"TOP-ILP",
}

MARKERS = {
    "greedy+gcp":       "o",
    "greedy+annealing": "s",
    "greedy+ilp":       "D",
    "top_annealing":    "^",
    "top_PathCover":    "v",
    "top_ILP":          "P",
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
    """Load the benchmark CSV, keeping both Budget and nPath."""
    keep = ["Method", "Dataset", "Budget", "nPath", "mapTiles",
            "SumProb", "SumProbBound", "TilingTime", "PlanningTime", "IsValid"]

    df = pd.read_csv(csv_path, usecols=keep, engine="python",
                     on_bad_lines="skip")

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

    return agg


# --------------------------------------------------------------------------- #
# Dataset helpers                                                             #
# --------------------------------------------------------------------------- #

def _sorted_datasets(df: pd.DataFrame) -> list[str]:
    """Datasets sorted ascending by mapTiles (problem size proxy)."""
    return (df.groupby("Dataset")["mapTiles"].first()
              .sort_values()
              .index.tolist())


def short_dataset(name: str, tiles: int | None = None) -> str:
    """'GW191219_163120_117.txt' -> '191219(117)' when tiles given."""
    s = name.replace("GW", "").replace(".txt", "")
    head = s.split("_")[0]
    return f"{head}({tiles})" if tiles is not None else head


def _resolve_dataset(df: pd.DataFrame, name: str) -> str:
    """Find the Dataset value matching `name`, tolerant of the .txt suffix."""
    available = list(df["Dataset"].unique())
    if name in available:
        return name

    def _norm(s: str) -> str:
        return s.replace("GW", "").replace(".txt", "").strip()

    target = _norm(name)
    matches = [d for d in available if _norm(d) == target]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(
            f"Map {name!r} not found. Available: "
            + ", ".join(sorted(available)))
    raise ValueError(f"Map {name!r} is ambiguous: {matches}")


def plot_detprob_vs_npath(df: pd.DataFrame,
                          dataset: str,
                          budget: int,
                          outdir: Path) -> plt.Figure:
    """Detection probability vs nPath for one map at fixed Budget.

    Solid line per method, with each ILP method additionally getting a
    dashed line (same colour, no marker) for its proven upper bound
    SumProbBound. Linear y.
    """
    dataset = _resolve_dataset(df, dataset)
    sub = df[(df["Dataset"] == dataset) & (df["Budget"] == budget)].copy()
    if sub.empty:
        raise ValueError(f"No rows for dataset={dataset}, Budget={budget}")

    tiles = int(sub["mapTiles"].iloc[0])
    tag = short_dataset(dataset)
    event = short_dataset(dataset, tiles)
    out = outdir / f"{tag}_detprob_vs_npath.png"

    x_vals = sorted(sub["nPath"].unique())

    fig, ax = plt.subplots(figsize=(COL_W, 2.6))

    for m in METHOD_ORDER:
        sm = sub[sub["Method"] == m]
        xs, ys, bs = [], [], []
        for xv in x_vals:
            row = sm[sm["nPath"] == xv]
            if not row.empty:
                xs.append(xv)
                ys.append(float(row["SumProb"].iloc[0]))
                bs.append(float(row["SumProbBound"].iloc[0]))
        if not xs:
            continue

        ax.plot(xs, ys, marker=MARKERS[m], color=PALETTE[m],
                label=METHOD_LABEL[m], linewidth=1.2, markersize=4.5)

        if m in ILP_METHODS:
            ax.plot(xs, bs, color=PALETTE[m], linestyle="--",
                    linewidth=1.0, marker="", alpha=0.9,
                    label=f"{METHOD_LABEL[m]} bound")

    if x_vals:
        ax.set_xticks(x_vals)
    ax.set_xlabel(r"Number of telescopes $n_{\mathrm{path}}$")
    ax.set_ylabel("Detection probability")
    # ax.set_title(fr"{event}, Budget = {budget}")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    ax.legend(loc="best", fontsize=7.5,
              frameon=True, framealpha=0.9, edgecolor="0.7",
              handlelength=1.8, handletextpad=0.5)

    fig.savefig(out)
    return fig


def plot_runtime_vs_npath(df: pd.DataFrame,
                          dataset: str,
                          budget: int,
                          outdir: Path) -> plt.Figure:
    """Planning time vs nPath for one map at fixed Budget.

    Solid line per method, log y. No bound lines.
    """
    dataset = _resolve_dataset(df, dataset)
    sub = df[(df["Dataset"] == dataset) & (df["Budget"] == budget)].copy()
    if sub.empty:
        raise ValueError(f"No rows for dataset={dataset}, Budget={budget}")

    tiles = int(sub["mapTiles"].iloc[0])
    tag = short_dataset(dataset)
    event = short_dataset(dataset, tiles)
    out = outdir / f"{tag}_runtime_vs_npath.png"

    x_vals = sorted(sub["nPath"].unique())

    fig, ax = plt.subplots(figsize=(COL_W, 2.6))

    for m in METHOD_ORDER:
        sm = sub[sub["Method"] == m]
        xs, ys = [], []
        for xv in x_vals:
            row = sm[sm["nPath"] == xv]
            if not row.empty:
                xs.append(xv)
                ys.append(float(row["PlanningTime"].iloc[0]))
        if not xs:
            continue

        ax.plot(xs, ys, marker=MARKERS[m], color=PALETTE[m],
                label=METHOD_LABEL[m], linewidth=1.2, markersize=4.5)

    ax.set_yscale("log")
    if x_vals:
        ax.set_xticks(x_vals)
    ax.set_xlabel(r"Number of telescopes $n_{\mathrm{path}}$")
    ax.set_ylabel("Planning time (s)")
    # ax.set_title(fr"{event}, Budget = {budget}")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    ax.legend(loc="best", fontsize=7.5,
              frameon=True, framealpha=0.9, edgecolor="0.7",
              handlelength=1.8, handletextpad=0.5)

    fig.savefig(out)
    return fig

if __name__ == "__main__":
    csv_path = Path("maxp_large_path.csv")
    outdir = Path("figs")

    fixed_budget = 100

    map = "GW191113_071753_952.txt"
    # map = None

    set_rtss_style()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load(csv_path)

    maps = [map] if map is not None else _sorted_datasets(df)
    figs: list[plt.Figure] = []
    for dataset in maps:
        figs.append(plot_detprob_vs_npath(df, dataset, fixed_budget, outdir))
        figs.append(plot_runtime_vs_npath(df, dataset, fixed_budget, outdir))

    print(f"Wrote {len(figs)} figures to {outdir}/")
    plt.show()