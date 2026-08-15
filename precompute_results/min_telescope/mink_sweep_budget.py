from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


COL_W  = 3.5
DCOL_W = 7.16

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
      - ILP rows where IsValid=false
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


# --------------------------------------------------------------------------- #
# Dataset helpers                                                             #
# --------------------------------------------------------------------------- #

_TILE_RE = re.compile(r"_(\d+)\.txt$|_(\d+)$")


def _tile_count(name: str) -> int:
    """Pull the trailing tile-count int out of a dataset filename."""
    m = _TILE_RE.search(name)
    if not m:
        return 0
    return int(m.group(1) or m.group(2))


def short_dataset(name: str) -> str:
    """Compact label, e.g. 'GW191105_143521_59.txt' -> '191105(59)'."""
    s = name.replace("GW", "").replace(".txt", "")
    parts = s.split("_")
    if len(parts) >= 3:
        return f"{parts[0]}({parts[-1]})"
    return s[:12]


def _sorted_datasets(df: pd.DataFrame) -> list[str]:
    return sorted(df["Dataset"].unique(), key=_tile_count)


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


def plot_npath_vs_budget(df: pd.DataFrame,
                         dataset: str,
                         outdir: Path) -> plt.Figure:
    """nPath vs Budget for one map.

    Solid line per method, with each ILP method additionally getting a
    dashed line (same colour, no marker) for its proven lower bound
    nPathBound. Integer y-axis. Linear y.
    """
    dataset = _resolve_dataset(df, dataset)
    sub = df[df["Dataset"] == dataset].copy()
    if sub.empty:
        raise ValueError(f"No rows for dataset={dataset}")

    tag = short_dataset(dataset)
    out = outdir / f"{tag}_npath_vs_budget.png"

    x_vals = sorted(sub["Budget"].unique())

    fig, ax = plt.subplots(figsize=(COL_W, 2.6))

    for m in METHOD_ORDER:
        sm = sub[sub["Method"] == m]
        xs, ys, bs = [], [], []
        for xv in x_vals:
            row = sm[sm["Budget"] == xv]
            if not row.empty:
                xs.append(xv)
                ys.append(float(row["nPath"].iloc[0]))
                bs.append(float(row["nPathBound"].iloc[0]))
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
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_xlabel("Budget")
    ax.set_ylabel(r"$n_{\mathrm{path}}$")
    # ax.set_title(tag)
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    ax.legend(loc="best", fontsize=7.5,
              frameon=True, framealpha=0.9, edgecolor="0.7",
              handlelength=1.8, handletextpad=0.5)

    fig.savefig(out)
    return fig


def plot_runtime_vs_budget(df: pd.DataFrame,
                           dataset: str,
                           outdir: Path) -> plt.Figure:
    """Planning time vs Budget for one map.

    Solid line per method, log y. No bound lines.
    """
    dataset = _resolve_dataset(df, dataset)
    sub = df[df["Dataset"] == dataset].copy()
    if sub.empty:
        raise ValueError(f"No rows for dataset={dataset}")

    tag = short_dataset(dataset)
    out = outdir / f"{tag}_runtime_vs_budget.png"

    x_vals = sorted(sub["Budget"].unique())

    fig, ax = plt.subplots(figsize=(COL_W , 2.6))

    for m in METHOD_ORDER:
        sm = sub[sub["Method"] == m]
        xs, ys = [], []
        for xv in x_vals:
            row = sm[sm["Budget"] == xv]
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
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_xlabel("Budget")
    ax.set_ylabel("Planning time (s)")
    # ax.set_title(tag)
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)

    ax.legend(loc="best", fontsize=7.5,
              frameon=True, framealpha=0.9, edgecolor="0.7",
              handlelength=1.8, handletextpad=0.5)

    fig.savefig(out)
    return fig



if __name__ == "__main__":
    csv_path = Path("mink_large_budget.csv")
    outdir   = Path("figs")

    map = "GW191113_071753_952.txt"
    # map = None

    set_rtss_style()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load(csv_path)

    maps = [map] if map is not None else _sorted_datasets(df)
    figs: list[plt.Figure] = []
    for dataset in maps:
        figs.append(plot_npath_vs_budget  (df, dataset, outdir))
        figs.append(plot_runtime_vs_budget(df, dataset, outdir))

    print(f"Wrote {len(figs)} figures to {outdir}/")
    plt.show()