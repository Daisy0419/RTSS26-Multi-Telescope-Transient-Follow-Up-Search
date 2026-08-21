#!/usr/bin/env python3
"""Sweep the deadline for the paper's representative large-map instance.

Place this script under ``results/``. Edit only the configuration block below for simplicity.

Output:
    results/min_telescope/mink_large_budget.csv
"""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration: edit values here
# ---------------------------------------------------------------------------

OVERWRITE_RESULTS = False

# Figure 11 of the paper uses these five deadlines on 191113(952). Edit this
# list to run a subset, for example [100, 400].
DEADLINES_SECONDS = [50, 100, 200, 400, 800]

W_MAX = 10.0
W_ACC = 10.0
DWELL_ZENITH_SECONDS = 1.0
IS_DEEPSLOW = False

# This is a large-map experiment, so R-ILP receives two hours per deadline.
# For a reduced evaluator run, set this to 600-1800 seconds (10-30 minutes).
ILP_TIME_LIMIT_SECONDS = 7200


def find_repo_root(start: Path) -> Path:
    """Find the repository containing the build/ and data/ directories."""
    for candidate in (start.resolve(), *start.resolve().parents):
        if (candidate / "data").is_dir() and (candidate / "build").is_dir():
            return candidate
    raise FileNotFoundError(
        "Could not find the repository root containing data/ and build/."
    )


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def prepare_output(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        return
    if not OVERWRITE_RESULTS:
        raise FileExistsError(
            f"Output already exists: {path}\n"
            "Move it, remove it, or set OVERWRITE_RESULTS = True."
        )
    path.unlink()


def run_case(
    executable: Path,
    map_file: Path,
    tiling_file: Path,
    output_file: Path,
    deadline_seconds: int,
) -> None:
    # ts_mink MAP TILE DEADLINE W_MAX W_ACC DWELL IS_DEEPSLOW OUT TIME_LIMIT
    command = [
        str(executable),
        str(map_file),
        str(tiling_file),
        str(deadline_seconds),
        str(W_MAX),
        str(W_ACC),
        str(DWELL_ZENITH_SECONDS),
        "1" if IS_DEEPSLOW else "0",
        str(output_file),
        str(ILP_TIME_LIMIT_SECONDS),
    ]
    print("\n$", shlex.join(command), flush=True)
    subprocess.run(command, cwd=executable.parent, check=True)


def main() -> None:
    repo_root = find_repo_root(Path(__file__).parent)
    executable = repo_root / "build" / "ts_mink"
    data_dir = repo_root / "data"

    map_file = (
        data_dir
        / "large_maps_4.0x2.0_tiling"
        / "GW191113_071753_952.txt"
    )
    tiling_file = data_dir / "tilings" / "4.0x2.0_tiling.csv"
    output_file = (
        repo_root
        / "results"
        / "min_telescope"
        / "mink_large_budget.csv"
    )

    require_file(executable, "ts_mink executable")
    require_file(map_file, "representative 191113(952) map")
    require_file(tiling_file, "4.0x2.0 large-map tiling")
    if not DEADLINES_SECONDS:
        raise ValueError("DEADLINES_SECONDS cannot be empty")
    if any(deadline <= 0 for deadline in DEADLINES_SECONDS):
        raise ValueError("Every deadline must be positive")
    prepare_output(output_file)

    for index, deadline in enumerate(DEADLINES_SECONDS, start=1):
        print(
            f"\nDeadline case {index}/{len(DEADLINES_SECONDS)}: "
            f"D={deadline} seconds"
        )
        run_case(
            executable,
            map_file,
            tiling_file,
            output_file,
            deadline,
        )

    print(f"\nCompleted {len(DEADLINES_SECONDS)} deadline runs.")
    print(f"Result: {output_file}")


if __name__ == "__main__":
    main()
