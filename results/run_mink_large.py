#!/usr/bin/env python3
"""Run the minimum-telescope experiment over all large-map instances.

Place this script under ``results/``. Edit only the configuration block below;
the script intentionally has no argparse/command-line parameter parser.

Output:
    results/min_telescope/mink_large.csv
"""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration: edit values here
# ---------------------------------------------------------------------------

OVERWRITE_RESULTS = False

# Set to a positive integer to run only the first N sorted large maps. For
# example, use 3 for a short representative evaluation.
MAP_LIMIT: int | None = None

DEADLINE_SECONDS = 100
W_MAX = 10.0
W_ACC = 10.0
DWELL_ZENITH_SECONDS = 1.0
IS_DEEPSLOW = False

# The paper limits R-ILP to two hours on large-map instances. For a reduced
# evaluator run, set this to 600-1800 seconds (10-30 minutes).
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


def list_maps(folder: Path, limit: int | None = None) -> list[Path]:
    if not folder.is_dir():
        raise FileNotFoundError(f"Large-map directory not found: {folder}")
    maps = sorted(folder.glob("*.txt"))
    if not maps:
        raise FileNotFoundError(f"No .txt maps found in {folder}")
    if limit is not None:
        if limit <= 0:
            raise ValueError("MAP_LIMIT must be a positive integer or None.")
        maps = maps[:limit]
    return maps


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
) -> None:
    # ts_mink MAP TILE DEADLINE W_MAX W_ACC DWELL IS_DEEPSLOW OUT TIME_LIMIT
    command = [
        str(executable),
        str(map_file),
        str(tiling_file),
        str(DEADLINE_SECONDS),
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

    map_dir = data_dir / "large_maps_4.0x2.0_tiling"
    tiling_file = data_dir / "tilings" / "4.0x2.0_tiling.csv"
    output_file = (
        repo_root / "results" / "min_telescope" / "mink_large.csv"
    )

    require_file(executable, "ts_mink executable")
    require_file(tiling_file, "4.0x2.0 large-map tiling")
    maps = list_maps(map_dir, MAP_LIMIT)
    prepare_output(output_file)

    for index, map_file in enumerate(maps, start=1):
        print(f"\nLarge-map case {index}/{len(maps)}: {map_file.name}")
        run_case(executable, map_file, tiling_file, output_file)

    print(f"\nCompleted {len(maps)} large-map runs.")
    print(f"Result: {output_file}")


if __name__ == "__main__":
    main()
