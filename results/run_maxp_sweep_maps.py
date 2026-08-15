#!/usr/bin/env python3
"""Run the fixed-map-set experiment for Figures 2 and 3.

Place this script under ``results/``. Edit only the configuration block below;
the script intentionally has no argparse/command-line parameter parser.
"""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration: edit values here
# ---------------------------------------------------------------------------

RUN_SMALL_MAPS = True
RUN_LARGE_MAPS = True
OVERWRITE_RESULTS = False

# Set these to a positive integer to run only the first N sorted maps in each
# class (for example, set both to 3 for a short six-instance evaluation).
SMALL_MAP_LIMIT: int | None = None
LARGE_MAP_LIMIT: int | None = None

DEADLINE_SECONDS = 100
N_PATHS = 4
W_MAX = 10.0
W_ACC = 10.0
DWELL_ZENITH_SECONDS = 1.0
IS_DEEPSLOW = False

# The paper uses a one-hour ILP limit for the cross-instance comparison. For a
# reduced evaluator run, set this to 600-1800 seconds (10-30 minutes).
ILP_TIME_LIMIT_SECONDS = 3600


def find_repo_root(start: Path) -> Path:
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
        raise FileNotFoundError(f"Map directory not found: {folder}")
    maps = sorted(folder.glob("*.txt"))
    if not maps:
        raise FileNotFoundError(f"No .txt maps found in {folder}")
    if limit is not None:
        if limit <= 0:
            raise ValueError("Map limits must be positive integers or None.")
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
    command = [
        str(executable),
        str(map_file),
        str(tiling_file),
        str(DEADLINE_SECONDS),
        str(W_MAX),
        str(W_ACC),
        str(DWELL_ZENITH_SECONDS),
        "1" if IS_DEEPSLOW else "0",
        str(N_PATHS),
        str(output_file),
        str(ILP_TIME_LIMIT_SECONDS),
    ]
    print("\n$", shlex.join(command), flush=True)
    subprocess.run(command, cwd=executable.parent, check=True)


def main() -> None:
    repo_root = find_repo_root(Path(__file__).parent)
    executable = repo_root / "build" / "ts_maxp"
    data_dir = repo_root / "data"
    output_dir = repo_root / "results" / "max_probability"

    # These map directories and tilings are paired. Do not mix their tile IDs.
    small_map_dir = data_dir / "small_maps_6.9x6.9_tiling"
    small_tiling = data_dir / "tilings" / "6.9x6.9_tiling.csv"
    small_output = output_dir / "maxp_small.csv"

    large_map_dir = data_dir / "large_maps_4.0x2.0_tiling"
    large_tiling = data_dir / "tilings" / "4.0x2.0_tiling.csv"
    large_output = output_dir / "maxp_large.csv"

    require_file(executable, "ts_maxp executable")

    jobs: list[tuple[list[Path], Path, Path]] = []
    if RUN_SMALL_MAPS:
        require_file(small_tiling, "small-map tiling")
        jobs.append((
            list_maps(small_map_dir, SMALL_MAP_LIMIT),
            small_tiling,
            small_output,
        ))
    if RUN_LARGE_MAPS:
        require_file(large_tiling, "large-map tiling")
        jobs.append((
            list_maps(large_map_dir, LARGE_MAP_LIMIT),
            large_tiling,
            large_output,
        ))
    if not jobs:
        raise RuntimeError("Both RUN_SMALL_MAPS and RUN_LARGE_MAPS are False.")

    # Validate/prepare every output before starting any long-running job.
    for _, _, output_file in jobs:
        prepare_output(output_file)

    run_count = 0
    for maps, tiling_file, output_file in jobs:
        for map_file in maps:
            run_case(executable, map_file, tiling_file, output_file)
            run_count += 1

    print(f"\nCompleted {run_count} map runs.")
    for _, _, output_file in jobs:
        print(f"Result: {output_file}")


if __name__ == "__main__":
    main()
