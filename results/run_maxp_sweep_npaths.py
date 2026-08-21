#!/usr/bin/env python3
"""Run the telescope-count sweep for Figure 7.

Place this script under ``results/``. Edit only the configuration block below for simplicity.
"""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration: edit values here
# ---------------------------------------------------------------------------

OVERWRITE_RESULTS = False

# Edit this list to run a subset, for example [1, 4].
N_PATH_VALUES = [1, 2, 4, 8]
DEADLINE_SECONDS = 100
W_MAX = 10.0
W_ACC = 10.0
DWELL_ZENITH_SECONDS = 1.0
IS_DEEPSLOW = False

# The paper uses a one-hour limit for the telescope-count sweep. For a reduced
# evaluator run, set this to 600-1800 seconds (10-30 minutes).
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
    n_paths: int,
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
        str(n_paths),
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

    map_file = (
        data_dir / "large_maps_4.0x2.0_tiling"
        / "GW191113_071753_952.txt"
    )
    tiling_file = data_dir / "tilings" / "4.0x2.0_tiling.csv"
    output_file = output_dir / "maxp_large_map_path.csv"

    require_file(executable, "ts_maxp executable")
    require_file(map_file, "representative large map")
    require_file(tiling_file, "large-map tiling")
    if not N_PATH_VALUES or any(value <= 0 for value in N_PATH_VALUES):
        raise ValueError("N_PATH_VALUES must contain positive integers.")
    prepare_output(output_file)

    for n_paths in N_PATH_VALUES:
        run_case(executable, map_file, tiling_file, output_file, n_paths)

    print(f"\nCompleted {len(N_PATH_VALUES)} telescope-count runs.")
    print(f"Result: {output_file}")


if __name__ == "__main__":
    main()
