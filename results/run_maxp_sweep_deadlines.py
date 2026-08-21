#!/usr/bin/env python3
"""Run the deadline sweeps for Figures 5 and 6.

Place this script under ``results/``. Edit only the configuration block below;
the script intentionally has no argparse/command-line parameter parser for simplicity.
"""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration: edit values here
# ---------------------------------------------------------------------------

RUN_LARGE_MAP_SWEEP = True
RUN_VERY_LARGE_MAP_SWEEP = True
OVERWRITE_RESULTS = False

N_PATHS = 4
W_MAX = 10.0
W_ACC = 10.0
DWELL_ZENITH_SECONDS = 1.0
IS_DEEPSLOW = False

# Figure 5: representative 952-tile large map. Edit either deadline list to
# run a subset, for example [100, 400].
LARGE_MAP_DEADLINES = [10, 50, 100, 200, 400, 800]
LARGE_MAP_ILP_LIMIT_BY_DEADLINE = {
    10: 3600,
    50: 3600,
    100: 3600,
    200: 3600,
    400: 7200,    # two hours in the paper
    800: 36000,   # ten hours in the paper
}

# Figure 6: representative 11,678-tile very-large map.
VERY_LARGE_MAP_DEADLINES = [200, 400, 800, 1600, 3200, 6400]

# Figure 6 reports only the two SA methods, so the paper specifies no ILP limit
# for this experiment. The current main_maxp entry point nevertheless launches
# its ILP routines; this is an explicit safety cap for those additional runs.
# If your production driver can select SA-only methods, use that driver instead.
VERY_LARGE_EXTRA_ILP_TIME_LIMIT_SECONDS = 36000

# Set to 600-1800 for a uniform 10-30 minute ILP cap in a reduced evaluator
# run. Leave as None to use the paper limits above.
ILP_TIME_LIMIT_OVERRIDE_SECONDS: int | None = None


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
    deadline_seconds: int,
    time_limit_seconds: int,
) -> None:
    if ILP_TIME_LIMIT_OVERRIDE_SECONDS is not None:
        time_limit_seconds = ILP_TIME_LIMIT_OVERRIDE_SECONDS
    command = [
        str(executable),
        str(map_file),
        str(tiling_file),
        str(deadline_seconds),
        str(W_MAX),
        str(W_ACC),
        str(DWELL_ZENITH_SECONDS),
        "1" if IS_DEEPSLOW else "0",
        str(N_PATHS),
        str(output_file),
        str(time_limit_seconds),
    ]
    print("\n$", shlex.join(command), flush=True)
    subprocess.run(command, cwd=executable.parent, check=True)


def main() -> None:
    repo_root = find_repo_root(Path(__file__).parent)
    executable = repo_root / "build" / "ts_maxp"
    data_dir = repo_root / "data"
    output_dir = repo_root / "results" / "max_probability"

    # Figure 5 inputs.
    large_map = (
        data_dir / "large_maps_4.0x2.0_tiling"
        / "GW191113_071753_952.txt"
    )
    large_tiling = data_dir / "tilings" / "4.0x2.0_tiling.csv"
    large_output = output_dir / "maxp_large_map_budget.csv"

    # Figure 6 inputs. The 11,678-tile map was generated with the 1.34x0.9
    # tiling; using a 6.9x6.9 tiling here would make map tile IDs incompatible.
    very_large_map = (
        data_dir / "very_large_maps_1.34x0.9_tiling"
        / "GW200105_162426_11678.txt"
    )
    very_large_tiling = data_dir / "tilings" / "1.34x0.9_tiling.csv"
    very_large_output = output_dir / "maxp_very_large_map_budget.csv"

    require_file(executable, "ts_maxp executable")
    if (ILP_TIME_LIMIT_OVERRIDE_SECONDS is not None
            and ILP_TIME_LIMIT_OVERRIDE_SECONDS <= 0):
        raise ValueError(
            "ILP_TIME_LIMIT_OVERRIDE_SECONDS must be positive or None."
        )
    jobs: list[tuple[Path, Path, Path, list[tuple[int, int]]]] = []

    if RUN_LARGE_MAP_SWEEP:
        require_file(large_map, "representative large map")
        require_file(large_tiling, "large-map tiling")
        missing_limits = set(LARGE_MAP_DEADLINES) - set(
            LARGE_MAP_ILP_LIMIT_BY_DEADLINE
        )
        if missing_limits:
            raise ValueError(
                f"Missing large-map ILP limits for deadlines: {sorted(missing_limits)}"
            )
        large_cases = [
            (deadline, LARGE_MAP_ILP_LIMIT_BY_DEADLINE[deadline])
            for deadline in LARGE_MAP_DEADLINES
        ]
        jobs.append((large_map, large_tiling, large_output, large_cases))

    if RUN_VERY_LARGE_MAP_SWEEP:
        require_file(very_large_map, "representative very-large map")
        require_file(very_large_tiling, "very-large-map tiling")
        very_large_cases = [
            (deadline, VERY_LARGE_EXTRA_ILP_TIME_LIMIT_SECONDS)
            for deadline in VERY_LARGE_MAP_DEADLINES
        ]
        jobs.append(
            (very_large_map, very_large_tiling, very_large_output,
             very_large_cases)
        )

    if not jobs:
        raise RuntimeError("Both deadline-sweep switches are False.")

    for _, _, output_file, _ in jobs:
        prepare_output(output_file)

    run_count = 0
    for map_file, tiling_file, output_file, cases in jobs:
        for deadline, time_limit in cases:
            run_case(
                executable, map_file, tiling_file, output_file,
                deadline, time_limit,
            )
            run_count += 1

    print(f"\nCompleted {run_count} deadline runs.")
    for _, _, output_file, _ in jobs:
        print(f"Result: {output_file}")


if __name__ == "__main__":
    main()
