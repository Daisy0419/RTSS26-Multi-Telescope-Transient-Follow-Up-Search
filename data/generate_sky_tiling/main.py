import math
import csv
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
from collections import defaultdict
from astropy_healpix import HEALPix
from regions import (
    CircleSkyRegion,
    PolygonSkyRegion,
    RectangleSkyRegion,
    Region,
    Regions,
)

from M4OPT_tiling import (
    footprint_healpix,
)

from shpere_tiling import (
    generate_rectangular_tiling,
    generate_rectangular_tiling_overlap,
    generate_rectangular_tiling_equator_edge
)

from visualize_tiling import (
    plot_coverage_with_boundaries_plotly_globe,
    plot_coverage_with_boundaries_plotly,
    plot_coverage_healpix_with_boundaries_plotly
)


# Binary mapping (store arrays of pixel IDs per tile)
def build_binary_mapping(tile_indices, tile_coords, fov_region, hpx):
    """dict[tile_id] -> np.ndarray of HEALPix pixels covered (binary)."""
    tile_to_healpix_map = {}
    for tile_id, ctr in zip(tile_indices, tile_coords):
        ipix = footprint_healpix(hpx, fov_region, ctr, rotation=0*u.deg)
        tile_to_healpix_map[int(tile_id)] = np.asarray(ipix, dtype=np.int64)
    return tile_to_healpix_map

# check Multiplicity
def multiplicity_map(tile_to_healpix_map, nside):
    m = np.zeros(hp.nside2npix(nside), dtype=np.int32)
    for arr in tile_to_healpix_map.values():
        m[arr] += 1
    return m

# Invert mapping (pixel -> list of tiles that cover it)
def invert_mapping(tile_to_healpix_map):
    pix_to_fovs = defaultdict(list)
    for fov_id, arr in tile_to_healpix_map.items():
        for ipix in np.asarray(arr, dtype=np.int64):
            pix_to_fovs[int(ipix)].append(int(fov_id))
    return pix_to_fovs

def invert_mapping(tile_to_healpix_map: dict[int | str, object]) -> dict[int, list]:
    """Return pixel -> list of candidate tiles."""
    pix_to_fovs: dict[int, list] = {}
    for fov_id, val in tile_to_healpix_map.items():
        if isinstance(val, dict):
            itr = val.keys()
        else:
            itr = val
        for ipix in itr:
            ip = int(ipix)
            if ip not in pix_to_fovs:
                pix_to_fovs[ip] = []
            pix_to_fovs[ip].append(fov_id)
    return pix_to_fovs

# Nearest-center helpers (use dict[int->SkyCoord] and a vectorized center array)
def assign_pixel_owner_nearest_center(ipix, hpx, candidate_fovs, center_array, center_ids):
    """Choose the candidate tile whose center is closest to the pixel center."""
    pcoord = hpx.healpix_to_skycoord(np.array([ipix]))[0]
    # Map candidate IDs to index in center_array
    cand_idx = np.searchsorted(center_ids, np.asarray(candidate_fovs, dtype=int))
    seps = pcoord.separation(center_array[cand_idx]).to_value(u.rad)
    return int(candidate_fovs[int(np.argmin(seps))])

def fill_uncovered_pixel_nearest_center(ipix, hpx, center_array, center_ids):
    pcoord = hpx.healpix_to_skycoord(np.array([ipix]))[0]
    seps = pcoord.separation(center_array).to_value(u.rad)
    return int(center_ids[int(np.argmin(seps))])

# Make disjoint assignment (each pixel owned by exactly one tile)
def make_disjoint_assignment(hpx: HEALPix,
                             tile_to_healpix_map: dict[int, np.ndarray],
                             tile_indices: list[int],
                             tile_coords: list[SkyCoord]):
    """
    Policy:
      - If a pixel is in multiple tiles, keep the tile with the nearest center
        (good proxy for largest fractional area).
      - If a pixel is in no tiles, assign it to the globally nearest center
        (guarantees every pixel is owned by exactly one tile).
    Returns: dict[tile_id] -> np.ndarray of owned pixel IDs (disjoint cover).
    """
    # Precompute center array and an aligned ID array for fast lookups
    center_ids = np.asarray(sorted(tile_indices), dtype=int)
    id_to_coord = {int(tid): coord for tid, coord in zip(tile_indices, tile_coords)}
    center_array = SkyCoord([id_to_coord[i].ra for i in center_ids],
                            [id_to_coord[i].dec for i in center_ids], frame=ICRS())

    pix_to_fovs = invert_mapping(tile_to_healpix_map)

    owner = np.empty(hp.nside2npix(hpx.nside), dtype=int)

    # Iterate all pixels to also fill gaps
    for ipix in range(owner.size):
        cands = pix_to_fovs.get(ipix, [])
        if len(cands) == 0:
            tid = fill_uncovered_pixel_nearest_center(ipix, hpx, center_array, center_ids)
        elif len(cands) == 1:
            tid = int(cands[0])
        else:
            tid = assign_pixel_owner_nearest_center(ipix, hpx, cands, center_array, center_ids)
        owner[ipix] = tid

    # Bucketize by tile
    result = {int(tid): [] for tid in tile_indices}
    for ipix, tid in enumerate(owner):
        result[int(tid)].append(int(ipix))
    for tid in result:
        result[tid] = np.asarray(result[tid], dtype=np.int64)
    return result


def write_tile_csv(
    filename: str,
    disjoint_mapping: dict,          # {tile_id -> iterable of ipix}
    tile_ids,                        # iterable of tile_ids (str or int)
    tile_coords,                     # iterable of SkyCoord (same order as tile_ids)
):
    """
    Write a CSV with columns: ID, RA, DEC, HEALPixels

    - RA/DEC are tile center coordinates (deg)
    - HEALPixels is a space-separated list of integer ipix owned by that tile
    - `tile_ids` and `tile_coords` must be aligned in order.

    Example row:
      FOV_00012, 123.45678900, -45.67890000, 10452 10453 10454
    """

    # Build a lookup: tile_id -> SkyCoord
    # (fix: use zip(...), not zip{...})
    centers = {tid: ctr for tid, ctr in zip(tile_ids, tile_coords)}

    # Sort helper: if id looks like "FOV_00012" pull out the integer part;
    # otherwise fall back to string/int comparison
    def sort_key(tid):
        if isinstance(tid, (int, np.integer)):
            return (0, int(tid))
        s = str(tid)
        digits = "".join(ch for ch in s if ch.isdigit())
        return (1, int(digits)) if digits else (2, s)

    with open(filename, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "RA", "DEC", "HEALPixels"])
        for tid in sorted(disjoint_mapping.keys(), key=sort_key):
            ctr = centers[tid]  # assumes every tid in mapping exists in centers
            # Normalize ipix list to a 1-D int array
            ipix = np.asarray(disjoint_mapping[tid], dtype=np.int64).ravel()
            ipix_str = " ".join(str(int(x)) for x in ipix)
            w.writerow([tid, f"{ctr.ra.deg:.8f}", f"{ctr.dec.deg:.8f}", ipix_str])


def verify_disjoint_mapping(hpx: HEALPix, disjoint_mapping: dict, preview: int = 20):
    """
    Verify that each HEALPix pixel is assigned to exactly one tile.
    
    Parameters
    ----------
    hpx : HEALPix
        The HEALPix geometry used for the mapping (defines npix).
    disjoint_mapping : dict
        {tile_id -> iterable of ipix} after your make_disjoint_assignment().
    preview : int
        Max number of example indices to show for each problem list.

    Returns
    -------
    report : dict
        {
          "ok": bool,
          "npix": int,
          "total_assigned": int,
          "unique_assigned": int,
          "missing_count": int,
          "duplicate_count": int,
          "out_of_range_count": int,
          "missing_examples": list[int],
          "duplicate_examples": list[int],
          "out_of_range_examples": list[int],
        }
    """
    npix = hpx.npix

    # Flatten all assigned pixels (accept list/ndarray; ignore empty tiles)
    parts = []
    out_of_range = []
    for arr in disjoint_mapping.values():
        a = np.asarray(arr, dtype=np.int64).ravel()
        if a.size == 0:
            continue
        # catch out-of-range right away
        bad = (a < 0) | (a >= npix)
        if np.any(bad):
            out_of_range.append(a[bad])
        parts.append(a)

    if len(parts) == 0:
        flat = np.empty((0,), dtype=np.int64)
    else:
        flat = np.concatenate(parts)

    total_assigned = flat.size

    # Fast unique/duplicate/missing checks
    uniq, counts = np.unique(flat, return_counts=True)
    unique_assigned = uniq.size
    duplicate_mask = counts > 1
    duplicate_pixels = uniq[duplicate_mask]
    missing_pixels = np.setdiff1d(np.arange(npix, dtype=np.int64), uniq, assume_unique=False)

    out_of_range_arr = np.concatenate(out_of_range) if len(out_of_range) else np.empty((0,), dtype=np.int64)

    ok = (
        (out_of_range_arr.size == 0) and
        (missing_pixels.size == 0) and
        (duplicate_pixels.size == 0) and
        (unique_assigned == npix)
    )

    print(f"""ok: {bool(ok)},
    npix: {int(npix)},
    total_assigned: {int(total_assigned)},
    unique_assigned: {int(unique_assigned)},
    missing_count: {int(missing_pixels.size)},
    duplicate_count: {int(duplicate_pixels.size)},
    out_of_range_count: {int(out_of_range_arr.size)},
    missing_examples: {missing_pixels[:preview].tolist()},
    duplicate_examples: {duplicate_pixels[:preview].tolist()},
    out_of_range_examples: {out_of_range_arr[:preview].tolist()}""")


def assert_disjoint_mapping(hpx: HEALPix, disjoint_mapping: dict):
    """Raise AssertionError with a helpful message if the mapping is not disjoint & complete."""
    rep = verify_disjoint_mapping(hpx, disjoint_mapping, preview=10)
    if rep["ok"]:
        return
    msgs = []
    if rep["out_of_range_count"]:
        msgs.append(f"out-of-range: {rep['out_of_range_count']} (e.g. {rep['out_of_range_examples']})")
    if rep["missing_count"]:
        msgs.append(f"missing: {rep['missing_count']} (e.g. {rep['missing_examples']})")
    if rep["duplicate_count"]:
        msgs.append(f"duplicates: {rep['duplicate_count']} (e.g. {rep['duplicate_examples']})")
    raise AssertionError(
        "Disjoint mapping failed: " + "; ".join(msgs) +
        f". unique_assigned={rep['unique_assigned']}, npix={rep['npix']}, total_assigned={rep['total_assigned']}"
    )

if __name__ == "__main__":

    nside = 64
    order ='nested'
    width = 10 # fov witdth in degree
    height = 6 # fov height in degree

    # Define HEALPix and base region
    hpx = HEALPix(nside=nside, order=order, frame=ICRS())
    base_region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), width=width * u.deg, height=height * u.deg, angle=0 * u.deg)

    # Generate tiling
    tile_indices, tile_coords = generate_rectangular_tiling(width=width * u.deg, height=height * u.deg)
    print(f"Generated {len(tile_indices)} tiles.")

    tile_to_healpix_map = {}
    for tile_id, coord in zip(tile_indices, tile_coords):
        pixels = footprint_healpix(hpx, base_region, coord)
        tile_to_healpix_map[tile_id] = {p: 1 for p in pixels}

    disjoint_mapping = make_disjoint_assignment(hpx, tile_to_healpix_map, tile_indices, tile_coords)
    print(verify_disjoint_mapping(hpx, disjoint_mapping))
    write_tile_csv('tiling.csv', disjoint_mapping, tile_indices, tile_coords)

    # # visualize tilings
    plot_coverage_with_boundaries_plotly_globe(hpx, tile_to_healpix_map, tile_coords, base_region)

    
