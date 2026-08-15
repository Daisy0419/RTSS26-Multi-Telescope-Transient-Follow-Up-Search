import math
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.coordinates import ICRS
from astropy_healpix import HEALPix
from regions import (
    CircleSkyRegion,
    PolygonSkyRegion,
    RectangleSkyRegion,
    Region,
    Regions,
)

def generate_rectangular_tiling_equator_edge(
    width: u.Quantity, height: u.Quantity,
    stagger: bool = True, pole_eps: float = 1e-4
):
    if not width.unit.is_equivalent(u.deg) or not height.unit.is_equivalent(u.deg):
        raise ValueError("width/height must be angular (e.g., 10*u.deg).")

    w_step = width
    h_step = height

    eps = pole_eps * u.deg
    dec_min_center = (-90*u.deg) + (height/2) - eps
    dec_max_center = ( +90*u.deg) - (height/2) + eps

    dec_centers = []
    d = dec_min_center
    while d <= dec_max_center + (1e-12 * u.deg):
        dec_centers.append(d)
        d = d + h_step

    tile_indices, tile_coords = [], []
    tid = 0

    for j, dec in enumerate(dec_centers):
        # equator-ward edge latitude for this row
        if dec.to_value(u.deg) >= 0:
            dec_edge = dec - height/2
        else:
            dec_edge = dec + height/2

        # keep away from exact poles for numerical stability
        dec_edge = np.clip(dec_edge.to_value(u.deg), -90 + pole_eps, 90 - pole_eps) * u.deg

        cosd = max(1e-9, math.cos(dec_edge.to_value(u.rad)))
        ra_step = (w_step / cosd)

        ra_offset = (0.5 * ra_step) if (stagger and (j % 2 == 1)) else (0.0 * u.deg)

        if ra_step.to_value(u.deg) >= 360.0 - 1e-9:
            ra_centers = [ra_offset % (360.0 * u.deg)]
        else:
            num_ra = int(np.ceil((360.0 * u.deg / ra_step).to_value(u.dimensionless_unscaled)))
            base = np.linspace(0.0, 360.0, num_ra, endpoint=False) * u.deg
            ra_centers = (base + ra_offset) % (360.0 * u.deg)

        for ra in ra_centers:
            tile_indices.append(tid)
            tile_coords.append(SkyCoord(ra, dec, frame=ICRS()))
            tid += 1

    return tile_indices, tile_coords


# Function to generate tiling covers sphere
def generate_rectangular_tiling(width: u.Quantity, height: u.Quantity,
                                stagger: bool = True, pole_eps: float = 1e-4):
    """
    Non-overlapping (butt-to-butt) RA/Dec-aligned centers.
    """
    if not width.unit.is_equivalent(u.deg) or not height.unit.is_equivalent(u.deg):
        raise ValueError("width/height must be angular (e.g., 10*u.deg).")

    w_step = width
    h_step = height

    eps = pole_eps * u.deg
    dec_min_center = (-90*u.deg) + (height/2) - eps
    dec_max_center = ( +90*u.deg) - (height/2) + eps

    # Build Dec centers robustly
    dec_centers = []
    d = dec_min_center
    while d <= dec_max_center + (1e-12 * u.deg):
        dec_centers.append(d)
        d = d + h_step

    tile_indices, tile_coords = [], []
    tid = 0
    for j, dec in enumerate(dec_centers):
        cosd = max(1e-9, math.cos(dec.to_value(u.rad)))
        ra_step = (w_step / cosd)  # Quantity in deg
        ra_offset = (0.5*ra_step) if (stagger and (j % 2 == 1)) else (0.0*u.deg)

        # one-tile row if step ≥ 360°
        if ra_step.to_value(u.deg) >= 360.0 - 1e-9:
            ra_centers = [ra_offset % (360.0*u.deg)]
        else:
            num_ra = int(np.ceil((360.0*u.deg / ra_step).to_value(u.dimensionless_unscaled)))
            base = np.linspace(0.0, 360.0, num_ra, endpoint=False) * u.deg
            ra_centers = (base + ra_offset) % (360.0*u.deg)

        for ra in ra_centers:
            tile_indices.append(tid)
            tile_coords.append(SkyCoord(ra, dec, frame=ICRS()))
            tid += 1

    return tile_indices, tile_coords


def generate_rectangular_tiling_overlap(
    width: u.Quantity,
    height: u.Quantity,
    overlap: float = 0.12,
    stagger: bool = True,
    pole_eps: float = 1e-5,  
):
    """
    Create RA/Dec grid of centers for RA/Dec-aligned rectangular FoVs (north-up).

    Parameters
    ----------
    width, height : astropy.units.Quantity
        FoV angular size (e.g., 10*u.deg, 6*u.deg).
    overlap : float
        Fractional overlap between neighboring FoVs in both RA and Dec (0..1).
    stagger : bool
        If True, apply brick pattern (alternate rows shifted by half an RA step).
    pole_eps : float
        Extra margin (deg) to keep centers away from the exact poles.

    Returns
    -------
    tile_indices : list[int]
    tile_coords  : list[SkyCoord]  (ICRS)
    """
    if not width.unit.is_equivalent(u.deg) or not height.unit.is_equivalent(u.deg):
        raise ValueError("width and height must be angle quantities (e.g., u.deg).")
    if not (0.0 <= overlap < 1.0):
        raise ValueError("overlap must be in [0, 1).")

    # Effective step (butt-to-butt minus overlap)
    w_step = width  * (1.0 - overlap)
    h_step = height * (1.0 - overlap)

    # Keep full tiles within [-90°, +90°]
    eps = pole_eps * u.deg
    dec_min_center = (-90*u.deg) + (height/2) - eps
    dec_max_center = ( +90*u.deg) - (height/2) + eps

    # Build Dec centers robustly
    dec_centers = []
    d = dec_min_center
    while d <= dec_max_center + (1e-12 * u.deg):
        dec_centers.append(d)
        d = d + h_step

    tile_indices, tile_coords = [], []
    tid = 0

    for j, dec in enumerate(dec_centers):
        cosd = max(1e-6, math.cos(dec.to_value(u.rad)))
        ra_step = (w_step / cosd) 
        # brick staggering
        ra_offset = (0.5 * ra_step) if (stagger and (j % 2 == 1)) else (0.0 * u.deg)

        if ra_step.to_value(u.deg) >= 360.0 - 1e-9:
            ra_centers = [ra_offset % (360.0*u.deg)]
        else:
            # Number of RA tiles so that spacing <= ra_step
            num_ra = int(np.ceil((360.0 * u.deg / ra_step).to_value(u.dimensionless_unscaled)))
            base = np.linspace(0.0, 360.0, num_ra, endpoint=False) * u.deg
            ra_centers = (base + ra_offset) % (360.0 * u.deg)

        for ra in ra_centers:
            tile_indices.append(tid)
            tile_coords.append(SkyCoord(ra, dec, frame='icrs'))
            tid += 1

    return tile_indices, tile_coords



def generate_healpix_pixel_centers(
    nside: int = 64,
    order: str = "nested",
):
    """
    Return all HEALPix pixel centers at the given NSIDE.

    Parameters
    ----------
    nside : int
        HEALPix NSIDE.
    order : {"nested","ring"}
        HEALPix ordering.

    Returns
    -------
    tile_ids : list[int]
        HEALPix pixel indices (0..12*nside^2-1) in the chosen ordering.
    tile_coords_list : list[SkyCoord]
        Pixel-center coordinates (ICRS).
    """
    if order not in ("nested", "ring"):
        raise ValueError("`order` must be 'nested' or 'ring'.")
    if not (isinstance(nside, int) and nside > 0):
        raise ValueError("`nside` must be a positive integer.")

    hpx = HEALPix(nside=nside, order=order, frame=ICRS())
    npix = hpx.npix

    ipix = np.arange(npix, dtype=np.int64)
    coords = hpx.healpix_to_skycoord(ipix)  # SkyCoord array

    tile_ids = ipix.tolist()
    tile_coords_list = list(coords)  # list[SkyCoord], one per pixel

    return tile_ids, tile_coords_list
