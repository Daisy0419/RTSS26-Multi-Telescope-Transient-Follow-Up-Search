import math
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
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
    footprint,
    footprint_healpix_inner,
    rectangle_to_polygon,
    skycoord_to_healpy_vec,
    skycoord_to_offset,
)


def plot_coverage_with_boundaries_plotly(
    hpx: HEALPix,
    tile_to_healpix_map: dict,
    tile_coords: list[SkyCoord],
    base_region: Region,
    output_html: str = 'coverage_map_with_boundaries.html'
):
    """
    Generate an interactive Plotly Mollweide projection plot of the HEALPix coverage map
    with tile boundaries overlaid.
    
    Parameters
    ----------
    hpx : HEALPix
        The HEALPix object.
    tile_to_healpix_map : dict
        Mapping of tile_id to {pixel: fraction}.
    tile_coords : list[SkyCoord]
        List of tile central coordinates.
    base_region : Region
        The base FOV region (e.g., RectangleSkyRegion).
    output_html : str, optional
        Path to save the HTML file (default: 'coverage_map_with_boundaries.html').
    """
    import healpy as hp
    import numpy as np
    import plotly.graph_objects as go
    from astropy.coordinates import SkyCoord, CartesianRepresentation
    from astropy import units as u
    
    # Create coverage map (sum of fractions)
    npix = hp.nside2npix(hpx.nside)
    coverage_map = np.zeros(npix)
    for tile_dict in tile_to_healpix_map.values():
        for pix, frac in tile_dict.items():
            coverage_map[pix] += frac
    
    # Function to create skew-symmetric matrix
    def skew(axis):
        ax, ay, az = axis
        return np.array([
            [0, -az, ay],
            [az, 0, -ax],
            [-ay, ax, 0]
        ])
    
    # Function to interpolate points along great circle
    def interpolate_great_circle(coord1, coord2, num_points=20):
        vec1 = coord1.cartesian.xyz.value
        vec2 = coord2.cartesian.xyz.value
        dot = np.dot(vec1, vec2)
        angle = np.arccos(np.clip(dot, -1, 1))
        if angle < 1e-6:
            return [coord1]
        axis = np.cross(vec1, vec2)
        norm = np.linalg.norm(axis)
        if norm < 1e-10:
            return [coord1]
        axis /= norm
        thetas = np.linspace(0, angle, num_points)
        points = []
        for theta in thetas:
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)
            rot_matrix = cos_theta * np.eye(3) + sin_theta * skew(axis) + (1 - cos_theta) * np.outer(axis, axis)
            vec = np.dot(rot_matrix, vec1)
            points.append(SkyCoord(CartesianRepresentation(vec * u.dimensionless_unscaled), frame='icrs'))
        return points
    
    # Prepare HEALPix pixel data for Plotly
    lon_list, lat_list = [], []
    for ipix in range(npix):
        # if coverage_map[ipix] > 0:
        ra, dec = hp.pix2ang(hpx.nside, ipix, nest=(hpx.order == 'nested'), lonlat=True)
        lon_list.append(ra - 180)
        lat_list.append(dec)
    
    # Create Plotly figure
    fig = go.Figure()
    
    # Add HEALPix pixels as scatter
    fig.add_trace(go.Scattergeo(
        lon=lon_list,
        lat=lat_list,
        mode='markers',
        marker=dict(
            size=3, 
            color='green', 
            showscale=False
        ),
        name='HEALpixels'
    ))
    
    # Add tile boundaries (no skip to show all tiles; interactive plot allows zoom)
    for tile_id, coord in enumerate(tile_coords):
        tile_footprint = footprint(base_region, coord)
        if isinstance(tile_footprint, PolygonSkyRegion):
            vertices = tile_footprint.vertices
            # Close loop
            vertices = SkyCoord(np.append(vertices.ra, vertices.ra[0]), np.append(vertices.dec, vertices.dec[0]))
            
            lon_all, lat_all = [], []
            for i in range(len(vertices) - 1):
                interp_coords = interpolate_great_circle(vertices[i], vertices[i+1])
                for c in interp_coords:
                    ra_deg = c.ra.deg
                    dec_deg = c.dec.deg
                    lon_all.append(ra_deg - 180)
                    lat_all.append(dec_deg)
            
            fig.add_trace(go.Scattergeo(
                lon=lon_all,
                lat=lat_all,
                mode='lines',
                line=dict(color='red', width=1),
                name=f'Tile {tile_id} Boundary'
            ))
    
    # Update layout for Mollweide projection
    fig.update_layout(
        geo=dict(
            projection_type='mollweide',
            showland=False,
            showcountries=False,
            showocean=False,
            bgcolor='rgba(0,0,0,0)',
            lataxis_showgrid=True,
            lonaxis_showgrid=True,
            lataxis_gridcolor='gray',
            lonaxis_gridcolor='gray'
        ),
        title='Tile Coverage Map with Boundaries',
        height=800,
        width=1500
    )
    
    fig.write_html(output_html)
    print(f'Interactive plot saved as {output_html}')
    fig.show()


def plot_coverage_with_boundaries_plotly_globe(
    hpx: HEALPix,
    tile_to_healpix_map: dict,
    tile_coords: list[SkyCoord],
    base_region: Region,
    output_html: str = 'coverage_globe.html',
    pixel_size: int = 2,
):
    import healpy as hp
    import numpy as np
    import plotly.graph_objects as go
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from regions import PolygonSkyRegion

    # --- coverage ---
    npix = hp.nside2npix(hpx.nside)
    cov = np.zeros(npix)
    for tile_dict in tile_to_healpix_map.values():
        for ip, frac in tile_dict.items():
            cov[ip] += frac

    # covered pixel centers → 3D
    nest = (hpx.order == 'nested')
    covered_idx = np.nonzero(cov > 0)[0]
    if covered_idx.size:
        theta, phi = hp.pix2ang(hpx.nside, covered_idx, nest=nest)  # theta=colat, phi=lon (rad)
        # Convert to Cartesian on unit sphere (HEALPix uses theta=0 at North pole)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
    else:
        x = y = z = np.array([])

    fig = go.Figure()

    # sphere wireframe (for context)
    longitudes = np.linspace(-np.pi, np.pi, 181)
    latitudes  = np.linspace(-np.pi/2, np.pi/2, 91)
    for lat in latitudes:
        xs = np.cos(longitudes) * np.cos(lat)
        ys = np.sin(longitudes) * np.cos(lat)
        zs = np.full_like(longitudes, np.sin(lat))
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                                   line=dict(color='lightgray', width=1),
                                   hoverinfo='skip', showlegend=False))
    for lon in longitudes[::15]:
        xs = np.cos(lon) * np.cos(latitudes)
        ys = np.sin(lon) * np.cos(latitudes)
        zs = np.sin(latitudes)
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                                   line=dict(color='lightgray', width=1),
                                   hoverinfo='skip', showlegend=False))

    # covered pixels as points
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(size=pixel_size, color='green', opacity=0.9),
        name='Covered HEALPix pixels',
        hoverinfo='skip'
    ))

    # helper: great-circle interpolation in 3D
    def interp_gc_3d(c1: SkyCoord, c2: SkyCoord, n=24):
        v1 = c1.cartesian.xyz.value
        v2 = c2.cartesian.xyz.value
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)
        dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
        ang = np.arccos(dot)
        if ang < 1e-8:
            return v1[None, :]
        # slerp
        t = np.linspace(0, 1, n)[:, None]
        s1 = np.sin((1-t)*ang) / np.sin(ang)
        s2 = np.sin(t*ang) / np.sin(ang)
        v = s1 * v1 + s2 * v2
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        return v  # (n,3)

    # tile boundaries as 3D lines
    for i, center in enumerate(tile_coords):
        poly = footprint(base_region, center)
        if isinstance(poly, PolygonSkyRegion):
            verts = poly.vertices
            verts = SkyCoord(
                ra=np.append(verts.ra.deg, verts.ra.deg[0])*u.deg,
                dec=np.append(verts.dec.deg, verts.dec.deg[0])*u.deg,
                frame="icrs"
            )
            X, Y, Z = [], [], []
            for k in range(len(verts)-1):
                pts = interp_gc_3d(verts[k], verts[k+1], n=24)
                X.extend(pts[:,0]); Y.extend(pts[:,1]); Z.extend(pts[:,2])
                # segment break (None) isn’t supported in 3D; just continues
            fig.add_trace(go.Scatter3d(
                x=np.array(X), y=np.array(Y), z=np.array(Z),
                mode='lines',
                line=dict(color='red', width=4),
                name=f'Tile {i} boundary',
                hoverinfo='skip',
                showlegend=False
            ))

    fig.update_layout(
        title='Coverage with tile boundaries (rotatable globe)',
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode='data',  # true sphere
            bgcolor='white'
        ),
        height=900, width=900,
        showlegend=True
    )
    fig.write_html(output_html)
    print(f"Saved: {output_html}")
    fig.show()


def plot_coverage_healpix_with_boundaries_plotly(
    hpx: HEALPix,
    tile_to_healpix_map: dict,
    tile_coords: list[SkyCoord],
    base_region: Region,
    output_html: str = 'coverage_map_with_boundaries.html'
):
    """
    Generate an interactive Plotly Mollweide projection plot of the HEALPix coverage map
    with tile boundaries overlaid.
    
    Parameters
    ----------
    hpx : HEALPix
        The HEALPix object.
    tile_to_healpix_map : dict
        Mapping of tile_id to {pixel: fraction}.
    tile_coords : list[SkyCoord]
        List of tile central coordinates.
    base_region : Region
        The base FOV region (e.g., RectangleSkyRegion).
    output_html : str, optional
        Path to save the HTML file (default: 'coverage_map_with_boundaries.html').
    """
    import healpy as hp
    import numpy as np
    import plotly.graph_objects as go
    from astropy.coordinates import SkyCoord, CartesianRepresentation
    from astropy import units as u
    
    # Create coverage map (sum of fractions)
    npix = hp.nside2npix(hpx.nside)
    coverage_map = np.zeros(npix)
    for tile_dict in tile_to_healpix_map.values():
        for pix, frac in tile_dict.items():
            coverage_map[pix] += frac
    
    # Function to create skew-symmetric matrix
    def skew(axis):
        ax, ay, az = axis
        return np.array([
            [0, -az, ay],
            [az, 0, -ax],
            [-ay, ax, 0]
        ])
    
    # Function to interpolate points along great circle
    def interpolate_great_circle(coord1, coord2, num_points=20):
        vec1 = coord1.cartesian.xyz.value
        vec2 = coord2.cartesian.xyz.value
        dot = np.dot(vec1, vec2)
        angle = np.arccos(np.clip(dot, -1, 1))
        if angle < 1e-6:
            return [coord1]
        axis = np.cross(vec1, vec2)
        norm = np.linalg.norm(axis)
        if norm < 1e-10:
            return [coord1]
        axis /= norm
        thetas = np.linspace(0, angle, num_points)
        points = []
        for theta in thetas:
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)
            rot_matrix = cos_theta * np.eye(3) + sin_theta * skew(axis) + (1 - cos_theta) * np.outer(axis, axis)
            vec = np.dot(rot_matrix, vec1)
            points.append(SkyCoord(CartesianRepresentation(vec * u.dimensionless_unscaled), frame='icrs'))
        return points
    
    # Prepare HEALPix pixel data for Plotly (only covered pixels, uniform color)
    lon_list, lat_list = [], []
    for ipix in range(npix):
        if coverage_map[ipix] > 0:
            ra, dec = hp.pix2ang(hpx.nside, ipix, nest=(hpx.order == 'nested'), lonlat=True)
            lon_list.append(ra - 180)
            lat_list.append(dec)
    
    # Create Plotly figure
    fig = go.Figure()
    
    # Add HEALPix pixel centers as dots (uniform color, no colorbar)
    fig.add_trace(go.Scattergeo(
        lon=lon_list,
        lat=lat_list,
        mode='markers',
        marker=dict(
            size=3,  # Adjust for density; larger for lower nside
            color='green',  # Uniform color
            showscale=False
        ),
        name='HEALpixel Centers'
    ))
    
    # Add HEALPix pixel boundaries (for all covered pixels)
    for ipix in range(npix):
        if coverage_map[ipix] > 0:
            # Get pixel boundaries in cartesian
            vectors = hp.boundaries(hpx.nside, ipix, nest=(hpx.order == 'nested'), step=1)
            # Convert to theta, phi
            theta, phi = hp.vec2ang(vectors.T)
            # Convert to SkyCoord (theta colat, phi lon)
            corners = SkyCoord(phi * u.rad, (np.pi/2 - theta) * u.rad, frame='icrs')
            # Close loop
            corners = SkyCoord(np.append(corners.ra, corners.ra[0]), np.append(corners.dec, corners.dec[0]))
            
            lon_all, lat_all = [], []
            for i in range(len(corners) - 1):
                interp_coords = interpolate_great_circle(corners[i], corners[i+1], num_points=5)  # Fewer points for boundaries
                for c in interp_coords:
                    ra_deg = c.ra.deg
                    dec_deg = c.dec.deg
                    lon_all.append(ra_deg - 180)
                    lat_all.append(dec_deg)
            
            fig.add_trace(go.Scattergeo(
                lon=lon_all,
                lat=lat_all,
                mode='lines',
                line=dict(color='black', width=0.5),
                name=f'Pixel {ipix} Boundary',
                showlegend=False
            ))
    
    # Add tile boundaries (no skip to show all tiles; interactive plot allows zoom)
    for tile_id, coord in enumerate(tile_coords):
        tile_footprint = footprint(base_region, coord)
        if isinstance(tile_footprint, PolygonSkyRegion):
            vertices = tile_footprint.vertices
            # Close loop
            vertices = SkyCoord(np.append(vertices.ra, vertices.ra[0]), np.append(vertices.dec, vertices.dec[0]))
            
            lon_all, lat_all = [], []
            for i in range(len(vertices) - 1):
                interp_coords = interpolate_great_circle(vertices[i], vertices[i+1])
                for c in interp_coords:
                    ra_deg = c.ra.deg
                    dec_deg = c.dec.deg
                    lon_all.append(ra_deg - 180)
                    lat_all.append(dec_deg)
            
            fig.add_trace(go.Scattergeo(
                lon=lon_all,
                lat=lat_all,
                mode='lines',
                line=dict(color='red', width=1),
                name=f'Tile {tile_id} Boundary'
            ))
    
    # Update layout for Mollweide projection
    fig.update_layout(
        geo=dict(
            projection_type='mollweide',
            showland=False,
            showcountries=False,
            showocean=False,
            bgcolor='rgba(0,0,0,0)',
            lataxis_showgrid=True,
            lonaxis_showgrid=True,
            lataxis_gridcolor='gray',
            lonaxis_gridcolor='gray'
        ),
        title='Tile Coverage Map with Boundaries',
        height=600,
        width=1200
    )
    
    fig.write_html(output_html)
    print(f'Interactive plot saved as {output_html}')
    fig.show()


def plot_coverage_healpix_with_boundaries_plotly_globe(
    hpx: HEALPix,
    tile_to_healpix_map: dict[int, dict[int, float]] | dict[str, dict[int, float]],
    tile_coords: list[SkyCoord],
    base_region,
    output_html: str = 'coverage_globe.html',
    pixel_point_size: int = 2,
    draw_pixel_boundaries: bool = True,
    max_boundaries: int = 5000,
):
    """
    Interactive 3D globe (Plotly Scatter3d).
    - Green points: centers of covered HEALPix pixels
    - Light gray polylines: pixel boundaries (optional)
    - Red polylines: tile footprints
    """
    import healpy as hp
    import numpy as np
    import plotly.graph_objects as go
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from regions import PolygonSkyRegion
    # `footprint` must be in scope

    npix = hp.nside2npix(hpx.nside)
    nest = (hpx.order == 'nested')

    cov = np.zeros(npix, dtype=float)
    for tile_dict in tile_to_healpix_map.values():
        for ipix, frac in tile_dict.items():
            cov[int(ipix)] += float(frac)
    covered = np.nonzero(cov > 0)[0]

    # covered pixel centers to Cartesian
    if covered.size:
        theta, phi = hp.pix2ang(hpx.nside, covered, nest=nest)  # theta=colatitude, phi=lon
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
    else:
        x = y = z = np.array([])

    fig = go.Figure()

    # add a sphere grid for context
    longs = np.linspace(-np.pi, np.pi, 181)
    lats  = np.linspace(-np.pi/2, np.pi/2, 91)
    for lat in lats:
        xs = np.cos(longs) * np.cos(lat)
        ys = np.sin(longs) * np.cos(lat)
        zs = np.full_like(longs, np.sin(lat))
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                                   line=dict(color='lightgray', width=1),
                                   hoverinfo='skip', showlegend=False))
    for lon in longs[::15]:
        xs = np.cos(lon) * np.cos(lats)
        ys = np.sin(lon) * np.cos(lats)
        zs = np.sin(lats)
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                                   line=dict(color='lightgray', width=1),
                                   hoverinfo='skip', showlegend=False))

    # covered centers
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z, mode='markers',
        marker=dict(size=pixel_point_size, color='green', opacity=0.9),
        name='Covered pixels', hoverinfo='skip'
    ))

    # helper: spherical linear interpolation (slerp) between two SkyCoord
    def slerp_unit(c1: SkyCoord, c2: SkyCoord, n=24):
        v1 = c1.cartesian.xyz.value; v1 /= np.linalg.norm(v1)
        v2 = c2.cartesian.xyz.value; v2 /= np.linalg.norm(v2)
        dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
        ang = np.arccos(dot)
        if ang < 1e-12:
            return v1[None, :]
        t = np.linspace(0, 1, n)[:, None]
        s1 = np.sin((1-t)*ang) / np.sin(ang)
        s2 = np.sin(t*ang) / np.sin(ang)
        v = s1*v1 + s2*v2
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        return v

    # draw pixel boundaries (optional; cap for speed)
    if draw_pixel_boundaries:
        idx = covered
        if len(idx) > max_boundaries:
            idx = covered[:max_boundaries]
        for ip in idx:
            b = hp.boundaries(hpx.nside, int(ip), step=1, nest=nest)  # (3, nv)
            v = b.T  # (nv, 3)
            # close loop
            v = np.vstack([v, v[0]])
            fig.add_trace(go.Scatter3d(
                x=v[:,0], y=v[:,1], z=v[:,2],
                mode='lines',
                line=dict(color='lightgray', width=2),
                hoverinfo='skip',
                name='Pixel boundary',
                showlegend=False
            ))

    # tile outlines in red
    for center in tile_coords:
        poly = footprint(base_region, center)
        if isinstance(poly, PolygonSkyRegion):
            v = poly.vertices
            v = SkyCoord(
                ra=np.append(v.ra.deg, v.ra.deg[0]) * u.deg,
                dec=np.append(v.dec.deg, v.dec.deg[0]) * u.deg,
                frame='icrs'
            )
            X, Y, Z = [], [], []
            for k in range(len(v) - 1):
                pts = slerp_unit(v[k], v[k+1], n=32)  # (n,3)
                X.extend(pts[:,0]); Y.extend(pts[:,1]); Z.extend(pts[:,2])
            fig.add_trace(go.Scatter3d(
                x=np.array(X), y=np.array(Y), z=np.array(Z),
                mode='lines',
                line=dict(color='red', width=6),
                name='Tile boundary',
                hoverinfo='skip',
                showlegend=False
            ))

    fig.update_layout(
        title='Coverage with tile boundaries (3D globe)',
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode='data', bgcolor='white'
        ),
        height=900, width=900,
        showlegend=True
    )
    fig.write_html(output_html)
    print(f"Saved: {output_html}")
    fig.show()


def plot_coverage_with_boundaries_plotly(
    hpx: HEALPix,
    tile_coords: list[SkyCoord] | None = None,
    base_region: "Region | None" = None,      # e.g., RectangleSkyRegion
    output_html: str = "coverage_all_pixels.html",
):
    import healpy as hp
    import numpy as np
    import plotly.graph_objects as go
    from astropy.coordinates import SkyCoord, CartesianRepresentation
    from astropy import units as u

    nside = hpx.nside
    npix  = hp.nside2npix(nside)

    # --- all pixels (no coverage mask) ---
    ra_all, dec_all = hp.pix2ang(
        nside, np.arange(npix), nest=(hpx.order == "nested"), lonlat=True
    )
    lon_all = ra_all - 180.0  # Plotly expects [-180, +180]
    lat_all = dec_all

    fig = go.Figure()

    # Plot every pixel as a light gray dot
    fig.add_trace(go.Scattergeo(
        lon=lon_all, lat=lat_all, mode="markers",
        marker=dict(size=3, color="green"),
        name=f"All pixels (NSIDE={nside})"
    ))

    # --- optional tile boundaries overlay ---
    if (tile_coords is not None) and (base_region is not None):
        # helpers: great-circle interpolation between two coords
        def _skew(axis):
            ax, ay, az = axis
            return np.array([[0, -az, ay],[az, 0, -ax],[-ay, ax, 0]])

        def _interp_gc(c1, c2, n=24):
            v1 = c1.cartesian.xyz.value
            v2 = c2.cartesian.xyz.value
            dot = np.dot(v1, v2)
            ang = np.arccos(np.clip(dot, -1.0, 1.0))
            if ang < 1e-12:
                return [c1]
            axis = np.cross(v1, v2)
            nm = np.linalg.norm(axis)
            if nm < 1e-14:
                return [c1]
            axis /= nm
            thetas = np.linspace(0, ang, n)
            pts = []
            for th in thetas:
                ct, st = np.cos(th), np.sin(th)
                R = ct*np.eye(3) + st*_skew(axis) + (1-ct)*np.outer(axis, axis)
                vec = R @ v1
                pts.append(SkyCoord(CartesianRepresentation(vec*u.dimensionless_unscaled), frame="icrs"))
            return pts

        for tid, ctr in enumerate(tile_coords):
            fp = footprint(base_region, ctr)
            if isinstance(fp, PolygonSkyRegion):
                verts = fp.vertices
                verts = SkyCoord(
                    np.append(verts.ra,  verts.ra[0]),
                    np.append(verts.dec, verts.dec[0])
                )
                lon_b, lat_b = [], []
                for i in range(len(verts)-1):
                    seg = _interp_gc(verts[i], verts[i+1], n=24)
                    for c in seg:
                        lon_b.append(c.ra.deg - 180.0)
                        lat_b.append(c.dec.deg)
                fig.add_trace(go.Scattergeo(
                    lon=lon_b, lat=lat_b, mode="lines",
                    line=dict(color="crimson", width=1),
                    name=f"Tile {tid} boundary",
                    showlegend=False
                ))

    fig.update_layout(
        title=f"All HEALPix pixels (NSIDE={nside}) with optional tile boundaries",
        geo=dict(
            projection_type="mollweide",
            showcoastlines=False, showland=False, showocean=False,
            lataxis_showgrid=True, lonaxis_showgrid=True,
            lataxis_gridcolor="gray", lonaxis_gridcolor="gray",
            bgcolor="rgba(0,0,0,0)",
        ),
        height=800, width=1500,
    )
    fig.write_html(output_html)
    print(f"Interactive plot saved as {output_html}")
    fig.show()
