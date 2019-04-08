import astropy.coordinates
import astropy.units as u
import numpy as np
import sunpy.coordinates.frames as frames
import pfsspy
import pfsspy.coords
import sunpy.io.fits


def compute_pfss(gong_map):
    [[br, header]] = sunpy.io.fits.read(gong_map)
    br = br - np.mean(br)
    br = np.roll(br, header['CRVAL1'] + 180, axis=1)
    nrho = 60
    rss = 2.5

    input = pfsspy.Input(br, nrho, rss)
    ssfile = gong_map.with_suffix('.ssmap')
    if ssfile.exists():
        ssmap = np.loadtxt(ssfile)
        output = None
    else:
        # Compute PFSS solution and source surface map
        output = pfsspy.pfss(input)
        br, _, _ = output.bg
        ssmap = br[:, :, -1].T
        np.savetxt(ssfile, ssmap)

    return input, ssmap, header


def trace(map_file, psp_coord, pfss_input, retrace=False):
    # Load field line
    fline_file = map_file.with_suffix('.fline')
    if fline_file.exists() and not retrace:
        fline = load_fline(fline_file)
    else:
        # Calculate field line
        psp_coord.representation_type = 'cartesian'
        coord = np.array((psp_coord.x.value, psp_coord.y.value, psp_coord.z.value))[:, 0]
        coord = coord * (2.5 - 0.001) / np.linalg.norm(coord)

        output = pfsspy.pfss(pfss_input)
        fline = output.trace(coord)
        if not retrace:
            fline_xyz = np.array([fline.x / u.m, fline.y / u.m, fline.z / u.m])
            np.savetxt(fline_file, fline_xyz)

    fline.representation_type = 'spherical'
    lon = fline.lon.value
    sinlat = np.sin(fline.lat).value
    lon, sinlat = insert_nans(lon, sinlat)

    return lon, sinlat


def load_fline(f):
    """
    Load a field line file.

    Parameters
    ----------
    f : str
        Field line file location
    """
    fline = np.loadtxt(f).T
    fline = astropy.coordinates.SkyCoord(x=fline[:, 0] * u.m,
                                         y=fline[:, 1] * u.m,
                                         z=fline[:, 2] * u.m,
                                         frame=frames.HeliographicCarrington,
                                         representation_type='cartesian')
    return fline


def insert_nans(lon, sinlat):
    insert = np.nonzero(np.abs(np.diff(lon)) > 180)[0]
    lon = np.insert(lon, insert + 1, np.nan)
    sinlat = np.insert(sinlat, insert + 1, np.nan)
    return lon, sinlat
