import astropy.coordinates
from astropy.time import Time
import astropy.units as u
import numpy as np

import sunpy.coordinates.frames as frames
import sunpy.io.fits
import sunpy.map

import pfsspy
import pfsspy.coords


def compute_pfss(gong_fname, dtime):
    gong_map = sunpy.map.Map(gong_fname)
    br = gong_map.data
    header = gong_map.meta

    br = br - np.mean(br)
    br = np.roll(br, header['CRVAL1'] + 180, axis=1)
    header['CRVAL1'] = 180
    header['DATE_ORI'] = header['DATE']
    header['date-obs'] = Time(dtime).isot

    gong_map = sunpy.map.Map(br, header)
    nrho = 60
    rss = 2.5

    input = pfsspy.Input(gong_map, nrho, rss)
    ssfile = gong_fname.with_suffix('.ssmap')

    if ssfile.exists():
        ssmap = np.loadtxt(ssfile)
        ssmap = sunpy.map.Map(ssmap, header)
    else:
        print('Calculating PFSS solution')
        # Compute PFSS solution and source surface map
        output = pfsspy.pfss(input)
        ssdata = output.source_surface_br.data
        np.savetxt(ssfile, ssdata)
        ssmap = output.source_surface_br

    return input, ssmap


def trace(map_file, psp_coord, pfss_input, retrace=False):
    # Load field line
    fline_file = map_file.with_suffix('.fline')
    if fline_file.exists() and not retrace:
        fline = load_fline(fline_file)
    else:
        print('Tracing field line')
        # Calculate field line

        output = pfsspy.pfss(pfss_input)
        tracer = pfsspy.tracing.PythonTracer()
        fline = tracer.trace(psp_coord, output)[0].coords
        if not retrace:
            fline_xyz = np.array([fline.x / u.m, fline.y / u.m, fline.z / u.m])
            np.savetxt(fline_file, fline_xyz)

    fline.representation_type = 'spherical'
    lon = fline.lon.to_value(u.deg)
    lat = fline.lat.to_value(u.deg)
    r = fline.radius.to_value(u.m)
    lon, lat, r = insert_nans(lon, lat, r)

    fline = astropy.coordinates.SkyCoord(lon * u.deg, lat * u.deg, r * u.m,
                                         frame='heliographic_carrington',
                                         obstime=fline.obstime)

    return fline


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


def insert_nans(lon, sinlat, r):
    insert = np.nonzero(np.abs(np.diff(lon)) > 180)[0]
    lon = np.insert(lon, insert + 1, np.nan)
    sinlat = np.insert(sinlat, insert + 1, np.nan)
    r = np.insert(r, insert + 1, np.nan)
    return lon, sinlat, r
