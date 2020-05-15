from datetime import datetime, timedelta
import pathlib
import urllib.request
import urllib.error

from astropy.coordinates import SkyCoord, Longitude
from astropy.time import Time
import astropy.units as u
from astropy.wcs import WCS
from reproject import reproject_interp

import numpy as np
import matplotlib.pyplot as plt

from sunpy.net import vso
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map, make_fitswcs_header
from sunpy.coordinates import get_earth
import sunpy.sun.constants

from time_helpers import start_of_day
from synoptic_helpers import synop_header, synop_weights


map_dir = pathlib.Path('/Users/dstansby/Data/aia')
map_dir.mkdir(exist_ok=True, parents=True)


def map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_193_{datestr}.fits'


def synoptic_map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'aia_193_synoptic_{datestr}.fits'


def download_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    # This is broken for now, see https://github.com/sunpy/sunpy/issues/4159
    print(f'Fetching AIA map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(193 * u.Angstrom))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No AIA map available for {dtime}')
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))
    """
    import parfive
    dl = parfive.Downloader(max_conn=1)
    url = (f"http://jsoc2.stanford.edu/data/aia/synoptic/nrt/"
           f"{dtime.year}/{dtime.month:02}/{dtime.day:02}/"
           f"H0000/AIA{dtime.year}{dtime.month:02}{dtime.day:02}_"
           f"000000_0193.fits")
    dl.enqueue_file(url, filename=map_path(dtime))
    res = dl.download()
    if len(res.errors):
        print(res.errors)
        raise RuntimeError('Download failed')'''
    """


def load_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    print(f'Loading AIA map for {dtime}')
    try:
        ret = Map(str(mappath))
        ret.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
        return ret
    except OSError as e:
        raise RuntimeError(f'No AIA map available for {dtime}') from e


def synop_reproject(dtime, shape_out):
    synop_map_path = synoptic_map_path(dtime)
    if not synop_map_path.exists():
        m = load_start_of_day_map(dtime)
        # Reproject
        print(f'Reprojecting {synop_map_path}')
        header = synop_header(shape_out, m.date)
        with np.errstate(invalid='ignore'):
            array, footprint = reproject_interp(m, WCS(header),
                                                shape_out=shape_out)
        # Save some memory
        array = np.int16(array)
        for key in m.meta:
            if key not in header:
                header[key] = m.meta[key]
        new_map = Map((array, header))
        new_map.save(str(synop_map_path))

    print(f'Loading {synop_map_path}')
    new_map = Map(synop_map_path)
    return new_map


def create_synoptic_map(endtime, aia_maps={}):
    """
    Create an AIA synoptic map, using 25 daily AIA 193 maps ending on the
    endtime given. Note that the maps are taken from the start of each day.

    Parameters
    ----------
    endtime :
    aia_maps : dict
        A mapping of `datetime.date` to `sunpy.map.GenericMap`.

    Returns
    -------
    sunpy.map.Map : synoptic map
    """
    if endtime > datetime.now():
        endtime = datetime.now()
    shape = [720, 1440]
    data = np.zeros(shape)
    weight_sum = np.zeros(shape)
    nmaps = 23

    dtimes = [endtime - timedelta(days=i) for i in range(nmaps)[::-1]]
    # Fill up aia_maps
    for dtime in dtimes:
        if dtime.date() in aia_maps:
            continue
        aia_maps[dtime.date()] = synop_reproject(dtime, shape)

    # Add up all the reprojected maps
    for dtime in dtimes:
        aia_synop_map = aia_maps[dtime.date()]
        weights = synop_weights(aia_synop_map)

        aia_data = aia_synop_map.data * weights
        weights[np.isnan(aia_data)] = 0
        aia_data[np.isnan(aia_data)] = 0
        data += aia_data
        weight_sum += weights

    weight_sum[weight_sum == 0] = np.nan
    data /= weight_sum

    meta = aia_synop_map.meta
    meta['date-obs'] = dtime.strftime('%Y-%m-%dT%H:%M:%S')
    data = np.roll(data, data.shape[1] // 2, axis=1)
    meta['crval1'] = 180
    meta['telescop'] = 'sdo'
    meta['instrume'] = 'AIA'
    meta['detector'] = 'AIA'
    meta['waveunit'] = 'angstrom'
    meta['wavelnth'] = 193

    synop_map = Map((data, meta))
    synop_map.plot_settings = aia_synop_map.plot_settings
    # synop_map.meta['crln_new'] = aia_map.meta['crln_obs']
    return synop_map


def aia_fov(dtime):
    l0 = sunpy.coordinates.sun.L0(dtime)
    bounds = Longitude([l0 - 90 * u.deg, l0 + 90 * u.deg])
    return bounds


if __name__ == '__main__':
    dtime = datetime.now()
    dtime = datetime(2020, 2, 1)
    map = create_synoptic_map(dtime)
    # Norm the data
    data = map.data
    data = map.plot_settings['norm'](data)

    map = Map((data, map.meta))
    datestr = dtime.strftime('%Y%m%d')
    map.save(f'aia193_synoptic_latest_{datestr}.fits')
