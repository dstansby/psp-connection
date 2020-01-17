from datetime import datetime, timedelta
import pathlib


import numpy as np

import astropy.units as u
from astropy.wcs import WCS
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map
import sunpy.sun.constants
from reproject import reproject_interp


from time_helpers import start_of_day
from synoptic_helpers import synop_header


map_dir = pathlib.Path('/Users/dstansby/Data/euvi')
map_dir.mkdir(exist_ok=True, parents=True)


def map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'euvi_195_{datestr}.fits'


def synoptic_map_path(dtime):
    datestr = dtime.strftime('%Y%m%d')
    return map_dir / f'euvi_195_synoptic_{datestr}.fits'


def download_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    print(f'Fetching EUVI map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('EUVI'))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No EUVI map available for {dtime}')
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def load_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    print(f'Loading EUVI map for {dtime}')
    return Map(str(mappath))


def synop_reproject(m, shape_out):
    """
    Reproject a helioprojective map into a synoptic map.
    """
    synop_map_path = synoptic_map_path(m.date.to_datetime())
    if not synop_map_path.exists():
        m.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
        header = synop_header(shape_out, m.date)
        array, footprint = reproject_interp(m, WCS(header),
                                            shape_out=shape_out)
        new_map = Map((array, header))
        new_map.save(str(synop_map_path))

    print(f'Loading {synop_map_path}')
    new_map = Map(synop_map_path)
    new_map.plot_settings = m.plot_settings
    return new_map


def create_synoptic_map(endtime):
    """
    Create an AIA synoptic map, using 27 daily AIA 193 maps ending on the
    endtime given. Note that the maps are taken from the start of each day.

    Returns
    -------
    sunpy.map.Map : synoptic map
    """
    shape = [720, 1440]
    data = np.zeros(shape)
    weight_sum = np.zeros(shape)
    nmaps = 27
    for i in range(nmaps):
        dtime = endtime - timedelta(days=i)
        try:
            aia_map = load_start_of_day_map(dtime)
        except RuntimeError:
            print(f'Failed to load map for {dtime}')
            continue

        aia_synop_map = synop_reproject(aia_map, shape)

        # Create weights
        coord = sunpy.map.all_coordinates_from_map(aia_synop_map)
        longs = coord.lon.to(u.deg).value
        l0 = sunpy.coordinates.sun.L0(dtime).to(u.deg).value - 78
        dcenterlong = (longs - l0 + 180) % 360 - 180
        weights = np.exp(-(dcenterlong / 10)**2)
        weights[weights < 0] = 0

        aia_data = aia_synop_map.data
        aia_data[np.isnan(aia_data)] = 0
        data += (aia_data * weights)
        weight_sum += weights

    weight_sum[weight_sum == 0] = np.nan
    data /= weight_sum
    synop_map = Map((data, aia_synop_map.meta))
    synop_map.plot_settings = aia_synop_map.plot_settings
    return synop_map


if __name__ == '__main__':
    map = create_synoptic_map(datetime.now() - timedelta(days=4))
    data = map.plot_settings['norm'](map.data)
    map = Map((data, map.meta))
    datestr = datetime.now().strftime('%Y%m%d')
    map.save(f'euvi195_synoptic_latest_{datestr}.fits')
