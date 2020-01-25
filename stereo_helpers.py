from datetime import datetime, timedelta
import pathlib


import numpy as np

import astropy.units as u
from astropy.wcs import WCS
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import AsinhStretch
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map
import sunpy.sun.constants
from reproject import reproject_interp


from time_helpers import start_of_day
from synoptic_helpers import synop_header, synop_weights


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
    nmaps = 23
    for i in range(nmaps)[::-1]:
        dtime = endtime - timedelta(days=i)
        try:
            euvi_map = load_start_of_day_map(dtime)
        except RuntimeError:
            print(f'Failed to load map for {dtime}')
            continue

        aia_synop_map = synop_reproject(euvi_map, shape)
        weights = synop_weights(aia_synop_map, euvi_map.meta['crln_obs'] * u.deg)

        aia_data = aia_synop_map.data
        aia_data[np.isnan(aia_data)] = 0
        data += (aia_data * weights)
        weight_sum += weights

    weight_sum[weight_sum == 0] = np.nan
    data /= weight_sum
    meta = aia_synop_map.meta
    meta['date-obs'] = dtime.strftime('%Y-%m-%dT%H:%M:%S')

    synop_map = Map((data, meta))
    synop_map.plot_settings = aia_synop_map.plot_settings
    synop_map.meta['crln_new'] = euvi_map.meta['crln_obs']
    return synop_map


if __name__ == '__main__':
    map = create_synoptic_map(datetime.now() - timedelta(days=4))
    # Add a fudge factor so we get the same ballpark as AIA data
    data = map.data / 3 + (71 - 959 / 3)
    data[data <= 0] = 0
    # Apply the same normalisation as an AIA map
    data = ImageNormalize(stretch=AsinhStretch(0.01), clip=False)(data)
    map = Map((data, map.meta))
    datestr = datetime.now().strftime('%Y%m%d')
    map.save(f'euvi195_synoptic_latest_{datestr}.fits')
