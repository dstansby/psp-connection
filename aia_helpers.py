from datetime import datetime, timedelta
import pathlib
import urllib.request
import urllib.error

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from sunpy.net import vso
from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map, make_fitswcs_header
import sunpy.sun.constants
from reproject import reproject_interp

map_dir = pathlib.Path('/Users/dstansby/Data/aia')
map_dir.mkdir(exist_ok=True, parents=True)


def map_path(dtime):
    return map_dir / f'aia_193_{dtime.year}{dtime.month}{dtime.day}.fits'


def synoptic_map_path(dtime):
    return map_dir / f'aia_193_synoptic_{dtime.year}{dtime.month}{dtime.day}.fits'


def start_of_day(dtime):
    return datetime(dtime.year, dtime.month, dtime.day)


def download_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    print(f'Fetching map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(193 * u.Angstrom))
    result = Fido.search(*query)
    try:
        mappath = Fido.fetch(result[0, 0])[0]
    except IndexError as e:
        raise RuntimeError(f'No map available for {dtime}')
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def load_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    print(f'Loading AIA map for {dtime}')
    return Map(str(mappath))


def synop_header(shape_out, dtime):
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime)
    header = make_fitswcs_header(
        shape_out, frame_out,
        scale=[180 / shape_out[0],
               360 / shape_out[1]] * u.deg / u.pix,
        projection_code="CAR")
    return header


def synop_reproject(m, shape_out):
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
    """
    shape = [360, 720]
    data = np.zeros(shape)
    weight_sum = np.zeros(shape)
    nmaps = 27
    nskip = 1
    for i in range(nmaps // nskip + 1):
        dtime = endtime - timedelta(days=i * nskip)
        aia_map = load_start_of_day_map(dtime)
        aia_synop_map = synop_reproject(aia_map, shape)

        # Create weights
        coord = sunpy.map.all_coordinates_from_map(aia_synop_map)
        longs = coord.lon.to(u.deg).value
        l0 = sunpy.coordinates.sun.L0(dtime).to(u.deg).value
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
    dtime = datetime(2019, 11, 1)
    create_synoptic_map(dtime)
