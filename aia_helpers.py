from datetime import datetime, timedelta
import pathlib
import urllib.request
import urllib.error

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.wcs import WCS
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


def start_of_day(dtime):
    return datetime(dtime.year, dtime.month, dtime.day)


def download_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    print(f'Fetching map for {dtime}')
    query = (a.Time(dtime, dtime + timedelta(days=1), dtime),
             a.Instrument('AIA'),
             a.Wavelength(193 * u.Angstrom))
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, 0])[0]
    mappath = pathlib.Path(mappath)
    mappath.replace(map_path(dtime))


def load_start_of_day_map(dtime):
    dtime = start_of_day(dtime)
    mappath = map_path(dtime)
    if not mappath.exists():
        download_start_of_day_map(dtime)

    print(f'Loading map for {dtime}')
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


def synop_reproject(m):
    m.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
    shape_out = [360, 720]
    header = synop_header(shape_out, m.date)
    array, footprint = reproject_interp(m, WCS(header),
                                        shape_out=shape_out)
    new_map = Map((array, header))
    new_map.plot_settings = m.plot_settings
    return new_map


def create_synoptic_map(endtime):
    """
    Create an AIA synoptic map, using 27 daily AIA 193 maps ending on the
    endtime given. Note that the maps are taken from the start of each day.
    """
    maps = []
    for i in range(2):
        dtime = endtime - timedelta(days=i)
        aia_map = load_start_of_day_map(dtime)
        aia_synop_map = synop_reproject(aia_map)
        aia_synop_map.peek()
        exit()


if __name__ == '__main__':
    dtime = datetime(2019, 11, 1)
    create_synoptic_map(dtime)
