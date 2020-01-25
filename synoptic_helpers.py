"""
Functions to help creation of synoptic maps.
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.map import make_fitswcs_header
import sunpy.coordinates.sun


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


def synop_weights(synop_map, crln_obs):
    """
    Get a set of weights to apply to a synoptic map.
    """
    coord = sunpy.map.all_coordinates_from_map(synop_map)
    longs = coord.lon.to_value(u.deg)
    crln_obs = crln_obs.to_value(u.deg)
    dcenterlong = (longs - crln_obs + 180) % 360 - 180
    weights = np.exp(-(dcenterlong / 10)**2)
    return weights
