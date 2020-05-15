"""
Functions to help creation of synoptic maps.
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.map import make_fitswcs_header
import sunpy.coordinates.sun
import matplotlib.pyplot as plt


def synop_header(shape_out, dtime):
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime)
    ref_pix = [(shape_out[1] - 1) // 2, (shape_out[0] - 1) // 2] * u.pix
    header = make_fitswcs_header(
        shape_out, frame_out,
        reference_pixel=ref_pix,
        scale=[360 / shape_out[1],
               180 / shape_out[0]] * u.deg / u.pix,
        projection_code="CAR")
    return header


def synop_weights(synop_map):
    """
    Get a set of weights to apply to a synoptic map.
    """
    crln_obs = synop_map.meta['crln_obs'] * u.deg
    coord = sunpy.map.all_coordinates_from_map(synop_map)
    longs = coord.lon.to_value(u.deg)
    crln_obs = crln_obs.to_value(u.deg)
    dcenterlong = longs - crln_obs
    with np.errstate(invalid='ignore'):
        dcenterlong[dcenterlong > 180] -= 360
        dcenterlong[dcenterlong < -180] += 360
    weights = np.exp(-(dcenterlong / 10)**2)
    return weights
