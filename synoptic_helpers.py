"""
Functions to help creation of synoptic maps.
"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.map import make_fitswcs_header


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
