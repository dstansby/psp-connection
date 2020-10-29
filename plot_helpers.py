import aia_helpers

import matplotlib.patches as mpatch
from astropy.coordinates import Longitude
import astropy.units as u
import numpy as np

map_width = 1440


def add_fov(ax, dtime):
    aia_fov = aia_helpers.aia_fov(dtime)
    stereo_fov = Longitude(aia_fov - 78 * u.deg)

    height = 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    y0 = ax.get_ylim()[1]

    def add_lon_fov(fov, y0, color):
        kwargs = dict(clip_on=False, color=color, alpha=0.7)
        x = np.array([lon.to(u.deg).value for lon in fov]) * map_width / 360
        # If FOV crosses zero longitude
        bar_height = 0.8 * height
        if x[0] > x[1]:
            rects = [mpatch.Rectangle((0, y0), x[1] - 0, bar_height, **kwargs),
                     mpatch.Rectangle((x[0], y0), map_width - x[0], bar_height, **kwargs)]
        else:
            rects = [mpatch.Rectangle((x[0], y0), x[1] - x[0], bar_height, **kwargs)]
        for rect in rects:
            ax.add_patch(rect)

    # Add SDO then STEREO longs
    add_lon_fov(aia_fov, y0 + 0.5 * height, 'C2')
    add_lon_fov(stereo_fov, y0 + 1.5 * height, 'C3')
