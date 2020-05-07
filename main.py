"""
Main code.
"""
import argparse
import pathlib
from datetime import datetime, timedelta

from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.patches as mpatch
import numpy as np
import pandas as pd
import pfsspy
import sunpy.map

import aia_helpers
import gong_helpers
import pfss_helpers
import psp_helpers
import plot_helpers
import website_helpers


def create_figure(dtime, aia_map):
    """
    dtime
        Datetime at which the PSP orbit is sampled
    aia_map : GenericMap
        A synoptic AIA map to show on the bottom panel.
    """
    dtime_str = Time(dtime).isot
    # Get PFSS/GONG data
    gong_map, infuture = gong_helpers.get_closest_map(dtime)
    input, ssmap, header = pfss_helpers.compute_pfss(gong_map, dtime)

    # Get PSP location data
    psp_loc = psp_helpers.solo_loc(dtime)
    psp_loc_ss = psp_helpers.spiral_correction(psp_loc, 350 * u.km / u.s)

    # Trace magnetic field line
    fline = pfss_helpers.trace(gong_map, psp_loc_ss, input, retrace=True)

    # Plot everything
    fig = plt.figure(figsize=(6, 8))

    dobs = psp_loc.obstime.isot[0]
    # Magnetogram
    gong_map = input._map_in
    gong_date = gong_map.meta['DATE_ORI']

    ax = fig.add_subplot(3, 1, 1, projection=gong_map)
    gong_map.plot(axes=ax, cmap='RdBu',
                  norm=mcolor.SymLogNorm(linthresh=5, vmin=-100, vmax=100, base=10))
    ax.set_title('Input GONG magnetogram')
    ax.text(0.01, 0.01, (f'Last updated {gong_date}'), color='black', fontsize=8, transform=ax.transAxes)
    for coord in ax.coords:
        coord.set_axislabel(' ')
    ax.plot_coord(fline, lw=1, color='k')
    ax.plot_coord(psp_loc_ss, color='black', marker='o', ms=5)

    # Source surface Br
    ax = fig.add_subplot(3, 1, 2, projection=ssmap)
    ssmap.plot(axes=ax, cmap='RdBu')
    ax.set_title('Source surface magnetic field')
    ax.contour(ssmap.data, levels=[0], colors='black', linewidths=0.5)
    for coord in ax.coords:
        coord.set_axislabel(' ')
    ax.text(0.01, 0.01, (f'Orbiter r = {psp_loc.radius[0].to_value(u.au):.03} AU, '
                       f't = {dtime}'),
            color='white', fontsize=8, transform=ax.transAxes)
    ax.plot_coord(psp_loc_ss, color='black', marker='o', ms=5)

    aia_map.meta['date-obs'] = dtime_str
    # AIA synoptic map
    ax = fig.add_subplot(3, 1, 3, projection=aia_map)
    aia_map.plot(axes=ax)
    ax.plot_coord(fline, lw=1, color='w')
    ax.plot_coord(psp_loc_ss, color='w', marker='o', ms=5)
    ax.set_title('AIA 193 synoptic map')
    # plot_helpers.add_fov(ax, dtime)

    fig.subplots_adjust(hspace=0.3)
    return fig


if __name__ == '__main__':
    # Set start date and end date
    #
    # For PSP, see peri_dates.csv for a list
    sdate = datetime.now()
    sdate = datetime(sdate.year, sdate.month, sdate.day, 0)
    edate = datetime.now() + timedelta(days=7)
    print(sdate, edate)

    # Get an AIA synoptic map
    if edate > datetime.now():
        aia_map = aia_helpers.create_synoptic_map(datetime.now())
    else:
        aia_map = aia_helpers.create_synoptic_map(edate)

    savedir = pathlib.Path('figures')
    savedir.mkdir(exist_ok=True, parents=True)
    # Loop through each day
    while sdate < edate:
        fname = savedir / f'{sdate.year}{sdate.month:02}{sdate.day:02}.png'
        # Check if we already have a file
        if not fname.exists() or sdate > datetime.now():
            print(f"Creating figure for {sdate}")
            fig = create_figure(sdate + timedelta(hours=12), aia_map)
            fig.savefig(savedir / f'{sdate.year}{sdate.month:02}{sdate.day:02}.png',
                        bbox_inches='tight', dpi=150)
            plt.close(fig)
        else:
            print(f'Figure already present for {sdate}')
        sdate += timedelta(days=1)

    # website_helpers.gen_html(peri_n)
