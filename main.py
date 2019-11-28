"""
Main code.
"""
import argparse
import pathlib
from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import numpy as np
import pandas as pd
import pfsspy

import aia_helpers
import gong_helpers
import pfss_helpers
import psp_helpers
import website_helpers


def create_figure(dtime, aia_map):
    # Get PFSS/GONG data
    gong_map, infuture = gong_helpers.get_closest_map(dtime)
    input, ssmap, header = pfss_helpers.compute_pfss(gong_map)
    gong_date = header['DATE']

    # Get PSP location data
    psp_loc = psp_helpers.psp_loc(dtime)
    psp_loc = psp_helpers.spiral_correction(psp_loc, 400 * u.km / u.s)

    # Trace magnetic field line
    lon, sinlat = pfss_helpers.trace(gong_map, psp_loc, input, retrace=True)

    def add_fline(ax):
        ax.plot(lon, sinlat, lw=1, color='k')
        psp_loc.representation_type = 'spherical'
        ax.scatter(psp_loc.lon / u.deg, np.sin(psp_loc.lat), color='black', s=5)

    # Plot everything
    fig = plt.figure(figsize=(6, 8))

    ax = fig.add_subplot(311)
    input.plot_input(ax, norm=mcolor.SymLogNorm(linthresh=5, vmin=-100, vmax=100))
    ax.set_title(f'Input GONG map ({gong_date})')
    ax.set_xlabel('')
    add_fline(ax)

    ax = fig.add_subplot(312)
    pfsspy.plot.radial_cut(input.grid.pg, input.grid.sg, ssmap, ax)
    phi, theta = np.meshgrid(input.grid.pg, input.grid.sg)
    ax.contour(np.rad2deg(phi), theta, ssmap, levels=[0])
    ax.set_title('Source surface magnetic field')
    add_fline(ax)
    ax.text(5, 0.85, (f'PSP r = {psp_loc.radius[0].value:.03} AU, '
                      f't = {dtime}'),
            color='white', fontsize=8)

    ax = fig.add_subplot(313)
    data = aia_map.data
    data = np.roll(data, data.shape[1] // 2, axis=1)
    ax.imshow(data, extent=(0, 360, -90, 90), origin='lower',
              cmap=aia_map.plot_settings['cmap'],
              norm=aia_map.plot_settings['norm'])
    ax.plot(lon, np.rad2deg(np.arcsin(sinlat)), lw=1, color='w')
    psp_loc.representation_type = 'spherical'
    ax.scatter(psp_loc.lon / u.deg, psp_loc.lat / u.deg, color='w', s=5)

    fig.subplots_adjust(hspace=0.3)
    return fig


if __name__ == '__main__':
    # Get perihelion number from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('peri_n', metavar='n', type=int,
                        help='perihelion number')
    args = parser.parse_args()
    peri_n = args.peri_n

    # Read in list of perihelion dates
    peri_dates = pd.read_csv(
        'peri_dates.csv', parse_dates=[1, 2]).set_index('n')
    peri_dates = peri_dates.loc[peri_n]
    sdate = peri_dates['start']
    edate = peri_dates['end']
    print(sdate, edate)

    # Get an AIA synoptic map
    if edate > datetime.now():
        aia_map = aia_helpers.create_synoptic_map(datetime.now())
    else:
        aia_map = aia_helpers.create_synoptic_map(edate)

    # Loop through each day
    while sdate < edate:
        print(f"Creating figure for {sdate}")
        fig = create_figure(sdate + np.timedelta64(12, 'h'), aia_map)
        savedir = pathlib.Path('figures') / str(peri_n)
        savedir.mkdir(exist_ok=True, parents=True)
        fig.savefig(savedir / f'{sdate.year}{sdate.month:02}{sdate.day:02}.png',
                    bbox_inches='tight', dpi=150)
        plt.close(fig)
        sdate += np.timedelta64(1, 'D')

    website_helpers.gen_html(peri_n)
