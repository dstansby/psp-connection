"""
Main code.
"""
import argparse

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


def create_figure(dtime):
    # Get PFSS/GONG data
    gong_map = gong_helpers.get_closest_map(dtime)
    input, ssmap, header = pfss_helpers.compute_pfss(gong_map)
    gong_date = header['DATE']

    # Get PSP location data
    psp_loc = psp_helpers.psp_loc(dtime)

    # Trace magnetic field line

    # Plot everything
    fig, axs = plt.subplots(nrows=2)
    fig.subplots_adjust(hspace=0.3)
    ax = axs[0]
    input.plot_input(ax, norm=mcolor.SymLogNorm(linthresh=5, vmin=-100, vmax=100))
    ax.set_title(f'Input GONG map ({gong_date})')
    ax.set_xlabel('')

    ax = axs[1]
    pfsspy.plot.radial_cut(input.grid.pg, input.grid.sg, ssmap, ax)
    phi, theta = np.meshgrid(input.grid.pg, input.grid.sg)
    ax.contour(np.rad2deg(phi), theta, ssmap, levels=[0])
    ax.set_title('Source surface magnetic field')

    ax.scatter(psp_loc.lon / u.deg, np.sin(psp_loc.lat), color='black', s=5)
    ax.text(5, 0.85, f'PSP r = {psp_loc.radius[0].value:.02} AU', color='white', fontsize=8)

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
    edate = peri_dates['end'] + np.timedelta64(1, 'D')
    print(sdate, edate)

    # Loop through each day
    while sdate < edate:
        print(sdate)
        fig = create_figure(sdate + np.timedelta64(12, 'h'))
        fig.savefig(f'figures/{peri_n}/{sdate.year}{sdate.month:02}{sdate.day:02}.png',
                    bbox_inches='tight', dpi=150)
        sdate += np.timedelta64(1, 'D')
