"""
Main code.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import aia_helpers
import gong_helpers
import pfss_helpers


def create_figure(dtime):
    gong_map = gong_helpers.get_closest_map(dtime)
    input, output, header = pfss_helpers.compute_pfss(gong_map)

    fig, axs = plt.subplots(nrows=2)
    ax = axs[0]
    input.plot_input(ax)
    ax.set_title('Input GONG map')

    ax = axs[1]
    mesh = output.plot_source_surface(ax)
    output.plot_pil(ax)
    ax.set_title('Source surface magnetic field')

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
        fig.savefig(f'{sdate.day}.png', bbox_inches='tight')
        sdate += np.timedelta64(1, 'D')
