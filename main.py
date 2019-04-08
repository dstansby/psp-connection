"""
Main code.
"""
import argparse

import numpy as np
import pandas as pd

import gong_helpers
import aia_helpers


def create_figure(dtime):
    gong_map = gong_helpers.get_closest_map(dtime)
    print(gong_map)


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
        create_figure(sdate + np.timedelta64(12, 'h'))
        sdate += np.timedelta64(1, 'D')
