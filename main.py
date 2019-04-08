"""
Main code.
"""
import argparse

import pandas as pd

import gong_helpers
import aia_helpers

# Get perihelion number from the command line
parser = argparse.ArgumentParser()
parser.add_argument('peri_n', metavar='n', type=int, help='perihelion number')
args = parser.parse_args()
peri_n = args.peri_n

# Read in list of perihelion dates
peri_dates = pd.read_csv('peri_dates.csv', parse_dates=[1, 2]).set_index('n')
peri_dates = peri_dates.loc[peri_n]
print(peri_dates)
