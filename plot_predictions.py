import datetime
import pathlib

from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd

import aia_helpers


def plot_predictions(footpoint_file):
    """
    Plot a set of magnetic field predictions.

    Paramters
    ---------
    footpoint_file: pathlike
        File in which the fooptoint predictions are stored.
        The following columns must be present:
            - 'date': dates (including times)
            - 'carlon': carrington longitude in degrees
              (must be between 0 and 360)
            - 'carlat': carrington latitutde in degrees
              (must be between -90 and 90)
    """
    # Read in footpoints
    footpoint_file = pathlib.Path(footpoint_file)
    footpoints_df = pd.read_csv(footpoint_file)

    footpoint_dates = pd.to_datetime(footpoints_df['date'])
    carlon = footpoints_df['carlon'].values * u.deg
    carlat = footpoints_df['carlat'].values * u.deg

    synoptic_map = aia_helpers.create_synoptic_map(datetime.datetime.now())

    footpoint_coords = SkyCoord(carlon, carlat, obstime=synoptic_map.date,
                                frame='heliographic_carrington')

    for date, coord in zip(footpoint_dates, footpoint_coords):
        print(f'Plotting {date}...')
        fig = plt.figure(figsize=(6, 8))
        # synoptic map
        ax = fig.add_subplot(211, projection=synoptic_map)
        synoptic_map.plot(axes=ax)
        ax.plot_coord(coord, marker='+', ms=10, color='white', lw=0)
        ax.set_title(str(date))
        ax.set_xlim(left=-0.49)
        ax.set_ylim(bottom=-0.49)

        # reproject map
        earth_map = aia_helpers.stonyhurst_reproject(synoptic_map, date)

        coord = SkyCoord(coord.lon, coord.lat, obstime=date,
                         frame='heliographic_carrington')
        coord = coord.transform_to(earth_map.coordinate_frame)
        # Earth view
        ax = fig.add_subplot(212, projection=earth_map)
        earth_map.draw_limb(axes=ax)
        earth_map.plot(axes=ax, alpha=0)

        # Plot footpoint
        if coord.distance < coord.observer.radius:
            ax.plot_coord(coord, marker='+', ms=10, color='white', lw=0)
        else:
            ax.text(0.5, 0.5, 'Footpoint behind Sun',
                    transform=ax.transAxes, color='white',
                    horizontalalignment='center',
                    verticalalignment='center')

        ax.set_xlabel('Solar-X (arcsec)')
        ax.set_ylabel('Solar-X (arcsec)')
        ax.set_title('Earth view')
        ax.set_facecolor("k")

        fig.savefig(f'figures/aia/{date}.png', bbox_inches='tight')
        plt.close('all')


if __name__ == '__main__':
    plot_predictions('UCB_psp_e4_20200111a_ADAPTEnsemble2pt5Rs.csv')
