from datetime import datetime
import numpy as np

from sunpy.map import Map
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm
from sunpy.coordinates.frames import HeliographicCarrington
from sunpy.coordinates.ephemeris import get_earth, get_horizons_coord
from astropy.coordinates import SkyCoord
import astropy.units as u

aia = Map('aia*.fits')[-1]
euvi = Map('euvi*.fits')[-1]


def update_line(m, offset):
    coord = SkyCoord(m.meta['crln_new'] * u.deg,
                     0 * u.deg, obstime=m.date,
                     frame='heliographic_carrington')
    pix = m.wcs.world_to_pixel(coord)[0]
    return pix - ((offset / 360) * m.data.shape[0])


earth_pix = update_line(aia, 60)
stereo_pix = update_line(euvi, 60)


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(211, projection=aia)
aia.plot(axes=ax, cmap='sdoaia193')
ax.set_title('AIA 193')
ax.axvline(earth_pix, color='w', linewidth=1, linestyle='--')
ax.text(earth_pix + 10, aia.data.shape[0] + 15, 'New', color='black')
ax.text(earth_pix - 80, aia.data.shape[0] + 15, 'Old', color='black')
ax.set_xlabel(' ')

ax = fig.add_subplot(212, projection=euvi)
euvi.plot(axes=ax, cmap='sdoaia193')
ax.set_title('EUVI 195')
ax.axvline(stereo_pix, color='w', linewidth=1, linestyle='--')
ax.text(stereo_pix + 10, aia.data.shape[0] + 15, 'New', color='black')
ax.text(stereo_pix - 80, aia.data.shape[0] + 15, 'Old', color='black')

fig.subplots_adjust(hspace=0.3)
datestr = datetime.now().strftime('%Y%m%d')
fig.savefig(f'AIA_EUVI_synoptic_latest_{datestr}.png', bbox_inches='tight')

'''
fig, ax = plt.subplots()
ax.hist(aia.data.ravel(), bins='auto', histtype='step', label='AIA')
ax.hist(euvi.data.ravel(), bins='auto', histtype='step', label='EUVI')
ax.legend()
'''

plt.show()
