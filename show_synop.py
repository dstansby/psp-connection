from datetime import datetime
import numpy as np

from sunpy.map import Map
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm
from sunpy.coordinates.ephemeris import get_earth, get_horizons_coord
import astropy.units as u

aia = Map('aia*.fits')[-1]
euvi = Map('euvi*.fits')[-1]

sta_coord = get_horizons_coord('STEREO-A', time=euvi.date)
earth_coord = get_earth(time=aia.date)


def update_line(m, coord, offset):
    # Get earth update line
    coord = coord.transform_to(m.coordinate_frame)
    pix = m.wcs.world_to_pixel(coord)[0]
    return pix - ((offset / 360) * m.data.shape[0])


earth_pix = update_line(aia, earth_coord, 65)
stereo_pix = update_line(euvi, sta_coord, 65)


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
