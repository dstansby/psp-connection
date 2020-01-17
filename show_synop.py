from sunpy.map import Map
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm
from sunpy.coordinates.ephemeris import get_earth, get_horizons_coord
import astropy.units as u

aia = Map('aia*.fits')
euvi = Map('euvi*.fits')


sta_coord = get_horizons_coord('STEREO-A', time=euvi.date)
earth_coord = get_earth(time=aia.date)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(211, projection=aia)
aia.plot(axes=ax, cmap='sdoaia193')
ax.set_title('AIA 193')
ax.plot_coord(earth_coord.transform_to(aia.coordinate_frame), marker='o', ms=5)

ax = fig.add_subplot(212, projection=euvi)
euvi.plot(axes=ax, cmap='sdoaia193')
ax.set_title('EUVI 195')
ax.plot_coord(sta_coord.transform_to(euvi.coordinate_frame), marker='o', ms=5)
fig.savefig('AIA_EUVI_synoptic.png', bbox_inches='tight')

'''
fig, ax = plt.subplots()
ax.hist(aia.data.ravel(), bins='auto', histtype='step', label='AIA')
ax.hist(euvi.data.ravel(), bins='auto', histtype='step', label='EUVI')
ax.legend()
'''



plt.show()
