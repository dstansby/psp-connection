"""
Script to generate and print perihelion dates for PSP.
"""
import heliopy.data.spice as spicedata
import heliopy.spice as spice
from datetime import datetime, timedelta
import astropy.units as u
import numpy as np

kernels = spicedata.get_kernel('psp')
kernels += spicedata.get_kernel('psp_pred')
spice.furnish(kernels)
psp = spice.Trajectory('SPP')

starttime = datetime(2018, 8, 14)
endtime = starttime + timedelta(days=1000)
times = []
while starttime < endtime:
    times.append(starttime)
    starttime += timedelta(days=1)

psp.generate_positions(times, 'Sun', 'ECLIPJ2000')
psp.change_units(u.au)

psp_times = np.array(psp.times, dtype=np.datetime64)
psp_r = psp.r
keep = psp.r < 0.3 * u.au

psp_times = psp_times[keep]
peri_dates = np.split(psp_times, np.where(np.diff(psp_times) > np.timedelta64(1, 'D'))[0] + 1)
for i, peri in enumerate(peri_dates):
    print(f'{i + 1}, {peri[0]}, {peri[-1]}')
