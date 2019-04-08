import astropy.units as u
import astropy.coordinates as coord
import sunpy.coordinates.frames as frames
import heliopy.data.spice as spicedata
from heliopy import spice

for kernel in ['psp', 'planet_trajectories', 'planet_orientations', 'psp_pred']:
    k = spicedata.get_kernel(kernel)
    spice.furnish(k)


def psp_loc(dtime):
    psp = spice.Trajectory('SPP')
    psp.generate_positions([dtime], 'Sun', 'IAU_SUN')
    psp.change_units(u.au)
    psp_coord = coord.SkyCoord(x=psp.x, y=psp.y, z=psp.z,
                               frame=frames.HeliographicCarrington,
                               representation_type='cartesian')
    psp_coord.representation_type = 'spherical'
    return psp_coord
