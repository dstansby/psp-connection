import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
import sunpy.coordinates.frames as frames
import heliopy.data.spice as spicedata
from heliopy import spice

for kernel in ['psp', 'planet_trajectories',
               'planet_orientations', 'psp_pred', 'solo']:
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


def solo_loc(dtime):
    solo = spice.Trajectory("Solar Orbiter")
    solo.generate_positions([dtime], 'Sun', 'IAU_SUN')
    solo.change_units(u.au)
    solo_coord = coord.SkyCoord(x=psp.x, y=psp.y, z=psp.z,
                                frame=frames.HeliographicCarrington,
                                representation_type='cartesian')
    solo_coord.representation_type = 'spherical'
    return solo_coord


def spiral_correction(psp_coord, vsw):
    omega_sun = 14.713 * u.deg / u.d

    def delta_long(r):
        return omega_sun * (r - 2.5 * const.R_sun) / vsw

    psp_solar_lon = psp_coord.lon + delta_long(psp_coord.radius)
    psp_solar_surface = coord.SkyCoord(radius=psp_coord.radius,
                                       lat=psp_coord.lat,
                                       lon=psp_solar_lon,
                                       frame=frames.HeliographicCarrington,
                                       representation_type='spherical')
    return psp_solar_surface


def psp_xyz(dtime):
    psp = spice.Trajectory('SPP')
    psp.generate_positions([dtime], 'Sun', 'IAU_SUN')
    psp.change_units(u.au)
    return np.array([psp.x.value, psp.y.value, psp.z.value])[:, 0]
