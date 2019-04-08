import numpy as np
import pfsspy
import sunpy.io.fits


def compute_pfss(gong_map):
    [[br, header]] = sunpy.io.fits.read(gong_map)
    br = br - np.mean(br)
    br = np.roll(br, header['CRVAL1'] + 180, axis=1)
    nrho = 60
    rss = 2.5

    input = pfsspy.Input(br, nrho, rss)
    output = pfsspy.pfss(input)
    return input, output, header
