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
    ssfile = gong_map.with_suffix('.ssmap')
    if ssfile.exists():
        ssmap = np.loadtxt(ssfile)
        output = None
    else:
        # Compute PFSS solution and source surface map
        output = pfsspy.pfss(input)
        br, _, _ = output.bg
        ssmap = br[:, :, -1].T
        np.savetxt(ssfile, ssmap)

    return input, ssmap, header
