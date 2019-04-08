from datetime import datetime
import pathlib
import sunpy.io.fits

import numpy as np

map_dir = pathlib.Path('/Users/dstansby/Data/gong/daily_synoptic')


def sync_gong(year=datetime.now().year,
              month=datetime.now().month):
    from ftplib import FTP
    import parfive
    year = str(year)
    month = f'{month:02}'

    local = map_dir / (year + month)

    ftp = FTP('gong2.nso.edu')
    ftp.login()
    ftp.cwd(f"/oQR/zqs/{year}{month}/")
    daydirs = ftp.nlst()

    # List to store files that need to be downloaded
    dl = parfive.Downloader()
    for daydir in daydirs:
        ftp.cwd(f"/oQR/zqs/{year}{month}/{daydir}")
        files = ftp.nlst()
        for file in files:
            local_file = local / daydir / file
            if not local_file.exists():
                dl.enqueue_file(f"ftp://gong2.nso.edu/oQR/zqs/{year}{month}/{daydir}/{file}", local / daydir, file)

    if dl.queued_downloads > 0:
        print(f'Downloading {dl.queued_downloads} files')
        dl.download()


def get_closest_map(dtime):
    latest_gong_map()


def gong_daily_files(year, month):
    '''
    Parameters
    ----------
    year : str
        Year number
    month : str
        Month number

    Returns
    -------
    files : list
        List of file locations in given month
    '''
    month_dir = map_dir / f'{year}{month}'
    files = []
    for d in month_dir.iterdir():
        if d.name == '.DS_Store':
            continue
        today_files = [str(f) for f in d.iterdir() if f.suffix == '.fits']
        if not len(today_files):
            continue
        today_files.sort()
        files += today_files
    files.sort()
    return [pathlib.Path(file) for file in files]


def extract_br(filepath):
    """
    Extract Br from a GONG daily map.

    Paramters
    ---------
    m : sunpy.map.Map

    Returns
    -------
    br : 2D array
    """
    [[data, header]] = sunpy.io.fits.read(filepath)
    br = data
    br = br - np.nanmean(br)
    # GONG maps have their LH edge at -180deg, so roll to get it at 0deg
    br = np.roll(br, header['CRVAL1'] + 180, axis=1)
    return br


def extract_date(filepath):
    """
    Get date of a GONG daily map.
    """
    [[data, header]] = sunpy.io.fits.read(filepath)
    return datetime.strptime(header['DATE'], '%Y-%m-%dT%H:%M:%S')
