from datetime import datetime
import gzip
import os
import pathlib
import warnings

import sunpy.io.fits
import numpy as np

map_dir = pathlib.Path('/Users/dstansby/Data/gong/daily_synoptic')


def gong_dir(year, month, day):
    year = str(year)
    month = f'{month:02}'
    day = f'{day:02}'
    return map_dir / (year + month) / f"mrzqs{year[2:]}{month}{day}"


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
    dl = parfive.Downloader(max_conn=2)
    for daydir in daydirs:
        ftp.cwd(f"/oQR/zqs/{year}{month}/{daydir}")
        files = ftp.nlst()
        file = files[0]
        local_file = local / daydir / file
        if not local_file.exists():
            dl.enqueue_file(f"ftp://gong2.nso.edu/oQR/zqs/{year}{month}/{daydir}/{file}", local / daydir, file)

    if dl.queued_downloads > 0:
        print(f'Downloading {dl.queued_downloads} files')
        res = dl.download()
        if len(res.errors):
            print(res.errors)
        unzip_gong()


def unzip_gong():
    """
    Unzip all GONG .fits.gz files
    """
    for month_dir in map_dir.iterdir():
        if month_dir.is_dir():
            for day_dir in month_dir.iterdir():
                if day_dir.is_dir():
                    for file in day_dir.iterdir():
                        if file.suffix == '.gz':
                            if not file.with_suffix('').exists():
                                print(f'Unzipping {file}')
                                with gzip.open(file, 'rb') as gzf:
                                    with open(file.with_suffix(''), 'wb') as g:
                                        g.write(gzf.read())


def get_closest_map(dtime):
    if dtime > datetime.now():
        dtime = datetime.now()
        infuture = True
    else:
        infuture = False
    dir = gong_dir(dtime.year, dtime.month, dtime.day)
    files = [x for x in dir.iterdir() if x.suffix == '.fits']
    if len(files) == 0:
        raise RuntimeError(f'Could not find any .fits files in {dir}')
    hours = [int(x.stem[12:14]) for x in files]
    fidx = np.argmin(np.abs(np.array(hours) - dtime.hour))
    return files[fidx], infuture


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


if __name__ == '__main__':
    sync_gong()
