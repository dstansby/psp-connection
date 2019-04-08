import pathlib
import urllib.request
import urllib.error

local_aia_dir = pathlib.Path('/Users/dstansby/Data/sdo/aia/synoptic')
remote_dir = 'http://jsoc2.stanford.edu/data/aia/synoptic'


def closest_aia_map(dtime):
    """
    Returns filename of the
    """
    pass


def get_aia_synoptic(wavelength, year, month, day, hour, minute):
    """
    wavelength : string or int
        Either '193',
    year, month, day, hour, minute : int
    """
    if minute % 2 != 0:
        raise ValueError('minute must be a multiple of 2')

    year = str(year)
    month = f'{month:02d}'
    day = f'{day:02d}'
    hour = f'{hour:02d}'
    minute = f'{minute:02d}'
    wavelength = f'{wavelength:04d}'

    fname = f'AIA{year}{month}{day}_{hour}{minute}_{wavelength}.fits'
    local_dir = local_aia_dir / year / month / day
    if not local_dir.exists():
        local_dir.mkdir(exist_ok=True, parents=True)

    local_file = local_dir / fname
    if not local_file.exists():
        remote_file = remote_dir + f'/{year}/{month}/{day}/H{hour}00/{fname}'
        print(f'Downloading {remote_file}')
        try:
            urllib.request.urlretrieve(remote_file, local_file)
        except urllib.error.HTTPError:
            return

    return local_file
