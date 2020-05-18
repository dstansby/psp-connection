import pathlib
import shutil
import subprocess

website_folder = pathlib.Path('/Users/dstansby/github/dstansby.github.io')
code_folder = pathlib.Path('/Users/dstansby/github/psp-connection')

width = 796
height = 1142


def html_img(fname):
    return '\n      <img class="mySlides" src="{}" style="width:100%">'.format(fname)


def copy_images():
    # Copy images across
    img_folder = code_folder / 'figures'
    web_img_folder = website_folder / 'images' / 'solo'
    web_img_folder.mkdir(parents=True, exist_ok=True)
    # Collect all files to copy
    imgs = [x for x in img_folder.iterdir() if x.suffix == '.png']
    # Copy files and create list of new images
    new_imgs = [shutil.copy(img, web_img_folder) for img in imgs]
    new_imgs = [pathlib.Path(img) for img in new_imgs]
    new_imgs.sort()
    return new_imgs


def gen_html():
    """
    Generate slideshow HTML.
    """
    print(f'Generating HTML')
    html = f'''<html>
<body>
<center>

 <video width="{width * 3 // 4}" height="{height * 3 // 4}" controls>
  <source src="/movies/solo_psp_connection/solo_psp_connection.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>

</center>
</body>
</html>
    '''
    with open(website_folder / '_includes' / f'solo_slideshow.html', "w") as f:
        f.write(html)


def gen_movie():
    fname = 'solo_psp_connection.mp4'
    ffmpeg_cmd = f"ffmpeg -framerate 5 -pattern_type glob -i 'figures/*.png' -pix_fmt yuv420p -vf scale=796:1142 {fname}"
    subprocess.check_output(ffmpeg_cmd, shell=True)
    shutil.copy(fname, website_folder / 'movies' / 'solo_psp_connection' / fname)


if __name__ == '__main__':
    gen_movie()
    gen_html()
