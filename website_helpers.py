import pathlib
import shutil

website_folder = pathlib.Path('/Users/dstansby/github/dstansby.github.io')
code_folder = pathlib.Path('/Users/dstansby/github/psp-connection')


def html_img(fname):
    return '\n      <img class="mySlides" src="{}" style="width:100%">'.format(fname)


def copy_images(n_peri):
    # Copy images across
    img_folder = code_folder / 'figures' / str(n_peri)
    web_img_folder = website_folder / 'images' / 'psp' / str(n_peri)
    web_img_folder.mkdir(parents=True, exist_ok=True)
    # Collect all files to copy
    imgs = [x for x in img_folder.iterdir() if x.suffix == '.png']
    # Copy files and create list of new images
    new_imgs = [shutil.copy(img, web_img_folder) for img in imgs]
    new_imgs = [pathlib.Path(img) for img in new_imgs]
    new_imgs.sort()
    return new_imgs


def gen_html(n_peri):
    """
    Generate slideshow HTML.
    """
    print(f'Generating HTML for perihelion {n_peri}')
    html_start = '''
<html>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://www.w3schools.com/w3css/4/w3.css">
<style>
.mySlides {display:none}
</style>
<body>

<div class="w3-content" style="max-width:800px">'''

    html_end = '''
</div>

<div class="w3-center">
  <div class="w3-section">
    <button class="w3-button w3-light-grey" onclick="plusDivs(-1)">❮ Prev</button>
    <button class="w3-button w3-light-grey" onclick="plusDivs(1)">Next ❯</button>
  </div>
</div>

<script>
var slideIndex = 1;
showDivs(slideIndex);

function plusDivs(n) {
  showDivs(slideIndex += n);
}

function currentDiv(n) {
  showDivs(slideIndex = n);
}

function showDivs(n) {
  var i;
  var x = document.getElementsByClassName("mySlides");
  var dots = document.getElementsByClassName("demo");
  if (n > x.length) {slideIndex = x.length}
  if (n < 1) {slideIndex = 1}
  for (i = 0; i < x.length; i++) {
    x[i].style.display = "none";
  }
  for (i = 0; i < dots.length; i++) {
    dots[i].className = dots[i].className.replace(" w3-red", "");
  }
  x[slideIndex-1].style.display = "block";
  dots[slideIndex-1].className += " w3-red";
}
</script>

</body>
</html>
    '''
    new_imgs = copy_images(n_peri)
    img_html = ''
    for img in new_imgs:
        img_html += html_img(f'/images/psp/{n_peri}/{img.name}')
    html = html_start + img_html + html_end
    with open(website_folder / '_includes' / 'psp_slideshow.html', "w") as f:
        f.write(html)
