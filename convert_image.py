"""
This script will convert individual .fits images into .png type
If you are looking to convert a cube, use the convert_cube.py script
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits
from scipy.misc import imsave
import warnings
warnings.simplefilter("ignore")


image_file = input('Please enter the name of the image (of the form .fits):')

# Sanity check to ensure the image inputted actually exists.
# If not, it will continue to loop until the correct name is inputted
while (not os.path.isfile(image_file)):
    print("Image '%s' not found" % image_file)
    image_file = input('Please enter the name of the image (of the form .fits):')

name = input('Please enter the new name of the image (of the form png):')

hdu_list = fits.open(image_file, ignore_missing_end=True)
image=hdu_list[0].data
print(image.shape)

flip = np.flipud(image)

imsave(name, flip)

print('\nEnd of script')
