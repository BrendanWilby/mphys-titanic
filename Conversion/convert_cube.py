"""
This script will split .fits cubes into individual images
It will then convert a frame of the user's choice into a .png image
If you are looking to convert an individual image:
use the convert_image.py script
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import math
from astropy.io import fits
from scipy.misc import imsave
import warnings
warnings.simplefilter("ignore")


image_file = input('Please enter the name of the cube (of the form .fits):')

# Sanity check to ensure the image inputted actually exists.
# If not, it will continue to loop until the correct name is inputted
while (not os.path.isfile(image_file)):
    print("Image '%s' not found" % image_file)
    image_file = input('Please enter the name of the image (of the form .fits):')

number = int(input('Please enter the frame number (as it appears on ds9) you wish to use:'))
name = input('Please enter the name of the new image (of the form png):')
num = number-1

hdu_list = fits.open(image_file, ignore_missing_end=True)
image=hdu_list[0].data
print(image.shape)

use = image[num]

flip = np.flipud(use)

imsave(name, flip)

print('\nEnd of script')
