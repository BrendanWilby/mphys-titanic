#This script will take in names and epochs for stars in DangerZone
#from a file called: star_data.txt
#It will read the exposure time from the first image for each epoch 
#of each star. It will use this and Alexis' correction factor f to
#compute a correction to calibrate the contrast curves.
#
#The outputted calibration factor should be multiplied with the
#contrast curves:
#5sig(new) = 5sig(old) * cal

import numpy as np
import os
from astropy.io import fits
import warnings

warnings.simplefilter("ignore")

data = np.loadtxt('star_data.txt', dtype=str)
star_names = data[:,0]
star_epochs = data[:,1]
f = data[:,2]
f_new = f.astype(np.float)
Number_stars = len(star_names)

os.chdir('../DangerZone')

exp_times = np.arange(0, Number_stars, dtype=float)

underscore = [None]*Number_stars

for i in range(0, Number_stars):
	underscore[i] = star_names[i].replace('_','')
	os.chdir('{name}/{date}'.format(name=underscore[i], date=star_epochs[i]))
	
	image_filename = np.loadtxt('{name}_filenames.txt'.format(name=star_names[i]), dtype=str, max_rows=2)
	
	hdulist = fits.open(image_filename[0], ignore_missing_end=True, verbose=False)
	coadd_cube = hdulist[0].header['COADDS']
	int_time_cube = hdulist[0].header['ITIME']
	exp_times[i] = coadd_cube * int_time_cube

	os.chdir('../..')


os.chdir('../calibration')

prelim = np.multiply(exp_times,f_new)
cal = np.divide(5,prelim)


np.savetxt('stars_calibration.txt', np.column_stack([star_names,star_epochs,cal]), fmt='%s')


print('End of script')
