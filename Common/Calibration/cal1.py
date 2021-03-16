import matplotlib.pyplot as plt
import numpy as np
import warnings
import math

warnings.simplefilter("ignore")

data = np.loadtxt('star_names.txt', dtype=str)
star_name = data[:,0]
star_epoch = data[:,1]
star_mag = data[:,2]
psf_mag = data[:,3]
Number_stars = len(star_name)

star_mag_new = star_mag.astype(np.float)
psf_mag_new = psf_mag.astype(np.float)
a = np.zeros(Number_stars)
b = np.zeros(Number_stars)
f = np.zeros(Number_stars)

for i in range(0, Number_stars):
	a[i]=star_mag_new[i]-psf_mag_new[i]
	b[i]=a[i]/-2.5
	f[i] = pow(10, b[i])
	
np.savetxt('star_data.txt', np.column_stack([star_name, star_epoch, f]), fmt='%s')

print('\nEnd of script')
