m = 0.1
m_psf = 7.051
psf_xy =[95,111]			#Centre of psf.
reference_xy=[0.1, 0.1]		#Reference bright point defined as a decimal so it loops through the decimal check.
star_xy=[0.1, 0.1] 			#Centre of the star co-ordinates, defined as a decimal so it loops through the decimal check.
averageflux=0				#Initialises the background flux.				#Just a range of flux's for the first guess algorithm .
flvlup = 15000
flvldown = 1000				#to loop over.
Loop_count=0				
pxscale_keck=0.00953 		#Pixel scale for keck 2.
sigma = 5					#Sets the value of sigma for the contrast curves.
sub_count = 0
im = 0.1

"""
End of variables 
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from math import cos
from math import sin
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import figure
import os 
import argparse
from argparse import ArgumentParser
import glob
import warnings

warnings.simplefilter("ignore")

	
def Contrastcurvedata(pxscale_keck, sigma, name, epoch, cal):
	#RDI contrast curve:
	#(the important one!)

	from astropy.io import fits
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	from pylab import figure
	
	underscore = name.replace('_','')

	os.chdir('{u}/RDI'.format(u=underscore))
	
	for file in glob.glob('*%s*'%('contrast_curve*.txt')):
		RDI_data = np.loadtxt(file)	
		length = len(RDI_data)
		
		#Initialises the RDI variables (length being the size).
		RDI_Distance=np.zeros(length)
		RDI_Sensitivity_Gauss=np.zeros(length)
		RDI_Sensitivity_Student=np.zeros(length)
		outputs=[None]*length

	
		#Loop loads the RDI outputs into the RDI variables.
		for i in range (0,length):
			RDI_Distance[i] =RDI_data[i,4]
			RDI_Sensitivity_Gauss[i]=RDI_data[i,0] * cal
			RDI_Sensitivity_Student[i]= RDI_data[i,1] * cal
			outputs[i]=(RDI_Distance[i], RDI_Sensitivity_Gauss[i], RDI_Sensitivity_Student[i])	
		remove = file.replace('_curve','')
		np.savetxt('calibrated_{name}'.format(name = remove), outputs)
	
		plt.plot(RDI_Distance, RDI_Sensitivity_Gauss, 'r-', label='Sensitivity (Gaussian)')
		plt.plot(RDI_Distance, RDI_Sensitivity_Student, 'b-', label='Sensitivity (Student)')
		plt.legend(loc='upper right')
		plt.xlabel('Angular separation (arcsec)')
		plt.ylabel('5 sigma contrast')
		plt.yscale('log')
		under = remove.replace('.txt','')
		plt.grid()
		plt.savefig('calibrated_{insert}.png'.format(insert=under))
		plt.close()

	os.chdir('../..')



##################################################################


data = np.loadtxt('stars_calibration.txt', dtype=str)
star_names = data[:,0]
star_epochs = data[:,1]
cal = data[:,2]
cal_new = cal.astype(np.float)
Number_stars = len(star_names)

os.chdir('../DangerZone')
	
for b in range (0,Number_stars):
	Contrastcurvedata(pxscale_keck, sigma, name=star_names[b], epoch=star_epochs[b], cal=cal_new[b])

os.chdir('../calibration')

print('\nEnd of program')
