"""
**************************************NOTES****************************************
**************************************INSTRUCTIONS**********************************
To use the script, you can put it anywhere you like on analysis. The default wdir is '/data/shinkley/Keck_Data/completed_stars/DangerZone/'.
An example output and files required has been uploaded on github under 'Common/Calibration'

Usage: 
	python cal.py <cal_name> <cal_date> <cal_name_file> <mag> <output>


	File you need: 
		cal_name [.txt] : any text file contains the name of the directory (e.g. HIP16095)
		cal_name_file [.txt] : any text file contains the name of the contrast output (e.g. HIP_16095)
		cal_date [.txt] : any text file contains the date of the dir (e.g. 2017sep01) 
		mag [.txt]: magnitude of the targets and the PSF 

	Ouput: 
		output: name.csv file which contains the ['name', 'date', 'mag', 'mag_psf','factor']. 

Example: 
	python cal.py cal_name.txt cal_date.txt cal_name_file.txt mag.txt output 


"""
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
wdir = '/data/shinkley/Keck_Data/completed_stars/DangerZone/'

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


"""
Contrastcurvedata writes the contrast curve results into 
a text file, saves and then loads it back, plots the LLSG
and the PCA data on the same graph and saves it to a file
"""
"""
Note: This function is terribly inefficient, could be much 
better looped over but haven't had the time to do so.
"""
	
def Contrastcurvedata(pxscale_keck, sigma,a,f,name):
	
	#Saves the PCA,ADI and LLSG curve outputs and reads them back in.
	from astropy.io import fits
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	from pylab import figure
	
	PCA_data=np.loadtxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name))
	LLSG_data=np.loadtxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name))
	ADI_data=np.loadtxt('new_ADI_{star_name}_curve_outputs'.format(star_name=name))
	
	#Initialises the PCA variables (500 being the size).
	PCA_Distance=np.zeros(500)
	PCA_Sensitivity_Gauss=np.zeros(500)
	PCA_Sensitivity_Student=np.zeros(500)
	LLSG_Distance=np.zeros(500)
	LLSG_Sensitivity_Gauss=np.zeros(500)
	LLSG_Sensitivity_Student=np.zeros(500)
	ADI_Distance=np.zeros(500)
	ADI_Sensitivity_Gauss=np.zeros(500)
	ADI_Sensitivity_Student=np.zeros(500)
	outputs=[None]*500

	#Loop loads the PCA outputs into the PCA variables.
	for i in range (0,500):
		PCA_Distance[i] =PCA_data[i,0] * pxscale_keck
		PCA_Sensitivity_Gauss[i]=PCA_data[i,2] / f
		PCA_Sensitivity_Student[i]= PCA_data[i,3] / f
	#Loop loads the LLSG outputs into the LLSG variables.
		LLSG_Distance[i] =LLSG_data[i,0] * pxscale_keck
		LLSG_Sensitivity_Gauss[i]=LLSG_data[i,2] / f
		LLSG_Sensitivity_Student[i]= LLSG_data[i,3] / f
	
	#Loop loads the ADI outputs into the ADI variables.
		ADI_Distance[i] =ADI_data[i,0] * pxscale_keck
		ADI_Sensitivity_Gauss[i]=ADI_data[i,2] / f
		ADI_Sensitivity_Student[i]= ADI_data[i,3] / f
		outputs[i]=(ADI_Distance[i], PCA_Sensitivity_Gauss[i], ADI_Sensitivity_Gauss[i],LLSG_Sensitivity_Gauss[i])

	#save the output
	np.savetxt('cal_contrast_{star_name}'.format(star_name=name), outputs)	
	
	#Plotting all 3 on one curve: 
		
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	ADI_Legend = mpatches.Patch(color='black', label='ADI Contrast Curve')
	plt.legend(handles=[PCA_Legend,LLSG_Legend,ADI_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = [0,0,0]	#Initialises the Curve variable.
	Curve[2] = plt.plot(ADI_Distance, ADI_Sensitivity_Gauss, linewidth =2, color='black')
	Curve[1] = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	Curve[0] = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	#Saves the figure and then shows the plot. 
	savefig('cal_new_Contrast_curves_{star_name}.png'.format(star_name=name))
	
	#LLSG contrast curve 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	#Saves the figure and then shows the plot. 
	savefig('cal_new_LLSG_curves_{star_name}.png'.format(star_name=name))

	
	#The ADI curve:
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	ADI_Legend = mpatches.Patch(color='black', label='ADI Contrast Curve')
	plt.legend(handles=[ADI_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(ADI_Distance, ADI_Sensitivity_Gauss, linewidth =2, color='black')
	
	#Saves the figure and then shows the plot. 
	savefig('cal_new_ADI_curve_{star_name}.png'.format(star_name=name))

	#PCA contrast curve. 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	plt.legend(handles=[PCA_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	
	#Saves the figure and then shows the plot. 
	savefig('cal_new_PCA_curve_{star_name}.png'.format(star_name=name))


##################################################################

def main ():
	parser = argparse.ArgumentParser() 

	parser.add_argument('name', 
		help = 'the .txt file which contain targets'
		)
	parser.add_argument('date', 
		help = 'the .txt file which contain dates'
	)

	parser.add_argument('filename', 
		help = 'the .txt file which contain filename'
	)

	parser.add_argument('magnitude', 
		help = 'the .txt file which contain magnitude'
	)

	parser.add_argument('output', 
		help = 'the output filename'
	)

	arg = parser.parse_args()
	
	name = [line.rstrip('\n') for line in open(arg.name)]
	#print name

	date = [line.rstrip('\n') for line in open(arg.date)]
	#print date

	name_file = [line.rstrip('\n') for line in open(arg.filename)]

	mag = np.loadtxt(arg.magnitude)

	m = np.zeros(len(name))
	m_psf = np.zeros(len(name))
	data = []

	os.chdir(wdir)

	for b in range (0,len(name)):
		os.chdir('{star_name}/{epoch}'.format(star_name=name[b],epoch=date[b]))
		# run the script in here
		m[b] = mag[b,0]
		m_psf[b]= mag[b,1]

		#need to have apparant magnitude of all stars then open it as a array, make m_psf as a array as well
		#m = input ('What is the apparant magnitude of {name}?'.format(name=name_input))

		a = (m[b] - m_psf[b])/(-2.5)
		a = np.asscalar(a)
		f = math.pow(10,a)
		ls = [name[b], date[b], m[b], m_psf[b],f]

		data.append(ls)

		#Contrastcurvedata(pxscale_keck, sigma,a,f,name=name[b])

		print(a)
		os.chdir('../..')
		b = b + 1

	
	keys = ['name', 'date', 'mag', 'mag_psf','factor']
	df = pd.DataFrame(data, columns = keys)
	df.to_csv(arg.output)
	print("The output file has been save to %s %s.csv" %(os.getcwd(), arg.output))

if __name__ == "__main__": 
	main()
	print("\nEnd of programme")




