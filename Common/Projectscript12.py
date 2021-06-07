""" 
**************************************NOTES****************************************

PROJECT SCRIPT 12

Script last updated 21/01/2020
Works for VIP 0.9.11 in Python3 environment (current version on analysis)

For any help, contact:
	Jorge Fernández
	MPhys student with Sasha, 2019/2020
	jorgefz.fernandez@gmail.com

**************************************INSTRUCTIONS**********************************

Place script on the folder where the science images are.
Run with python.
Enter the name of the star to run the script, e.g. HIP_99542. 


Ignore this for the moment:
	--Please upload contrcurve.py, cosmetics.py, fakecomp.py and shapes.py when running this script (avaliable in DangerZone)



**************************************UPDATES**************************************

01/06/17
	1. Option to break from loop when looping over different centre postions without
	   having to wait till all positions tried
	2. Stops users enter more six planets at planet injection stage, recomends use of 6
	
05/06/17
	1. Includes allowing for coordinate re-entry after looping over centre positions
	2. allow user to enter what multiple of background flux is used for planet brightness

06/06/17
	1. Clearer instructions and error checking when entering planet flux
	2. Background flux check now takes four 100x100 pixel grids from each image in the 
	   data cube. If the average flux in any of each of this grids is greater than 1 stdev
	   away from the population mean that grid is rejected and not used to build a mean
	   background flux
	3. Minor adjustment to Checkfuntion() allowing for more input options
	4. Changed Checkfunction to allow for exit to be called.	

07/06/17
	1. Addition of loop to allow planet injection phase to run mutliple times.
	2. Loop to redo entirly the re-centring if ds9 inspection shows problem
	3. Made if more obvious when you need to enter reference coordinates
	4. fixed bug in coordinate re-entry
	5. Made Image removal function more efficient and removed the necessity to know how many
	   images were going to be removed before entering numbers.
	6. Background flux now takes grids within 0.5 stdev of mean
08/06/17
	1. Removed Planet injection loop due to bug
	2. The .png files of each reduction were outputing blank files. No way was found to 
	   fix this therefore the code was removed. still saved as .fits but no .png.
	   This appears to be beacuse it is saving wrong figure but this is due to a VIP
	   hence no solution was found.
14/06/17
	1. Fixed the ininate loop when asking to re enter coodinates
	2. Fixed the re entered coordinates not being correctly over writtten

28/05/18
	1. Update VIP from 0.5.4 to 0.9.4
	2. Note that any script older than Projectscript 7 will not work now as structure of packages are different 
	3. Update the script due to the difference in packages

28/06/18 
	1. Add planet subtraction as a loop
	2. Edited the script hence user can run the script by entering the name

28/05/19 
	1. Avoid the crash of contrast curve by doing : 
		import fakecomp.py cosmetics.py shapes.py contrcurve.py (the script I have uploaded on analysis)
		from contrcurve import *
		change "vip.metrics.contrcurve.contrast_curve" to "contrast_curve"

29/05/19 
	1. Fixed bugs due to the updates
	2. Disable Planet subtraction (there should a bug in matplotlib)


----------------------------------------

DEVELOPMENT CONTINUED BY JORGE FERNANDEZ

Contact: jorgefz.fernandez@gmail.com


02/10/19
	1. Porting to Python3
		Changed 'print' to 'print()'
		Changed 'raw_input()' to 'input()'
	
03/10/19
	1. Updating VIP from 0.9.8 to 0.9.11
		pp_subplots no longer supported. Use HCIplot module
		Previous scripts will no longer work, Python 3 required!!
		Function open_adicube() no longer in VIP. Use fits.open_fits()
		Function vip.var.fit_2dgaussian() returns pandas dataframe, 
			output columns in different order now.
		
	2. Porting to Python 3
		Reading filenames changed from np.loadtxt() to standard Python open()
		because Python 3 treats strings and bytes differently
	
	3. Changed requiring-int-loops to simply reading input and using int()
	
04/10/19
	1. Cleaning code
		backgroundflux function is repeated!
		Adding spaces and indentation
		
	2. Porting to VIP 0.9.11
		Function vip.pca.pca_optimize_snr() doesn't work anymore.
		
07/10/19
	Conflict of names: HIP_number and HIPnumber
	Created new function to fit 2d gaussian

09/10/19
	Fixed some Contrast curves problems 
		Still needs optimal number of PCs
	PCA function with full output doesn't have the expected number of outputs
		-Need residual cube for STIM maps
		
15/11/19
	Changed name to Projectscript12
	Removed nested loops that would try different values to recenter star
	Added option to choose between using 2D Gaussian fit or manual fit when recentering
	
18/11/19
	Implementing Andromeda

19/11/19
	Reordered post-processing algorithms into functions
	Obtained optimal PCs from FFPCA output
	Implemented Annular PCA
	Implemented STIM map using residual cube from Annular PCA
	Implemented STIm map using residual cube from LLSG

20/01/2020
	Updated Andromeda Oversampling Factor

21/01/2020
	Fixed Contrast Curves with James Davies' work.

23/03/2020
	Reworked Recentering function
	Added more recentering algorithms

		
		
		
		

**************************************************************************************   
	
"""

"""
Script designed to build a data cube from fits files with the aid of the VIP module
in python, reduce the images using the PCA and LLSG method of imaging reduction. 
Script also then calculates the flux from the exposure times of the fits files and 
then generates contrast curves based on the data cube produced and this value of flux.


STUCTURE

	1. Packages
	
	2. Global variables
	
	3. Functions

	4. Main Functions

	5. Image Processing


""" 


"""
============================	PACKAGES	=================================
"""
import matplotlib.pyplot as plt
import numpy as np
import vip_hci as vip
from vip_hci.preproc import cube_recenter_2dfit, cube_recenter_dft_upsampling, cosmetics, badpixremoval, recentering
from vip_hci.negfc import firstguess
from vip_hci.negfc import cube_planet_free
from vip_hci.negfc import show_corner_plot, show_walk_plot
from vip_hci.negfc import mcmc_negfc_sampling

# VIP 0.9.11 does not use pp_subplots anymore. Use matplotlib instead.
#plots = vip.var.pp_subplots
import hciplot

import math
import pandas as pd
from math import cos
from math import sin
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import figure
import warnings
warnings.simplefilter("ignore")

from vip_hci.metrics.contrcurve import *
"""
End of package imports 
"""





"""
============================	GLOBAL VARIABLES		=================================
"""
psf_xy = [95,111]			#Centre of psf.
reference_xy = [0.1, 0.1]	#Reference bright point defined as a decimal so it loops through the decimal check.
star_xy = [0.1, 0.1] 		#Centre of the star co-ordinates, defined as a decimal so it loops through the decimal check.

pxscale_keck = 0.00953 		#Pixel scale for keck 2 [arcsec/pixel]


#Initialises the background flux.
#Just a range of flux's for the first guess algorithm .
averageflux = 0

#Initialises the loop counter, loop counter is used to prevent error warning on first trial
#when looping through integer check functions.
Loop_count = 0				
seperation = 60 			#Radial seperation between synthetic planets in terms of pixels.
starphot = 0				#Initialises the Starphot parametre.
sigma = 5					#Sets the value of sigma for the contrast curves.
sub_count = 0
im = 0.1
""" End of variables """ 






"""
============================	FUNCTION DEFINITIONS		=================================
"""


"""
Reading fits image names from text file
"""
def ReadImageFilenames(name_input):
	
	#Open file containing filenames and saving the contents
	f = open('{name}_filenames.txt'.format(name = name_input))
	file_names = f.readlines()
	f.close()
	
	# Reads filenames as "xxxx.fits\n" so delete 'newline'
	for i in range(len(file_names)):
		file_names[i] = file_names[i].rstrip()

	return file_names




"""
Starphot calculation calculates the flux from the psf, reads the number of
coadds and the integration time from the headers. The function then uses these
values to calculate the parametre starphot the flux from the images in the cube.
"""

def StarphotCalculation(Number_images, psf, fwhm):
	
	int_time_cube = 0
	coadd_cube = 0
	for i in range(0, Number_images):	
		hdulist = fits.open(file_names[i], ignore_missing_end = True, verbose = False)
			#Creates list of open files
		coadd_cube += hdulist[0].header['COADDS']
			#Reads the number of coadds from the first fits file
		int_time_cube += hdulist[0].header['ITIME']
			#Reads the integration time for each coadd from the first fits file

	hdulist2 = fits.open('psf.fits', ignore_missing_end = True, verbose = False)
	
	int_time_psf = hdulist2[0].header['ITIME']
		#Reads the integration time per coadd from the psf fits file
	coadd_psf = hdulist2[0].header['COADDS']
		#Reads the number of coadds from the psf fits file
	
	ycord = np.full((Number_images), psf.shape[1], dtype=int)
	xcord = np.full((Number_images), psf.shape[0], dtype=int)
	
	flux = vip.metrics.contrcurve.aperture_flux(
			psf, ycord, xcord, fwhm, verbose = False)
	#Returns the sum of pixel values in a circular aperture centered on the input coordinates
	
	psfflux = flux[1]
		
	starphot = psfflux * ((int_time_cube * coadd_cube)/(int_time_psf * coadd_psf))
	
	if (starphot != starphot):
		print(" Error: starphot is NaN")
		exit()

	return starphot



				
"""
Backgroundflux uses a 100 by 100 square away from the centre 
of the images to get a rough guess of the background noise. 
It takes four 100x100 squares for each image in the data cube 
at oposing points of the image. It then checks that each cube has an
average flux within 0.5 stdev of the population mean. If this criteria
is failed the square is thown out. The average flux is calcuated from 
the remaining squares
"""

def backgroundflux(cube,averageflux):

	number_images = len(cube)
		#calculates number of image in cube
		
	totalflux = np.zeros(shape = [4, number_images])
		#Initialises the total flux array variable.
		
	averageflux_array = np.zeros(shape = [4, number_images])
		#Initialises an array to contain average of each 100x100 grid
													
	for k in range (0, (number_images)):
		for i in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[0,k] = totalflux[0,k] + cube[k,i,j]
				
		for i in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[1,k] = totalflux[1,k] + cube[k,i,j]
			
		for i in range (800, 901): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[2,k] = totalflux[2,k] + cube[k,i,j]
			
		for i in range (900, 801): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[3,k] = totalflux[3,k] + cube[k,i,j]
		
		for i in range(0, 4):
			averageflux_array[i,k] = totalflux[i,k]/10000
			#Calculates the average flux in each grid by dividing by number of pixels (10000)
			
	upper_bound = np.mean(averageflux_array) + (np.std(averageflux_array)/2)
	lower_bound = np.mean(averageflux_array) - (np.std(averageflux_array)/2)
		#Creates upper and lower bounds of 0.5 stdev above and below population mean
	
	
			
	check_array = np.ones(shape = [4, number_images])
	#Initialises array of ones so Boolean logic will accept as 'true'. Same size and shape as averageflux_array
	
	for k in range (0, (number_images)):
		for i in range(0, 4):
			if averageflux_array[i,k] > upper_bound or averageflux_array[i,k] < lower_bound: 
				check_array[i,k] = 0
				#If flux out of bounds changes check_array data point to zero so boolean logic would give false
				
	Loop_count = 1 #counts how many images are used
	totalflux = 0 #Initialises total flux variable that counts all flux
	
	for k in range (0, (number_images)): 
		for i in range(0, 4):
			if check_array[i,k]:
				#If check_array point [i,k] is true, it adds the corresponding averageflux_array point and adds to totalflux
				totalflux = totalflux + averageflux_array[i,k]
				Loop_count = Loop_count + 1
			
	averageflux = totalflux/Loop_count
	#Divides the totalflux by the number of images used.
	
	print( "Number of grids used = \n", Loop_count )
	print( "Average flux = \n", averageflux )
	print( "Standard Deviation = \n", np.std(averageflux_array) )
	Loop_count = 0
	
	return averageflux
	
	
	
"""
Readangles reads in the filenames and their parralactic 
angles from their headers and writes them to a text file.
"""

def Readangles(Number_images,name_input):
	
	file_names = ReadImageFilenames(name_input)
	#Loads the filenames into python from the textfile containing them.
	
	Number_images = len(file_names)			#Counts the length of the list.
	angles = np.zeros((1, Number_images))	#Initialises an array for the angles.
	
	#Loop creates a list of open files and reads the 'PARANG' parametre from the headers.
	for i in range(0, Number_images):
		hdulist = fits.open(file_names[i], ignore_missing_end = True, verbose = False)
		angles[0,i] = hdulist[0].header['PARANG']
	
	#Alters the angles from a 2D array into a 1D list for ease of use.
	
	np.savetxt('{name}_angles.txt'.format(name = name_input), angles)
	# Saves angles into txt
	angles = np.loadtxt('{name}_angles.txt'.format(name=name_input))
	#Loads the angles back into python.



"""	
Buildcube initialises an empty array and loads the images into 
this array, removing the badpixels and also subtracting the flatfield.
"""	

def Buildcube(Number_images,name_input):
	
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names)
	#Counts the length of the list.
	
	#Initialises two cubes to loop through.
	cube0 = np.zeros((Number_images, 1024, 1024))				
	cube = np.zeros((Number_images, 1024, 1024))				
	
	flatfield0 = './flat_Kp.fits'	#Loads in the flat-field.	
		
	flatfield = vip.fits.open_fits(
		flatfield0, n=0, header=False, 
		ignore_missing_end=True, verbose=True)
	#Opens the flat_field.
	
	Images_loaded =0	#Initialises a counter.
	
	#Loop opens the images and loads them into the data cube format required.
	#Loop also divides by the flatfield. 
	for i in range (0, Number_images):
	
		Images_loaded = Images_loaded + 1
			#Counter.
		print( "Processing {:d}/{:d}.".format(Images_loaded, Number_images) )
		
		cube0 = vip.fits.open_fits(
			file_names[i], n=0, header=False,
			ignore_missing_end=True, verbose=False)
			#Opens the fits file using VIP.
		
		for j in range(0, 1024):
			for k in range(0, 1024):
				cube[i][j][k] = cube0[j][k] / flatfield[j][k]
		
		if Images_loaded == Number_images:
			print( "\nImages loaded in.\n" )
	
	#Removes the 'bad pixels' in the image using a VIP function to do so.
	print( "Removing the badpixels in the Cube...\n" )
	
	cube = vip.preproc.badpixremoval.cube_fix_badpix_isolated(
		cube, bpm_mask=None, sigma_clip=3,
		num_neig=5, size=5, protect_mask=False,
		radius=30, verbose=True, debug=False)		
		
	#Writes the created cube to a fits file.	
	vip.fits.write_fits('Newcube.fits', cube,  verbose=True)	
	
	#Loads the angles and appends them to the cube, then overwrites the previous cube
	#with the appended one.
	   
	angles = np.loadtxt('{name}_angles.txt'.format(name=name_input))
	vip.fits.info_fits('Newcube.fits')
	
		


"""
Removeimage removes the bad images by erasing 
the name of the file from the list of filenames 
"""

def Removeimage(Number_images,name_input):

	#Loads in the filenames.
	file_names = ReadImageFilenames(name_input)
	user_input = input(
		"Type the file numbers (eg. 1,4,5.. ) of the images you would like to remove\n Seperated by single commas\n")
	
	Bad_images = [int(i) for i in user_input.split(',') if i.isdigit()]
	Bad_images = sorted(Bad_images, reverse =True)	#Sorts the numbers in descending order.
	
	#Loop removes the filenames corresponding to the numbers entered from the list of filenames.
	for i in range(0,len(Bad_images)):
		file_names = np.delete(file_names, (Bad_images[i]-1))
	
	return file_names

"""		
Checkfunction allows for yes no answers and prompts the 
user to re-enter if initially enetered incorrectly.

Yes = 0
No = 1

"""

def Checkfunction():

	check=2		
	while check == 2:	#Ensures that this will be looped over unless a yes or no is entered.
	
		yes=set(['yes', 'y', 'Y','Yes','YES'])	#Sets any of these entries to equal yes.
		no=set(['no', 'n', 'N','No','NO'])	#Sets any of these entries to equal no.
		exit=set(['EXIT','exit','exit()','Exit()','Exit','EXIT()'])
		question=input('Yes (y) or no (n)?\n')
		if question in yes:
			check=0
			
		if question in no:
			check=1	
		
		if question in exit:
			raise SystemExit()
		
		#An error if loop to prompt the user to change their response.	
		if question not in yes and question not in no:
			print( '\nUnrecognised answer\n' )
		
	return check

"""
Contrastcurvedata writes the contrast curve results into 
a text file, saves and then loads it back, plots the LLSG
and the PCA data on the same graph and saves it to a file
"""
"""
Note: This function is terribly inefficient, could be much 
better looped over but haven't had the time to do so.
"""
def Contrastcurvedata(contrast_curve, pxscale_keck, sigma, name_input, algo):

	#Saves the curve outputs
	np.savetxt('new_{algo}_{star_name}_curve_outputs'.format(algo=algo, star_name=name_input), contrast_curve)
	
	#Initialises the contrast curve variables
	size = 496

	distance = np.zeros(size)
	sensitivity_Gauss = np.zeros(size)
	sensitivity_Student = np.zeros(size)
	
	contrast_curve=np.array(contrast_curve) #changes from pandas dataframe to numpy array

	"""
    PANDAS DATAFRAME: (IF STUDENT == TRUE)
    0 - sensitivity_gaussian
    1 - sensitivity_student
    2 - throughput
    3 - distance
    4 - distance_arcsec
    5 - noise
    6 - sigma corr

    PANDAS DATAFRAME: (IF STUDENT == FALSE)
    0 - sensitivity_gaussian
    1 - throughput
    2 - distance
    3 - distance_arcsec
    4 - noise
    """

	#Loop loads the PCA outputs into the PCA variables.
	for a in range (0, size):
		distance[a] = contrast_curve[a,3] * pxscale_keck
		sensitivity_Gauss[a] = contrast_curve[a,0]
		sensitivity_Student[a] = contrast_curve[a,1]  

	#Plotting all 3 on one curve: 
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	legend = mpatches.Patch(color='#fcb141', label='%s Contrast Curve' % algo)

	plt.legend(handles=[legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.

	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	y_min = 1E-6
	y_max = 1E-2
	plt.ylim((y_min, y_max))
	
	#Creates the variables that will be plotted.
	Curve = plt.plot(distance, sensitivity_Gauss, linewidth =2, color='red')

	plt.show()
	



"""
Asks for positions of candidates/reference stars??
"""
#Input the position(x,y) and flux as array
def candidate(number):
	newcent = [512,512]
	ox = [newcent[0]]
	oy = [newcent[1]]
	cent_x = []
	cent_y = []
	position =[None]*num
	cent_x = ox*num
	cent_y = oy*num
	x = list()
	y = list()
	#f = list()
	m = None
	n = None 
	#l = None

	cent_x = ox*num
	cent_y = oy*num
	for i in range(int(num)):
		while True:
			m = input("Input the x co-ordinate of the reference/candidate position\n(Integer number)\n")

			while isinstance(m,int) == False:
				Loop_count = None
				print( "\nAt least one of these co-ordinates is not an integer.\nPlease enter an integer number\n" )
				m = input("Input the x co-ordinate of the reference/candidate position\n(Integer number)\n")

			if isinstance(m,int) == True: 
				break				

		while True:
			n = input("Input the y co-ordinate of the reference/candidate position\n(Integer number)\n")

			while isinstance(n,int) == False:
				print( "Please enter an integer number\n" )
				n = input("Input the y co-ordinate of the reference/candidate position\n(Integer number)\n")

			if isinstance(n,int) == True: 
				break	

		x.append(int(m))
		y.append(int(n))
		
	for i in range(int(num)):
		position [i] = (x[i],y[i])

	return position


"""
Converts from Cartesian to Polar coordinates
"""
def cart2pol (x,y,flux):

	rho = np.sqrt( x**2 + y**2)
	phi = np.degrees( np.arctan2(y,x) )
	phi = np.mod(phi,360)
	f_p = flux

	return rho,phi,f_p
	
	


"""
Calculates the FWHM using a Gaussian 2D fit.
"""
def Calculate_fwhm(psf):
	
	
	gauss = vip.var.fit_2dgaussian(
		psf, crop=True, cropsize=30, cent=(psf_xy[0], psf_xy[1]),
		full_output=True, debug=True
		)
	
	# 'gauss' is a pandas dataframe
	# Columns are: 
	#	0-centroid_y		1-centroid_x	2-fwhm_y	3-fwhm_x	4-amplitude		5-theta
	# Hence take columns 2 and 3 for FWHM y and x components
	fwhm_x = gauss.iloc[0,3]
	fwhm_y = gauss.iloc[0,2]
	
	fwhm = np.mean([fwhm_x, fwhm_y])
	print( "FWHM = ", fwhm )
	
	return fwhm


"""
Function takes the cube and the PSF (plus its centre)
and fits a 2d Gaussian into it.
Using that, it recenters the cube and plots the result.
Returns the re-centered cube, and the x/y shifts
"""

def Gaussian_2d_Fit( psf, cube_orig, centre):
	
	fwhm = Calculate_fwhm(psf)

	fwhm_int = math.ceil(fwhm)
	fwhm_int = int(fwhm_int)
	print( "Round up = ", fwhm_int )
	
	#Shy1 = Shift in y, shx1 = shift in x. Calculates the shifts needed to centre the image.
	#Negative = True (Bright spot at the centre) Negative = False (No bright spot at centre).

	cube1, shy1, shx1 = cube_recenter_2dfit(
		cube_orig, (centre[0], centre[1]), fwhm,
		model='gauss', nproc = 1, subi_size = fwhm_int, 
		negative=False, full_output=True, debug=False
		)
		
	#Plots the shifts in x and y used to centre the images.
	plt.figure(figsize = (8,3))								
	plt.plot(shy1, 'o-', label='shifts in y', alpha=0.5)	
	plt.plot(shx1, 'o-', label='shifts in x', alpha=0.5)
	plt.legend(loc = 'best')
	plt.grid('on')
	plt.suptitle('Positive 2d Gaussian fit')
	plt.show()
	
	return cube1, shy1, shx1, fwhm
	
""" 
End of defining functions 
"""


def CenteringLoop(cube_orig, angs, psf):

	
	#---------------------------------------------------------------
	print("Recentering Cube")
	print(" There are two ways to recenter the cube")
	print(" 1. Use a 2D Gaussian Fit")
	print(" 2. Recenter manually. (Use if Gaussian fit doesn't work. Works if all images are already aligned but not in the center of the image)")
	
	while True:
		print(" Choose 1 (Gaussian fit) or 2 (manual fit): ")
		fit_method = input()
		if fit_method == '1' or fit_method == '2':
			break
		else:
			print(" ->Option not recognised, please input '1' or '2'. ")
	
	# Initialise variable
	star_xy = [0.1, 0.1]
	
	print("Using DS9, open any image in Newcube.fits")
	star_xy[0] = input("Input the x coordinate of the central position of the star: ")
	star_xy[0] = int(star_xy[0])
		
	star_xy[1]=input("Input the x coordinate of the central position of the star: ")
	star_xy[1] = int(star_xy[1])
	
	
	#Opens the cube (using VIP) with cube_orig as HDU:0 and calls it cube_orig, 
	#and the parallactic angles from a text file.
	#Uses VIP to also open the point spread function previously loaded in.
		
	
	cube1 = cube_orig
	
	# Gaussian Fit
	if fit_method == '1':
		print(" --2D Gaussian Fit")
		print( "Fitting a 2D gaussian to centre the images..." )
		#Uses VIP's 2D gaussian fitting algorithm to centre the cube.
		cube1, shy1, shx1, fwhm = Gaussian_2d_Fit( psf, cube_orig, star_xy)
	
	# Manual Fit
	elif fit_method == '2':
		print(" --Manual Fit")
		# Calculate shifts here
		image_centre = [512, 512]
		print(" Image centre is at", image_centre)
		shift_x = image_centre[0] - star_xy[0]
		shift_y = image_centre[1] - star_xy[1]
		cube1 = vip.preproc.recentering.cube_shift(cube_orig, shift_y, shift_x)
		fwhm = Calculate_fwhm(psf)
		
	
	#Writes the values of the centered cube into a fits file.
	vip.fits.write_fits('centeredcube_{name}.fits'.format(name=name_input), cube1, verbose=True)	

	cube = cube1
	#Loads up the centered cube.
	#Plots the original cube vs the new centered cube.

	im1 = vip.preproc.cosmetics.frame_crop(cube_orig[0], 1000, verbose=False)
	im2 = vip.preproc.cosmetics.frame_crop(cube[0], 1000, verbose=False)

	hciplot.plot_frames( 
		(im1, im2), 
		label = ('Original first frame', 'First frame after recentering'), 
		grid = True, 
		size_factor = 4
		)
	
	print( "Open DS9 and look through the centred data cube at each frame making sure it is centred.")
	print( "If you're not happy with it, redo centring" )
	print( "Redo Centring?" )
	cent_loop = Checkfunction()

	return cube, cent_loop








""" 
************************************************************************************************
=============================	CODE STARTS HERE =============================
************************************************************************************************
"""



# ****************************************************************
#				MAIN FUNCTIONS
# ****************************************************************





"""
-----------------	FULL FRAME PCA	-----------------

Returns FFPCA frame and optimal number of principal components
"""
def FullFramePCA(cube, angs, fwhm, name_input):

	print( "\n --- FULL FRAME PCA")
	print( "Use DS9 to obtain the reference coordinates" )
	print( "A 'reference' is simply any kind of stable bright spot across the images, like a static speckle or a candidate.")

	reference_xy[0] = input("Input the x coordinate of the reference: ")
	reference_xy[0] = int(reference_xy[0])

	reference_xy[1] = input("Input the y coordinate of the reference: ")
	reference_xy[1] = int(reference_xy[1])


	print( 'Applying Full Frame PCA...' )

	ffpca_output = vip.pca.pca_fullfr.pca(
				cube, angs, fwhm = fwhm,
				source_xy = (reference_xy[0], reference_xy[1]),
				mask_center_px = None, ncomp = (1, len(cube) ,1),
				full_output=True
				)


	# vip_hci.pca.pca_fullfr.pca returns this outputs:
	# (if full_output=True and source_xy != None)
	# 0 - Final residuals cube (NOT TO USE IN STIM MAP)
	# 1 - PCA Frame
	# 2 - Pandas dataframe with PCs, S/Ns and fluxes
	#
	#<class 'numpy.ndarray'> 10,1024,1024
	#<class 'numpy.ndarray'> 1024,1024
	#<class 'pandas.core.frame.DataFrame'> 10 (columns = PCs, S/Ns, fluxes)

	# Save FFPCA output
	fr_pca1 = ffpca_output[1]

	#	----- Getting Optimal Number of PCs
	# Get pandas dataframe table output into a numpy array
	pca_data = ffpca_output[2].rename_axis('ID').values

	# Extract and save S/N ratio and PCs columns in arrays
	snr_data = []
	pcs_data = []

	for i in range(0, len(pca_data)):
		pcs_data.append(pca_data[i][0])
		snr_data.append(pca_data[i][1])

	# Get the index of the maximum value of the S/N column, 
	# and retrieve that same value in the PCs column
	# This will be the optimal number of principal components
	snr_max = np.argmax(snr_data, axis=0)
	optimal_pcs = int(pcs_data[snr_max])

	print("Optimal PCs", optimal_pcs)


	hciplot.plot_frames( 
			fr_pca1, 
			label = 'FFPCA reduction of {name}'.format(name=name_input), 
			grid = False, 
			size_factor = 5
			)


	#Loop asks user if they would like to save the image.
	print( 'Would you like to save the Full Frame PCA reduced images?' )
	questionsave = Checkfunction()

	if questionsave == 0:	
		vip.fits.write_fits('new_FFPCA_{name}.fits'.format(name=name_input), fr_pca1, verbose=True) 


	return fr_pca1, optimal_pcs








"""
-----------------	ANNULAR PCA	-----------------

Returns annular pcs frame, and residuals cube for STIM map
"""
def AnnularPCA(cube, angs, fwhm, name_input):

	print( "\n --- ANNULAR PCA")

	print("Performing Annular PCA...")

	ann_pca_output = vip.pca.pca_local.pca_annular(
							cube, angs, fwhm=fwhm,
							full_output = True, verbose = True
							)

	# Outputs three objects in this order:
	#	cube_out: cube of residuals
	#	cub_der: derotated cube of residuals
	#	frame: Annular PCA frame

	hciplot.plot_frames( 
			ann_pca_output[2], 
			label = 'AnnPCA reduction of {name}'.format(name=name_input), 
			grid = False, 
			size_factor = 5
			)

	print( 'Would you like to save the Annular PCA reduced images?' )
	questionsave = Checkfunction()

	if questionsave == 0:	
		vip.fits.write_fits('AnnPCA_{name}.fits'.format(name=name_input), ann_pca_output[2], verbose=True)
		vip.fits.write_fits('AnnPCA_residuals_{name}.fits'.format(name=name_input), ann_pca_output[1], verbose=True) 

	annpca_frame = ann_pca_output[2]
	residuals_stim_input = ann_pca_output[1]

	return annpca_frame, residuals_stim_input






"""
-----------------	LLSG	-----------------

Returns LLSG frame, and residuals cube for STIM map
"""

def LLSG(cube, angs, fwhm, name_input):

	print( "\n   --- LLSG")
	rank = input("What rank of llsg would you like to do? (default = 6)\n")
	rank = int(rank)
	
	#Reduces the image using the LLSG algorithm then plots the image.
	print(fwhm)
	llsg_output = vip.llsg.llsg(cube, angs, fwhm=fwhm, 
							rank = rank, thresh = 2, 
							max_iter = 20, random_seed = 10,
							full_output = True)


	# LLSG full output is:
	#0	list_l_array_der
	#1	list_s_array_der
	#2	list_g_array_der
    #3	frame_l
    #4	frame_s
    #5	frame_g
    # Each is the frame and residual cube for the three noise types:
    # l (low-rank), s (sparse), g (gaussian)
    # Companions are in the Sparse noise.
    # To get STIM map of LLSG, input 'list_s_array_der' into STIM Map

	llsg_frame = llsg_output[4]
	llsg_residual_sparse = np.asarray( llsg_output[1][0] )
		#list_s_array_der is a list of one member,
		# this one member is another list of the 24 residual images
		# So get that one element, and make it an array.

	hciplot.plot_frames(llsg_frame, label ='LLSG reduced image of {name}'.format(name=name_input))
	
	
	#Loop asks user if the would like to save the image.
	print( 'Would you like to save the LLSG reduced image?\n' )
	questionsave = Checkfunction()
	
	if questionsave == 0:
		vip.fits.write_fits('new_LLSG_{name}.fits'.format(name=name_input), llsg_frame, verbose=True)


	return llsg_frame, llsg_residual_sparse, rank





"""
-----------------	STIM MAP	-----------------
"""
def StimMap(residuals_cube, name_input, origin_algo):

	print("Computing STIM map using {algo} output...".format(algo = origin_algo))
	stim_map = vip.metrics.compute_stim_map(residuals_cube)


	hciplot.plot_frames( 
			stim_map, 
			label = 'STIM map of {name} using {algo}'.format(name=name_input, algo = origin_algo), 
			grid = False, 
			size_factor = 5
			)

	#Loop asks user if they would like to save the image.
	print( 'Would you like to save the image?' )
	questionsave = Checkfunction()

	if questionsave == 0:	
		vip.fits.write_fits('STIM_{algo}_{name}.fits'.format(algo=origin_algo, name=name_input), stim_map, verbose=True)







"""
-----------------	CONTRAST CURVES	-----------------
"""
# 'param' input refers to optimal PCs if using PCA, and the LLSG rank if using LLSG.
# 'algo' refers to the algorithm used: 'PCA', 'annPCA', or 'LLSG'.
def ContrastCurves(cube=None, angs=None, psf=None, fwhm=None, starname=None, param=None, algo=None):

	print(" Performing Starphot calculation...")
	#Starphot is the  total aperture flux from the star during
	#the integration time (if no coronograph).
	starphot = StarphotCalculation(cube.shape[0], psf, fwhm)
	starphot_array = np.zeros(500)

	starphot_array = np.full(500, starphot)
	np.savetxt('%s_starphot' % starname, starphot_array)


	print(" Building contrast curves...")

	psf_norm = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=101, threshold=None, mask_core=None)

	contrast_curve = None

	if (algo == 'PCA'):
		contrast_curve = vip.metrics.contrcurve.contrast_curve(
			cube, angs, psf_norm, fwhm, pxscale_keck, starphot,
			sigma=sigma, nbranch=1, algo=vip.pca.pca, 
			ncomp=param, debug=True, save_plot='PCA'
			)
		Contrastcurvedata(contrast_curve, pxscale_keck, sigma, name_input, 'PCA')
	
	elif (algo == 'annPCA'):
		contrast_curve = vip.metrics.contrcurve.contrast_curve(
	    	cube, angs, psf, fwhm, pxscale_keck, starphot,
	    	sigma=sigma, nbranch=1, algo=vip.pca.pca_local.pca_annular, 
	    	debug=True,save_plot='AnnPCA', verbose=2, student=True,
	    	ncomp=param
	    	)
		Contrastcurvedata(contrast_curve, pxscale_keck, sigma, name_input, 'annPCA')

	elif (algo == 'LLSG'):
		contrast_curve =  vip.metrics.contrcurve.contrast_curve(
			cube, angs, psf_norm, fwhm, pxscale_keck, starphot,
			sigma=sigma, nbranch=1, algo=vip.llsg.llsg, 
			debug=True,save_plot='LLSG', rank = param
			)
		Contrastcurvedata(contrast_curve, pxscale_keck, sigma, name_input, 'LLSG')









"""
-----------------	ANDROMEDA	-----------------
"""
def Andromeda(cube, angs, psf, fwhm):
	
	print(" Running ANDROMEDA algorithm...")
	print(" (Takes around 45 min to 1h 30min)")
	print(" (So sit back, relax, and have a cuppa)")
	
	# Normalises PSF (Done above too)
	psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=100,
										threshold=None, 
										mask_core=None,
										force_odd = False)

	print("Cube has shape", cube.shape)
	print("PSF has shape", psf.shape)

	# Andromeda requires even-sized images, hence "force_odd = False" is used.


	# In NIRC2, the Kp filter has central wavelength 2.124 um
	# Source: https://www2.keck.hawaii.edu/inst/nirc2/filters.html
	Kp_wavelength = 2.124E-6
	
	# Keck Telescope has an aperture diameter of 10 meters
	diam_tel = 10
	
	# Pixel scale in mas/pixel
	pxscale = pxscale_keck * 1E3
	
	# Nyquist pixel scale at Shannon wavelength in mas/pixel
	pxscale_nyq = 1/2 * (180*3600*1E3 / np.pi) * (Kp_wavelength / diam_tel)
	
	oversampling_factor = pxscale_nyq / pxscale

	# Set min_sep to 0.36 because it works.
	andromeda_output = vip.andromeda.andromeda(
						cube,
						oversampling_factor,
						angles = angs,
						psf = psf,
						iwa = 1,
						min_sep = 0.36,
						verbose = True
						)

	print("Andromeda output length:", len(andromeda_output))

	contrast_andr = andromeda_output[0]			# Contrast map / flux map
	snr_map_andr = andromeda_output[1]			# SNR map
	snr_norm_map_andr = andromeda_output[2]		# Normalised SNR Map
	stdcontrast_map_andr = andromeda_output[3]	# Standard Deviation Contrast map
	stdcontrast_norm_andr = andromeda_output[4]	# Normalised Standard Deviation Contrast map
	likelihood_andr = andromeda_output[5]		# Likelihood map
	ext_radius_andr = andromeda_output[6]		# External/cut-off processing radius

	vip.fits.write_fits('ANDR_snrmap_{name}.fits'.format(name=name_input), snr_map_andr, verbose=True)
	vip.fits.write_fits('ANDR_snrnorm_{name}.fits'.format(name=name_input), snr_norm_map_andr, verbose=True)
	vip.fits.write_fits('ANDR_contrast_{name}.fits'.format(name=name_input), contrast_andr, verbose=True)
	vip.fits.write_fits('ANDR_contrast_stddev_{name}.fits'.format(name=name_input), stdcontrast_map_andr, verbose=True)
	vip.fits.write_fits('ANDR_contrastnorm_stddev_{name}.fits'.format(name=name_input), stdcontrast_norm_andr, verbose=True)
	vip.fits.write_fits('ANDR_likelihood_{name}.fits'.format(name=name_input), likelihood_andr, verbose=True)

	print("ext_radius = ", ext_radius_andr)

		



#*******************************************************************
#*******************************************************************
#*******************************************************************
#*******************************************************************




# ****************************************************************
#						PRE PROCESSING
# ****************************************************************


print("\n\n\n ================= Project Script =================")
print("\n Version 12 \n\n")



	
"""
------ 1 - STAR NAME AND PRE-PROCESSING LOOP
"""
name_input = input('Name of the star (example = HIP_544):  ')

#Loads the filenames into python and counts the number of files.
file_names = ReadImageFilenames(name_input)
Number_images = len(file_names)


"""
	------ 2 - CREATE DATA CUBE
"""
print( "Would you like to create the data cube?")

#Loop runs the Readangles and Buildcube functions, building the data cube.
if Checkfunction() == 0:
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names)
	Readangles(Number_images,name_input)
	Buildcube(Number_images,name_input)


"""
------ 3 - REMOVE BAD IMAGES
"""
print( "Would you like to remove any bad frames?" )

#Loop runs the Removeimage function and then overwrites the new shortened list to
#the original text file. 
if Checkfunction() == 0:		
	file_names = Removeimage(Number_images,name_input)
	file_names = file_names.tofile('{name}_filenames.txt'.format(name=name_input), sep ="\n", format='%s')
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names)
	Readangles(Number_images,name_input)
	Buildcube(Number_images,name_input)


"""
------ 4 - RECENTERING
"""

# Retrieving angles and PSF
angs = np.loadtxt( '{name}_angles.txt'.format(name = name_input) )
psf = vip.fits.open_fits('./psf.fits', n=0, header=False, ignore_missing_end=True, verbose=False)

print(" Would you like to recenter the cube? ")
cent_loop = Checkfunction()

if cent_loop == 0:
	cube_orig = vip.fits.open_fits('./Newcube.fits')
	re_loop = 0
	while re_loop == 0:
		cube, re_loop = CenteringLoop(cube_orig, angs, psf)

elif cent_loop == 1:
	cube = vip.fits.open_fits('centeredcube_{name}.fits'.format(name=name_input))
	#Cropping cube to even size (1022x1022)
	cube = vip.preproc.cosmetics.cube_crop_frames(
									cube, 1022, verbose=False, force=True)

fwhm = Calculate_fwhm(psf)
star_xy = [512.0, 512.0]


print("Cube =", cube.shape )
print("PSF =", psf.shape )

	#---------------------------------------------------------------





# ****************************************************************
#						POST PROCESSING
# ****************************************************************




# Full Frame PCA
#ffpca_frame, optimal_pcs = FullFramePCA(cube, angs, fwhm, name_input)

# Annular PCA
#annpca_frame, pca_residuals = AnnularPCA(cube, angs, fwhm, name_input)



# LLSG
#llsg_frame, llsg_residuals, llsg_rank = LLSG(cube, angs, fwhm, name_input)

# STIM Map
#StimMap(pca_residuals, name_input, "PCA")
#StimMap(llsg_residuals, name_input, "LLSG")

# Contrast Curves
#ContrastCurves(cube=cube, angs=angs, psf=psf, fwhm=fwhm, starname=name_input, param=optimal_pcs, algo='PCA')
#ContrastCurves(cube=cube, angs=angs, psf=psf, fwhm=fwhm, starname=name_input, param=optimal_pcs, algo='annPCA')
llsg_rank = 6
ContrastCurves(cube=cube, angs=angs, psf=psf, fwhm=fwhm, starname=name_input, param=llsg_rank, algo='LLSG')

# Andromeda
#Andromeda(cube, angs, psf, fwhm)


print( "===========================	End of the program	===========================" )










# Some space here
