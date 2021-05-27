"""
    PROJECT NAME: CONTRAST CURVE
    
    DESCRIPTION: Turns contrast/radius curves to mass/radius
    Edited version of original mass_sensitivity.py so that it can be run on
    analysis.
    
    Input Files:
        text file containing Baraffe+ 2003 isochrones
        csv file containing the columns:
            star_name, epoch, best_Age(Myr), oldest_Age(Myr), youngest_age(Myr), 
            distance(pc), apparent magnitude (k-band)
    
    VERSION	3.1
    DATE 19/04/21
    
    """

import sys
import numpy as np
import math
import csv
from scipy import interpolate
from scipy.interpolate import griddata
import warnings
warnings.simplefilter("ignore")

"""
    DEFINED VARIABLES
    """         
INPUT_BARAFFE_DATA = 'baraffe_final.txt' #data from baraffe
INPUT_NAMES = 'star_names.csv' #list of stars with information specified above
METHOD = 'RDI' #method used to create that perticular contrast curve
seperation_column = 4                      #column of input data containg seperation in arcsec
contrast_column = 1                       #cloumn of input data containing contrast
ONE_PARSEC = 206268.0                    #AU scale to use on axis of output
JUP_MASS = 1047.34
PLOT_PRECISION = 100               #defines grid side for interpolation
SAVEFIG = True                    # Change to 'False' to show plots of mass/radius while program is running (does not save fig)
SHOW_SURFACE = False #Change to 'True' if plot of surface of baraffe data is required. requires mpl_toolkits to be available

"""
    FUNCTION DEFINITIONS
    """

#loads contrast data from file into array and returns orbital seperation and mag
def load_contrast_data(name,epoch):
    
    path= "inputs/" + name + "_" + METHOD + "_contrast_curve_" + epoch + "_new.txt" #format of contrast curve file name, change accordingly
    contrast_curve=np.zeros((497,2)) #each _new contrast curve file has 497 lines --> if the RDI script is changed this may be different so double check
    with open(path) as f:
        lines = f.readlines()
        contrast_curve[:,0] = [float(line.split()[seperation_column]) for line in lines]
        contrast_curve[:,1] = [float(line.split()[contrast_column]) for line in lines]
    
    return contrast_curve


#returns absolute magnitude of star
def absolute_mag(distance, mag):
    M = mag - ((5 * math.log10(distance)) - 5.0)
    return M

#returns the magnitude of the companion as array
#contrast_data is the sensitivity column of contrast output
#abs_mag is absolute magnitude of star
def companion_mag_curve(contrast_data, abs_mag):
    return (np.absolute(2.5 * np.log10(contrast_data)) + abs_mag)

#turns angular seperation into physical speration and returns
def orbital_sep(angular_sep, distance):
    rad_coef = math.pi /180.0 #convert from degrees to radians
    arc_sec_coef = 1 / 3600.0 #convert from arcseconds to degrees
    angular_sep = angular_sep * rad_coef * arc_sec_coef

    return distance * np.tan(angular_sep) * ONE_PARSEC

#plots or saves plot of mass against orbital seperation sensitivity curve
def plot_rad_mass(rad_mass,file_dest, name):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 7))
    plt.rc('font', family='serif', size='20')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='18')
    plt.rc('ytick', labelsize='18')
    plt.clf()
    best, = plt.plot(rad_mass[:,0], rad_mass[:,1],'r', label='Best age estimate')
    oldest, = plt.plot(rad_mass[:,0],rad_mass[:,2],'b', label='Oldest age estimate')
    youngest, = plt.plot(rad_mass[:,0],rad_mass[:,3],'g', label='Youngest age estimate')

    plt.legend(loc='upper right')
    plt.ylabel(r'Mass($M_{jup}$)')
    plt.xlabel('Orbital separation (AU)')
    plt.title('Mass - Orbital Separation Sensitivity \n' + name.replace("_"," "))
    if SAVEFIG:
        plt.savefig(file_dest,dpi=300)
    else:
        plt.show()


#finds nan's in interpolated grid and changes them to nearest neighbur value (flatterns area) but stops interpolation errors
def check_nan(grid):
    lum=0.0
    for i in range(0,PLOT_PRECISION):
        for j in range(0,PLOT_PRECISION):
            if np.isnan(grid[i,j]):
                grid[i,j]=lum
            else:
                lum = grid[i,j]
    return grid

#does chi suqered test of data to determine how good a fit function is to baraffe data
def chi_squared(points, lum, func):
    chi = 0.0
    for i in range(0,len(points)):
        chi += ((lum[i]-func(points[i,0], points[i,1]))**2)/func(points[i,0], points[i,1])
    return chi

#saves figures and output data to seperate files (creatng file if not present)
def save_data(name,epoch,rad_mass):
        
    path="outputs/" + METHOD + "_" + name + "_" + epoch + "_mass_sensitivity"
    
    
    data_file = path + "_data.txt"
    fig_file = path + ".png"

    plot_rad_mass(rad_mass, fig_file, name)
    np.savetxt(data_file, rad_mass)

def plot_surface_baraffe(ai_grid, mi_grid, lumi):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig=plt.figure()
    ax=fig.gca(projection='3d')
    ax.plot_surface(ai_grid, mi_grid,lumi)
    plt.show()

"""
    CODE IMPLENTATION
    """
if __name__ == "__main__":
    
    with open(INPUT_NAMES, 'r') as fd:
      name_data = list(csv.reader(fd))
    
    name_data=np.array(name_data) #converts list to array
    baraffe_data = np.loadtxt(INPUT_BARAFFE_DATA)
    
    points = np.zeros((len(baraffe_data),2))        #seperates data points from baraffe into points(age,mass) and magnitude list, lum
    points[:,0] = baraffe_data[:,0]
    points[:,1] = baraffe_data[:,1]
    lum = np.zeros((len(baraffe_data),1))
    lum = baraffe_data[:,2]
    
    ai = np.linspace(np.amin(points[:,0]),np.amax(points[:,0]),num = PLOT_PRECISION, endpoint=True) #Sets up mesh of values to interpolate surface from the scattered data points from baraffe
    mi = np.linspace(np.amin(points[:,1]),np.amax(points[:,1]),num = PLOT_PRECISION, endpoint=True) #the precision of the grid i set by PLOT_PRECISION (ai=age, mi=mass)
    ai_grid, mi_grid = np.meshgrid(ai, mi)
    
    lumi = griddata(points, lum, (ai_grid, mi_grid), method='cubic')  #interpolate from baraffe data points to surface grid lumi
    lumi = check_nan(lumi)                                            #remove any nan values form grid
    func_lum = interpolate.interp2d(ai, mi, lumi,kind='quintic')      #create function from grid
    
    print('Chi Squared test of function verse known data points =')   #print results of chi squared test
    print(chi_squared(points, lum, func_lum))
    
    if SHOW_SURFACE:
        plot_surface_baraffe(ai_grid, mi_grid, lumi)
    
    
    for i in range(0,len(name_data)):                   #loop over all stars in file outputing data for each star
        
        contrast_curve = load_contrast_data(name_data[i,0],name_data[i,1])
        
        mag_curve  = np.zeros((len(contrast_curve),2))      #creates array of orbital seperation and companion magnitude
        mag_curve[:,0] = orbital_sep(contrast_curve[:,0], float(name_data[i,5]))
        mag_curve[:,1] = companion_mag_curve(contrast_curve[:,1], absolute_mag(float(name_data[i,5]), float(name_data[i,6])))
        
        star_age = np.zeros((3))
        
        rad_mass = np.zeros((len(mag_curve),4)) #creates new array to contain mass and orbital seperation
        rad_mass[:,0] = mag_curve[:,0]
        
        for j in range(0,3):
            star_age[j] = np.log10(np.multiply(1.0e6,float(name_data[i,j+2]))) #loads star age
            massi = np.linspace(np.amin(points[:,1]),np.amax(points[:,1]),num=PLOT_PRECISION, endpoint=True) #creates new linespace of mass to create function of mass in terms of magnitude
            
            lumi = func_lum(star_age[j], massi)        #takes slice of baraffe data surface defined by star age
            check = False
            if lumi[0] > lumi[-1]:
                lumi = lumi[::-1]
                massi = massi[::-1]
                check=False
        
          
            func_mass = interpolate.interp1d(np.ravel(lumi), massi, kind='linear', fill_value='extrapolate') #interpolates function mass(magnitude) from baraffe data slice and new mass linespace
            
            mag = mag_curve[:,1]
       
            #if np.amin(mag) < np.amin(lumi):
             #   print('Warning data below Baraffe model bounds....Extrapolating')
                        
            #if np.amax(mag) > np.amax(mag):
             #   print('Warning data above Baraffe model bounds....Extrapolating')
    
            rad_mass[:,j+1] = np.multiply(JUP_MASS,func_mass(mag)) #uses function mass(magnitude) to output mass for a give magnitude for each point of the load ontrast curves from VIP
        
        save_data(name_data[i,0],name_data[i,1],rad_mass)#saves all relavent data
    
