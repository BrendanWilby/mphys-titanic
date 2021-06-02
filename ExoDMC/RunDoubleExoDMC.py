"""
DETECTION PROBABILITY MAP

DATE: 25/05/21

WRITTEN BY: Sophia Stasevic

Runs ExoDMC, which calculates mass and separation detection probability from the contrast curve 
output from ADI/RDI that has been converted to a mass sensitivity curve

Plots the 68% and 95% probability contours for the comparison of two different methods,
with the full map of method 1 being plot, and just the contours for method 2 overlaid

Script adapted from DMC_DImode.py provided by Mariangela Bonavita at https://github.com/mbonav/Exo_DMC

Input Files:
        
        text file containing output of mass_sensitivity.py
        csv file containing the columns:
            star_name, epoch, best_Age(Myr), oldest_Age(Myr), youngest_age(Myr), 
            distance(pc), apparent magnitude (k-band)
    
Change lines 189-191 to reflect the file name containing list of stars, and post-processing method used on the data:
    for a comparison of two methods, change lines 190 + 191 to the methods used 
    for a full map of a single method, change 'method_2' on line 191 to say 'SINGLE'
    
"""

from DMC import *
import DoubleExoDMC as DEDMC

import numpy as np
import pandas
import csv
import re
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import ticker
from matplotlib import colors as mcolors
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate
from numpy import random as rn
import scipy.ndimage as ndimage
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

NAME_COL = 0
EPOCH_COL = 1
DIST_COL = 5

def save_map_data(name,epoch,method,sma,M2,mapT):
    
    save_path= "outputs/" + method + "_" + name + "_" + epoch + "_det_sensitivity.txt"
    sma=sma.T #transposes array
    M2=M2.T
    axes=np.stack((sma,M2))
    save_data= np.concatenate((axes,mapT))
    
    np.savetxt(save_path,save_data)

def plot_map(name,epoch,dist,ID,method,sma,M2,mapT):
    
    fig = plt.figure(figsize=(10, 8))
    plt.rc('font', family='serif', size='20')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='18')
    plt.rc('ytick', labelsize='18')
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.7])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(" Mass (M$_{Jup}$)")
    ax.set_xlabel(" Semi major axis (AU) ")
    
    levels=[0,20,40,60,80,90,95,100]
    norm = mcolors.Normalize(0, 100)
    cont=plt.contour(sma, M2, mapT, levels, cmap='bone', zorder=10, linewidths=2)
    
    cf0=ax.contourf(sma, M2, mapT, norm=norm, levels=np.arange(0,100,0.1), extend='neither', cmap=plt.cm.Spectral_r, antialiased=False, zorder=0)
    CB = plt.colorbar(cf0,  extend='both', cmap=plt.cm.Spectral_r, ticks=levels)
    CB.add_lines(cont)
    CB.set_ticks(levels)
    CB.ax.set_yticklabels(["{:.0f}".format(i) for i in CB.get_ticks()]) # set ticks of your format   .format(i)+"%" for percentage in cb label
    CB.ax.set_ylabel('Detection Probability [%]', rotation=270,linespacing=5,fontsize=20,labelpad=25)               
    
    fmt={}
    string="%"
    strs=["{}{}".format(i,string) for i in levels]
    for l, s in zip(levels, strs):
        fmt[l] = s
                
    h1,_ = cont.legend_elements()
    plt.legend([h1[0]], [method],loc='lower left')
        
    plt.clabel(cont,levels,inline=True, manual=False, colors = 'bone', fmt=fmt, fontsize=18)	
    plt.title(name.replace("_"," ") + " " + method + ' Sensitivity Map ')
    plt.savefig("plots/" + name + "_" + epoch + "_detprob_" + method+ "_.png", dpi=300) 
    
def plot_comparison_map(name,epoch,dist,ID,method_1,method_2,sma_1,M2_1,mapT_1,sma_2,M2_2,mapT_2):
    
    fig = plt.figure(figsize=(10, 8))
    plt.rc('font', family='serif', size='20')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='18')
    plt.rc('ytick', labelsize='18')
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.7])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(" Mass (M$_{Jup}$)")
    ax.set_xlabel(" Semi major axis (AU) ")
    
    levels=[68,95]
    norm = mcolors.Normalize(0, 100)
    cont_1=plt.contour(sma_1, M2_1, mapT_1, levels, cmap='bone', zorder=10, linewidths=2)
    cont_2=plt.contour(sma_2, M2_2, mapT_2, levels, cmap='bone', zorder=10, linewidths=2, linestyles='--', alpha=0.7)
    
    cf0=ax.contourf(sma_1, M2_1, mapT_1, norm=norm, levels=np.arange(0,100,0.1), extend='neither', cmap=plt.cm.Spectral_r, antialiased=False, zorder=0)
    CB = plt.colorbar(cf0,  extend='both', cmap=plt.cm.Spectral_r, ticks=levels)
    CB.add_lines(cont_1)
    CB.set_ticks(levels)
    CB.ax.set_yticklabels(["{:.0f}".format(i) for i in CB.get_ticks()]) # set ticks of your format   .format(i)+"%" for percentage in cb label
    CB.ax.set_ylabel('Detection Probability [%]', rotation=270,linespacing=5,fontsize=20,labelpad=25)               
    
    fmt={}
    string="%"
    strs=["{}{}".format(i,string) for i in levels]
    for l, s in zip(levels, strs):
        fmt[l] = s
                
    h1,_ = cont_1.legend_elements()
    h2,_ = cont_2.legend_elements()
    plt.legend([h1[0], h2[0]], [method_1, method_2],loc='lower left')
        
    plt.title(name.replace("_"," ") + ' Sensitivity Map ')
    plt.savefig("plots/" + name + "_" + epoch + "_detprob_comp_" + method_1 + "_" + method_2 + ".png", dpi=300) 


def run_exo_dmc(name,epoch,dist,ID,method_1,method_2):
    
            
    map=DEDMC.exodmc(ID, dist)       
    # the set_grid method allows to change range or resolution of the grid 
    map.set_grid(x_min=1, x_max=1000, y_min=10, y_max=100, logx=True, logy=True)
    
    path="inputs/" + method_1 + "_" + name + "_" + epoch + "_mass_sensitivity_data.txt"       
    data=np.loadtxt(path)
    xlim=data[:,0] # separation (arcsec)
    ylim=data[:,2] # mass (mjup )    
            
    prob = map.DImode(xlim, ylim) #creates detection probability map
    sma_1=prob[0] #semi major axis
    M2_1=prob[1] #mass
    mapT_1=prob[2] #probability
    
    save_map_data(name,epoch,method_1,sma_1,M2_1,mapT_1)
    
    #if comparing detection probability for 2 different methods, this repeats the above
    #process for the second method
    if method_2 != 'SINGLE':
        
        #probably don't need to repeat this part of exodmc but I may be wrong so if it
        #stops working just uncomment the two lines below:
        #map=DEDMC.exodmc(ID, dist)
        #map.set_grid(x_min=1, x_max=1000, y_min=10, y_max=100, logx=True, logy=True)
                
        path="input/" + method_2 + "_" + name + "_" + epoch + "_mass_sensitivity_data.txt"       
        data=np.loadtxt(path)
        xlim=data[:,0] # separation (arcsec)
        ylim=data[:,2] # mass (mjup ) 
                
        prob = map.DImode(xlim, ylim)
        sma_2=prob[0]
        M2_2=prob[1]
        mapT_2=prob[2]
        
        save_map_data(name,epoch,method_2,sma_2,M2_2,mapT_2)
        plot_comparison_map(name,epoch,dist,ID,method_1,method_2,sma_1,M2_1,mapT_1,sma_2,M2_2,mapT_2)
        
    else:
        
        plot_map(name,epoch,dist,ID,method_1,sma_1,M2_1,mapT_1)
    

if __name__ == "__main__":
    
    input_names='star_names_trimmed.csv'
    method_1='RDI'
    method_2='ADI' # ='SINGLE' for only one method
    
    name_data=np.loadtxt(input_names,delimiter=',',dtype='str')
    name=name_data[:,NAME_COL]
    distance=name_data[:,DIST_COL]
    epoch=name_data[:,EPOCH_COL]
    
    for i in range(0,len(name)):
        
        temp = re.findall(r'\d+', name[i]) #finds just the number part of the star name
        res = list(map(int, temp))
        ID=np.array(res)
        
        run_exo_dmc(name[i],epoch[i],float(distance[i]),ID[0],method_1,method_2)
