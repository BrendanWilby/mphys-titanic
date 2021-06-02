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
