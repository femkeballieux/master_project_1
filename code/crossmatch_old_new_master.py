"""
Code written by Femke Ballieux, takes in the old master sample and new master sample
crossmatches it, to make sure all of Martjes PS sources were found
"""
#run on laptop
import os

#Isolation radius of 47 arcsec, this does mean 15 arcsec in the params parameter

os.system('java -jar topcat-full.jar -stilts \
          tmatchn join1=always join2=always matcher=sky multimode=pairs nin=2 params=15 \
    in1=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\fit_vals_power_law_NVSS_intflux.fits values1="RA Dec" \
    in2=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\data_before_new_LoLSS\\fit_vals_power_law_NVSS_intflux.fits values2="RA Dec" \
    out=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\crossmatch_old_new_master.fits')

print("yey we did it")