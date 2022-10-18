"""
Code written by Femke Ballieux, takes in the preliminary and the first LoLSS data release
crossmatches it to get a bit more insight in the discrepancy
"""
#run on laptop
import os

#Isolation radius of 47 arcsec, this does mean 15 arcsec in the params parameter

os.system('java -jar topcat-full.jar -stilts \
          tmatchn join1=match join2=match matcher=sky multimode=pairs nin=2 params=15 \
    in1=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\LoLSS_new_names.fits values1="RA DEC" \
    in2=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\LoLSS.fits values2="RA DEC" \
    out=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\compare_old_new_LoLSS\\crossmatch_old_new_LoLSS.fits')

print("yey we did it")