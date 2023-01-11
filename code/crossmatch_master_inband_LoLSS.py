"""
Code written by Femke Ballieux, takes in the master sample and crossmatches it with LoLSS inband spectra
"""
#run on laptop
import os


os.system('java -jar topcat-full.jar -stilts \
          tmatchn join1=always matcher=sky multimode=pairs nin=8 params=15 \
    in1=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\fit_vals_power_law_NVSS_intflux.fits values1="LoLSS_RA LoLSS_Dec" \
    in2=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_0_source.fits values2="RA DEC" \
    in3=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_1_source.fits values3="RA DEC" \
    in4=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_2_source.fits values4="RA DEC" \
    in5=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_3_source.fits values5="RA DEC"\
    in6=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_4_source.fits values6="RA DEC"\
    in7=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\LoLSS_inband\\channel_5_source.fits values7="RA DEC" \
    in8=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\WENSS.fits values8="_RAJ2000 _DEJ2000" \
    out=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\crossmatch_master_inbandLoLSS.fits')

print("yey we did it")