"""
Original code written by Joe Callingham and Martje Slob, adapted by Femke Ballieux
Crossmatching LoTTs unresolved isolated sources, with some other surveys
"""
#run on laptop
import os


os.system('java -jar topcat-full.jar -stilts \
          tmatchn matcher=sky multimode=pairs nin=2 params=15 \
    in1=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\unresolved_isolated_S_source.fits  values1="RA DEC" \
    in2=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\NVSS.fits values2="RAJ2000 DEJ2000" \
    out=C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\crossmatch_NVSS_LoTSS.fits')

print("yey we did it")