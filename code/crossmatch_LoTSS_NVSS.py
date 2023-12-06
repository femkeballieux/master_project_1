
#run on laptop
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import matplotlib.pyplot as plt


os.system('java -jar topcat-full.jar -stilts \
          tmatchn matcher=sky multimode=pairs nin=2 params=10\
    in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/unresolved_isolated_S_source_45.fits  values1="RA DEC" \
    in2=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/NVSS.fits values2="RAJ2000 DEJ2000" \
    out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_NVSS_LoTSS.fits')

print("yey we did it")




# hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_NVSS_LoTSS.fits') #For justifying the crossmatching radius
# tbdata = hdulist[1].data
# orig_cols = hdulist[1].columns
# hdulist.close()
# print(orig_cols)
# plt.figure(figsize=(8,8))
# print(len(tbdata['RA']))

# delta_RA =(tbdata['RA'] - tbdata['RAJ2000']) *np.cos(np.deg2rad(tbdata['DEC']))
# delta_Dec = tbdata['DEC'] - tbdata['DEJ2000'] 

# plt.scatter(delta_RA*3600, delta_Dec*3600, alpha=0.5, s=0.5)
# plt.xlabel('Delta RA')
# plt.ylabel('Delta Dec')
# plt.grid()

# r = [1,2,4,6,7,8,10,12,14,16, 18,20]
# N=np.array([26258,56105,96084,121268,130158,137100,146975,153210,157226,159905, 161727,163011])

# c = (2493565*1773484*5635 /(41253*0.82*5635*3600*3600)) * np.square(r) * np.pi
# plt.figure()
# plt.plot(r,N, label='N (total number of associations)')
# plt.plot(r,c, label='C (chance associations)')
# plt.plot(r,N-c, label='N-C')
# # plt.plot(r,c*100/N, label='C/N')
# # plt.hlines(1,0,20, label='1%')
# # plt.yscale('log')
# # plt.xscale('log')
# plt.legend()
# plt.ylabel('number of sources')
# plt.xlabel('r in arcsec')