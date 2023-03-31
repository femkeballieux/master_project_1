"""
Written by Femke Ballieux
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

# os.system('java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts \
#           tmatchn matcher=sky multimode=pairs nin=2 params=15 \
#     in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/official_VLASS_no_duplicates.fits values1="RA DEC" \
#     in2=/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/vlass_catalog_highest.fits values2="ra dec" \
#     out=/net/vdesk/data2/bach1/ballieux/master_project_1/compare_VLASS/crossmatch_VLASS_official_mine.fits')


#read in the data
path = '/net/vdesk/data2/bach1/ballieux/master_project_1/compare_VLASS/'
hdulist = fits.open(path + 'crossmatch_VLASS_official_mine.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()
print(orig_cols)


plt.figure(figsize=(10,10))
official_flux=tbdata['Total_flux']/1000
my_flux = tbdata['int_flux']
# mask=np.where(np.logical_or((tbdata['S_code']=='S'), (tbdata['S_code']=='M') ) )
mask=[tbdata['S_code']=='S']
plt.scatter(official_flux[mask], my_flux[mask], alpha=0.2, s=2, zorder=2)
plt.xlabel('official flux [Jy]')
plt.ylabel('mine flux [Jy]')
plt.plot(official_flux, official_flux, color='red', zorder=1)
plt.xscale('log')
plt.yscale('log')
plt.xlim(4e-4, 5)
plt.ylim(4e-4, 5)


plt.figure(figsize=(8,8))


SNR = tbdata['Peak_flux_1']/tbdata['isl_rms']

plt.scatter(SNR, official_flux/my_flux, alpha=0.1, s=2, zorder=2)
plt.xlabel('SNR')
plt.ylabel('ratio official/mine')
plt.hlines(1,5,4000, color='red', zorder=1)
#plt.vlines(7, 1e-4, 10, color='red', label='SNR=7', zorder=1)
plt.legend()
plt.xlim(xmin=4, xmax=5000)

plt.xscale('log')
plt.yscale('log')

plt.figure(figsize=(8,8))
SNR = tbdata['Total_flux']/tbdata['isl_rms']

plt.scatter(SNR, official_flux/my_flux, alpha=0.1, s=2, zorder=2)
plt.xlabel('SNR')
plt.ylabel('ratio official/mine')
plt.hlines(1,5,4000, color='red', zorder=1)
#plt.vlines(7, 1e-4, 10, color='red', label='SNR=7', zorder=1)
plt.legend()
plt.xlim(xmin=4, xmax=5000)

plt.xscale('log')
plt.yscale('log')


