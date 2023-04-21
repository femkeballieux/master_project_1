# !/bin/python
# Measuring the rms in Vlotss

import numpy as np
import os
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
import scipy.optimize as opt


fields = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/rms_filenames.txt', dtype='U')

# rms_meas_2deg = np.zeros(len(fields))
# centre_ra_store = np.zeros(len(fields))
# centre_dec_store = np.zeros(len(fields))
# i = 0
# for field in fields:

#     print(i, 'Processing the rms of field '+field[2:15],'...')

#     # Reading in Stokes V data
#     try:
#         hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/' + field[2:])
#     except:
#         print('Bad field '+field[2:])
#         continue
#     stokesv = hdulist[0].data
#     # stokesv = hdulist[0].data[0,1,:,:]
#     header = hdulist[0].header

#     centre_ra_store[i] = header['CRVAL1']
#     centre_dec_store[i] = header['CRVAL2']
#     centre_ra_pix = header['CRPIX1']
#     centre_dec_pix = header['CRPIX2']
#     num_pix_1 = header['NAXIS1']
#     num_pix_2 = header['NAXIS2']
    
#     if num_pix_1 != num_pix_2:
#         continue
#     # hdulist.close()

#     rms_meas_2deg[i] = np.median(stokesv)
#     #    rms_meas_2deg[i] = np.median(stokesv[int(centre_ra_pix)-(num_pix_1//2):int(centre_ra_pix)+(num_pix_1//2),\
#                                        #  int(centre_dec_pix)-(num_pix_2//2):int(centre_dec_pix)+(num_pix_2//2)]) #* 1e6 #in muJy
#     print(rms_meas_2deg[i])
#     i = i+1

# cols = fits.ColDefs([
# fits.Column(name='Field', format='11A',
# array=fields),
# fits.Column(name='ra_centre_field', format='D',
# array=centre_ra_store),
# fits.Column(name='dec_centre_field', format='D',
# array=centre_dec_store),
# fits.Column(name='median_rms', format='D',
# array=rms_meas_2deg)
#     ])

# hdu_save = fits.BinTableHDU.from_columns(cols)
# hdu_save.writeto('VLASS_median_rms.fits',overwrite=True)

#read in the file again
hdulist = fits.open('VLASS_median_rms.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

median_rms = tbdata['median_rms'] * 1e6 #muJy
#print('median', np.median(median_rms[median_rms>0.]))

#make the plot
fig = plt.figure(1,figsize=(12, 8))
gs = plt.GridSpec(1,1)
ax = plt.subplot(gs[0])
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
ax.tick_params(axis='both',which='major',length=8,width=1.5)
ax.hist(median_rms[median_rms>0.],bins=10000,histtype='step',lw=2,color='k')
ax.axvline(np.median(median_rms[median_rms>0.]),color='crimson',lw=2,ls='--', \
           label='median rms at {:.4} muJy'.format(np.median(median_rms[median_rms>0.])))
ax.set_ylabel('Number', fontsize = 25)
ax.set_xlabel(r'Median RMS Noise ($\mu$Jy) ', fontsize = 25)
ax.tick_params(axis='both',which='both',labelsize=17)
ax.set_xlim([90, 350])


#Here we find the missing fields, start by isolating the indexes
indexes = tbdata['Field']
indexes_1 = indexes.split('_')
indexes_2 = []
for i, j in enumerate(indexes_1):
    inx = (indexes_1[i])[0][2:]
    indexes_2.append(int(inx))
indexes_2_arr= np.array(indexes_2) #final array with the present indexes, contain duplicates
max_index = max(indexes_2_arr)
mask_arr = np.zeros(max_index+1) #empty array

for i in indexes_2_arr: #run over all present indexes
    mask_arr[i]=1
    
missing_fields = len(mask_arr[mask_arr==0]) / len(mask_arr)* 100
print('We find {:.3} % missing fields'.format(missing_fields))

#Here we compute the bad fields
#A field is bad when the faintest source in NVSS (0.0021Jy) is undetected in VLASS 
#For a spectral index of -0.8
#This would be any VLASS source with
alpha_mean=-1.1
VLASS_limit = 0.0021 * ((3000/1400)**(alpha_mean))  
#TODO: get a more robust number here. Cite something 
#A source is detected at 5 sigma, so a field is bad when the rms is higher than
VLASS_rms_limit = (VLASS_limit / 5) *1e6 #muJy
print('VLASS rms limit', VLASS_rms_limit)
ax.axvline(VLASS_rms_limit,color='blue',lw=2,ls='--', label='Bad fields limit at {:.4} muJy'.format(VLASS_rms_limit))
ax.legend()
plt.savefig('rms_analysis_histo.png',bbox_inches='tight')
#The amount of bad fields
num_too_faint = len(median_rms[median_rms > VLASS_rms_limit])


print('We find {:.3} % Bad fields'.format(num_too_faint/len(median_rms[median_rms>0])*100))


