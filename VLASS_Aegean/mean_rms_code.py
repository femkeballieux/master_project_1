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

#This is a file containing the filenames of the bane rms files
fields = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/rms_filenames.txt', dtype='U')

#In this code we compute the median of all these files and write them to a file, along with the original LoTSS index

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
#print(orig_cols)
hdulist.close()
median_rms = tbdata['median_rms'] * 1e6 #muJy


#make the plot of the distribution of the median rms per field, I did not filter on doubles here
fig = plt.figure(1,figsize=(12, 8))
gs = plt.GridSpec(1,1)
ax = plt.subplot(gs[0])
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
ax.tick_params(axis='both',which='major',length=8,width=1.5)
ax.hist(median_rms[median_rms>0.],bins=10000,histtype='step',lw=2,color='k')
ax.axvline(np.median(median_rms[median_rms>0.]),color='crimson',lw=2,ls='--', \
           label='Median of rms = {:.4} muJy'.format(np.median(median_rms[median_rms>0.])))
ax.set_ylabel('Number of sources', fontsize = 25)
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
indexes_2_arr= np.array(indexes_2) #final array with the present indexes, contain duplicates, unordered
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
alpha_mean=-0.8
VLASS_limit = 0.0021 * ((3000/1400)**(alpha_mean))  
#TODO: get a more robust number here. Cite something 
#A source is detected at 5 sigma, so a field is bad when the rms is higher than
VLASS_rms_limit = (VLASS_limit / 5) *1e6 #muJy
print('VLASS rms limit', VLASS_rms_limit)
ax.axvline(VLASS_rms_limit,color='blue',lw=2,ls='--', label='High-noise limit = {:.4} muJy'.format(VLASS_rms_limit))
ax.legend(fontsize=20)
plt.savefig('rms_analysis_histo.pdf',bbox_inches='tight')
#The amount of bad fields
num_too_faint = len(median_rms[median_rms > VLASS_rms_limit])


print('We find {:.3} % Bad fields'.format(num_too_faint/len(median_rms[median_rms>0])*100))


# hdulist2 = fits.open(\
#    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/Official_VLASS_no_dups/official_mega_master_clean_6.fits')

# tbdata2 = hdulist2[1].data
# orig_cols2 = hdulist2[1].columns
# # print(orig_cols2)
# hdulist2.close()


# LoTSS_int=tbdata2['LoTSS_flux']
# LoTSS_int_err = tbdata2['e_LoTSS_flux']
# LoTSS_peak=tbdata2['LoTSS_flux_peak']
# e_LoTSS_peak=tbdata2['stat_e_LoTSS_flux_peak']
# # rms=tbdata['LoTSS_rms']


# # ordered_median_rms = np.zeros(len(LoTSS_int))
# # for i, rms in enumerate(indexes):
# #     #i runs over the indexes of the unordered median_rms array. Clean_array runs over the ordered LoTSS array
# #     dirty_index=(indexes[i]).split('_')[0]
# #     clean_index=int(dirty_index[2:])
# #     ordered_median_rms[clean_index] = median_rms[i] 
 

# ordered_median_rms = np.zeros(len(LoTSS_int))
# for i, LoTSS_index in enumerate(indexes_2):
#     #We have already computed before that indexes_2 is the array of indexes 
#     ordered_median_rms[LoTSS_index] = median_rms[i] 
 

# # mask = np.where((ordered_median_rms<=500.))
# mask = np.where((tbdata2['VLASS_flux']==0.)&(ordered_median_rms<=228.3)) 
# #Since we only need those with missing VLASS flux, and we dont want to include those we can already explain by bad fields

# #There are a bunch of negative median-rms values. Check with Joe if this is weird
# R=np.log(LoTSS_int/LoTSS_peak)
# SNR=(LoTSS_int/LoTSS_int_err)
# plt.figure(figsize=(10,8))
# plt.semilogx()

# #plt.scatter(SNR,R, alpha=1, s=1, label='all', color='crimson')
# plt.scatter(SNR[mask],R[mask],c=ordered_median_rms[mask], alpha=1, s=0.8, label='missing VLASS')
# plt.colorbar(label='VLASS rms in mJy')

# x_array=np.linspace(1,10000,1000)

# def R_99(SNR):
#     return 0.42 + (1.08 / (1+(SNR/96.57)**2.49))
# plt.plot(x_array, R_99(x_array), color='black', label='Shimwell')
# plt.legend()
# plt.xlim(1, 10000)
# plt.xlabel('SNR')
# plt.ylabel('$\ln(S_I/S_P)$')


# hdulist3 = fits.open(\
#     '/net/vdesk/data2/bach1/ballieux/master_project_1/data/Official_VLASS_no_dups/official_mega_master_intermediate_crossmatchtest_6.fits')

# tbdata3 = hdulist3[1].data
# orig_cols3 = hdulist3[1].columns
# print(orig_cols3)
# hdulist3.close()


# LoTSS_maj=tbdata3['Maj_1']

# plt.figure(figsize=(10,8))
# # plt.scatter(SNR[mask],LoTSS_maj[mask],c=ordered_median_rms[mask], alpha=1, s=0.8, label='missing VLASS')
# plt.scatter(SNR,LoTSS_maj,c=ordered_median_rms, alpha=1, s=0.8, label='missing VLASS')
# plt.colorbar(label='VLASS rms in mJy')
# plt.semilogx()
# plt.semilogy()
# plt.xlabel('SNR')
# plt.ylabel('LoTSS major axis')
# plt.hlines(6,0.5,5000, label='6 arcseconds')



