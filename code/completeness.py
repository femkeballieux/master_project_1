"""
Written by Femke Ballieux, used to analyse completeness properties of our sample
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import iqr

# plt.style.use('style.mplstyle')

# Load in the table
hdulist = fits.open(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/official_mega_master_clean.fits')
sample = 'VLASS'  # This van be either VLASS for the sample LoTSS NVSS VLASS, or LoLSS for LoLSS LoTSS NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
#print(orig_cols)
hdulist.close()

#Define the samples
flux_sample = tbdata[sample + '_flux']
flux_NVSS= tbdata['NVSS_flux']
flux_LoTSS= tbdata['LoTSS_flux']

alpha_low=tbdata['alpha_low_'+sample]
alpha_high=tbdata['alpha_high_'+sample]
alpha_low_e=tbdata['e_alpha_low_'+sample]
alpha_high_e=tbdata['e_alpha_high_'+sample]

#We want to exclude all those without a detection in VLASS/LoLSS. This is a base mask that always counts.
mask = np.where((flux_sample!=0.) & (alpha_high_e != 0.) & (alpha_high !=0.)\
                & (alpha_low_e != 0.) & (alpha_low !=0.))

print('with a VLASS detection',len(alpha_low[mask]))


NVSS_cut = 0.1
above_NVSS = len(alpha_low[flux_NVSS>=NVSS_cut])
print("")
print("NVSS cut", NVSS_cut * 1000, 'mJy')
print('Above NVSS cut', above_NVSS )

completeness_mask =  np.where((flux_sample!=0.) & (alpha_high_e != 0.) & (alpha_high !=0.)\
                & (alpha_low_e != 0.) & (alpha_low !=0.)\
                    & (flux_NVSS>=NVSS_cut))
with_VLASS_above_NVSS = len(alpha_low[completeness_mask])

print('With VLASS and above NVSS cut', with_VLASS_above_NVSS)

complete_PS_mask =  np.where((flux_sample!=0.) & (alpha_high_e != 0.) & (alpha_high !=0.)\
                & (alpha_low_e != 0.) & (alpha_low !=0.)\
                    & (flux_NVSS>=NVSS_cut) & (alpha_low>=np.median(alpha_low_e))& (alpha_high<=-np.median(alpha_high_e)))
    
with_VLASS_above_NVSS_PS = len(alpha_low[complete_PS_mask])

print('PS with VLASS and above NVSS cut', with_VLASS_above_NVSS_PS)

print('completeness', 100 * with_VLASS_above_NVSS / above_NVSS, '%')
print('percentage PS', 100 * with_VLASS_above_NVSS_PS / with_VLASS_above_NVSS , '%')