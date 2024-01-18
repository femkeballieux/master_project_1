"""
Code for studying and adapting the uncertainty in alpha_low for the 2 samples
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import iqr
import scipy.stats as st

#Open the clean master sample
hdulist = fits.open(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_clean.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#Arrays we need
LoLSS_flux = tbdata['LoLSS_flux']
LoLSS_err = tbdata['e_LoLSS_flux']
LoLSS_rms = tbdata['LoLSS_rms']
LoTSS_flux = tbdata['LoTSS_flux']
LoTSS_err = tbdata['e_LoTSS_flux']
LoTSS_rms = tbdata['LoTSS_rms']
VLASS_flux = tbdata['VLASS_flux']
VLASS_err = tbdata['e_VLASS_flux']


#finding the non-zero sources for the LF sample
source_bright_ind_LF = np.where((LoLSS_flux > 0.))
alpha_low_LF = tbdata['alpha_low_LoLSS'][source_bright_ind_LF]
err_alpha_low_LF = tbdata['e_alpha_low_LoLSS'][source_bright_ind_LF]
LoLSS_flux_LF = LoLSS_flux[source_bright_ind_LF]
LoLSS_err_LF = LoLSS_err[source_bright_ind_LF]
LoLSS_rms_LF = LoLSS_rms[source_bright_ind_LF]

#Study how the LoLSS SNR is divided
SNR_LF_LoLSS = LoLSS_flux_LF / LoLSS_rms_LF

#The HF sample
source_bright_ind_HF = np.where((VLASS_flux > 0.))
LoTSS_flux_HF = LoTSS_flux[source_bright_ind_HF]
LoTSS_err_HF = LoTSS_err[source_bright_ind_HF]
LoTSS_rms_HF = LoTSS_rms[source_bright_ind_HF]

alpha_low_HF = tbdata['alpha_low_VLASS'][source_bright_ind_HF]
err_alpha_low_HF = tbdata['e_alpha_low_VLASS'][source_bright_ind_HF]

SNR_HF_LoTSS = LoTSS_flux_HF / LoTSS_rms_HF

#Make a new array with the lentgh of our HF sample, to store the new artificial uncertainties in the HF alpha_low
new_err_alpha_low_HF = np.zeros(len(LoTSS_flux_HF))


#break up the LF sample into LoLSS SNR bins
percentages = [0.,20., 40., 60., 80., 100.] #The boundaries for the SNR bins
for i, percentual in enumerate(percentages[1:]):
    # print(np.percentile(SNR_LF_LoLSS, percentages[i]))
    SNR_boundary_LF_low = np.percentile(SNR_LF_LoLSS, percentages[i])
    SNR_boundary_LF_high = np.percentile(SNR_LF_LoLSS, percentages[i+1])

    # #The mask for the SNR bins in the LF sample
    mask_SNRbin_LF = np.where((SNR_LF_LoLSS >= SNR_boundary_LF_low) & (SNR_LF_LoLSS<SNR_boundary_LF_high))
    print('The amount of LF sources in SNR bin', SNR_boundary_LF_low, '-', SNR_boundary_LF_high,':', len(SNR_LF_LoLSS[mask_SNRbin_LF]))

    # #The mask for the SNR bins in the HF sample
    # print(np.percentile(SNR_HF_LoTSS, percentages[i]))
    SNR_boundary_HF_low = np.percentile(SNR_HF_LoTSS, percentages[i])
    SNR_boundary_HF_high = np.percentile(SNR_HF_LoTSS, percentages[i+1])
    mask_SNRbin_HF = np.where((SNR_HF_LoTSS >= SNR_boundary_HF_low) & (SNR_HF_LoTSS<SNR_boundary_HF_high))
    print('The amount of HF sources in SNR bin',SNR_boundary_HF_low, '-',SNR_boundary_HF_high,':',len(SNR_HF_LoTSS[mask_SNRbin_HF]))

    # #The amount of HF sources in each SNR bin
    length_HF_in_bin = len(SNR_HF_LoTSS[mask_SNRbin_HF])

    # #This is where the magic happens, a new value for each HF datapoint in the bin
    new_err_alpha_low_HF[mask_SNRbin_HF] = np.random.choice(err_alpha_low_LF[mask_SNRbin_LF], size=length_HF_in_bin)

    # print(np.median(err_alpha_low_LF[mask_SNRbin_LF]))
    # print(np.median(np.random.choice(err_alpha_low_LF[mask_SNRbin_LF], size=length_HF_in_bin)))

print(np.median(new_err_alpha_low_HF), 'new HF ')
print(np.median(err_alpha_low_LF), 'LF')
print(np.median(err_alpha_low_HF), 'old HF')

"""
Without the SNR bins, it is just a random distribution
"""

col = fits.Column(name='new_e_alpha_low_HF', format = 'E', array = new_err_alpha_low_HF)

cols = fits.ColDefs([col])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  


tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/new_e_alpha_low_HF.fits', overwrite = True)

# def plot_e_alpha(e_alpha_low, name, title, xmin=0):
#     plt.figure()
#     plt.hist(e_alpha_low, bins=50)
#     plt.title(title)
#     plt.xlabel(r'Uncertainty in $ \alpha_{low}$')
#     plt.ylabel('Number of sources')
#     plt.xlim(xmin=xmin, xmax = np.max(e_alpha_low)+0.01)
#     plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/'+ name)

# plot_e_alpha(new_err_alpha_low_HF,'new_e_alpha_HF', 'HF sample new distribution', xmin=0.09)
# plot_e_alpha(err_alpha_low_HF, 'old_e_alpha_HF', 'HF sample old distribution')
# plot_e_alpha(err_alpha_low_LF, 'e_alpha_LF',  'LF sample',  xmin=0.09)


# import matplotlib.pyplot as plt
# import numpy as np

def plot_e_alpha(ax, e_alpha_low, title):
    ax.hist(e_alpha_low, bins=50)
    ax.set_title(title)
    ax.set_xlabel(r'Uncertainty in $ \alpha_{low}$')
    ax.set_ylabel('Number of sources')
    ax.set_xlim(xmin=0.085, xmax=0.47)

# Create a single figure and three subplots arranged vertically
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 12))


# Call the function for each subplot
plot_e_alpha(ax1, err_alpha_low_LF, 'LF sample')
plot_e_alpha(ax3, new_err_alpha_low_HF, 'HF sample new distribution')
plot_e_alpha(ax2, err_alpha_low_HF, 'HF sample old distribution')


# Adjust layout for better spacing
plt.tight_layout()

# Save the final plot
plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/alpha_distributions.pdf', dpi=100)
