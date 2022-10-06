# -*- coding: utf-8 -*-
"""
Written by Femke Ballieux
takes in master sample, and produces some plots, so for individual sources as well
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#read in the data
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/fit_vals_power_law_NVSS_intflux.fits')
high_survey = 'NVSS' #could also be first, but we use NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#print the columnnames
#print(orig_cols)

#make lists of the most used columns
alpha_low = tbdata['alpha_low']
alpha_high = tbdata['alpha_high']
name_list = tbdata['LoTSS_name']

flux_LoTSS = tbdata['LoTSS_flux']
flux_LoLSS = tbdata['LoLSS_flux']
flux_NVSS = tbdata['NVSS_flux']
flux_TGSS = tbdata['TGSS_flux']
flux_VLSSr = tbdata['VLSSr_flux']
flux_FIRST = tbdata['FIRST_flux']
flux_inband_low = tbdata['S_inband_low']
flux_inband_mid = tbdata['S_inband_mid']
flux_inband_high = tbdata['S_inband_high']



#list of frequencies in MHz
freq_list =     144., 1400., 150., 74., 54., 1400.001, 128., 144.00001, 160.
label_list=['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high']    

def find_index(galaxy_name, name_list = name_list):
    """
    Gives the index of a galaxy when the name is entered
    """
    return int(np.where(name_list == galaxy_name )[0])


def make_sed_singular(galaxy_name, alpha_low=alpha_low):
    """
    returns the sed for a singular galaxy
    """
    index = find_index(galaxy_name)
    flux_list = [flux_LoTSS[index],flux_NVSS[index],flux_TGSS[index], flux_VLSSr[index],flux_LoLSS[index],   \
                 flux_FIRST[index], flux_inband_low[index], flux_inband_mid[index], flux_inband_high[index] ]
    
    plt.figure(figsize=(10,8))
    for s in range(len(freq_list)):
        plt.scatter(freq_list[s], flux_list[s], label=label_list[s])
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

    


name_of_interest = '104702+474945' 
make_sed_singular(name_of_interest)
