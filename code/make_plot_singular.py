# -*- coding: utf-8 -*-
"""
Written by Femke Ballieux
takes in master sample, and produces some plots, so for individual sources as well
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.optimize as opt
import gpscssmodels


path_laptop = 'C:/Users/Femke/Documents/GitHub/master_project_1/data'
path_vdesk= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'

#read in the data
hdulist = fits.open(path_laptop + '/fit_vals_power_law_NVSS_intflux.fits')
high_survey = 'NVSS' #could also be first, but we use NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#print the columnnames
#print(orig_cols)

#make lists of the most used columns
name_list = tbdata['LoTSS_name']

alpha_low = tbdata['alpha_low']
alpha_high = tbdata['alpha_high']
a_low = tbdata['a_low']
a_high = tbdata['a_high']

flux_LoTSS = tbdata['LoTSS_flux']
flux_LoLSS = tbdata['LoLSS_flux']
flux_NVSS = tbdata['NVSS_flux']
flux_TGSS = tbdata['TGSS_flux']
flux_VLSSr = tbdata['VLSSr_flux']
flux_FIRST = tbdata['FIRST_flux']
flux_inband_low = tbdata['S_inband_low']
flux_inband_mid = tbdata['S_inband_mid']
flux_inband_high = tbdata['S_inband_high']

error_LoTSS = tbdata['e_LoTSS_flux']
error_LoLSS = tbdata['e_LoLSS_flux']
error_NVSS = tbdata['e_NVSS_flux']
error_TGSS = tbdata['e_TGSS_flux']
error_VLSSr = tbdata['e_VLSSr_flux']
error_FIRST = tbdata['e_FIRST_flux']
error_inband_low = tbdata['e_S_inband_low']
error_inband_mid = tbdata['e_S_inband_mid']
error_inband_high = tbdata['e_S_inband_high']


#list of frequencies in MHz
freq_list = 144., 1400., 150., 74., 54., 1400.001, 128., 144.00001, 160.
freq_array = np.array(freq_list)

#labels in order they are used
label_list = ['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high']    

#used for plotting
x_range_low = np.linspace(50, 144, 1000)
x_range_high = np.linspace(144, 1400, 1000)
x_range_full= np.linspace(50, 1400, 1000)

def find_index(galaxy_name, name_list = name_list):
    """
    Gives the index of a galaxy when the name is entered
    """
    return int(np.where(name_list == galaxy_name )[0])

def curve(freq, freq_peak, flux_peak, alphathick, alphathin): 
    """
    Model taken from Tschager et al. 2003. General fit not based on any physics.
    """
    return flux_peak/(1 -np.exp(-1))*((freq/freq_peak)**alphathick)*(1-np.exp(-(freq/freq_peak)**(alphathin-alphathick)))


def spectral_index_eval_curve(freq, flux, flux_err):
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data
    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'
    
    fluxplot = flux[(flux > 0)]
    freqplot = freq[flux > 0]
    flux_errplot = flux_err[(flux > 0)]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0gen = [peakfreq[0], peakflux, 0.8, -0.7]

    try:
        poptgen, pcovgen = opt.curve_fit(curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    return (fluxplot, freqplot, flux_errplot, poptgen, pcovgen)
    


def make_sed_singular(galaxy_name, use_index=False):
    """
    returns the sed for a singular galaxy
    """
    #extract the data
    if use_index:
        index = galaxy_name
        #TODO: make sure the name goes right when you do this
    else:
        index = find_index(galaxy_name)
        
    #for this particular galaxy the fluxes    
    flux_list = [flux_LoTSS[index],flux_NVSS[index],flux_TGSS[index], flux_VLSSr[index],flux_LoLSS[index],   \
                 flux_FIRST[index], flux_inband_low[index], flux_inband_mid[index], flux_inband_high[index] ]
    error_list = [error_LoTSS[index],error_NVSS[index],error_TGSS[index], error_VLSSr[index],error_LoLSS[index],   \
                 error_FIRST[index], error_inband_low[index], error_inband_mid[index], error_inband_high[index] ]
    flux_array = np.array(flux_list)
    error_array = np.array(error_list)
    
    #general plotting settings
    plt.figure(figsize=(10,8))
    plt.xlabel('frequency [MHz]', fontsize=14)
    plt.ylabel('flux [Jy]',fontsize=14)
    plt.title(galaxy_name, fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    
    #Plot the data
    for s in range(len(freq_list)):
        if flux_list[s] != 0:
            plt.errorbar(freq_list[s], flux_list[s], yerr=error_list[s],\
                         label=label_list[s], zorder=10, fmt='o')
    
    #plot the power law fits  
    plt.plot(x_range_low, a_low[index] * (x_range_low ** alpha_low[index]),\
             color='green', label='powerlaw low') 
    plt.plot(x_range_high, a_high[index] * (x_range_high ** alpha_high[index]),\
             color='black', label='powerlaw high') 
    
    #fit the curved models
    fluxplot, freqplot, flux_errplot, poptgen, pcovgen = \
                spectral_index_eval_curve(freq_array, flux_array, error_array)
    print('parameters for the curve are', poptgen )
    plt.plot(x_range_full, curve(x_range_full, poptgen[0],poptgen[1] , poptgen[2], poptgen[3] ), color='red', label='curve')
    
    plt.legend()
    plt.show()

#104732+472532 no longer PS, was in Martjes paper
#cannot find '110038+523620 from Martjes paper'
#'104713+470331' example PS
#'104657+482724' example non-PS
name_of_interest = '104713+470331'
make_sed_singular(name_of_interest, use_index=False)
