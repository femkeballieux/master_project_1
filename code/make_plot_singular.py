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
#hdulist = fits.open(path_vdesk + 'PS_with_vlass.fits')
hdulist = fits.open(path_vdesk + 'mega_master_clean.fits')
high_survey = 'NVSS' #could also be first, but we use NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#print the columnnames
print(orig_cols)

#make lists of the most used columns
name_list = tbdata['LoTSS_name']

alpha_low_LoLSS = tbdata['alpha_low_LoLSS']
alpha_high_LoLSS = tbdata['alpha_high_LoLSS']
a_low_LoLSS = tbdata['a_low_LoLSS']
a_high_LoLSS = tbdata['a_high_LoLSS']
alpha_low_VLASS = tbdata['alpha_low_VLASS']
alpha_high_VLASS = tbdata['alpha_high_VLASS']
a_low_VLASS = tbdata['a_low_VLASS']
a_high_VLASS = tbdata['a_high_VLASS']
flux_LoTSS = tbdata['LoTSS_flux']
flux_LoLSS = tbdata['LoLSS_flux']
flux_NVSS = tbdata['NVSS_flux']
flux_TGSS = tbdata['TGSS_flux']
flux_VLSSr = tbdata['VLSSr_flux']
flux_FIRST = tbdata['FIRST_flux']
flux_inband_low = tbdata['S_inband_low']
flux_inband_mid = tbdata['S_inband_mid']
flux_inband_high = tbdata['S_inband_high']
flux_channel0 = tbdata['channel0_flux']
flux_channel1 = tbdata['channel1_flux']
flux_channel2 = tbdata['channel2_flux']
flux_channel3 = tbdata['channel3_flux']
flux_channel4 = tbdata['channel4_flux']
flux_channel5 = tbdata['channel5_flux']

error_LoTSS = tbdata['e_LoTSS_flux']
error_LoLSS = tbdata['e_LoLSS_flux']
error_NVSS = tbdata['e_NVSS_flux']
error_TGSS = tbdata['e_TGSS_flux']
error_VLSSr = tbdata['e_VLSSr_flux']
error_FIRST = tbdata['e_FIRST_flux']
error_inband_low = tbdata['e_S_inband_low']
error_inband_mid = tbdata['e_S_inband_mid']
error_inband_high = tbdata['e_S_inband_high']
error_channel0 = tbdata['e_channel0_flux']
error_channel1 = tbdata['e_channel1_flux']
error_channel2 = tbdata['e_channel2_flux']
error_channel3 = tbdata['e_channel3_flux']
error_channel4 = tbdata['e_channel4_flux']
error_channel5 = tbdata['e_channel5_flux']

flux_vlass = tbdata['VLASS_flux']
e_flux_vlass = tbdata['e_VLASS_flux']


#frequencies all in Mhz
freq_LoTSS= 144.
freq_NVSS = 1400.
freq_TGSS = 150.
freq_VLSSr = 74.
freq_LoLSS = 54.
freq_FIRST = 1400.001
freq_inband_low = 128. 
freq_inband_mid =  144.00001
freq_inband_high = 160.
freq_LoLLS_ch0 = 44.
freq_LoLLS_ch1 = 48.
freq_LoLLS_ch2 = 52.
freq_LoLLS_ch3 = 56.
freq_LoLLS_ch4 = 60.
freq_LoLLS_ch5 = 64.
freq_vlass = 3000.

#list of frequencies in MHz
freq_array = np.array([freq_LoTSS, freq_NVSS, freq_TGSS, freq_VLSSr, freq_LoLSS, freq_FIRST , freq_inband_low,\
            freq_inband_mid, freq_inband_high,freq_LoLLS_ch0, freq_LoLLS_ch1, freq_LoLLS_ch2 , freq_LoLLS_ch3 \
                , freq_LoLLS_ch4 , freq_LoLLS_ch5, freq_vlass ])


#labels in order they are used
label_array = np.array(['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high', \
              'LoLLS_ch0', 'LoLLS_ch1', 'LoLLS_ch2', 'LoLLS_ch3', 'LoLLS_ch4', 'LoLLS_ch5', 'VLASS' ] )   

#used for plotting the spectral indices
#TODO: add wenss
x_range_low = np.linspace(50, 144, 1000)
x_range_mid = np.linspace(144, 1400, 1000)
x_range_high = np.linspace(1400, 3000, 1000)
x_range_full= np.linspace(50, 3000, 1000)

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

    did_it_fit = True
    try:
        poptgen, pcovgen = opt.curve_fit(curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        print('Curve_fit could not fit curve model.')
        did_it_fit = False
        return (fluxplot,freqplot,flux_errplot,0,0, did_it_fit)
        
        print('fitting went right')
    return (fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit)
    


def make_sed_singular(galaxy_name, use_index=False, save_fig=False):
    """
    returns the sed for a singular galaxy
    """
    index = find_index(galaxy_name)
        
    #for this particular galaxy the fluxes    
    flux_array = np.array([flux_LoTSS[index], flux_NVSS[index], flux_TGSS[index], flux_VLSSr[index], flux_LoLSS[index],\
                 flux_FIRST[index], flux_inband_low[index], flux_inband_mid[index], flux_inband_high[index], flux_channel0[index],\
                     flux_channel1[index], flux_channel2[index], flux_channel3[index], flux_channel4[index], flux_channel5[index],\
                         flux_vlass[index]])
    error_array = np.array([error_LoTSS[index], error_NVSS[index], error_TGSS[index], error_VLSSr[index], error_LoLSS[index],\
                 error_FIRST[index], error_inband_low[index], error_inband_mid[index], error_inband_high[index], \
                     error_channel0[index], error_channel1[index], error_channel2[index], error_channel3[index],\
                         error_channel4[index], error_channel5[index], e_flux_vlass[index]])
    
    #general plotting settings
    plt.figure(figsize=(10,8))
    plt.xlabel('frequency [MHz]', fontsize=14)
    plt.ylabel('flux [Jy]',fontsize=14)
    plt.title(galaxy_name, fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    
    #Plot the data
    inband_lolss = 0 #fixing the legends
    inband_lotss = 0
    for s in range(len(freq_array)):
        
        if flux_array[s] != 0.:
            
            if 'LoLLS_ch' in str(label_array[s]):
                if inband_lolss==0:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                             label='LoLSS inband', zorder=10, fmt='o', color='m', alpha=0.4)
                    inband_lolss=1   
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                         zorder=10, fmt='o', color='m', alpha=0.4)
                        
            elif 'inband_' in str(label_array[s]):
                if inband_lotss == 0:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                         label='LoTSS inband', zorder=10, fmt='o', color='darkgreen', alpha=0.35)
                    inband_lotss = 1
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                         zorder=10, fmt='o', color='darkgreen', alpha=0.35)
                    
            elif label_array[s]=='VLASS':
                if error_array[s]==-1.:
                    plt.plot(freq_array[s], flux_array[s],\
                           label=label_array[s], zorder=10, marker='v')  
                    print('UPPERLIMIT')
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                     label=label_array[s], zorder=10, fmt='o')
            else:    
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                     label=label_array[s], zorder=10, fmt='o')
    
    #plot the power law fits for LoLSS
    if alpha_low_LoLSS[index]!=0.:
        plt.plot(x_range_low, a_low_LoLSS[index] * (x_range_low ** alpha_low_LoLSS[index]),\
                 color='green', label='powerlaw low LoLSS') 
        plt.plot(x_range_mid, a_high_LoLSS[index] * (x_range_mid ** alpha_high_LoLSS[index]),\
                 color='black', label='powerlaw high LoLSS') 
    #Plot the powerlaw fits for VLASS        
    if alpha_low_VLASS[index]!=0.:
        
        plt.plot(x_range_mid, a_low_VLASS[index] * (x_range_mid ** alpha_low_VLASS[index]),\
                 color='indigo', label='powerlaw low VLASS') 
        plt.plot(x_range_high, a_high_VLASS[index] * (x_range_high ** alpha_high_VLASS[index]),\
                 color='blue', label='powerlaw high VLASS') 
    else:
        print('no spectral index fit')
    
    #fit the curved models
    
    fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit = \
                spectral_index_eval_curve(freq_array, flux_array, error_array)
    print('parameters for the curve are', poptgen )
    # if did_it_fit:
    #     plt.plot(x_range_full, curve(x_range_full, poptgen[0], poptgen[1] , poptgen[2], poptgen[3] ), color='red', label='curve')
    
    plt.legend()
    if save_fig:
        plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/'+ galaxy_name + '.pdf')
        print('figure saved for', galaxy_name)
    plt.show()


make_sed_singular(tbdata['LoTSS_name'][46178], save_fig=False)
# counter=0
# for i, name in enumerate(name_list):
#     #if (alpha_low[i] >= 0.1) & (alpha_high[i]<=0): #select when a source is PS
#     make_sed_singular(name, save_fig=True)
#     counter +=1
#     print(counter, '/767')

        
