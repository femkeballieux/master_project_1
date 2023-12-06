# -*- coding: utf-8 -*-
"""
Written by Femke Ballieux
takes in master sample, and produces some plots, so for individual sources as well
I sincerely apologize for how ugly this code is... it works though.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.optimize as opt
import gpscssmodels
from matplotlib.ticker import ScalarFormatter


#read in the data
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_clean.fits') #for checking with whole master sample
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

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
e_flux_vlass = 0.1 * tbdata['VLASS_flux']


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
freq_LoLSS_ch0 = 44.
freq_LoLSS_ch1 = 48.
freq_LoLSS_ch2 = 52.
freq_LoLSS_ch3 = 56.
freq_LoLSS_ch4 = 60.
freq_LoLSS_ch5 = 64.
freq_vlass = 3000.
 
#fixing the legends
inband_lolss = 0 
inband_lotss = 0
lotss = 0 
lolss = 0
nvss = 0
tgss = 0
first = 0
vlass = 0
vla = 0
green = 0

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
    


def make_sed_singular(galaxy_name, use_index=False, save_fig=False, ax=None):
    """
    returns the sed for a singular galaxy
    """
    global lotss, inband_lolss, inband_lotss, lolss, nvss, tgss, first, vlass, vla, green
    index = find_index(galaxy_name)
    x_range_low = np.linspace(50, 144, 1000)
    x_range_mid = np.linspace(144, 1400, 1000)
    x_range_high = np.linspace(1400, 3000, 1000)
    x_range_full= np.linspace(50, 3000, 1000)
    
    #for this particular galaxy the fluxes    
    flux_array = np.array([flux_LoTSS[index], flux_NVSS[index], flux_TGSS[index], flux_VLSSr[index], flux_LoLSS[index],\
                 flux_FIRST[index], flux_inband_low[index], flux_inband_mid[index], flux_inband_high[index], flux_channel0[index],\
                     flux_channel1[index], flux_channel2[index], flux_channel3[index], flux_channel4[index], flux_channel5[index],\
                         flux_vlass[index]])

        
    error_array = np.array([error_LoTSS[index], error_NVSS[index], error_TGSS[index], error_VLSSr[index], error_LoLSS[index],\
                 error_FIRST[index], error_inband_low[index], error_inband_mid[index], error_inband_high[index], error_channel0[index],\
                     error_channel1[index], error_channel2[index], error_channel3[index], error_channel4[index], error_channel5[index],\
                         e_flux_vlass[index]])
    #list of frequencies in MHz
    freq_array = np.array([freq_LoTSS, freq_NVSS, freq_TGSS, freq_VLSSr, freq_LoLSS, freq_FIRST , freq_inband_low,\
                freq_inband_mid, freq_inband_high,freq_LoLSS_ch0, freq_LoLSS_ch1, freq_LoLSS_ch2 , freq_LoLSS_ch3 \
                    , freq_LoLSS_ch4 , freq_LoLSS_ch5, freq_vlass])
    
    
    #labels in order they are used
    label_array = np.array(['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high', \
                  'LoLSS_ch0', 'LoLSS_ch1', 'LoLSS_ch2', 'LoLSS_ch3', 'LoLSS_ch4', 'LoLSS_ch5', 'VLASS' ] )  


    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='both', labelsize=15, size=5)
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title(galaxy_name, fontsize=20)
    
    """
    Comment back in when plotting these specific galaxies
    """
    if galaxy_name == tbdata['LoTSS_name'][39286]:
        ax.set_title('87GB B0903+6656', fontsize=20)
        ax.set_xlabel('Frequency [GHz]', fontsize=16)
        ax.set_ylabel('$S_{total}$ [mJy]',fontsize=16)

        flux_array = np.array([flux_LoTSS[index]*1000, flux_NVSS[index]*1000,  74.9, 195 ,flux_vlass[index]*1000]) #mJy
        error_array = np.array([error_LoTSS[index]*1000, error_NVSS[index]*1000, 11.235, 29.25, e_flux_vlass[index]*1000])
        freq_array = np.array([freq_LoTSS/1000, freq_NVSS/1000, 8.400 , 4.850 , freq_vlass/1000]) #GHz
        
        label_array = np.array(['LoTSS', 'NVSS',  'VLA', 'NRAO Green Bank' ,'VLASS' ] ) 
        x_range_full= np.linspace(0.050, 9.000, 1000)
        x_range_low = np.linspace(0.050, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)
        
    elif galaxy_name in [tbdata['LoTSS_name'][103164], tbdata['LoTSS_name'][98418] , tbdata['LoTSS_name'][83458]]:
        x_range_full= np.linspace(0.050, 9.000, 1000)
        x_range_low = np.linspace(0.050, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)
        
        ax.set_xlabel('Frequency [GHz]', fontsize=16)
        ax.set_ylabel('$S_{total}$ [mJy]',fontsize=16)
        
        flux_array = flux_array* 1000
        error_array = error_array * 1000
        freq_array = freq_array / 1000




    markers=['v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4']
    colors=['darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'deeppink', \
            'mediumorchid', 'darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'blue', 'mediumorchid','darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'green', \
                    'mediumorchid', 'darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'green', 'mediumorchid',]
    
    #This is a very ugly way of getting the legend right
    legend_handles = []
    legend_labels = []
    for s in range(len(freq_array)):
        if flux_array[s] != 0.:
            
            if 'LoLSS_ch' in str(label_array[s]):
                if inband_lolss==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='LoLSS inband', zorder=1, fmt='o', color='grey', alpha=0.6,markersize=8)
                    inband_lolss=1   
                    legend_handles.append(handle)
                    legend_labels.append('LoLSS\ninband')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=1, fmt='o', color='grey', alpha=0.6,markersize=8)
                    

            elif 'inband_' in str(label_array[s]):
                if inband_lotss == 0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          label='LoTSS\ninband', zorder=1, fmt='o', color='deeppink', alpha=0.5, markersize=8)
                    inband_lotss = 1
                    legend_handles.append(handle)
                    legend_labels.append('LoTSS\ninband')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=1, fmt='o', color='deeppink', alpha=0.5 ,markersize=8)
                                        
            elif 'LoTSS' in str(label_array[s]):
                if lotss==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='LoTSS', zorder=10, fmt='v', color='darkgreen', markersize=8)
                    lotss=1   
                    legend_handles.append(handle)
                    legend_labels.append('LoTSS')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='v', color='darkgreen', markersize=8)            
            elif 'LoLSS' in str(label_array[s]):
                if lolss==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='LoLSS', zorder=10, fmt='^', color='black', markersize=8)
                    lolss=1   
                    legend_handles.append(handle)
                    legend_labels.append('LoLSS')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='^', color='black', markersize=8)                         
            elif 'NVSS' in str(label_array[s]):
                if nvss==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='NVSS', zorder=11, fmt='>', color='red', markersize=8)
                    nvss=1   
                    legend_handles.append(handle)
                    legend_labels.append('NVSS')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=11, fmt='>', color='red', markersize=8)   
            elif 'TGSS' in str(label_array[s]):
                if tgss==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='TGSS', zorder=10, fmt='<', color='blue', markersize=8)
                    tgss=1   
                    legend_handles.append(handle)
                    legend_labels.append('TGSS')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='<', color='blue', markersize=8)     
            elif 'FIRST' in str(label_array[s]):
                if first==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='FIRST', zorder=10, fmt='s', color='darkorange', markersize=8)
                    first=1   
                    legend_handles.append(handle)
                    legend_labels.append('FIRST')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='s', color='darkorange', markersize=8) 
            elif 'VLASS' in str(label_array[s]):
                if vlass==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='VLASS', zorder=10, fmt='d', color='darkmagenta', markersize=8)
                    vlass=1   
                    legend_handles.append(handle)
                    legend_labels.append('VLASS')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='d', color='darkmagenta', markersize=8) 
                    
            elif 'VLA' in str(label_array[s]):
                if vla==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='VLA', zorder=10, fmt='X', color='magenta', markersize=8)
                    vla=1   
                    legend_handles.append(handle)
                    legend_labels.append('VLA')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='X', color='purple', markersize=8) 
            elif 'NRAO' in str(label_array[s]):
                if green==0:
                    handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='NRAO Green Bank', zorder=10, fmt='p', color='green', markersize=8)
                    green=1   
                    legend_handles.append(handle)
                    legend_labels.append('NRAO\nGreen Bank')
                else:
                    ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='p', color='green', markersize=8)                     


            else:    
                handle=ax.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                      label=label_array[s], zorder=10, fmt='o', marker=markers[s], markersize=8, color=colors[s])
                legend_handles.append(handle)
                legend_labels.append(label_array[s])
                print('there are unspecified datapoints')


    #fit the curved models
    fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit = \
                spectral_index_eval_curve(freq_array, flux_array, error_array)
    # print('parameters for the curve are', poptgen )
    if did_it_fit:
        ax.plot(x_range_full, curve(x_range_full, poptgen[0], poptgen[1] , poptgen[2], poptgen[3] ), color='black')


    return legend_handles, legend_labels


def make_sed_subplot(galaxy_names, save_fig=False):
    """
    Plots the SEDs of multiple galaxies in a subplot with a singular legend.
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=False)
    handles, labels = [],[]

    # Iterate over the provided galaxy names
    for i, galaxy_name in enumerate(galaxy_names):
        ax = axs[i // 2, i % 2]

        # Call your make_sed_singular function with the current galaxy name and axis
        handles_sub, labels_sub = make_sed_singular(galaxy_name, ax=ax)
        handles.extend(handles_sub)
        labels.extend(labels_sub)

    plt.tight_layout()

    # Add a common legend
    plt.subplots_adjust(left=0.19)
    fig.legend(handles, labels, bbox_to_anchor=(0.17,0.7), fontsize=13.5)

    # Save figure if required
    if save_fig:
        plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/PS_extreme/combinedfigure.pdf')
        print('Figure saved.')



# List of galaxy names to include in the subplot
galaxies_to_plot = [tbdata['LoTSS_name'][98418], tbdata['LoTSS_name'][103164],
                    tbdata['LoTSS_name'][83458], tbdata['LoTSS_name'][39286]]

# Call the make_sed_subplot function with the list of galaxy names
make_sed_subplot(galaxies_to_plot, save_fig=True)

        
