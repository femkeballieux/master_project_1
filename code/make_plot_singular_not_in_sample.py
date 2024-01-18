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
from matplotlib.ticker import ScalarFormatter



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

# def spectral_index_eval_curve_convex(freq, flux, flux_err):
#     fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data

#     if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
#         return 'No flux?'
    
#     fluxplot = flux[(flux > 0)]
#     freqplot = freq[flux > 0]
#     flux_errplot = flux_err[(flux > 0)]

#     peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
#     peakflux = max(fluxplot)
#     p0gen = [peakfreq[0], peakflux, 0.8, -0.7, 200, -0.8] #TODO: change initial values

#     did_it_fit = True
#     try:
#         poptgen, pcovgen = opt.curve_fit(curve_convex, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
#     except (RuntimeError, TypeError, ValueError):
#         print('Curve_fit could not fit curve model.')
#         did_it_fit = False
#         return (fluxplot,freqplot,flux_errplot,0,0, did_it_fit)
        
#         print('fitting went right')
#     return (fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit)
    


def make_sed_singular(galaxy_name, use_index=False, save_fig=False):
    """
    returns the sed for a singular galaxy
    """
    fig, ax = plt.subplots(figsize=(10,6))

    plt.yscale('log')
    plt.xscale('log')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='both', labelsize=15, size=5)
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        
    if galaxy_name == 'B1315+415':
        # plt.title('B1315+415', fontsize=20)
        plt.xlabel('Frequency [GHz]', fontsize=20)
        plt.ylabel('Flux Density [mJy]',fontsize=20)

        flux_array = np.array([140, 266.5, 270, 235, 207.3, 171, 153.7]) #mJy
        error_array = np.array([20, 8, 27, 30, 20, 17, 7.4])
        freq_array = np.array([0.408, 1.4, 2.3, 4.85, 5, 8.4, 10.45]) #GHz
        label_array = np.array(['Northern Cross Radiotelescope', 'NVSS', 'VCS', 'NRAO Green Bank', 'VIPS', 'VLA', 'Effelsberg'] ) 

        x_range_full= np.linspace(0.3, 13.000, 1000)
        x_range_low = np.linspace(0.3, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)
        

    inband_lolss = 0 #fixing the legends
    inband_lotss = 0
    VLBI=0
    VLBA=0
    markers=['v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4', 'v', '>', '<', '^', 'o', 's', 'D', '1', '2', '3', '4']
    colors=['darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'deeppink', \
            'mediumorchid', 'darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'blue', 'mediumorchid','darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'green', \
                    'mediumorchid', 'darkorange', 'red', 'blue', 'darkmagenta', 'black', 'magenta', 'green', 'mediumorchid',]
    for s in range(len(freq_array)):
        if flux_array[s] != 0.:
            
            if 'LoLLS_ch' in str(label_array[s]):
                if inband_lolss==0:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='LoLSS inband', zorder=1, fmt='o', color='grey', alpha=0.6,markersize=8)
                    inband_lolss=1   
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=1, fmt='o', color='grey', alpha=0.6,markersize=8)
                        
            elif 'inband_' in str(label_array[s]):
                if inband_lotss == 0:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          label='LoTSS inband', zorder=10, fmt='o', color='darkgreen', alpha=0.35,markersize=8)
                    inband_lotss = 1
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=1, fmt='o', color='darkgreen', alpha=0.35 ,markersize=8)
                    
            elif 'LoTSS' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='LoTSS', zorder=10, fmt='v', color='darkgreen', markersize=8)
                    
            elif 'LoLSS' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                            label='LoLSS', zorder=10, fmt='^', color='black', markersize=8)
                
            elif 'NVSS' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='NVSS', zorder=11, fmt='>', color='red', markersize=8)
                    
            elif label_array[s]=='VLASS':
                if error_array[s]==-1.:
                    plt.plot(freq_array[s], flux_array[s],\
                            label=label_array[s], zorder=10, marker='d', color=colors[s] ,markersize=8)  
                    print('UPPERLIMIT')
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                      label=label_array[s], zorder=10, fmt='d', color=colors[s],markersize=8)
                    
            elif 'FIRST' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='FIRST', zorder=10, fmt='s', color='darkorange', markersize=8)

            elif 'Effelsberg' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='Effelsberg', zorder=10, fmt='^', color='blue', markersize=8)                

            elif 'Green Bank' in str(label_array[s]):
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                              label='NRAO Green Bank', zorder=10, fmt='p', color='green', markersize=8)   
                
            elif 'VLBI' in str(label_array[s]):
                if VLBI == 0:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          label='VLBI', zorder=10, fmt='X', color='magenta', markersize=8)
                    VLBI = 1
                else:
                    plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                          zorder=10, fmt='X', color='magenta', markersize=8)


            else:    
                plt.errorbar(freq_array[s], flux_array[s], yerr=error_array[s],\
                      label=label_array[s], zorder=10, fmt='o', marker=markers[s], markersize=8, color=colors[s])


    if galaxy_name == 'B1315+415':
        fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit = \
                    spectral_index_eval_curve(freq_array, flux_array, error_array)
        # print('parameters for the curve are', poptgen )
        if did_it_fit:
            plt.plot(x_range_full, curve(x_range_full, poptgen[0], poptgen[1] , poptgen[2], poptgen[3] ), color='black', linestyle='solid')
    
    #make legend
    plt.legend(bbox_to_anchor=(0.7, 0.5), fontsize=15)
    plt.xlim(.3,13)
    plt.ylim(90,325)
    plt.tight_layout()
    
    #save figure
    if save_fig:
        plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/+ ' + galaxy_name + '.pdf') #edit when doing random PS sources
        # plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/known_PS/'+str(index)+'_' + galaxy_name + '.pdf') #edit when doing random PS sources

        print('figure saved for', galaxy_name)


#convex source
make_sed_singular('B1315+415', save_fig=True) 


