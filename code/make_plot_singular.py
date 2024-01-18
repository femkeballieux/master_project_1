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




#read in the data
# hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/match_knownPS_mine_10arcsec_old.fits') #For checking with known PS
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_clean.fits') #for checking with whole master sample

# high_survey = 'NVSS' #could also be first, but we use NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#print the columnnames
# print(orig_cols)

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
freq_LoLLS_ch0 = 44.
freq_LoLLS_ch1 = 48.
freq_LoLLS_ch2 = 52.
freq_LoLLS_ch3 = 56.
freq_LoLLS_ch4 = 60.
freq_LoLLS_ch5 = 64.
freq_vlass = 3000.
 


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


def curve_convex(freq, freq_peak, flux_peak, alphathick, alphathin, a, alpha_low): 
    """
    Model taken from Tschager et al. 2003. General fit not based on any physics. Adapted for an uptick
    """
    return flux_peak/(1 -np.exp(-1))*((freq/freq_peak)**alphathick)*(1-np.exp(-(freq/freq_peak)**(alphathin-alphathick))) + a * freq ** alpha_low





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

def spectral_index_eval_curve_convex(freq, flux, flux_err):
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data

    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'
    
    fluxplot = flux[(flux > 0)]
    freqplot = freq[flux > 0]
    flux_errplot = flux_err[(flux > 0)]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0gen = [peakfreq[0], peakflux, 0.8, -0.7, 200, -0.8] #TODO: change initial values

    did_it_fit = True
    try:
        poptgen, pcovgen = opt.curve_fit(curve_convex, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
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
                freq_inband_mid, freq_inband_high,freq_LoLLS_ch0, freq_LoLLS_ch1, freq_LoLLS_ch2 , freq_LoLLS_ch3 \
                    , freq_LoLLS_ch4 , freq_LoLLS_ch5, freq_vlass])
    
    
    #labels in order they are used
    label_array = np.array(['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high', \
                  'LoLLS_ch0', 'LoLLS_ch1', 'LoLLS_ch2', 'LoLLS_ch3', 'LoLLS_ch4', 'LoLLS_ch5', 'VLASS' ] )  

    fig, ax = plt.subplots(figsize=(10,8))
    plt.xlabel('Frequency [MHz]', fontsize=20)
    plt.ylabel('Flux Density [Jy]',fontsize=20)
    # plt.title(str(index)+'_'+galaxy_name, fontsize=20)
    plt.title(str(index) + galaxy_name, fontsize=20)
    plt.yscale('log')
    plt.xscale('log')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='both', labelsize=15, size=5)
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        
    
    """
    Comment back in when plotting these specific galaxies
    """
    if galaxy_name == tbdata['LoTSS_name'][39286]:
        plt.title('87GB B0903+6656', fontsize=20)
        plt.xlabel('Frequency [GHz]', fontsize=20)
        plt.ylabel('Flux Density [mJy]',fontsize=20)

        flux_array = np.array([flux_LoTSS[index]*1000, flux_NVSS[index]*1000,  74.9, 195 ,flux_vlass[index]*1000]) #mJy
        error_array = np.array([error_LoTSS[index]*1000, error_NVSS[index]*1000, 11.235, 29.25, e_flux_vlass[index]*1000])
        freq_array = np.array([freq_LoTSS/1000, freq_NVSS/1000, 8.400 , 4.850 , freq_vlass/1000]) #GHz
        
        label_array = np.array(['LoTSS', 'NVSS',  'VLA', 'NRAO Green Bank' ,'VLASS' ] ) 
        x_range_full= np.linspace(0.050, 9.000, 1000)
        x_range_low = np.linspace(0.050, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)
        
    elif galaxy_name == tbdata['LoTSS_name'][76243]:
        flux_array = np.array([flux_LoTSS[index], flux_NVSS[index], flux_TGSS[index], flux_VLSSr[index], flux_LoLSS[index],\
                     flux_FIRST[index], flux_inband_low[index], flux_inband_mid[index], flux_inband_high[index], flux_channel0[index],\
                         flux_channel1[index], flux_channel2[index], flux_channel3[index], flux_channel4[index], flux_channel5[index],\
                             flux_vlass[index], 0.78, 0.71, 0.305, 0.3718, 0.4829, 0.318, 0.14])
        error_array = np.array([error_LoTSS[index], error_NVSS[index], error_TGSS[index], error_VLSSr[index], error_LoLSS[index],\
                     error_FIRST[index], error_inband_low[index], error_inband_mid[index], error_inband_high[index], error_channel0[index],\
                         error_channel1[index], error_channel2[index], error_channel3[index], error_channel4[index], error_channel5[index],\
                             e_flux_vlass[index], 0.1, 0.058, 0.045, 0.03718,0.04829, 0.002, 0.014])
        #list of frequencies in MHz
        freq_array = np.array([freq_LoTSS, freq_NVSS, freq_TGSS, freq_VLSSr, freq_LoLSS, freq_FIRST , freq_inband_low,\
                    freq_inband_mid, freq_inband_high,freq_LoLLS_ch0, freq_LoLLS_ch1, freq_LoLLS_ch2 , freq_LoLLS_ch3 \
                        , freq_LoLLS_ch4 , freq_LoLLS_ch5, freq_vlass, 74, 365, 4850, 5000, 8400, 15000, 22000])
            
        label_array = np.array(['LoTSS', 'NVSS', 'TGSS', 'VLSSr', 'LoLSS', 'FIRST', 'inband_low', 'inband_mid', 'inband_high', \
                      'LoLLS_ch0', 'LoLLS_ch1', 'LoLLS_ch2', 'LoLLS_ch3', 'LoLLS_ch4', 'LoLLS_ch5', 'VLASS' , 'VLA', 'Texas', '', 'VLBA', 'VLA', 'OVRO', 'VERA'] )  
        
        plt.title('SDSS J121711.02+583526.6', fontsize=20)
        plt.xlabel('Frequency [MHz]', fontsize=20)
        plt.ylabel('$S_{total}$ [Jy]',fontsize=20)
         
        # x_range_full= np.linspace(0.050, 9.000, 1000)
        # x_range_low = np.linspace(0.050, 0.144, 1000)
        # x_range_mid = np.linspace(0.144, 1.400, 1000)
        # x_range_high = np.linspace(1.400, 3.000, 1000)

    elif galaxy_name == tbdata['LoTSS_name'][66495]: #the convex SED
        plt.title('NGC 3894', fontsize=20)
        plt.xlabel('Frequency [GHz]', fontsize=20)
        plt.ylabel('Flux Density [mJy]',fontsize=20)

        flux_array = np.array([flux_LoTSS[index]*1000, flux_NVSS[index]*1000, flux_TGSS[index]*1000, flux_vlass[index]*1000, flux_LoLSS[index]*1000, flux_FIRST[index]*1000, \
                               flux_channel0[index]*1000, flux_channel1[index]*1000, flux_channel2[index]*1000, flux_channel3[index]*1000, flux_channel4[index]*1000, flux_channel5[index]*1000,\
                                524 , 385, 450, 279]) #mJy
        error_array = np.array([error_LoTSS[index]*1000, error_NVSS[index]*1000, error_TGSS[index]*1000, e_flux_vlass[index]*1000, error_LoLSS[index]*1000, error_FIRST[index]*1000, \
                                error_channel0[index]*1000, error_channel1[index]*1000, error_channel2[index]*1000, error_channel3[index]*1000, error_channel4[index]*1000, error_channel5[index]*1000,\
                                   52, 38.5, 10, 27.9 ]) #For the 5GHz, 15GHz, took 10% uncertainty as no reported
        freq_array = np.array([freq_LoTSS/1000, freq_NVSS/1000, freq_TGSS/1000, freq_vlass/1000, freq_LoLSS/1000, freq_FIRST/1000, freq_LoLLS_ch0/1000, freq_LoLLS_ch1/1000, \
                               freq_LoLLS_ch2/1000 , freq_LoLLS_ch3 /1000, freq_LoLLS_ch4/1000 , freq_LoLLS_ch5/1000, 5, 8, 10.45, 15]) #GHz
        
        label_array = np.array(['LoTSS', 'NVSS',  'TGSS', 'VLASS' ,'LoLSS', 'FIRST', 'LoLLS_ch0', 'LoLLS_ch1', 'LoLLS_ch2', 'LoLLS_ch3', 'LoLLS_ch4', 'LoLLS_ch5', 'NRAO', 'VLBI','Effelsberg', 'VLBI' ] ) 
        x_range_full= np.linspace(0.040, 18.000, 3000)
        x_range_low = np.linspace(0.050, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)    
        
    elif galaxy_name in [tbdata['LoTSS_name'][103164], tbdata['LoTSS_name'][98418] , tbdata['LoTSS_name'][83458]]:
        x_range_full= np.linspace(0.050, 9.000, 1000)
        x_range_low = np.linspace(0.050, 0.144, 1000)
        x_range_mid = np.linspace(0.144, 1.400, 1000)
        x_range_high = np.linspace(1.400, 3.000, 1000)
        
        plt.xlabel('Frequency [GHz]', fontsize=20)
        plt.ylabel('$S_{total}$ [mJy]',fontsize=20)
        
        flux_array = flux_array* 1000
        error_array = error_array * 1000
        freq_array = freq_array / 1000


    print(label_array)
    inband_lolss = 0 #fixing the legends
    inband_lotss = 0
    VLBI=0
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

    #Plot the powerlaw fits for VLASS        
    # if alpha_low_VLASS[index]!=0.:
    #     #special sources
    #     if galaxy_name==tbdata['LoTSS_name'][39286]:
    #         plt.plot(x_range_mid, a_low_VLASS[index] *1000 * (1000** (alpha_low_VLASS[index])) * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='indigo', linewidth=2) #, label='powerlaw low GPS') 
    #         plt.plot(x_range_high, a_high_VLASS[index]  *1000 * (1000** (alpha_high_VLASS[index])) * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='darkorange', linewidth=2) #, label='powerlaw high GPS') 
                
    #     elif galaxy_name == tbdata['LoTSS_name'][103164]:
    #         plt.plot(x_range_mid, a_low_VLASS[index] *1000 * (1000** (alpha_low_VLASS[index])) * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='indigo', linewidth=2) #, label='powerlaw low GPS') 
    #         plt.plot(x_range_high, a_high_VLASS[index]  *1000 * (1000** (alpha_high_VLASS[index])) * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='darkorange', linewidth=2) #, label='powerlaw high GPS')   
    #         if alpha_low_LoLSS[index]!=0.:
    #             plt.plot(x_range_low, a_low_LoLSS[index] * 1000 * (1000**alpha_low_LoLSS[index])* (x_range_low ** alpha_low_LoLSS[index]),\
    #                       color='green', linewidth=2) 
                
    #     elif galaxy_name == tbdata['LoTSS_name'][98418]:
    #         plt.plot(x_range_mid, a_low_VLASS[index] *1000 * (1000** (alpha_low_VLASS[index])) * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='indigo', linewidth=2) #, label='powerlaw low GPS') 
    #         plt.plot(x_range_high, a_high_VLASS[index]  *1000 * (1000** (alpha_high_VLASS[index])) * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='darkorange', linewidth=2.5) #, label='powerlaw high GPS') 
    #         if alpha_low_LoLSS[index]!=0.:
    #             plt.plot(x_range_low, a_low_LoLSS[index] * 1000 * (1000**alpha_low_LoLSS[index])* (x_range_low ** alpha_low_LoLSS[index]),\
    #                       color='green', linewidth=2) 
                
    #     elif galaxy_name == tbdata['LoTSS_name'][83458]:
    #         plt.plot(x_range_mid, a_low_VLASS[index] *1000 * (1000** (alpha_low_VLASS[index])) * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='indigo', linewidth=2) #, label='powerlaw low GPS') 
    #         plt.plot(x_range_high, a_high_VLASS[index]  *1000 * (1000** (alpha_high_VLASS[index])) * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='darkorange', linewidth=2) #, label='powerlaw high GPS') 
    #         if alpha_low_LoLSS[index]!=0.:
    #             plt.plot(x_range_low, a_low_LoLSS[index] * 1000 * (1000**alpha_low_LoLSS[index])* (x_range_low ** alpha_low_LoLSS[index]),\
    #                       color='green', linewidth=2) 
                
    #     elif galaxy_name == tbdata['LoTSS_name'][66495]:
    #         plt.plot(x_range_mid, a_low_VLASS[index] *1000 * (1000** (alpha_low_VLASS[index])) * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='indigo', linewidth=2) #, label='powerlaw low GPS') 
    #         plt.plot(x_range_high, a_high_VLASS[index]  *1000 * (1000** (alpha_high_VLASS[index])) * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='darkorange', linewidth=2) #, label='powerlaw high GPS') 
    #         if alpha_low_LoLSS[index]!=0.:
    #             plt.plot(x_range_low, a_low_LoLSS[index] * 1000 * (1000**alpha_low_LoLSS[index])* (x_range_low ** alpha_low_LoLSS[index]),\
    #                       color='green', linewidth=2) 
          
    # #If VLASS exists, plot the VLASS high power law, and the LoLSS low powerlaw   
    #     else:      
    #         plt.plot(x_range_mid, a_low_VLASS[index] * (x_range_mid ** alpha_low_VLASS[index]),\
    #                   color='green', label='powerlaw low VLASS', linewidth=2) 
    #         plt.plot(x_range_high, a_high_VLASS[index] * (x_range_high ** alpha_high_VLASS[index]),\
    #                   color='black', label='powerlaw high VLASS', linewidth=2)
                              
                        
    #         if alpha_low_LoLSS[index]!=0.:
    #             plt.plot(x_range_low, a_low_LoLSS[index] * (x_range_low ** alpha_low_LoLSS[index]),\
    #                   color='magenta', label='powerlaw low LoLSS', linewidth=2) 
    
                    
        
    #         elif alpha_low_LoLSS[index]!=0.: #if VLASS does not exist, plot the LoLSS powerlaw
            
    #             plt.plot(x_range_low, a_low_LoLSS[index] * (x_range_low ** alpha_low_LoLSS[index]),\
    #                       color='magenta', label='powerlaw low LoLSS', linewidth=2) 
    #             plt.plot(x_range_mid, a_high_LoLSS[index] * (x_range_mid ** alpha_high_LoLSS[index]),\
    #                       color='green', label='powerlaw high LoLSS', linewidth=2)
                

    # else:
    #     print('no spectral index fit')
    
    #fit the curved models

    if galaxy_name == tbdata['LoTSS_name'][66495]:
        fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit = \
                    spectral_index_eval_curve_convex(freq_array, flux_array, error_array)
        plt.ylim(ymin=200, ymax=600)
        # print('parameters for the curve are', poptgen )
        if did_it_fit:
            plt.plot(x_range_full, curve_convex(x_range_full, poptgen[0], poptgen[1] , poptgen[2], poptgen[3], poptgen[4], poptgen[5] ), color='black')
    else:
        fluxplot, freqplot, flux_errplot, poptgen, pcovgen, did_it_fit = \
                    spectral_index_eval_curve(freq_array, flux_array, error_array)
        # print('parameters for the curve are', poptgen )
        if did_it_fit:
            plt.plot(x_range_full, curve(x_range_full, poptgen[0], poptgen[1] , poptgen[2], poptgen[3] ), color='black', linestyle='dotted')
    
    #make legend
    plt.legend(bbox_to_anchor=(0.8, 0.54), fontsize=15)
    
    #save figure
    if save_fig:
        plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/known_PS/'+ galaxy_name + '.pdf') #edit when doing random PS sources
        # plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/code/PS_seds/known_PS/'+str(index)+'_' + galaxy_name + '.pdf') #edit when doing random PS sources

        print('figure saved for', galaxy_name)
    # plt.show()


# For checking with known PS, make sure to select the right file at the top:
    #in HF not PS: 4, 15,21, 67, 90,93, 94,102,110,129,136,162,165,169,171,184,190,200,204,217,221,279,280,284,295,296,299, 315, 316
# for i in [31,35,48,50,53,62,63,93,103,106,122,138,146,148, 181,202,205,216,219,234,237,238,243,245,248,249,250,255,260,261,263,267]:
#     make_sed_singular(tbdata['LoTSS_name'][i], save_fig=True)


# These 4 are for the extreme spectral index plots, from master clean
# make_sed_singular(tbdata['LoTSS_name'][98418], save_fig=True)
# make_sed_singular(tbdata['LoTSS_name'][103164], save_fig=True)
# make_sed_singular(tbdata['LoTSS_name'][83458], save_fig=True)
# make_sed_singular(tbdata['LoTSS_name'][39286], save_fig=True)

#convex source
make_sed_singular(tbdata['LoTSS_name'][66495], save_fig=True) 


#Make all PS seds
# counter=0
# for i, name in enumerate(name_list):
#     #if (alpha_low[i] >= 0.1) & (alpha_high[i]<=0): #select when a source is PS
#     make_sed_singular(name, save_fig=True)
#     counter +=1
#     print(counter, '/767')

        
