"""
Code written by Femke Ballieux, takes in the master sample and crossmatches it with vlass
"""
#run on laptop
import os
import numpy as np
from astropy.io import fits
import time
from tqdm import tqdm
import scipy.optimize as opt
import gpscssmodels


def crossmatching1():
    """
    Crossmatching is done seperatel for LoTSS-NVSS with VLASS, since VLASS has a high resolution.
    params=3 arcsec
    """
    os.system('java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts \
              tmatchn join1=always matcher=sky multimode=pairs nin=2 params=3 \
        in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_NVSS_LoTSS.fits values1="RA DEC" \
        in2=/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/vlass_catalog_10000.fits values2="ra dec" \
        out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/mega_master_10000_intermediate.fits')
    
    print("crossmatching 1 done")
    
def crossmatching2():
    
    """
    Crossmatching for the LoTSS NVSS VLASS sample with all the other radio surveys.
    params=15 arcsec
    """
    os.system('java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts \
              tmatchn join1=always matcher=sky multimode=pairs nin=13 params=15 \
        in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/mega_master_10000_intermediate.fits values1="RA_1 DEC_1" \
        in2=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/TGSS.fits values2="RAJ2000 DEJ2000" \
        in3=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/VLSSr.fits values3="RAJ2000 DEJ2000" \
        in4=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/LoLLS_new_names.fits values4="RA DEC"\
        in5=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/first_14dec17.fits values5="RA DEC"\
        in6=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/inband_spec_LoTSS.fits values6="RA DEC" \
        in7=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband/channel_0_source.fits values7="RA DEC" \
        in8=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband/channel_1_source.fits values8="RA DEC" \
        in9=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband/channel_2_source.fits values9="RA DEC" \
        in10=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband/channel_3_source.fits values10="RA DEC"\
        in11=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband//channel_4_source.fits values11="RA DEC"\
        in12=/net/vdesk/data2/bach1/ballieux/master_project_1/data/LoLSS_inband/channel_5_source.fits values12="RA DEC" \
        in13=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/WENSS.fits values13="_RAJ2000 _DEJ2000" \
        out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/mega_master_10000.fits')
    
    print("crossmatching 2 done")
# crossmatching1()
# crossmatching2()

"""
We now have a crossmatched table, this needs to be cleaned, put in the same units,
and for some the errors need to be convoluted with a percentage of the flux.
Below the spectral indices will be calculated as well 
"""

#TODO: select only isoladed LoLSS sources? we do something somewhere with LoLSS rms? See make_tab_PS

#Load in the data
hdulist = fits.open("/net/vdesk/data2/bach1/ballieux/master_project_1/data/mega_master_10000.fits")
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()
# print(orig_cols)

#For all the surveys, define properly the arrays such that they are in the right units, etc
#LoTSS
Name = tbdata['Source_Name_1'].rstrip()
RA = tbdata['RA_1']
DEC = tbdata['DEC_1']

S_LoTSS, stat_e_S_LoTSS = np.array(tbdata['Total_flux_1'])/1000., np.array(tbdata['E_Total_flux_1'])/1000. # Jy
e_S_LoTSS = np.sqrt(stat_e_S_LoTSS**2 + (0.1*S_LoTSS)**2) # Combine the statistical LoTSS error and overall error on LoTSS

#NVSS
S_NVSS, e_S_NVSS = np.array(tbdata['S1_4'])/1000., 0.1 * np.array(tbdata['S1_4'])/1000. # Jy
S_NVSS, e_S_NVSS = np.where(np.isnan(S_NVSS), 0., S_NVSS), np.where(np.isnan(S_NVSS), 0., e_S_NVSS) #deal with any non-detections

#VLASS
S_VLASS, stat_e_S_VLASS = np.array(tbdata['int_flux']), np.array(tbdata['err_int_flux']) # Jy
S_VLASS, stat_e_S_VLASS = np.where(np.isnan(S_VLASS), 0., S_VLASS), np.where(np.isnan(S_VLASS), 0., stat_e_S_VLASS) #deal with any non-detections
badfit_mask = (stat_e_S_VLASS== -1. ) #anything with a bad fit has this flag, then no error can be attributed

#make sure any non-detections become upper limits with error -1, value equal to sensitivity of vlass
# vlass_sensitivity = 120e-6 #Jy
# S_VLASS, stat_e_S_VLASS = np.where(np.isnan(S_VLASS),\
#             vlass_sensitivity, S_VLASS), np.where(np.isnan(stat_e_S_VLASS), -1., stat_e_S_VLASS)
#propagate the error with 15 percent
e_S_VLASS= np.sqrt(stat_e_S_VLASS**2 + (0.15 * S_VLASS)**2)
#e_S_VLASS= np.where(e_S_VLASS==1., -1., e_S_VLASS) #non-detections have error -1
e_S_VLASS[badfit_mask] = 0.25 * S_VLASS[badfit_mask] #bad fits have error of 25% instead of 15% convolution
#TODO: maybe just change this to anything with a flag

#TGSS
S_TGSS, e_S_TGSS = np.array(tbdata['Stotal'])/1000., np.array(tbdata['e_Stotal'])/1000. # Jy/beam
S_TGSS, e_S_TGSS = np.where(np.isnan(S_TGSS), 0, S_TGSS), np.where(np.isnan(S_TGSS), 0, e_S_TGSS)

#VLSSr
S_VLSSr, e_S_VLSSr = np.array(tbdata['Sp']), 0.1 * np.array(tbdata['Sp']) # Jy
S_VLSSr, e_S_VLSSr = np.where(np.isnan(S_VLSSr), 0, S_VLSSr), np.where(np.isnan(S_VLSSr), 0, e_S_VLSSr)

#LoLSS
S_LoLSS, e_S_LoLSS = np.array(tbdata['Total_flux_LoLLS'])/1000., np.array(tbdata['E_Total_flux_LoLLS'])/1000. # Jy/beam
S_LoLSS, e_S_LoLSS = np.where(np.isnan(S_LoLSS), 0, S_LoLSS), np.where(np.isnan(S_LoLSS), 0, e_S_LoLSS)

#FIRST
S_FIRST, e_S_FIRST = np.array(tbdata['FINT'])/1000, 0.1*(np.array(tbdata['FINT'])/1000) #Jy, set error to 10% of flux
S_FIRST, e_S_FIRST = np.where(np.isnan(S_FIRST), 0, S_FIRST), np.where(np.isnan(S_FIRST), 0, e_S_FIRST)

#inband LoTSS
S_inband_low, e_S_inband_low = np.array(tbdata['M1_L1flux']), 0.1 * np.array(tbdata['M1_L1flux']) #mJy
S_inband_low, e_S_inband_low = np.where(np.isnan(S_inband_low), 0, S_inband_low), np.where(np.isnan(S_inband_low), 0, e_S_inband_low)

S_inband_mid, e_S_inband_mid = np.array(tbdata['M1_L2flux']), 0.1 * np.array(tbdata['M1_L2flux']) #mJy
S_inband_mid, e_S_inband_mid = np.where(np.isnan(S_inband_mid), 0, S_inband_mid), np.where(np.isnan(S_inband_mid), 0, e_S_inband_mid)

S_inband_high, e_S_inband_high = np.array(tbdata['M1_L3flux']), 0.1 * np.array(tbdata['M1_L3flux']) #mJy
S_inband_high, e_S_inband_high = np.where(np.isnan(S_inband_high), 0, S_inband_high), np.where(np.isnan(S_inband_high), 0, e_S_inband_high)

#inband_LoLSS
S_channel0, e_S_channel0 = tbdata['Total_flux_7']/1000., tbdata['E_Total_flux_7']/1000. # Jy
S_channel0, e_S_channel0 = np.where(np.isnan(S_channel0), 0, S_channel0), np.where(np.isnan(S_channel0), 0, e_S_channel0)

S_channel1, e_S_channel1 = tbdata['Total_flux_8']/1000., tbdata['E_Total_flux_8']/1000. # Jy
S_channel1, e_S_channel1 = np.where(np.isnan(S_channel1), 0, S_channel1), np.where(np.isnan(S_channel1), 0, e_S_channel1)

S_channel2, e_S_channel2 = tbdata['Total_flux_9']/1000., tbdata['E_Total_flux_9']/1000. # Jy
S_channel2, e_S_channel2 = np.where(np.isnan(S_channel2), 0, S_channel2), np.where(np.isnan(S_channel2), 0, e_S_channel2)

S_channel3, e_S_channel3 = tbdata['Total_flux_10']/1000., tbdata['E_Total_flux_10']/1000. # Jy
S_channel3, e_S_channel3 = np.where(np.isnan(S_channel3), 0, S_channel3), np.where(np.isnan(S_channel3), 0, e_S_channel3)

S_channel4, e_S_channel4 = tbdata['Total_flux_11']/1000., tbdata['E_Total_flux_11']/1000. # Jy
S_channel4, e_S_channel4 = np.where(np.isnan(S_channel4), 0, S_channel4), np.where(np.isnan(S_channel4), 0, e_S_channel4)

S_channel5, e_S_channel5 = tbdata['Total_flux_12']/1000., tbdata['E_Total_flux_12']/1000. # Jy
S_channel5, e_S_channel5 = np.where(np.isnan(S_channel5), 0, S_channel5), np.where(np.isnan(S_channel5), 0, e_S_channel5)

#WENSS
WENSS_intflux= tbdata['Sint']/1000. #Jy
WENSS_mask = WENSS_intflux < 0.
WENSS_intflux[WENSS_mask] = 0.
WENSS_intflux[~WENSS_mask] = 0.9 * WENSS_intflux[~WENSS_mask] #overestimation correction
WENSS_intflux_error = 0.15 * WENSS_intflux #15 percent of flux density

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
freq_VLASS = 3000.


def spectral_index_eval_norms(freq,flux,flux_err, alpha_in=0.7):
    """
    This is a function that fits the power law to an array of frequencies, fluxes and flux errors
    that are entered as input. First it is checked whether these are positive, if there are no
    positive fluxes then no fit is done
    """
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        print('Curve_fit could not fit powerlaw.')
        return 
    
    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

  #  peakfreq = freqplot[np.where(fluxplot == max(fluxplot))] #Frequency where flux is highest
    peakflux = max(fluxplot) #highest flux in the array
    p0pow = [peakflux, alpha_in]#[peakflux,0.7] #initial values for the fit, a=peakflux, alpha=0.7. Not sure why a=peak_flux as initial value?

    try:
        poptpowlaw, pcovpowlaw = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError): # If curve fit can not find a good fit, skip the source.    
        return 'Curve_fit could not fit powerlaw.'

    return (poptpowlaw, pcovpowlaw, fluxplot, freqplot, flux_errplot)

def spectral_index_eval_curve(freq,flux,flux_err):
    """
    This is a function that fits a curve to an array of frequencies, fluxes and flux errors
    that are entered as input. First it is checked whether these are positive, if there are no
    positive fluxes then no fit is done
    """
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data
    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0gen = [peakfreq, peakflux, 0.8, -0.7] #initial values, here I changed some things so check if it goes allright
    #TODO: 

    try:
        poptgen, pcovgen = opt.curve_fit(gpscssmodels.curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    return (fluxplot, freqplot, flux_errplot, poptgen, pcovgen)
    


def make_fluxes_LoLSS():
    """
    For the enitre catalog, makes arrays with all non-zero fluxes(freq)
    corresponding frequencies(flux), corresponding errors.
    Then an array with high arrays of freq and flux for LoTSs en NVSS, and low for loLSS en LoTSS is made.
    """
    S_match = np.where(((S_NVSS > 0.) & (S_LoLSS > 0.)))

    flux = np.vstack((S_LoTSS, S_LoLSS, S_NVSS))
    flux = flux[:,S_match[0]]
    flux_err = np.vstack((e_S_LoTSS, e_S_LoLSS, e_S_NVSS))
    flux_err = flux_err[:,S_match[0]] 

    freq = np.array([freq_LoTSS, freq_LoLSS, freq_NVSS])
    freq_high = np.array([freq_LoTSS, freq_NVSS])
    freq_low = np.array([freq_LoLSS, freq_LoTSS])

    flux_high = np.vstack((S_LoTSS, S_NVSS))
    flux_high = flux_high[:,S_match[0]]
    flux_low = np.vstack((S_LoLSS, S_LoTSS))
    flux_low = flux_low[:,S_match[0]]

    flux_err_high = np.vstack((e_S_LoTSS, e_S_NVSS))
    flux_err_high = flux_err_high[:,S_match[0]]
    flux_err_low = np.vstack((e_S_LoLSS, e_S_LoTSS))
    flux_err_low = flux_err_low[:,S_match[0]]
    
    return freq, freq_high, freq_low, flux, flux_err, flux_high, flux_low, flux_err_high, flux_err_low, S_match

def make_fluxes_VLASS():
    """
    For the catalog, makes arrays with all non-zero fluxes(freq)
    corresponding frequencies(flux), corresponding errors.
    Then an array with high arrays of freq and flux for NVSS and VLASS, en low for LoTSS and NVSS
    I still need to make another version of this function that does it for VLASS NVSS LoTSS
    """
    S_match = np.where(((S_NVSS > 0.) & (S_LoTSS > 0.) & (S_VLASS > 0.)) )

    flux = np.vstack((S_LoTSS, S_NVSS, S_VLASS))
    flux = flux[:,S_match[0]]
    flux_err = np.vstack((e_S_LoTSS, e_S_NVSS, e_S_VLASS))
    flux_err = flux_err[:,S_match[0]] 

    freq = np.array([freq_LoTSS, freq_NVSS, freq_VLASS])
    freq_high = np.array([freq_NVSS, freq_VLASS])
    freq_low = np.array([freq_LoTSS, freq_NVSS])

    flux_high = np.vstack((S_NVSS, S_VLASS))
    flux_high = flux_high[:,S_match[0]]
    flux_low = np.vstack((S_LoTSS, S_NVSS))
    flux_low = flux_low[:,S_match[0]]

    flux_err_high = np.vstack((e_S_NVSS, e_S_VLASS))
    flux_err_high = flux_err_high[:,S_match[0]]
    flux_err_low = np.vstack((e_S_LoTSS, e_S_NVSS))
    flux_err_low = flux_err_low[:,S_match[0]]
    
    return freq, freq_high, freq_low, flux, flux_err, flux_high, flux_low, flux_err_high, flux_err_low, S_match



"""
Here we do the fitting for LoLSs LoTSS NVSS, further in the code we do the same for LoTSS, NVSS, VLASS
"""
freq, freq_high, freq_low, flux, flux_err, flux_high, flux_low, flux_err_high, flux_err_low, S_match = make_fluxes_LoLSS()

#These are empty arrays
alpha_low, alpha_err_low, alpha_high, alpha_err_high, alpha_tot, alpha_err_tot, peak_flux_man, peak_freq_man, alphathick_err,\
      alphathin_err, freq_peak_gen_err, S_max_gen_err, norm_high, curve_fit_qual, norm_low = [np.zeros(np.shape(flux_low[1])[0]) for dummy in range(15)]

alpha_err_high_b = np.zeros((np.shape(flux_low[1])[0],2))
alpha_err_low_b = np.zeros((np.shape(flux_low[1])[0],2))


poptgen_store = np.zeros((4,np.shape(flux[1])[0]))

Name_source = Name[S_match]
Name_source_list = []
bad_sources_name = []
error = np.zeros((4,np.shape(flux[1])[0]))
PS_count = 0

#Here we do the actual fitting, to LoLSS LoTSS NVSS
start = time.time()
for i in tqdm(range(np.shape(flux_high[1])[0])):
    # print('Processing '+ Name_source[i])
    # Fit LoTSS to NVSS
    if spectral_index_eval_norms(freq_high,flux_high[:,i],flux_err_high[:,i]) == 'Curve_fit could not fit powerlaw.':
        print('Could not fit error in high powerlaw to '+ Name_source[i])
        continue
    else:
        poptpowlaw_high, pcovpowlaw_high, fluxplot_high, freqplot_high, flux_errplot_high = spectral_index_eval_norms(freq_high,flux_high[:,i],flux_err_high[:,i])
        alpha_high[i]=poptpowlaw_high[1]
        norm_high[i]=poptpowlaw_high[0]

        # Calculate uncertainty on alpha_high
        # find min and max lotss flux, min and max nvss flux
        lotss_min = flux[0,i] - flux_err[0,i]
        lotss_max = flux[0,i] + flux_err[0,i]
        nvss_min = flux[2,i] - flux_err[2,i]
        nvss_max = flux[2,i] + flux_err[2,i]

        # find corresponding a_high, a_high
        p0pow = [144.,alpha_high[i]]
        popt_err_high_l, pcov_err_high_l = opt.curve_fit(gpscssmodels.powlaw, [144.,1400.], [lotss_min, nvss_max], p0 = p0pow, maxfev = 10000)
        popt_err_high_up, pcov_err_high_up = opt.curve_fit(gpscssmodels.powlaw, [144.,1400.], [lotss_max, nvss_min], p0 = p0pow, maxfev = 10000)
        alpha_err_high_b[i] = [alpha_high[i] - popt_err_high_l[1], popt_err_high_up[1] - alpha_high[i]]
        alpha_err_high[i] = np.abs(np.median(alpha_err_high_b[i]))
        
    # Fit LoTSS to LoLSS  
    if spectral_index_eval_norms(freq_low,flux_low[:,i],flux_err_low[:,i]) == 'Curve_fit could not fit powerlaw.':
        print('Could not fit low powerlaw to '+ Name_source[i])
        continue
    else:
        poptpowlaw_low, pcovpowlaw_low, fluxplot_low, freqplot_low, flux_errplot_low = spectral_index_eval_norms(freq_low,flux_low[:,i],flux_err_low[:,i])
        alpha_low[i] = poptpowlaw_low[1]
        norm_low[i] = poptpowlaw_low[0]

        # Calculate uncertainty on alpha_low
        # find min and max lotss flux, min and max nvss flux
        lotss_min = flux[0,i] - flux_err[0,i]
        lotss_max = flux[0,i] + flux_err[0,i]
        lolss_min = flux[1,i] - flux_err[1,i]
        lolss_max = flux[1,i] + flux_err[1,i]

        # find corresponding a_low, a_high
        p0pow = [144.,alpha_low[i]]
        popt_err_low_l, pcov_err_low_l = opt.curve_fit(gpscssmodels.powlaw, [54.,144.], [lolss_min, lotss_max], p0 = p0pow, maxfev = 10000)
        popt_err_low_up, pcov_err_low_up = opt.curve_fit(gpscssmodels.powlaw, [54.,144.], [lolss_max, lotss_min], p0 = p0pow, maxfev = 10000)
        alpha_err_low_b[i] = [alpha_low[i] - popt_err_low_l[1], popt_err_low_up[1] - alpha_low[i]]
        alpha_err_low[i] = np.abs(np.median(alpha_err_low_b[i]))

    if -alpha_low[i] >= 0.1 and -alpha_high[i] <= 0.:
        PS_count += 1

    
print(PS_count, " PS sources for LoLSS LoTSS NVSS")
#need to turn the spectral indices into an array the same shape as the table
alpha_low_array = np.zeros(len(S_LoTSS))
alpha_low_array[S_match] = -alpha_low

alpha_err_low_array = np.zeros(len(S_LoTSS))
alpha_err_low_array[S_match] = alpha_err_low

alpha_high_array=np.zeros(len(S_LoTSS))
alpha_high_array[S_match] = -alpha_high

alpha_err_high_array = np.zeros(len(S_LoTSS))
alpha_err_high_array[S_match] = alpha_err_high

norm_low_array=np.zeros(len(S_LoTSS))
norm_low_array[S_match]= norm_low

norm_high_array=np.zeros(len(S_LoTSS))
norm_high_array[S_match]= norm_high

"""
Here we do the fitting for LoTSS NVSS VLASS. I have added a '2' to all variables since this helps keep track of them.
Not pretty but it works

"""
freq2, freq_high2, freq_low2, flux2, flux_err2, flux_high2, flux_low2, flux_err_high2, flux_err_low2, S_match2 = make_fluxes_VLASS()

#These are empty arrays
alpha_low2, alpha_err_low2, alpha_high2, alpha_err_high2, alpha_tot2, alpha_err_tot2, peak_flux_man2, peak_freq_man2, alphathick_err2,\
     alphathin_err2, freq_peak_gen_err2, S_max_gen_err2, norm_high2, curve_fit_qual2, norm_low2 = [np.zeros(np.shape(flux_low2[1])[0]) for dummy in range(15)]

alpha_err_high_b2 = np.zeros((np.shape(flux_low2[1])[0],2))
alpha_err_low_b2 = np.zeros((np.shape(flux_low2[1])[0],2))


poptgen_store2 = np.zeros((4,np.shape(flux2[1])[0]))

Name_source2 = Name[S_match2]
Name_source_list2 = []
bad_sources_name2 = []
error2 = np.zeros((4,np.shape(flux2[1])[0]))
PS_count2 = 0

#Here we do the actual fitting, to LoLSS LoTSS NVSS
start2 = time.time()
for i in tqdm(range(np.shape(flux_high2[1])[0])):
    # print('Processing '+ Name_source[i])
    # Fit NVSS to VLASS
    if spectral_index_eval_norms(freq_high2, flux_high2[:,i], flux_err_high2[:,i]) == 'Curve_fit could not fit powerlaw.':
        print('Could not fit high powerlaw to '+ Name_source2[i])
        continue
    else:
        poptpowlaw_high2, pcovpowlaw_high2, fluxplot_high2, freqplot_high2, flux_errplot_high2 = spectral_index_eval_norms(freq_high2,flux_high2[:,i],flux_err_high2[:,i])
        alpha_high2[i]=poptpowlaw_high2[1]
        norm_high2[i]=poptpowlaw_high2[0]

        # Calculate uncertainty on alpha_high, by looking at the boundary cases, plus and minus 1 sigma
        # find min and max lotss flux, min and max nvss flux
        NVSS_min = flux2[1,i] - flux_err2[1,i]
        NVSS_max = flux2[1,i] + flux_err2[1,i]
        VLASS_min = flux2[2,i] - flux_err2[2,i]
        VLASS_max = flux2[2,i] + flux_err2[2,i]

        # find corresponding a_high, a_high
        p0pow2 = [1400.,alpha_high2[i]]
        popt_err_high_l2, pcov_err_high_l2 = opt.curve_fit(gpscssmodels.powlaw, [1400., 3000.], [NVSS_min, VLASS_max], p0 = p0pow2, maxfev = 10000)
        try:
            popt_err_high_up2, pcov_err_high_up2 = opt.curve_fit(gpscssmodels.powlaw, [1400., 3000.], [NVSS_max, VLASS_min], p0 = p0pow2, maxfev = 100000)
        except:
            print('Could not fit error high powerlaw to '+ Name_source2[i])
            continue
        #TODO: find out why I need to do an except here, better initial values?

            
        alpha_err_high_b2[i] = [alpha_high2[i] - popt_err_high_l2[1], popt_err_high_up2[1] - alpha_high2[i]]
        alpha_err_high2[i] = np.abs(np.median(alpha_err_high_b2[i]))
        
    # Fit LoTSS to NVSS
    if spectral_index_eval_norms(freq_low2,flux_low2[:,i],flux_err_low2[:,i]) == 'Curve_fit could not fit powerlaw.':
        print('Could not fit low powerlaw to '+ Name_source2[i])
        continue
    else:
        poptpowlaw_low2, pcovpowlaw_low2, fluxplot_low2, freqplot_low2, flux_errplot_low2 = spectral_index_eval_norms(freq_low2,flux_low2[:,i],flux_err_low2[:,i])
        alpha_low2[i] = poptpowlaw_low2[1]
        norm_low2[i] = poptpowlaw_low2[0]

        # Calculate uncertainty on alpha_low
        # find min and max lotss flux, min and max nvss flux
        lotss_min = flux2[0,i] - flux_err2[0,i]
        lotss_max = flux2[0,i] + flux_err2[0,i]
        NVSS_min = flux2[1,i] - flux_err2[1,i]
        NVSS_max = flux2[1,i] + flux_err2[1,i]

        # find corresponding a_low, a_high
        p0pow2 = [1400.,alpha_low2[i]]
        popt_err_low_l2, pcov_err_low_l2 = opt.curve_fit(gpscssmodels.powlaw, [144., 1400.], [lotss_min, NVSS_max], p0 = p0pow2, maxfev = 10000)
        popt_err_low_up2, pcov_err_low_up2 = opt.curve_fit(gpscssmodels.powlaw, [144., 1400.], [lotss_max, NVSS_min], p0 = p0pow2, maxfev = 10000)
        alpha_err_low_b2[i] = [alpha_low2[i] - popt_err_low_l2[1], popt_err_low_up2[1] - alpha_low2[i]]
        alpha_err_low2[i] = np.abs(np.median(alpha_err_low_b2[i]))

    if -alpha_low2[i] >= 0. and -alpha_high2[i] <= 0.:
        PS_count2 += 1
        #TODO: fix these values

    
print(PS_count2, " PS sources for LoTSS NVSS VLASS")
#IMPORTANT, from here on out we use alpha high and low as defined in paper by Slob (positive)

alpha_low_array2 = np.zeros(len(S_LoTSS))
alpha_low_array2[S_match2] = -alpha_low2

alpha_err_low_array2 = np.zeros(len(S_LoTSS))
alpha_err_low_array2[S_match2] = alpha_err_low2

alpha_high_array2=np.zeros(len(S_LoTSS))
alpha_high_array2[S_match2] = -alpha_high2

alpha_err_high_array2 = np.zeros(len(S_LoTSS))
alpha_err_high_array2[S_match2] = alpha_err_high2

norm_low_array2=np.zeros(len(S_LoTSS))
norm_low_array2[S_match2]= norm_low2

norm_high_array2=np.zeros(len(S_LoTSS))
norm_high_array2[S_match2]= norm_high2



col1 = fits.Column(name='LoTSS_name', format = '34A', array = Name)
col2 = fits.Column(name='RA', format = 'E', array = RA)
col3 = fits.Column(name='Dec', format = 'E', array = DEC)
col4 = fits.Column(name='LoTSS_flux', format = 'E', array = S_LoTSS)
col5 = fits.Column(name='e_LoTSS_flux', format = 'E', array = e_S_LoTSS)

col6 = fits.Column(name='NVSS_RA', format = 'E', array = tbdata['RAJ2000_1'])
col7 = fits.Column(name='NVSS_Dec', format = 'E', array = tbdata['DEJ2000_1'])
col8 = fits.Column(name='NVSS_flux', format = 'E', array = S_NVSS)
col9 = fits.Column(name='e_NVSS_flux', format = 'E', array = e_S_NVSS)

col10 = fits.Column(name='VLASS_RA', format = 'E', array = tbdata['ra_2'])
col11 = fits.Column(name='VLASS_Dec', format = 'E', array = tbdata['dec_2'])
col12 = fits.Column(name='VLASS_flux', format = 'E', array = S_VLASS)
col13 = fits.Column(name='e_VLASS_flux', format = 'E', array = e_S_VLASS)
col14 = fits.Column(name='VLASS_index', format = '8A', array = tbdata['index'])
col43 = fits.Column(name='VLASS_flags', format = '8A', array = tbdata['flags'])

col15 = fits.Column(name='TGSS_RA', format = 'E', array = tbdata['RAJ2000_2'])
col16 = fits.Column(name='TGSS_Dec', format = 'E', array = tbdata['DEJ2000_2'])
col17 = fits.Column(name='TGSS_flux', format = 'E', array = S_TGSS)
col18 = fits.Column(name='e_TGSS_flux', format = 'E', array = e_S_TGSS)

col19 = fits.Column(name='VLSSr_RA', format = 'E', array = tbdata['RAJ2000_3'])
col20 = fits.Column(name='VLSSr_Dec', format = 'E', array = tbdata['DEJ2000_3'])
col21 = fits.Column(name='VLSSr_flux', format = 'E', array = S_VLSSr)
col22 = fits.Column(name='e_VLSSr_flux', format = 'E', array = e_S_VLSSr)

col23 = fits.Column(name='LoLSS_RA', format = 'E', array = tbdata['RA_4'])
col24 = fits.Column(name='LoLSS_Dec', format = 'E', array = tbdata['DEC_4'])
col25 = fits.Column(name='LoLSS_flux', format = 'E', array = S_LoLSS)
col26 = fits.Column(name='e_LoLSS_flux', format = 'E', array = e_S_LoLSS)

col27 = fits.Column(name='FIRST_RA', format = 'E', array = tbdata['RA_5'])
col28 = fits.Column(name='FIRST_Dec', format = 'E', array = tbdata['DEC_5'])
col29 = fits.Column(name='FIRST_flux', format = 'E', array = S_FIRST)
col30 = fits.Column(name='e_FIRST_flux', format = 'E', array = e_S_FIRST) 

col31 = fits.Column(name='inband_RA', format = 'E', array = tbdata['RA_6'])
col32 = fits.Column(name='inband_Dec', format = 'E', array = tbdata['DEC_6'])
col33 = fits.Column(name='S_inband_low', format = 'E', array = S_inband_low)
col34 = fits.Column(name='e_S_inband_low', format = 'E', array = e_S_inband_low)
col35 = fits.Column(name='S_inband_mid', format = 'E', array = S_inband_mid)
col36 = fits.Column(name='e_S_inband_mid', format = 'E', array = e_S_inband_mid)
col37 = fits.Column(name='S_inband_high', format = 'E', array = S_inband_high)
col38 = fits.Column(name='e_S_inband_high', format = 'E', array = e_S_inband_high)

col39 = fits.Column(name='channel0_RA', format = 'E', array= tbdata['RA_7'])
col40 = fits.Column(name='channel0_Dec', format = 'E', array = tbdata['DEC_7'])
col41 = fits.Column(name='channel0_flux', format = 'E', array = S_channel0) #Jansky
col42 = fits.Column(name='e_channel0_flux', format = 'E', array = e_S_channel0)

col44 = fits.Column(name='channel1_RA', format = 'E', array= tbdata['RA_8'])
col45 = fits.Column(name='channel1_Dec', format = 'E', array = tbdata['DEC_8'])
col46 = fits.Column(name='channel1_flux', format = 'E', array = S_channel1)
col47 = fits.Column(name='e_channel1_flux', format = 'E', array = e_S_channel1)

col48 = fits.Column(name='channel2_RA', format = 'E', array= tbdata['RA_9'])
col49 = fits.Column(name='channel2_Dec', format = 'E', array = tbdata['DEC_9'])
col50 = fits.Column(name='channel2_flux', format = 'E', array = S_channel2)
col51 = fits.Column(name='e_channel2_flux', format = 'E', array = e_S_channel2)

col52 = fits.Column(name='channel3_RA', format = 'E', array= tbdata['RA_10'])
col53 = fits.Column(name='channel3_Dec', format = 'E', array = tbdata['DEC_10'])
col54 = fits.Column(name='channel3_flux', format = 'E', array = S_channel3)
col55 = fits.Column(name='e_channel3_flux', format = 'E', array = e_S_channel3)

col56 = fits.Column(name='channel4_RA', format = 'E', array= tbdata['RA_11'])
col57 = fits.Column(name='channel4_Dec', format = 'E', array = tbdata['DEC_11'])
col58 = fits.Column(name='channel4_flux', format = 'E', array = S_channel4)
col59 = fits.Column(name='e_channel4_flux', format = 'E', array = e_S_channel4)

col60 = fits.Column(name='channel5_RA', format = 'E', array= tbdata['RA_12'])
col61 = fits.Column(name='channel5_Dec', format = 'E', array = tbdata['DEC_12'])
col62 = fits.Column(name='channel5_flux', format = 'E', array = S_channel5)
col63 = fits.Column(name='e_channel5_flux', format = 'E', array = e_S_channel5)

col64 = fits.Column(name='WENSS_RA', format = 'E', array= tbdata['_RAJ2000'])
col65 = fits.Column(name='WENSS_DEC', format = 'E', array = tbdata['_DEJ2000'])
col66 = fits.Column(name='WENSS_intflux', format = 'E', array = WENSS_intflux)
col67 = fits.Column(name='WENSS_intflux_error', format = 'E', array = WENSS_intflux_error)

#LoLSS, LoTSS, NVSS
col68 = fits.Column(name='a_low_LoLSS', format = 'E', array = norm_low_array)
col69 = fits.Column(name='alpha_low_LoLSS', format = 'E', array = alpha_low_array)
col70 = fits.Column(name='e_alpha_low_LoLSS', format = 'E', array = alpha_err_low_array)
col71 = fits.Column(name='a_high_LoLSS', format = 'E', array = norm_high_array)
col72 = fits.Column(name='alpha_high_LoLSS', format = 'E', array = alpha_high_array)
col73 = fits.Column(name='e_alpha_high_LoLSS', format = 'E', array = alpha_err_high_array)

#LoTSS, NVSS, VLASS
col74 = fits.Column(name='a_low_VLASS', format = 'E', array = norm_low_array2)
col75 = fits.Column(name='alpha_low_VLASS', format = 'E', array = alpha_low_array2)
col76 = fits.Column(name='e_alpha_low_VLASS', format = 'E', array = alpha_err_low_array2)
col77 = fits.Column(name='a_high_VLASS', format = 'E', array = norm_high_array2)
col78 = fits.Column(name='alpha_high_VLASS', format = 'E', array = alpha_high_array2)
col79 = fits.Column(name='e_alpha_high_VLASS', format = 'E', array = alpha_err_high_array2)

cols = fits.ColDefs([col68, col69, col70, col71,col72, col73, col74, col75, col76, col77, col78, col79, col1, col2, col3, col4, col5, col6, col7,\
                     col8, col9, col10, col11, col12, col13, col14, col43, col15, col16, col17, col18, col19,\
    col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, \
        col40, col41, col42,col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55\
            ,col56, col57, col58, col59, col60, col61, col62, col63, col64, col65, col66, col67])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/mega_master_10000_clean.fits', overwrite = True)