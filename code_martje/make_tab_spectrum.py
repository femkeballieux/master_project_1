import numpy as np
from astropy.io import fits
import gpscssmodels
import seds_plot_func_extreme
from tqdm import tqdm
import scipy.optimize as opt
import warnings
import time
warnings.filterwarnings("ignore")

#Run on vdesk

# fit power-law without rms cuts
def spectral_index_eval_norms(freq,flux,flux_err):
    fluxposind = np.where((flux > 0))
    if len(fluxposind[0]) <= 1: # some sources are all as not detected in subbands? 
        print('Curve_fit could not fit powerlaw.')
        return 
    
    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0pow = [peakflux,0.7]

    try:
        poptpowlaw, pcovpowlaw = opt.curve_fit(gpscssmodels.powlaw, freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError): # If curve fit can not find a good fit, skip the source.    
        return 'Curve_fit could not fit powerlaw.'

    return (poptpowlaw, pcovpowlaw, fluxplot, freqplot, flux_errplot)

# fit curvature models
def spectral_index_eval_curve(freq,flux,flux_err):
    fluxposind = np.where((flux > 0)) # need to concatenate and insert to make sure you still select TGSS, MRC, SUMSS and NVSS data
    if len(fluxposind[0]) == 0: # some sources are all as not detected in subbands? 
        return 'No flux?'

    fluxplot = flux[fluxposind]
    freqplot = freq[fluxposind]
    flux_errplot = flux_err[fluxposind]

    peakfreq = freqplot[np.where(fluxplot == max(fluxplot))]
    peakflux = max(fluxplot)
    p0gen = [peakflux,300.,0.8,-0.7]

    try:
        poptgen, pcovgen = opt.curve_fit(gpscssmodels.curve, freqplot, fluxplot, p0 = p0gen, sigma = flux_errplot, maxfev = 10000)
    except (RuntimeError, TypeError, ValueError):
        return 'Curve_fit could not fit curve model.'
    return (fluxplot, freqplot, flux_errplot, poptgen, pcovgen)
    

def make_fluxes(catalogue, e_catalogue, freq_catalogue):
    S_match = np.where(((catalogue > 0.) & (S_LoLSS > 0.)))

    flux = np.vstack((S_LoTSS, S_LoLSS, S_VLSSr, S_TGSS, S_NVSS))
    flux = flux[:,S_match[0]]
    flux_err = np.vstack((e_S_LoTSS, e_S_LoLSS, e_S_VLSSr, e_S_TGSS, e_S_NVSS))
    flux_err = flux_err[:,S_match[0]] 

    freq = np.array([freq_LoTSS, freq_LoLSS, freq_VLSSr, freq_TGSS, freq_NVSS])
    freq_high = np.array([freq_LoTSS, freq_catalogue])
    freq_low = np.array([freq_LoLSS, freq_LoTSS])
    # freq = np.array(freq_LoLSS, freq_LoTSS, freq_catalogue)

    flux_high = np.vstack((S_LoTSS, catalogue))
    flux_high = flux_high[:,S_match[0]]
    flux_low = np.vstack((S_LoLSS, S_LoTSS))
    flux_low = flux_low[:,S_match[0]]

    flux_err_high = np.vstack((e_S_LoTSS, e_catalogue))
    flux_err_high = flux_err_high[:,S_match[0]]
    flux_err_low = np.vstack((e_S_LoLSS, e_S_LoTSS))
    flux_err_low = flux_err_low[:,S_match[0]]
    
    return freq, freq_high, freq_low, flux, flux_err, flux_high, flux_low, flux_err_high, flux_err_low, S_match

hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_NVSS_TGSS_VLSSr_LoLSS_inband.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()# fit curvature models with rms cuts

S_LoTSS, stat_e_S_LoTSS = np.array(tbdata['Peak_flux'])/1000., np.array(tbdata['E_Peak_flux'])/1000. # Jy/beam
e_S_LoTSS = np.sqrt(stat_e_S_LoTSS**2 + (0.1*S_LoTSS)**2) # Combine the statistical LoTSS error and overall error on LoTSS

S_NVSS, e_S_NVSS = np.array(tbdata['S1_4'])/1000., 0.1 * np.array(tbdata['e_S1_4'])/1000. # Jy
S_NVSS, e_S_NVSS = np.where(np.isnan(S_NVSS), 0, S_NVSS), np.where(np.isnan(S_NVSS), 0, e_S_NVSS)

S_TGSS, e_S_TGSS = np.array(tbdata['Stotal'])/1000., np.array(tbdata['e_Stotal'])/1000. # Jy/beam
S_TGSS, e_S_TGSS = np.where(np.isnan(S_TGSS), 0, S_TGSS), np.where(np.isnan(S_TGSS), 0, e_S_TGSS)

S_VLSSr, e_S_VLSSr = np.array(tbdata['Sp']), 0.1 * np.array(tbdata['Sp']) # Jy
S_VLSSr, e_S_VLSSr = np.where(np.isnan(S_VLSSr), 0, S_VLSSr), np.where(np.isnan(S_VLSSr), 0, e_S_VLSSr)

S_LoLSS, e_S_LoLSS = np.array(tbdata['Total_flux_LoLSS'])/1000., np.array(tbdata['E_Total_flux_LoLSS'])/1000. # Jy/beam
S_LoLSS, e_S_LoLSS = np.where(np.isnan(S_LoLSS), 0, S_LoLSS), np.where(np.isnan(S_LoLSS), 0, e_S_LoLSS)

S_FIRST, e_S_FIRST = np.array(tbdata['FPEAK'])/1000, 0.1*(np.array(tbdata['FPEAK'])/1000) #Jy, set error to 10% of flux
S_FIRST, e_S_FIRST = np.where(np.isnan(S_FIRST), 0, S_FIRST), np.where(np.isnan(S_FIRST), 0, e_S_FIRST)

S_inband_low, e_S_inband_low = np.array(tbdata['M1_L1flux']), 0.1 * np.array(tbdata['M1_L1flux']) #mJy
S_inband_low, e_S_inband_low = np.where(np.isnan(S_inband_low), 0, S_inband_low), np.where(np.isnan(S_inband_low), 0, e_S_inband_low)

S_inband_mid, e_S_inband_mid = np.array(tbdata['M1_L2flux']), 0.1 * np.array(tbdata['M1_L2flux']) #mJy
S_inband_mid, e_S_inband_mid = np.where(np.isnan(S_inband_mid), 0, S_inband_mid), np.where(np.isnan(S_inband_mid), 0, e_S_inband_mid)

S_inband_high, e_S_inband_high = np.array(tbdata['M1_L3flux']), 0.1 * np.array(tbdata['M1_L3flux']) #mJy
S_inband_high, e_S_inband_high = np.where(np.isnan(S_inband_high), 0, S_inband_high), np.where(np.isnan(S_inband_high), 0, e_S_inband_high)

freq_LoTSS, freq_NVSS, freq_TGSS, freq_VLSSr, freq_LoLSS, freq_FIRST, freq_inband_low, freq_inband_mid, freq_inband_high = 144., 1400., 150., \
                                                                                                                            74., 54., 1400.001,\
                                                                                                                            128., 144.00001, 160. # MHz

Name = tbdata['Source'].rstrip()
RA = tbdata['RAJ2000_1']
DEC = tbdata['DEJ2000_1']

NVSS = True

freq, freq_high, freq_low, flux, flux_err, flux_high, flux_low, flux_err_high, flux_err_low, S_match = make_fluxes(S_NVSS, e_S_NVSS, freq_NVSS)

alpha_low, alpha_err_low, alpha_high, alpha_err_high, alpha_tot, alpha_err_tot, peak_flux_man, peak_freq_man, alphathick_err,\
     alphathin_err, freq_peak_gen_err, S_max_gen_err, norm_high, curve_fit_qual, norm_low = [np.zeros(np.shape(flux_low[1])[0]) for dummy in range(15)]

alpha_err_high_b = np.zeros((np.shape(flux_low[1])[0],2))
alpha_err_low_b = np.zeros((np.shape(flux_low[1])[0],2))

freq_cont = np.linspace(1,6000,20000)

poptgen_store = np.zeros((4,np.shape(flux[1])[0]))

Name_source = Name[S_match]
Name_source_list = []
bad_sources_name = []
error = np.zeros((4,np.shape(flux[1])[0]))
PS_count = 0

start = time.time()
for i in tqdm(range(np.shape(flux_high[1])[0])):
    # print('Processing '+ Name_source[i])
    # Fit LoTSS to high source (NVSS or FIRST)
    if spectral_index_eval_norms(freq_high,flux_high[:,i],flux_err_high[:,i]) == 'Curve_fit could not fit powerlaw.':
        print('Could not fit high powerlaw to '+ Name_source[i])
        continue
    else:
        poptpowlaw_high, pcovpowlaw_high, fluxplot_high, freqplot_high, flux_errplot_high = spectral_index_eval_norms(freq_high,flux_high[:,i],flux_err_high[:,i])
        alpha_high[i]=poptpowlaw_high[1]
        norm_high[i]=poptpowlaw_high[0]

        # Calculate uncertainty on alpha_high
        # find min and max lotss flux, min and max nvss flux
        lotss_min = flux[0,i] - flux_err[0,i]
        lotss_max = flux[0,i] + flux_err[0,i]
        nvss_min = flux[4,i] - flux_err[4,i]
        nvss_max = flux[4,i] + flux_err[4,i]

        # find corresponding a_low, a_high
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
        # fitting curved model to entire data apart from inband spectra
        if spectral_index_eval_curve(freq,flux[:,i],flux_err[:,i]) == 'Curve_fit could not fit curve model.' or spectral_index_eval_curve(freq,flux[:,i],flux_err[:,i]) == 'No flux?':
            poptcurve = np.array([np.nan,np.nan,np.nan,np.nan])
            poptgen = np.array([np.nan, np.nan, np.nan, np.nan])
        else:    
            fluxplot, freqplot, flux_errplot, poptgen, pcovgen = spectral_index_eval_curve(freq,flux[:,i],flux_err[:,i])
    else:
        poptcurve = np.array([np.nan,np.nan,np.nan,np.nan])
        poptgen = np.array([np.nan, np.nan, np.nan, np.nan])


    if np.isnan(poptgen[0]):
        poptgen_store[:,i] = np.nan
        peak_flux_man[i] = np.nan
        peak_freq_man[i] = np.nan
        freq_peak_gen_err[i] = np.nan
        S_max_gen_err[i] = np.nan
        alphathick_err[i] = np.nan
        alphathin_err[i] = np.nan
        curve_fit_qual[i] = 0
    else:
        poptgen_store[:,i] = poptgen
        S_max_gen_err[i] = np.sqrt(pcovgen[0][0])
        freq_peak_gen_err[i] = np.sqrt(pcovgen[1][1])
        alphathin_err[i] = np.sqrt(pcovgen[2][2])
        alphathick_err[i] = np.sqrt(pcovgen[3][3])
        error[:,i] = [S_max_gen_err[i], freq_peak_gen_err[i], alphathin_err[i], alphathick_err[i]]
        peak_flux_man[i] = max(gpscssmodels.curve(freq_cont,*poptgen_store[:,i]))
        try:
            peak_freq_man[i] = freq_cont[np.where((gpscssmodels.curve(freq_cont,*poptgen_store[:,i]) == peak_flux_man[i]))]
            curve_fit_qual[i] = 1
        except(ValueError):
            peak_freq_man[i] = np.nan
            bad_sources_name.append(Name_source[i])
            curve_fit_qual[i] = 0
        
    # making sure S_TGSS is plotted 
    sedplot_freq = np.concatenate((freqplot_low,freqplot_high))
    sedplot_flux = np.concatenate((fluxplot_low,fluxplot_high))
    sedplot_flux_err = np.concatenate((flux_errplot_low,flux_errplot_high))
    
    if S_TGSS[S_match][i] > 0.:
        sedplot_freq = np.insert(sedplot_freq,1,freq_TGSS)
        sedplot_flux = np.insert(sedplot_flux,1,S_TGSS[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_TGSS[S_match][i])
    if S_VLSSr[S_match][i] > 0.:
        sedplot_freq = np.insert(sedplot_freq,1,freq_VLSSr)
        sedplot_flux = np.insert(sedplot_flux,1,S_VLSSr[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_VLSSr[S_match][i])
    if NVSS == True and S_FIRST[S_match][i] > 0.:
        sedplot_freq = np.insert(sedplot_freq,1,freq_FIRST)
        sedplot_flux = np.insert(sedplot_flux,1,S_FIRST[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_FIRST[S_match][i])
    if NVSS == False and S_NVSS[S_match][i] > 0.:
        sedplot_freq = np.insert(sedplot_freq,1,freq_NVSS)
        sedplot_flux = np.insert(sedplot_flux,1,S_NVSS[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_NVSS[S_match][i])
    if S_inband_low[S_match][i] > 0.:
        sedplot_freq = np.insert(sedplot_freq,1,freq_inband_low)
        sedplot_flux = np.insert(sedplot_flux,1,S_inband_low[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_inband_low[S_match][i])
        sedplot_freq = np.insert(sedplot_freq,1,freq_inband_mid)
        sedplot_flux = np.insert(sedplot_flux,1,S_inband_mid[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_inband_mid[S_match][i])
        sedplot_freq = np.insert(sedplot_freq,1,freq_inband_high)
        sedplot_flux = np.insert(sedplot_flux,1,S_inband_high[S_match][i])
        sedplot_flux_err = np.insert(sedplot_flux_err,1,e_S_inband_high[S_match][i])
    
    if -alpha_low[i] >= 0. and -alpha_high[i] >= 0.:
        Name_source_list.append('Q1/'+Name_source[i])
    elif -alpha_low[i] >= 0.1 and -alpha_high[i] <= 0.:
        Name_source_list.append('PS/'+Name_source[i])
    elif -alpha_low[i] >= 0. and -alpha_high[i] <= 0.:
        Name_source_list.append('Q2/'+Name_source[i])
    elif -alpha_low[i] <= 0. and -alpha_high[i] <= 0.:
        Name_source_list.append('Q3/'+Name_source[i])
        curve_fit_qual[i] = 0
    elif -alpha_low[i] <= 0. and -alpha_high[i] >= 0.:
        Name_source_list.append('Q4/'+Name_source[i])

    Name_source_source = 'NVSS/'+Name_source_list[i]+'_NVSS_total'
    
    seds_plot_func_extreme.sed([gpscssmodels.powlaw,gpscssmodels.powlaw, gpscssmodels.curve],[poptpowlaw_high,poptpowlaw_low, poptgen],\
        sedplot_freq,sedplot_flux,sedplot_flux_err, Name_source_source, Name_source[i], freq_labels = True, savefig = True, resid = False, error = error[:,i])


print(PS_count)

alpha_low = -alpha_low
alpha_high = -alpha_high

col1 = fits.Column(name='LoTSS_name', format = '34A', array = Name[S_match])
col2 = fits.Column(name='RA', format = 'E', array = RA[S_match])
col3 = fits.Column(name='Dec', format = 'E', array = DEC[S_match])
col4 = fits.Column(name='LoTSS_flux', format = 'E', array = S_LoTSS[S_match])
col5 = fits.Column(name='e_LoTSS_flux', format = 'E', array = e_S_LoTSS[S_match])
col6 = fits.Column(name='a_low', format = 'E', array = norm_low)
col7 = fits.Column(name='alpha_low', format = 'E', array = alpha_low)
col8 = fits.Column(name='e_alpha_low', format = 'E', array = alpha_err_low)
col9 = fits.Column(name='a_high', format = 'E', array = norm_high)
col10 = fits.Column(name='alpha_high', format = 'E', array = alpha_high)
col11 = fits.Column(name='e_alpha_high', format = 'E', array = alpha_err_high)
col12 = fits.Column(name='LoLSS_RA', format = 'E', array = tbdata['RA_5'][S_match])
col13 = fits.Column(name='LoLSS_Dec', format = 'E', array = tbdata['DEC_5'][S_match])
col14 = fits.Column(name='LoLSS_flux', format = 'E', array = S_LoLSS[S_match])
col15 = fits.Column(name='e_LoLSS_flux', format = 'E', array = e_S_LoLSS[S_match])
col16 = fits.Column(name='NVSS_RA', format = 'E', array = tbdata['RAJ2000_2'][S_match])
col17 = fits.Column(name='NVSS_Dec', format = 'E', array = tbdata['DEJ2000_2'][S_match])
col18 = fits.Column(name='NVSS_flux', format = 'E', array = S_NVSS[S_match])
col19 = fits.Column(name='e_NVSS_flux', format = 'E', array = e_S_NVSS[S_match])
col20 = fits.Column(name='TGSS_RA', format = 'E', array = tbdata['RAJ2000_3'][S_match])
col21 = fits.Column(name='TGSS_Dec', format = 'E', array = tbdata['DEJ2000_3'][S_match])
col22 = fits.Column(name='TGSS_flux', format = 'E', array = S_TGSS[S_match])
col23 = fits.Column(name='e_TGSS_flux', format = 'E', array = e_S_TGSS[S_match])
col24 = fits.Column(name='VLSSr_RA', format = 'E', array = tbdata['RAJ2000_4'][S_match])
col25 = fits.Column(name='VLSSr_Dec', format = 'E', array = tbdata['DEJ2000_4'][S_match])
col26 = fits.Column(name='VLSSr_flux', format = 'E', array = S_VLSSr[S_match])
col27 = fits.Column(name='e_VLSSr_flux', format = 'E', array = e_S_VLSSr[S_match])
col28 = fits.Column(name='FIRST_RA', format = 'E', array = tbdata['RA_6'][S_match])
col29 = fits.Column(name='FIRST_Dec', format = 'E', array = tbdata['DEC_6'][S_match])
col30 = fits.Column(name='FIRST_flux', format = 'E', array = S_FIRST[S_match])
col31 = fits.Column(name='e_FIRST_flux', format = 'E', array = e_S_FIRST[S_match]) 
col32 = fits.Column(name='inband_RA', format = 'E', array = tbdata['RA_7'][S_match])
col33 = fits.Column(name='inband_Dec', format = 'E', array = tbdata['DEC_7'][S_match])
col34 = fits.Column(name='S_inband_low', format = 'E', array = S_inband_low[S_match])
col35 = fits.Column(name='e_S_inband_low', format = 'E', array = e_S_inband_low[S_match])
col36 = fits.Column(name='S_inband_mid', format = 'E', array = S_inband_mid[S_match])
col37 = fits.Column(name='e_S_inband_mid', format = 'E', array = e_S_inband_mid[S_match])
col38 = fits.Column(name='S_inband_high', format = 'E', array = S_inband_high[S_match])
col39 = fits.Column(name='e_S_inband_high', format = 'E', array = e_S_inband_high[S_match])

cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19,\
    col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/fit_vals_power_law_NVSS_intflux.fits', overwrite = True)


