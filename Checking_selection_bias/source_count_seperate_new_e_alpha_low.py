# !/bin/python
# Script for computing PS source counts

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
from matplotlib import rc
from astropy import units as u
from matplotlib.ticker import ScalarFormatter
from scipy import stats
from astropy.io.votable import parse_single_table
from scipy.optimize import curve_fit
import pandas as pd
from scipy.odr import *


rc('text', usetex=False)
rc('font', family='sans-serif')
plt.style.use(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/code/style.mplstyle')

def poisson_errors(n):
    """
    Compute log poisson errors
    """
    n_err = np.zeros((2, np.shape(n)[0]))

    for i, n_i in enumerate(n):
        if n_i > 0:
            # Use 84% confidence upper and lower limits based on Poisson statistics, tabulated in Gehrels (1986)
            m = n_i+1
            lambda_u = m*(1-1/(9*m)+1/(3*np.sqrt(m)))**3
            lambda_l = n_i*(1-1/(9*n_i)-1/(3*np.sqrt(n_i)))**3
            err_high = n_i*lambda_u/n_i - n_i
            err_low = n_i - n_i*lambda_l/n_i
            n_err[:, i] = [err_low, err_high]

    return n_err


def source_counts(flux, bin_low, bin_high, num_bins, area_str):
    """
    Calculating number of sources in flux bin, rescaling with S^2.5
    Flux is the flux array entered, 
    bin_low corresponds to the extrapolated completeness limit
    bin_high corresponds to the highest flux in the sample
    num_bins is number of bins. How many depends on what looks good, you want 
    order of magnitude same amount of sources in each bin
    """
    bins = np.logspace(np.log10(bin_low), np.log10(
        bin_high), num=num_bins)[:-2]  # log bins
    bins = np.concatenate((bins, [(bin_high)]))
    n, bins = np.histogram(flux, bins=bins)
    print('# of sources per bin:', n, '\n Total sources sampled:', np.sum(n))
    n_err = poisson_errors(n)

    # dividing by area of the survey.
    n_area = n / area_str
    n_err_area = n_err / area_str

    # dividing by bin width (in log space)
    bins_width = np.diff(bins)
    n_area_bin_width_norm = n_area / bins_width
    n_err_area_bin_width_norm = n_err_area / bins_width

    # Normalising by Eucledian counts
    # centers of bins
    mean_flux_per_bin = stats.binned_statistic(
        flux, flux, statistic='mean', bins=bins)[0]
    dn_ds = n_area_bin_width_norm * mean_flux_per_bin**2.5
    dn_ds_err = n_err_area_bin_width_norm * mean_flux_per_bin**2.5

    # For plotting, first have to get central bin value
    central_bins_PS = np.array([(j+i)/2. for i, j in zip(bins[:-1], bins[1:])])

    for i, count in enumerate(dn_ds):
        print(np.round(central_bins_PS[i], 3), "& $", np.round(
            count, 2), "^{+", np.round(dn_ds_err[1, i], 2), "}_{-", np.round(dn_ds_err[0, i], 2), "}$")

    return (dn_ds, dn_ds_err, central_bins_PS, bins_width, n)


def scale_to_freq(dn_ds, dn_ds_err, bins, bin_width, freq_in, alpha_thick, freq_out=144, printing=True):
    """
    Scaling dn/ds and bins to another frequency, freq_out. S propto nu^alpha. We only do this for 
    Snellen sample at 2GHz.
    alpha_thick is the scaling factor used to shift from one sample to another
    """
    factor = (freq_out / freq_in) ** alpha_thick
    dn_ds_freq = dn_ds * factor**1.5
    dn_ds_err_freq = dn_ds_err * factor**1.5
    bins_freq = bins * factor
    bin_width_freq = bin_width * factor
    if printing:
        print("")
        for i, count in enumerate(dn_ds_freq):
            print(np.round(bins_freq[i], 3), "& & $", np.round(dn_ds_freq[i], 2), "^{+", np.round(
                dn_ds_err_freq[1, i], 2), "}_{-", np.round(dn_ds_err_freq[0, i], 2), "}$")

    return dn_ds_freq, dn_ds_err_freq, bins_freq, bin_width_freq


def hist(LoTSS_flux, sample, freq=str(144)):
    """
    Histogram of freqMHz flux distribution for PS sources
    """
    gs = plt.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    bins = np.logspace(np.log10(0.1), np.log10(1000), 50)
    hist_PS_flux, bin_edges_PS_flux = np.histogram(
        LoTSS_flux*1000, bins=bins)  # mJy
    ax.stairs(hist_PS_flux, bin_edges_PS_flux,
              edgecolor='crimson', lw=2, hatch='//')
    ax.set_xlabel('S$_{'+freq+'}$ (mJy)')
    ax.set_ylabel('Number of PS sources')
    ax.set_xscale('log')
    ax.set_title('Flux hist for ' + sample)
    ax.xaxis.set_major_formatter(ScalarFormatter())

    plt.savefig(
        "/net/vdesk/data2/bach1/ballieux/master_project_1/plots/peak_flux_dist_new" + sample + ".png")
    plt.show()


def source_counts_plot(sample):
    """
    Here we make the 144MHz source count plot
    """
    
    fig = plt.figure(2, figsize=(10, 16))
    gs = plt.GridSpec(3, 1, height_ratios=[2, 1, 1]) 
    ax = plt.subplot(gs[0])
    
    # Values from Mandal 2021 (LoTSS Deep Field) (At 144MHz)
    if sample == 'MPS':
        S_deep_field = np.array([0.22, 0.31, 0.43, 0.61, 0.86, 1.22, 1.73, 2.45,
                                3.46, 4.89, 6.92, 9.79, 16.5, 32.9, 65.8, 132., 372., 1489.])/1000
        N_final_deep_field = np.array([37.15, 42.65, 47.08, 50.55, 47.97, 43.17, 42.77,
                                      42.30, 49.81, 54.21, 79.79, 109.6, 145.9, 341.5, 607.1, 1041, 1739, 4921])
    
        S_deep_field, N_final_deep_field, N_deep_field_up_err, N_deep_field_low_err = np.loadtxt(
            '/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/Lofar_deep_fields.csv', dtype=float, delimiter=',', unpack=True, skiprows=1, usecols=[0, 1, 2, 3])
        S_deep_field /= 1000
        N_deep_field_err = np.array([N_deep_field_low_err, N_deep_field_up_err])
        ax.errorbar(S_deep_field, N_final_deep_field, yerr=N_deep_field_err, zorder=1,
                    marker="^", label='LoTSS Deep Field \n\
Mandal et al. (2021) AGN',  linestyle='none', capthick=2, color='skyblue')

    # Values from Snellen 1998, which were read off from the plot. The error propagation is weird here, do not do that anywhere else.
    #Rescaled to 1400 for GPS, it was originally at 2000MHz
    if sample == 'GPS':
        dn_ds_snellen = np.array([12.73, 5.589, 2.105])
        dn_ds_snellen_err = np.array([[7.593, 17.346], [2.927, 8.143], [1.331, 2.817]])
        bins_snellen = np.array([300., 150., 75.])/1000.
    
        upp_bin_snellen = np.array([400, 200, 100])/1000
        low_bin_snellen = np.array([200, 100, 50])/1000
        bin_width_snellen = upp_bin_snellen - low_bin_snellen
    
        alpha_thick_snellen = -0.8
        dn_ds_1400_snellen, dn_ds_err_1400_snellen, bins_1400_snellen, bin_width_1400_snellen = \
            scale_to_freq(dn_ds_snellen, dn_ds_snellen_err, bins_snellen, bin_width_snellen,
                         2000., alpha_thick_snellen, freq_out=1400, printing=False)
        dn_ds_err_1400_snellen = np.abs(np.transpose(
            dn_ds_err_1400_snellen) - dn_ds_1400_snellen)
        
        ax.errorbar(bins_1400_snellen, dn_ds_1400_snellen, yerr=(dn_ds_err_1400_snellen), zorder=1,
                    xerr=bin_width_1400_snellen/2, marker='o', linestyle='none', label='Snellen et al. (1998) GPS', capthick=2, color='red')

    # de Zotti/Masardi model (read in at 150MHz). This file came from Joe.
    dZm_150 = parse_single_table(
        "/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/de_Zotti_model.vot").array #150MHz 
    if sample == 'MPS': 
        scaling = 44
        ax.plot(10**(dZm_150["S_mid"]), 10**(dZm_150["scount"]),
                color='#999999', label="Massardi (2010) 150MHz model")
        ax.plot(10**(dZm_150["S_mid"]), 10**(dZm_150["scount"])/scaling, ls='--',
                color='#999999', label="Massardi (2010)  150MHz model * 1/44")
        

    #This file came from a digitized webplot, therefore I smoothe it out by fitting to it
    dZm_1400 = parse_single_table(
        "/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/Masardi_1400.vot").array #1400MHz
    polycoeffs = np.polyfit(dZm_1400["S_mid"],dZm_1400["scount"],10)
    fitmodel = np.poly1d(polycoeffs)
    if sample == 'GPS': 
        scaling = 29
        ax.plot((10**(dZm_1400["S_mid"])), 10**fitmodel((dZm_1400["S_mid"])),
                color='#999999', label="Massardi (2010) 1400MHz model")
        ax.plot((10**(dZm_1400["S_mid"])), (10**fitmodel(dZm_1400["S_mid"]))/scaling, ls='--',
                color='#999999', label="Massardi (2010) 1400MHz model * 1/29")

    # Callingham Source Counts (143 MHz)
    if sample == 'MPS':
        dn_ds_call, dn_ds_err_call_low, dn_ds_corr_call_up, bins_call, bin_width_call, n_call = \
            np.vsplit(np.load(
                '/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/call_source_counts.npy'), 6)
        dn_ds_err_call = np.vstack((dn_ds_err_call_low, dn_ds_corr_call_up))
        ax.errorbar(bins_call[0, :], dn_ds_call[0, :], yerr=dn_ds_err_call[0, :], zorder=1,
                    xerr=bin_width_call[0, :]/2., marker='D',  linestyle='none', label='Callingham et al. (2017)', capthick=2, color='red')

    # Slob sample
    if sample == 'MPS':
        data_slob = (np.loadtxt(
            '/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/Slob_source_counts.txt'))
        bins_slob = data_slob[:, 0]
        dn_ds_slob = data_slob[:, 1]
        dn_ds_err_slob = data_slob[:, 2:4]
        bins_err_slob = data_slob[:, 4]
        ax.errorbar(bins_slob, dn_ds_slob, xerr=bins_err_slob, yerr=dn_ds_err_slob.T, zorder=1,
                    marker='v',  linestyle='none', label='Slob et al. (2022)', capthick=2, capsize=3,  color='blue', markersize='6')
        
    if sample == 'GPS':
        data_heywood = pd.read_csv('/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/heywood.csv')
        bins_heywood = np.array(data_heywood['Bin_Mid'])/1000
        bins_error_heywood = 0.5 * (np.array(data_heywood['Bin_Upper'])-np.array(data_heywood['Bin_Lower']))/1000
        dn_ds_heywood = np.array(data_heywood['Corrected_dNdS_S2.5'])
        dn_ds_err_upp_heywood = np.array(data_heywood['Raw_dNdS_S2.5_Upper_Error'])
        dn_ds_err_low_heywood = np.array(data_heywood['Raw_dNdS_S2.5_Lower_Error'])        
        # ax.errorbar(bins_heywood, dn_ds_heywood, xerr=bins_error_heywood, yerr=[dn_ds_err_low_heywood,dn_ds_err_upp_heywood], zorder=1,
        #             marker='v',  linestyle='none', label='Heywood et al. (2020)', capthick=2, capsize=3,  color='blue', markersize='6')
    
    if sample == 'GPS':
        data_dezotti = pd.read_csv('/net/vdesk/data2/bach1/ballieux/master_project_1/source_counts/DeZotti.csv')
        bins_dezotti= np.array(data_dezotti['log_s_bins'])
        dn_ds_dezotti = np.array(data_dezotti['s_counts'])
        dn_ds_err_upp_dezotti = np.array(data_dezotti['upper_error'])
        dn_ds_err_low_dezotti = np.array(data_dezotti['lower_error'])        
        ax.errorbar(10**bins_dezotti, dn_ds_dezotti, yerr=[dn_ds_err_low_dezotti, dn_ds_err_upp_dezotti], zorder=1,
                    marker='v',  linestyle='none', label='De Zotti et al. (2010)', capthick=2, capsize=3,  color='green', markersize='6')

    # My samples
    if sample == 'GPS':
        ax.errorbar(central_bins_PS_GPS, dn_ds_PS_GPS, yerr=dn_ds_err_PS_GPS, xerr=bins_width_PS_GPS/2.,
                    marker='s', label='GPS sample', linestyle='none', capthick=2, zorder=10, markersize='6', color='black' )
        # ax.errorbar(central_bins_PS_HF, dn_ds_PS_HF, yerr=dn_ds_err_PS_HF, xerr=bins_width_PS_HF/2.,
        #             marker='s', label='HF sample at 1400MHz', linestyle='none', capthick=2, zorder=10, markersize='6' )
        ax.set_xlim([0.001, 10])
        ax.set_ylim([0.05, 1000])
    if sample == 'MPS':
        # ax.errorbar(central_bins_PS_MPS[1:], dn_ds_PS_MPS[1:], yerr=dn_ds_err_PS_MPS[:,1:], xerr=bins_width_PS_MPS[1:]/2.,
        #             marker='s', label='MPS sample', linestyle='none', capthick=2, zorder=10, markersize='6', color='black')
        ax.errorbar(central_bins_PS_MPS, dn_ds_PS_MPS, yerr=dn_ds_err_PS_MPS, xerr=bins_width_PS_MPS/2.,
                    marker='s', label='MPS sample', linestyle='none', capthick=2, zorder=10, markersize='6', color='black')
        # ax.errorbar(central_bins_PS_LF, dn_ds_PS_LF, yerr=dn_ds_err_PS_LF, xerr=bins_width_PS_LF/2.,
        #             marker='s', label='LF sample', linestyle='none', capthick=2, zorder=10, markersize='6')
        ax.set_xlim([0.001, 100])
        ax.set_ylim([0.05, 10000])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='both', top=True, right=True, bottom=True, labelbottom=False,  labeltop=True)
    ax.legend(loc='lower right', fontsize=12)  # , bbox_to_anchor=(1, 0.5))

    ax.set_ylabel('S$^{5/2}$ dN/dS (Jy$^{1.5}$ sr$^{-1}$ )')

    #Make lower plot, residuals
    ax2 = plt.subplot(gs[1])
    ax2.sharex(ax)

    # def Massardi_model_MPS(S, scale):
    #     #Interpolates the model
    #     return np.interp(S, 10**(dZm_150["S_mid"]), 10**(dZm_150["scount"])/scale)
    

    # def Massardi_model_GPS(S, scale):
    #     #Interpolates the model
    #     return np.interp(S, 10**(dZm_1400["S_mid"]), (10**fitmodel(dZm_1400["S_mid"]))/scale)    
    
    #Below we fit to the model of Massardi for MPS, GPS and Snellen
    # if sample == 'MPS':
    #     popt_MPS, pcov_MPS = curve_fit(Massardi_model_MPS, central_bins_PS_MPS[1:], dn_ds_PS_MPS[1:], sigma=np.median(dn_ds_err_PS_MPS, axis=0)[1:])
    #     scaling_MPS = popt_MPS[0]
    #     err_scaling_MPS = np.sqrt(np.diag(pcov_MPS))[0]
    #     print('Best scaling for MPS sample', scaling_MPS, '$\pm$', err_scaling_MPS)
    #     interpol_MPS = Massardi_model_MPS(central_bins_PS_MPS[1:], scaling_MPS) #interpolated Massardi at MPS flux
    # if sample == 'GPS':
    #     popt_GPS, pcov_GPS = curve_fit(Massardi_model_GPS, central_bins_PS_GPS, dn_ds_PS_GPS, sigma=np.median(dn_ds_err_PS_GPS, axis=0))
    #     scaling_GPS = popt_GPS[0]
    #     err_scaling_GPS = np.sqrt(np.diag(pcov_GPS))[0]
    #     print('Best scaling for VLASS/GPS sample', scaling_GPS, '$\pm$', err_scaling_GPS)
    #     interpol_GPS = Massardi_model_GPS(central_bins_PS_GPS, scale=scaling_GPS)#interpolated Massardi at GPS flux

    def ODR_massardi_model_MPS(p,x):
        """
        This is a model such that we can do the ODR fitting for the MPS sample, where we have interpolate the Masardi 150 Mhz
        model at the values of x, divided by a parameter p = scale.
        """
        scale=p #Parameterspace
        return np.interp(x, 10**(dZm_150["S_mid"]), 10**(dZm_150["scount"])/scale)
    
    if sample == 'MPS':
        #Do the ODR fitting
        ODR_model_MPS = Model(ODR_massardi_model_MPS)
        data_MPS = RealData(central_bins_PS_MPS, dn_ds_PS_MPS, sx=bins_width_PS_MPS/2., sy=np.median(dn_ds_err_PS_MPS, axis=0))
        # data_MPS = RealData(central_bins_PS_MPS[1:], dn_ds_PS_MPS[1:], sx=bins_width_PS_MPS[1:]/2., sy=np.median(dn_ds_err_PS_MPS, axis=0)[1:])
        odr_MPS = ODR(data_MPS, ODR_model_MPS, beta0=[45])
        out_MPS = odr_MPS.run()
        print('using ODR, this is the best scaling for MPS',out_MPS.beta[0], '+/-', out_MPS.sd_beta[0])
        scaling_MPS = out_MPS.beta[0]
        err_scaling_MPS = out_MPS.sd_beta[0]
        interpol_MPS = ODR_massardi_model_MPS(scaling_MPS, central_bins_PS_MPS)#interpolated Massardi at GPS flux
        

        """Here I test a method to fit the GPS model with uncertainties in both x and y using orthogonal distance regression"""
    def ODR_massardi_model_GPS(p,x):
        """
        This is a model such that we can do the ODR fitting for the GPS sample, where we interpolate the Masardi 1400 Mhz
        model at the values of x, divided by a parameter p = scale. Fitmodel is because we have previously smoothed this function since I 
        read it off from a picture.
        """
        scale=p #Parameterspace
        return np.interp(x, 10**(dZm_1400["S_mid"]), (10**fitmodel(dZm_1400["S_mid"]))/scale) 
    
    if sample == 'GPS':
        #TODO: remove here the incomplete bins
        ODR_model_GPS = Model(ODR_massardi_model_GPS)
        #Here we neglect the first 5 bins, which is kind of based on nothing 
        data_GPS = RealData(central_bins_PS_GPS[5:], dn_ds_PS_GPS[5:], sx=bins_width_PS_GPS[5:]/2., sy=np.median(dn_ds_err_PS_GPS, axis=0)[5:])
        odr_GPS = ODR(data_GPS, ODR_model_GPS, beta0=[30])
        out_GPS = odr_GPS.run()
        print('using ODR, this is the best scaling for GPS',out_GPS.beta[0], '+/-', out_GPS.sd_beta[0])
        scaling_GPS = out_GPS.beta[0]
        err_scaling_GPS = out_GPS.sd_beta[0]
        interpol_GPS = ODR_massardi_model_GPS(scaling_GPS, central_bins_PS_GPS)#interpolated Massardi at GPS flux

    if sample == 'GPS':
        data_snellen = RealData(bins_1400_snellen, dn_ds_1400_snellen, sx=bin_width_snellen/2., sy=np.median(dn_ds_err_1400_snellen, axis=0))
        odr_snellen = ODR(data_snellen, ODR_model_GPS, beta0=[50])
        out_snellen = odr_snellen.run()
        print('using ODR, this is the best scaling for Snellen',out_snellen.beta[0], '+/-', out_snellen.sd_beta[0])   


    # if sample == 'GPS':
    #     popt_Snellen, pcov_Snellen = curve_fit(Massardi_model_GPS, bins_1400_snellen, dn_ds_1400_snellen, sigma=np.median(dn_ds_err_1400_snellen, axis=0))
    #     scaling_snellen = popt_Snellen[0]
    #     err_scaling_snellen = np.sqrt(np.diag(pcov_Snellen))[0]
    #     print('Best scaling for Snellen sample at 1400MHz', scaling_snellen, '$\pm$', err_scaling_snellen) 
        
    # rescale_061=(150/1400)**(-0.61)
    # popt_check, pcov_check = curve_fit(Massardi_model_check, \
    #                 (10**(dZm_1400["S_mid"]))*rescale_061, 10**fitmodel((dZm_1400["S_mid"]))*rescale_061**1.5)
    # scaling_check = popt_check[0]
    # print('check', scaling_check) 

    #Residuals
    #These for chi^2
    if sample == 'GPS':
        ax2.errorbar(central_bins_PS_GPS, (dn_ds_PS_GPS - interpol_GPS) / np.median(dn_ds_err_PS_GPS, axis=0), yerr=1,  xerr=bins_width_PS_GPS/2.,
                      marker='s', label='GPS residual to Massardi (2010) 1400MHz model, for 29 $\pm$ 1 scaling,'.format(scaling_GPS, err_scaling_GPS), linestyle='none')
        ax2.set_ylim([-8.5, 5])
    if sample == 'MPS':
        ax2.errorbar(central_bins_PS_MPS, (dn_ds_PS_MPS - interpol_MPS) / np.median(dn_ds_err_PS_MPS, axis=0), yerr=1,  xerr=(bins_width_PS_MPS)/2.,
                      marker='s', label='MPS residual to Massardi (2010) 150MHz model, for 45 $\pm$ 2 scaling'.format(scaling_MPS, err_scaling_MPS), linestyle='none')

        ax2.set_ylim([-6, 3])
    ax2.set_ylabel('$\chi$')
    ax2.set_xscale('log')
    
    #These for absolute errors
    # ax2.errorbar(central_bins_PS_VLASS, (dn_ds_PS_VLASS-interpol_VLASS), yerr=dn_ds_err_PS_VLASS, xerr=bins_width_PS_VLASS/2.,
    #              linestyle='none', marker='s', label='residual GPS sample rescaled to 144MHz')
    # ax2.errorbar(central_bins_PS_LoLSS, (dn_ds_PS_LoLSS-interpol_LoLSS),  yerr=dn_ds_err_PS_LoLSS, xerr=bins_width_PS_LoLSS/2.,
    #              linestyle='none',marker='s', label='residual MPS sample')
    # ax2.set_ylabel('residual to Massandi (2010) /40')
    # ax2.set_ylim([-50, 50])


    ax2.tick_params(axis='both', which='both', top=True, right=True, labeltop=False)
    ax2.legend(loc='lower right', fontsize=12)  # , bbox_to_anchor=(1, 0.5))
    plt.subplots_adjust(hspace=0.)
    ax2.hlines(0, 1e-3, 100)
    if sample == 'MPS':
        ax2.set_xlabel('S$_{144\, \mathrm{MHz}}$(Jy)')
    if sample == 'GPS':
        ax2.set_xlabel('S$_{1400\, \mathrm{MHz}}$(Jy)')
    
    if sample == 'MPS':
        plt.savefig(
            "/net/vdesk/data2/bach1/ballieux/master_project_1/plots/source_counts_MPS.pdf")
    if sample == 'GPS':
        plt.savefig(
            "/net/vdesk/data2/bach1/ballieux/master_project_1/plots/source_counts_GPS_new_e_alpha_low.pdf")
    # plt.show()
    plt.close()
    
    # plt.figure(3, figsize=(10,8))
    # plt.semilogx()
    # plt.semilogy()
    # plt.plot(10**(dZm_150["S_mid"]), 10**(dZm_150["scount"]),
    #         color='#999999', label="Massardi 150MHz")
    # rescale_08=(150/1400)**(-0.8)
    # plt.plot((10**(dZm_1400["S_mid"]))*rescale_08, 10**fitmodel((dZm_1400["S_mid"]))*rescale_08**1.5,
    #         color='blue', label="Massardi 1400MHz rescaled -0.8")
    # rescale_061=(150/1400)**(-0.61)
    # plt.plot((10**(dZm_1400["S_mid"]))*rescale_061, 10**fitmodel((dZm_1400["S_mid"]))*rescale_061**1.5,
    #         color='red', label="Massardi 1400MHz rescaled -0.61")
    # plt.xlim(1e-3, 1e2)
    # plt.ylim(1e-1,1e4)
    # plt.legend()
    
    
    
    
    #This is for a seperate figure
    # plt.figure(4)
    # scal = scaling_LoLSS
    # plt.errorbar(central_bins_PS_LoLSS, dn_ds_PS_LoLSS, yerr=dn_ds_err_PS_LoLSS, xerr=bins_width_PS_LoLSS/2.,
    #             marker='s', label='MPS', linestyle='none', capthick=2)
    # # plt.errorbar(central_bins_PS_LoLSS, dn_ds_PS_LoLSS, yerr=dn_ds_err_PS_LoLSS, xerr=bins_width_PS_LoLSS/2.,
    # #             marker='s', label='MPS sample', linestyle='none', capthick=2)
    # plt.plot(10**(dZm["S_mid"]), 10**(dZm["scount"])/scal, ls='--',
    #         color='#999999', label="Massardi (2010) model * 1/{:.3}".format(scal))
    # plt.xlim([0.001, 100])
    # plt.ylim([0.05, 10000])
    # plt.semilogx()
    # plt.semilogy()
    # plt.legend()

    # plt.ylabel('S$^{5/2}$ dN/dS (Jy$^{1.5}$ sr$^{-1}$ )')
    # plt.xlabel('S$_{144\, \mathrm{MHz}}$ (Jy)')



    # return dn_ds_call, dn_ds_err_call, bins_call, bin_width_call


# Read in the master sample
hdulist = fits.open(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_clean.fits')
tbdata = hdulist[1].data
hdulist.close()

# These are still the full arrays
name = tbdata['LoTSS_name']
alpha_low_VLASS = tbdata['alpha_low_VLASS']
alpha_high_VLASS = tbdata['alpha_high_VLASS']
e_alpha_low_VLASS = tbdata['e_alpha_low_VLASS']
e_alpha_high_VLASS = tbdata['e_alpha_high_VLASS']

alpha_low_LoLSS = tbdata['alpha_low_LoLSS']
alpha_high_LoLSS = tbdata['alpha_high_LoLSS']
e_alpha_low_LoLSS = tbdata['e_alpha_low_LoLSS']
e_alpha_high_LoLSS = tbdata['e_alpha_high_LoLSS']

NVSS_flux = tbdata['NVSS_flux']
LoTSS_flux = tbdata['LoTSS_flux']
VLASS_flux = tbdata['VLASS_flux']
LoLSS_flux = tbdata['LoLSS_flux']

# Now we want to only select those in the right sample
non_zero_VLASS = np.where(VLASS_flux > 0.)
non_zero_LoLSS = np.where(LoLSS_flux > 0.)

# Already define the sample-specfic columns to be non-zero
alpha_low_VLASS = alpha_low_VLASS[non_zero_VLASS]
alpha_high_VLASS = alpha_high_VLASS[non_zero_VLASS]
e_alpha_low_VLASS = e_alpha_low_VLASS[non_zero_VLASS]
e_alpha_high_VLASS = e_alpha_high_VLASS[non_zero_VLASS]

alpha_low_LoLSS = alpha_low_LoLSS[non_zero_LoLSS]
alpha_high_LoLSS = alpha_high_LoLSS[non_zero_LoLSS]
e_alpha_low_LoLSS = e_alpha_low_LoLSS[non_zero_LoLSS]
e_alpha_high_LoLSS = e_alpha_high_LoLSS[non_zero_LoLSS]


"""
For the MPS sample: 
"""
# these are in MPS sample
ind_peaked_MPS = np.where((alpha_low_LoLSS >= (e_alpha_low_LoLSS))
                            & (alpha_high_LoLSS <= 0))

# The LoTSS flux that is non-zero for the MPS sample
LoTSS_flux_1 = LoTSS_flux[non_zero_LoLSS]
# The LoTSS flux that is PS for the LoLSS sample
LoTSS_flux_MPS = LoTSS_flux_1[ind_peaked_MPS]

hist(LoTSS_flux_MPS, 'MPS')
dA_deg_LoLSS = 650 * (u.deg)**2#740 * (u.deg)**2  # LoLSS survey area in deg2
dA_LoLSS = dA_deg_LoLSS.to(u.sr)  # survey area in sr
print("")
print("The source counts for the LoLSS sample")
# We take the median of the error of alpha_low, coming from the selection criterion
alpha_low_scaling_MPS = np.median(alpha_low_LoLSS[ind_peaked_MPS])
# alpha_low_scaling_MPS = np.median(e_alpha_low_LoLSS)
print(alpha_low_scaling_MPS)
lower_limit_PS_MPS = (11/1000) * (144/54) ** (alpha_low_scaling_MPS)

print('95% completeness limit in LoLSS in Jy', lower_limit_PS_MPS)
dn_ds_PS_MPS, dn_ds_err_PS_MPS, central_bins_PS_MPS, bins_width_PS_MPS, n_PS_MPS =\
    source_counts(LoTSS_flux_MPS, lower_limit_PS_MPS,
                  3.65, 10, dA_LoLSS.value)

"""
For the GPS sample
"""
hdulist2 = fits.open(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/new_e_alpha_low_HF.fits')
tbdata2 = hdulist2[1].data
orig_cols2 = hdulist[1].columns
hdulist2.close()
new_e_alpha_low_HF = tbdata2['new_e_alpha_low_HF']

# Which ones are in GPS sample
ind_peaked_GPS = np.where((alpha_low_VLASS >= (new_e_alpha_low_HF))
                            & (alpha_high_VLASS <= 0))

NVSS_flux_2 = NVSS_flux[non_zero_VLASS]  # NVSS flux non-zero for VLASS sample
NVSS_flux_GPS = NVSS_flux_2[ind_peaked_GPS]
hist(NVSS_flux_GPS, 'GPS', '1400')
dA_deg_VLASS = 5634 * (u.deg)**2  # loTSS survey area
dA_VLASS = dA_deg_VLASS.to(u.sr)  # survey area in sr

# We take the median of alpha_low as the scaling, not specifically for the PS sources
alpha_low_scaling_VLASS = np.median(alpha_low_VLASS)


print("")
print("The NVSS source counts for the GPS sample, not rescaled to 144MHZ")
dn_ds_PS_GPS, dn_ds_err_PS_GPS, central_bins_PS_GPS, bins_width_PS_GPS, n_PS_GPS \
    = source_counts(NVSS_flux_GPS, 3.24/1000, 2.4104, 18, dA_VLASS.value)  # in NVSS frame, lower limit is also in NVSS
    #TODO: check these limits are still correct

    

print("")
print('Minimum NVSS flux in GPS sample:', np.min(NVSS_flux_GPS))
print('Maximum NVSS in GPS sample:', np.max(NVSS_flux_GPS))
print('Minimum LoTSS flux in MPS sample:', np.min(LoTSS_flux_MPS))
print('Maximum LoTSS in MPS sample:', np.max(LoTSS_flux_MPS))
print("")

# make the source counts plot
source_counts_plot(sample='MPS')
source_counts_plot(sample='GPS')