# -*- coding: utf-8 -*-
"""
code written by Femke Ballieux, takes the crossmatch between master_sample and inband fluxes of LoLSS
and throws away all that is not in the master sample. Also renames the columns
run on computer
"""
import numpy as np
from astropy.io import fits
import gpscssmodels
import seds_plot_func_extreme
from tqdm import tqdm
import scipy.optimize as opt
import warnings
import time
warnings.filterwarnings("ignore")

#read in the data
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_master_inbandLoLSS.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()# fit curvature models with rms cuts

print(orig_cols)

S_channel0, e_S_channel0 = tbdata['Total_flux_2']/1000., tbdata['E_Total_flux_2']/1000. # Jy
S_channel0, e_S_channel0 = np.where(np.isnan(S_channel0), 0, S_channel0), np.where(np.isnan(S_channel0), 0, e_S_channel0)

S_channel1, e_S_channel1 = tbdata['Total_flux_3']/1000., tbdata['E_Total_flux_3']/1000. # Jy
S_channel1, e_S_channel1 = np.where(np.isnan(S_channel1), 0, S_channel1), np.where(np.isnan(S_channel1), 0, e_S_channel1)

S_channel2, e_S_channel2 = tbdata['Total_flux_4']/1000., tbdata['E_Total_flux_4']/1000. # Jy
S_channel2, e_S_channel2 = np.where(np.isnan(S_channel2), 0, S_channel2), np.where(np.isnan(S_channel2), 0, e_S_channel2)

S_channel3, e_S_channel3 = tbdata['Total_flux_5']/1000., tbdata['E_Total_flux_5']/1000. # Jy
S_channel3, e_S_channel3 = np.where(np.isnan(S_channel3), 0, S_channel3), np.where(np.isnan(S_channel3), 0, e_S_channel3)

S_channel4, e_S_channel4 = tbdata['Total_flux_6']/1000., tbdata['E_Total_flux_6']/1000. # Jy
S_channel4, e_S_channel4 = np.where(np.isnan(S_channel4), 0, S_channel4), np.where(np.isnan(S_channel4), 0, e_S_channel4)

S_channel5, e_S_channel5 = tbdata['Total_flux_7']/1000., tbdata['E_Total_flux_7']/1000. # Jy
S_channel5, e_S_channel5 = np.where(np.isnan(S_channel5), 0, S_channel5), np.where(np.isnan(S_channel5), 0, e_S_channel5)

LoLSS_rms = tbdata['LoLSS_rms']/1000 #Jy

col1 = fits.Column(name='LoTSS_name', format = '34A', array = tbdata['LoTSS_name'])
col2 = fits.Column(name='RA', format = 'E', array= tbdata['RA_1'])
col3 = fits.Column(name='Dec', format = 'E', array = tbdata['Dec_1'])
col4 = fits.Column(name='LoTSS_flux', format = 'E', array = tbdata['LoTSS_flux'])
col5 = fits.Column(name='e_LoTSS_flux', format = 'E', array = tbdata['e_LoTSS_flux'])
col6 = fits.Column(name='a_low', format = 'E', array = tbdata['a_low'])
col7 = fits.Column(name='alpha_low', format = 'E', array = tbdata['alpha_low'])
col8 = fits.Column(name='e_alpha_low', format = 'E', array = tbdata['e_alpha_low'])
col9 = fits.Column(name='a_high', format = 'E', array = tbdata['a_high'])
col10 = fits.Column(name='alpha_high', format = 'E', array = tbdata['alpha_high'])
col11 = fits.Column(name='e_alpha_high', format = 'E', array = tbdata['e_alpha_high'])
col12 = fits.Column(name='LoLSS_RA', format = 'E', array = tbdata['LoLSS_RA'])
col13 = fits.Column(name='LoLSS_Dec', format = 'E', array = tbdata['LoLSS_Dec'])
col14 = fits.Column(name='LoLSS_flux', format = 'E', array = tbdata['LoLSS_flux'])
col15 = fits.Column(name='e_LoLSS_flux', format = 'E', array = tbdata['e_LoLSS_flux'])
col16 = fits.Column(name='NVSS_RA', format = 'E', array = tbdata['NVSS_RA'])
col17 = fits.Column(name='NVSS_Dec', format = 'E', array = tbdata['NVSS_Dec'])
col18 = fits.Column(name='NVSS_flux', format = 'E', array = tbdata['NVSS_flux'])
col19 = fits.Column(name='e_NVSS_flux', format = 'E', array = tbdata['e_NVSS_flux'])
col20 = fits.Column(name='TGSS_RA', format = 'E', array = tbdata['TGSS_RA'])
col21 = fits.Column(name='TGSS_Dec', format = 'E', array = tbdata['TGSS_Dec'])
col22 = fits.Column(name='TGSS_flux', format = 'E', array = tbdata['TGSS_flux'])
col23 = fits.Column(name='e_TGSS_flux', format = 'E', array = tbdata['e_TGSS_flux'])
col24 = fits.Column(name='VLSSr_RA', format = 'E', array = tbdata['VLSSr_RA'])
col25 = fits.Column(name='VLSSr_Dec', format = 'E', array = tbdata['VLSSr_Dec'])
col26 = fits.Column(name='VLSSr_flux', format = 'E', array = tbdata['VLSSr_flux'])
col27 = fits.Column(name='e_VLSSr_flux', format = 'E', array = tbdata['e_VLSSr_flux'])
col28 = fits.Column(name='FIRST_RA', format = 'E', array = tbdata['FIRST_RA'])
col29 = fits.Column(name='FIRST_Dec', format = 'E', array = tbdata['FIRST_Dec'])
col30 = fits.Column(name='FIRST_flux', format = 'E', array = tbdata['FIRST_flux'])
col31 = fits.Column(name='e_FIRST_flux', format = 'E', array = tbdata['e_FIRST_flux']) 
col32 = fits.Column(name='inband_RA', format = 'E', array = tbdata['inband_RA'])
col33 = fits.Column(name='inband_Dec', format = 'E', array = tbdata['inband_Dec'])
col34 = fits.Column(name='S_inband_low', format = 'E', array = tbdata['S_inband_low'])
col35 = fits.Column(name='e_S_inband_low', format = 'E', array = tbdata['e_S_inband_low'])
col36 = fits.Column(name='S_inband_mid', format = 'E', array = tbdata['S_inband_mid'])
col37 = fits.Column(name='e_S_inband_mid', format = 'E', array = tbdata['e_S_inband_mid'])
col38 = fits.Column(name='S_inband_high', format = 'E', array = tbdata['S_inband_high'])
col39 = fits.Column(name='e_S_inband_high', format = 'E', array =tbdata['e_S_inband_high'])

col40 = fits.Column(name='channel0_RA', format = 'E', array= tbdata['RA_2'])
col41 = fits.Column(name='channel0_Dec', format = 'E', array = tbdata['Dec_2'])
col42 = fits.Column(name='channel0_flux', format = 'E', array = S_channel0) #Jansky
col43 = fits.Column(name='e_channel0_flux', format = 'E', array = e_S_channel0)

col44 = fits.Column(name='channel1_RA', format = 'E', array= tbdata['RA_3'])
col45 = fits.Column(name='channel1_Dec', format = 'E', array = tbdata['Dec_3'])
col46 = fits.Column(name='channel1_flux', format = 'E', array = S_channel1)
col47 = fits.Column(name='e_channel1_flux', format = 'E', array = e_S_channel1)

col48 = fits.Column(name='channel2_RA', format = 'E', array= tbdata['RA_4'])
col49 = fits.Column(name='channel2_Dec', format = 'E', array = tbdata['Dec_4'])
col50 = fits.Column(name='channel2_flux', format = 'E', array = S_channel2)
col51 = fits.Column(name='e_channel2_flux', format = 'E', array = e_S_channel2)

col52 = fits.Column(name='channel3_RA', format = 'E', array= tbdata['RA_5'])
col53 = fits.Column(name='channel3_Dec', format = 'E', array = tbdata['Dec_5'])
col54 = fits.Column(name='channel3_flux', format = 'E', array = S_channel3)
col55 = fits.Column(name='e_channel3_flux', format = 'E', array = e_S_channel3)

col56 = fits.Column(name='channel4_RA', format = 'E', array= tbdata['RA_6'])
col57 = fits.Column(name='channel4_Dec', format = 'E', array = tbdata['Dec_6'])
col58 = fits.Column(name='channel4_flux', format = 'E', array = S_channel4)
col59 = fits.Column(name='e_channel4_flux', format = 'E', array = e_S_channel4)

col60 = fits.Column(name='channel5_RA', format = 'E', array= tbdata['RA_7'])
col61 = fits.Column(name='channel5_Dec', format = 'E', array = tbdata['Dec_7'])
col62 = fits.Column(name='channel5_flux', format = 'E', array = S_channel5)
col63 = fits.Column(name='e_channel5_flux', format = 'E', array = e_S_channel5)

col64 = fits.Column(name='LoLSS_rms', format = 'E',  array = LoLSS_rms)
col65 = fits.Column(name='LoLSS_maj', format = 'E', array = tbdata['LoLSS_maj'])
col66 = fits.Column(name='LoLSS_S_code', format = '1A', array = tbdata['LoLSS_S_code'])

cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19,\
    col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,\
        col64, col65, col66, col41, col42 ,col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, \
                    col58, col59, col60, col61, col62, col63])
    
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_LoLSS_with_inband.fits', overwrite = True)