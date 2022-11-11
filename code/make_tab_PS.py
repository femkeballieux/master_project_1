#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code that takes in the master sample including LoLSS inband, and makes a PS sample
Only PS sources, with demands on maj, SNR and s_code from LoLSS
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

#code can be run from eother laptop or pc
path_laptop = 'C:/Users/Femke/Documents/GitHub/master_project_1/data'
path_vdesk= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'

#read in the data
hdulist = fits.open(path_vdesk + 'master_LoLSS_with_inband.fits')
tbdata = Table(hdulist[1].data)
hdulist.close()
orig_cols = hdulist[1].columns

print(orig_cols)

#get some useful arrays. Units should be allright
SNR = tbdata['LoLSS_flux']/tbdata['LoLSS_rms']
#These are all the PS sources with some demands 
mask = (tbdata['alpha_low'] >= 0.1) & (tbdata['alpha_high']<=0) #& (tbdata['LoLSS_S_code'] == 'S') & (tbdata['LoLSS_maj']<=30) & (SNR >=1.)

print(len(tbdata['RA'][mask])) 

"""
I seemed a little confused, because when I set SNR to 10 or something, way more sources fall away than
the plot from SNR, ratio would suggest. However, there is one key difference: This plot only shows the sources with a flux measurement in
both the PDR as DR1, so it is not representable for what I am doing here. Therefore make a SNR histogram to see if we are going right.
Changes in number of sources
#767 only PS
#734 only LoLSS S_code> apparantly this matters?
#732 when also maj<30
#482 when SNR>10

However, after meeting with Joe we decided to do VLASS atleast for the entire sample of PS sources, things will always change after some time
"""
plt.figure(figsize=(10,8))
plt.hist(SNR[SNR<=200], bins=100, zorder=1)
plt.vlines(10, 0, 700, zorder=10, color='black', label='SNR=10')
plt.vlines(20, 0, 700, zorder=10, color='red', label='SNR=20')
plt.legend()
plt.xlabel('SNR')
plt.show()


#PS_to_write = tbdata[:][mask] #This is how you get all columns, only masked rows
tbdata.remove_rows(np.where(~mask))
tbdata.write('/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_sample.fits', format='fits', overwrite = True)
tbdata.write('/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_sample.csv', format='csv', overwrite = True)
