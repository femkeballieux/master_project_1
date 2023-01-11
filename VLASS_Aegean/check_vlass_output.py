# -*- coding: utf-8 -*-
"""
Takes in the vlass files, turns it into a complete catalogue
Also checks the duplicate indexes
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u


#checking the individual output files
# path = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images/catalog_outputs/'
# batch = '402_501/'
# filename = '37_206.800+60.7957_2020-08-04_comp.fits'
# hdulist = fits.open(path + batch + filename)

# tbdata = hdulist[1].data
# orig_cols = hdulist[1].columns

# hdulist.close()

# print(orig_cols)
# print('fluxes: should be in Jy')
# print(tbdata['int_flux']) 
# print('coordinates')
# print(tbdata['ra'])
# print(tbdata['dec']) 

#checking the entire output catalog
path = '/net/vdesk/data2/bach1/ballieux/master_project_1/data/' 
filename = 'PS_with_vlass.fits'

PS_with_vlass = fits.open(path + filename)
t = PS_with_vlass[1].data
print(PS_with_vlass[1].columns)
PS_with_vlass.close()

#here we check everything without a vlass match

#these do not have a very clear blop (364, 466 is boundary)
too_faint_indexes = np.array([83, 99, 108, 120, 133, 269, 361, 364, 413, 457, 465, 623, 677, 742, 746, 755]) 
slice_indexes = np.array([168,310, 359]) #these do have a clear source, but no output from aegean. Due to slice?
bad_fit_indexes = np.array([174, 395, 410]) #these have a -1 error in vlass, convolute flux with 15%
wrong_coords_index=np.array([436]) #does have an aegean output but does not get crossmatched
outside_footprint_indexes = np.array([466, 492]) #just no flux data there
def extrapolated_flux(i, freq_VLASS=3000 ):
    """
    s_v = a * v ^ alpha
    """
    return (freq_VLASS** t['alpha_high'][i]) *  t['a_high'][i]

#print(extrapolated_flux(too_faint_indexes), 'are the expected fluxes in Jy of the sources that seem too faint')

# for i in too_faint_indexes:
#     print(t['RA_1'][i],',', t['Dec_1'][i])

#TODO: Give these sources an upper limit corresponding to their individual RMS



print(t['RA_1'][slice_indexes], t['Dec_1'][slice_indexes])



