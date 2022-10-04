# -*- coding: utf-8 -*-
"""
change a column name in a fits file
"""
from astropy.table import Table
from astropy.io import fits

file = 'C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\crossmatch_NVSS_TGSS_VLSSr_LoLSS_inband.fits'
survey = fits.open(file, memmap=True)
data = Table(survey[1].data)
survey.close()

print('hoi')
data['Total_flux_LoLLS'].name = 'Total_flux_LoLSS'
data['E_Total_flux_LoLLS'].name = 'E_Total_flux_LoLSS'
data['Peak_flux_LoLLS'].name = 'Peak_flux_LoLSS'
data['E_Peak_flux_LoLLS'].name = 'E_Peak_flux_LoLSS'

data.write('C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\crossmatch_NVSS_TGSS_VLSSr_LoLSS_inband.fits', format='fits', overwrite = True)