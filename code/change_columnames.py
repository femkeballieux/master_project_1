# -*- coding: utf-8 -*-
"""
change a column name in a fits file
"""
from astropy.table import Table
from astropy.io import fits

file = 'C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\LoLSS_new.fits'
survey = fits.open(file, memmap=True)
data = Table(survey[1].data)
survey.close()

print('hoi')
data['Total_flux'].name = 'Total_flux_LoLLS'
data['E_Total_flux'].name = 'E_Total_flux_LoLLS'
data['Peak_flux'].name = 'Peak_flux_LoLLS'
data['E_Peak_flux'].name = 'E_Peak_flux_LoLLS'

data.write('C:\\Users\\Femke\\Documents\\GitHub\\master_project_1\\data\\surveys\\LoLLS_new_names.fits', format='fits', overwrite = True)