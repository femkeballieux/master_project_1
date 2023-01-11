"""
Code written by Femke Ballieux, takes in the master sample and crossmatches it with vlass
"""
#run on laptop
import os
import numpy as np


os.system('java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts \
          tmatchn join1=always matcher=sky multimode=pairs nin=2 params=15 \
    in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_sample.fits values1="LoLSS_RA LoLSS_Dec" \
    in2=/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/vlass_catalog.fits values2="ra dec" \
    out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_with_vlass.fits')

print("yey we did it")

from astropy.table import Table
#change the columnnames
t = Table.read("/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_with_vlass.fits")
t['island'].name = 'island_VLASS'
t['source'].name = 'source_VLASS'
t['background'].name = 'background_VLASS'
t['local_rms'].name = 'local_rms_VLASS'
t['ra_2'].name = 'RA_VLASS'
t['err_ra'].name = 'E_RA_VLASS'
t['dec_2'].name = 'Dec_VLASS'
t['err_dec'].name = 'E_Dec_VLASS'
t['peak_flux'].name = 'peak_flux_VLASS'
t['err_peak_flux'].name = 'E_peak_flux_VLASS'
t['int_flux'].name = 'int_flux_VLASS'
t['err_int_flux'].name = 'E_int_flux_VLASS'

#All with error -1 from aegean have a bad fit
badfit_mask = (t['E_peak_flux_VLASS']== -1. )

#make sure that non-detections become upper limits (does not deal yet with failed aegean)
vlass_sensitivity = 120e-6 #jansky
t['peak_flux_VLASS'], t['E_peak_flux_VLASS'] = np.where(np.isnan(t['peak_flux_VLASS']),\
                                vlass_sensitivity, t['peak_flux_VLASS']), np.where(np.isnan(t['E_peak_flux_VLASS']), -1., t['E_peak_flux_VLASS'])

t['int_flux_VLASS'], t['E_int_flux_VLASS'] = np.where(np.isnan(t['int_flux_VLASS']),\
                                vlass_sensitivity, t['int_flux_VLASS']), np.where(np.isnan(t['E_int_flux_VLASS']), -1., t['E_int_flux_VLASS'])

#propagate 10 % flux error
t['E_int_flux_VLASS_full']= np.sqrt(t['E_int_flux_VLASS']**2 + (0.15 * t['int_flux_VLASS'])**2)
t['E_peak_flux_VLASS_full']= np.sqrt(t['E_peak_flux_VLASS']**2 + (0.15 * t['peak_flux_VLASS'])**2)

t['E_peak_flux_VLASS_full']=np.where(t['E_peak_flux_VLASS_full']==1., -1., t['E_peak_flux_VLASS_full'])
t['E_int_flux_VLASS_full']=np.where(t['E_int_flux_VLASS_full']==1., -1., t['E_int_flux_VLASS_full'])

#all with bad fit have error of 15 percent of flux
t['E_int_flux_VLASS_full'][badfit_mask] = 0.25 * t['int_flux_VLASS'][badfit_mask] 
t['E_peak_flux_VLASS_full'][badfit_mask] = 0.25 * t['peak_flux_VLASS'][badfit_mask]



t.write("/net/vdesk/data2/bach1/ballieux/master_project_1/data/PS_with_vlass.fits", overwrite=True)
