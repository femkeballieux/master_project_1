"""
Code to check the resolved issue, Shimwell formula. 
We use the 6'' crossmatched sample, of VLASS only duplicates cut out, nothing else. Then we take the sources that have NO VLASs

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdulist = fits.open(\
   '/net/vdesk/data2/bach1/ballieux/master_project_1/data/Official_VLASS_no_dups/official_mega_master_clean_6.fits')

tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
print(orig_cols)
hdulist.close()


LoTSS_int=tbdata['LoTSS_flux']
LoTSS_int_err = tbdata['e_LoTSS_flux']
LoTSS_peak=tbdata['LoTSS_flux_peak']
e_LoTSS_peak=tbdata['stat_e_LoTSS_flux_peak']
rms=tbdata['LoTSS_rms']

mask = np.where((tbdata['VLASS_flux']==0.)&(rms>0.0)&(rms<0.0005))

R=np.log(LoTSS_int/LoTSS_peak)
SNR=(LoTSS_int/LoTSS_int_err)
plt.figure(figsize=(10,8))
plt.semilogx()

#plt.scatter(SNR,R, alpha=1, s=1, label='all', color='crimson')
plt.scatter(SNR[mask],R[mask], c=rms[mask], alpha=0.5, s=1, label='missing VLASS')
plt.colorbar()

x_array=np.linspace(1,10000,1000)

def R_99(SNR):
    return 0.42 + (1.08 / (1+(SNR/96.57)**2.49))
plt.plot(x_array, R_99(x_array), color='black', label='Shimwell')
plt.legend()
plt.xlim(1, 10000)
plt.xlabel('SNR')
plt.ylabel('$\ln(S_I/S_P)$')