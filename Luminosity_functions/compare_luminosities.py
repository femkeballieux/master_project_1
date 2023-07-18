"""
Compare_luminosities
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.gridspec as gridspec
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy import stats

hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_sample_with_redshift.fits')
tbdata = hdulist[1].data
hdulist.close()

cosmo = FlatLambdaCDM(H0 = 70 * u.km / u.s / u.Mpc, Om0=0.27)

def L_z(z, alpha, flux):
    dl = cosmo.luminosity_distance(z).to(u.m).value
    k_corr = np.power((1 + z),-(1+alpha))
    lum = 4 * np.pi * dl**2 * flux * k_corr * 1e-26 # W Hz^-1 m^-2
    return lum

z = tbdata['z_best']

name = tbdata['LoTSS_name_1']
alpha_low_LF = tbdata['alpha_low_LoLSS']
e_alpha_low_LF = tbdata['e_alpha_low_LoLSS']
alpha_high_LF = tbdata['alpha_high_LoLSS']
e_alpha_high_LF = tbdata['e_alpha_high_LoLSS']
alpha_low_HF = tbdata['alpha_low_VLASS']
e_alpha_low_HF = tbdata['e_alpha_low_VLASS']
alpha_high_HF = tbdata['alpha_high_VLASS']
e_alpha_high_HF = tbdata['e_alpha_high_VLASS']
S_1400 = tbdata['NVSS_flux']
# norm_high = tbdata['norm_high']
S_144 = tbdata['LoTSS_flux']

ind_LF = ((z>0.) & (tbdata['LoLSS_flux']>0.))
ind_HF =((z>0.) & (tbdata['VLASS_flux']>0.))
ind_MPS = ((alpha_low_LF > e_alpha_low_LF)\
                & (alpha_high_LF < - e_alpha_high_LF))
ind_GPS = ((alpha_low_HF > e_alpha_low_HF)\
                & (alpha_high_HF < - e_alpha_high_HF))


scaling = -0.8#-0.6146366 #np.median(alpha_low_HF[np.where((tbdata['VLASS_flux']>0.))])
print(scaling)
L_LF=L_z(z[np.where(ind_LF)],alpha_high_LF[np.where(ind_LF)],S_144[np.where(ind_LF)])
L_HF=L_z(z[np.where(ind_HF)],alpha_high_LF[np.where(ind_HF)],S_1400[np.where(ind_HF)]*(144/1400)**scaling)

L_MPS=L_z(z[np.where(ind_LF & ind_MPS)],alpha_high_LF[np.where(ind_LF & ind_MPS)],S_144[np.where(ind_LF & ind_MPS)])
L_GPS=L_z(z[np.where(ind_HF & ind_GPS)],alpha_high_LF[np.where(ind_HF & ind_GPS)], S_1400[np.where(ind_HF & ind_GPS)]*(144/1400)**scaling)
print(np.mean(L_LF),np.mean(L_HF), np.mean(L_MPS), np.mean(L_GPS))

plt.figure(0,figsize=(10,8))
plt.hist(L_MPS, bins=50, log=True, histtype='stepfilled', label='MPS')
plt.hist(L_GPS, bins=50, log=True, histtype='stepfilled',alpha=0.5, label='GPS')
plt.xlabel('Luminosity [...]')
plt.legend()

plt.figure(1,figsize=(10,8))
plt.scatter(z[np.where(ind_LF & ind_MPS)], L_MPS, zorder=2, label='MPS')
plt.scatter(z[np.where(ind_HF & ind_GPS)], L_GPS, zorder=1, alpha=0.3, label='GPS')
plt.yscale('log')
plt.ylabel('Luminosity [...]')
plt.xlabel('z')
plt.legend()

plt.figure(2,figsize=(10,8))
plt.scatter(z[np.where(ind_LF)], L_LF, zorder=2, alpha=0.3, label='HF')
plt.scatter(z[np.where(ind_HF)], L_HF, zorder=1, alpha=0.3, label='LF')
plt.ylabel('Luminosity [...]')
plt.xlabel('z')
plt.yscale('log')
plt.legend()

"""
Need to rescale NVSS flux to LoTSS with 
"""









