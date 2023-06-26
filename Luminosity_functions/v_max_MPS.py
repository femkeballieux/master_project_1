"""
Code written by Martje Slob, edited by Femke Ballieux
Calculated the v_max for a sample of radio sources, that has redshifts and SDSS
photometry. See Femkes thesis for more details.
This piece of code was adapted for my MPS sample
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from matplotlib import rc
from tqdm import tqdm
from calc_kcor import calc_kcor

plt.style.use('/net/vdesk/data2/bach1/ballieux/master_project_1/code/style.mplstyle')
rc('text', usetex=False)
rc('font', family='sans-serif')

def L_z(z, alpha, flux_lim):
    dl = cosmo.luminosity_distance(z).to(u.m).value
    k_corr = np.power((1 + z),-(1+alpha))
    lum = 4 * np.pi * dl**2 * flux_lim * k_corr * 1e-26 # W Hz^-1 m^-2
    return lum

def z_max_radio(alpha, Power, flux_lim):
    # Calculate the maximum redshift at which a source can be detected
    z_max = np.zeros(np.shape(Power))
    z_range = np.arange(z_step,25,z_step)

    for index, L_source in enumerate(tqdm(Power)):
        alpha_source = alpha[index]

        # Calculate the luminosity corresponding to the flux limit as function of redshift
        L = L_z(z_range, alpha_source, flux_lim)

        # Find the z-indices where luminosity < luminosity of source
        L_max_ind = np.where((L > L_source))

        try:
            z_max[index] = np.min(z_range[L_max_ind])
        except:
            print(index)
            print(alpha_source)
    return z_max


def K_corr_opt(z):
    return (-2.58*z + 6.67*(z**2) - 5.73*(z**3) - 0.42*(z**4))/(1-2.36*z + 3.82*(z**2) -3.53*(z**3) + 3.35*(z**4))


def i_MAG_z(z, i_mag, g_i):
    dl_z = cosmo.luminosity_distance(z) * 1e6 / u.Mpc # in pc
    k_corr = calc_kcor('r',z,'g - r', g_i)
    i_MAG = i_mag - 5*np.log10(dl_z/10) - k_corr
    return i_MAG


def z_max_opt(mag_lim, i_MAG, g_i):
    # Calculate the maximum redshift at which a source can be detected
    z_max = np.zeros(np.shape(i_MAG))
    z_range = np.arange(z_step,5,z_step)
    bad_i = []
    for index, i_MAG_source in enumerate(tqdm(i_MAG)):
        # Calculate the luminosity corresponding to the flux limit as function of redshift
        i_MAG_test = i_MAG_z(z_range, mag_lim, g_i[index])
        
        # Find the z-indices where luminosity < luminosity of source
        i_max_ind = np.where((i_MAG_test < i_MAG_source))

        try:
            z_max[index] = np.max(z_range[i_max_ind])
        except:
            bad_i.append(index)
            z_max[index] = np.inf 
            continue
    print('# of bad z-max for optical:', len(name[bad_i]))
    return z_max


def local_correction(z):
    h=cosmo.H(0).value/100
    gamma_s=1.66
    s0=3.76/h
    dl_limit=20.0/h
    dl=cosmo.luminosity_distance(z).value
    corr = np.zeros(np.shape(dl))
    # print(corr[0])
    for i, dl_z in enumerate(dl):
        if dl_z<=dl_limit:
            xi=(s0/dl_z)**gamma_s
            corr[i]=1.0+3.0/(3.0-gamma_s)*xi
        else:
            corr[i]=1.0
    return corr


def calc_volume(z_max):
    d = cosmo.comoving_distance(z_max).value
    corr = local_correction(z_max)
    return 4/3 * np.pi * d**3 * corr # Mpc^3


def poisson_errors(lum_func, err_lum_func, lum_func_counts):
    # Compute log poisson errors
    log_err_lum_func = np.zeros((2,np.shape(lum_func)[0]))
    
    for i, n in enumerate(lum_func_counts):
        # Use 84% confidence upper and lower limits based on Poisson statistics, tabulated in Gehrels (1986)
        m=n+1
        lambda_u = m*(1-1/(9*m)+1/(3*np.sqrt(m)))**3
        lambda_l = n*(1-1/(9*n)-1/(3*np.sqrt(n)))**3
        err_high = lum_func[i]*lambda_u/n - lum_func[i]
        err_low = lum_func[i] - lum_func[i]*lambda_l/n
        log_err_high = (1/np.log(10)) * (err_high / lum_func[i])
        log_err_low = (1/np.log(10)) * (err_low / lum_func[i])
        log_err_lum_func[:,i] = [log_err_low, log_err_high]
    return log_err_lum_func


cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.27)
z_step = 0.0001

#Import the data with SDSS photometry and redshifts
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_sample_redshift_SDSS.fits')
tbdata = hdulist[1].data
hdulist.close()


z = tbdata['z_best']

#We only us sources in our LF sample (LoLSS flux) that have redshift available
ind_LF = np.where((z>0.)&(tbdata['LoLSS_flux']>0.))

name = tbdata['LoTSS_name_1'][ind_LF]
alpha_low = tbdata['alpha_low_LoLSS'][ind_LF]
e_alpha_low = tbdata['e_alpha_low_LoLSS'][ind_LF]
alpha_high = tbdata['alpha_high_LoLSS'][ind_LF]
e_alpha_high = tbdata['e_alpha_high_LoLSS'][ind_LF]
S_1400 = tbdata['NVSS_flux'][ind_LF]
# norm_high = tbdata['norm_high'][ind_LF]
S_144 = tbdata['LoTSS_flux'][ind_LF]
i_mag = tbdata['i'][ind_LF]
g_mag = tbdata['g'][ind_LF]
r_mag = tbdata['r'][ind_LF]
z = z[ind_LF]

# Compute the cosmic luminosity distance (in m)
dl = cosmo.luminosity_distance(z) # in Mpc

# Do the k-correction
k_corr = np.power((1 + z),-(1+alpha_high))

# Compute power at 144MHz to check what happens without extrapolating
Power_144 = (4*np.pi*pow(dl,2)*S_144 * k_corr).to(u.m**2) * (10**-26)

# Remove sources with nan or 0 values, and redshifts > 3
Power_144_ind = np.where((np.isfinite(Power_144) & (Power_144 != 0) & (z < 3.)\
                          & (i_mag > 0.) & (r_mag > 0.) & (g_mag > 0.)))
    
Power_144 = Power_144[Power_144_ind].value
z_144 = z[Power_144_ind]
alpha_low_144, alpha_high_144 = alpha_low[Power_144_ind], alpha_high[Power_144_ind]
e_alpha_low_144, e_alpha_high_144 = e_alpha_low[Power_144_ind], e_alpha_high[Power_144_ind]
i_mag, g_mag, r_mag = i_mag[Power_144_ind], g_mag[Power_144_ind], r_mag[Power_144_ind]

# The ordinary radio sources have a different completeness limit compared to PS sources
# Use flux density limit at 144MHz & 5GHz from the extrapolated LoLSS flux density limit
# Rescaled from 95%LoLSS limit at 11mJy with alpha_low=-1.04... as found from 95% limit to alpha_low, see source_counts code

alpha_low_95 = -1.0456912517547607
flux_lim_144 = (11/1000) * (144/54) ** (alpha_low_95)

#TODO: is this still up to date with out DR?
mag_lim_i = 21.3
#mag_lim_r = 17.77

# Calculate absolute i_band magnitude for all sources
g_i = g_mag - i_mag
g_r = g_mag - r_mag
i_MAG = i_MAG_z(z_144, i_mag, g_i)
#r_MAG = i_MAG_z(z_144, r_mag, g_r)

# Calculate maximum redshift at which a source can be detected
z_max_opt = z_max_opt(mag_lim_i, i_MAG, g_i)
z_max_144 = z_max_radio(alpha_high_144, Power_144, flux_lim_144)


# define PS sample, which will have a different limit
ind_peaked = np.where((alpha_low_144 > e_alpha_low_144)\
                & (alpha_high_144 < - e_alpha_high_144))

#We handle the PS sources in a different way than the normal sources
flux_lim_144_PS = 0.013 #Jy, extrapolated from 95% complete LoLSS at 11mJy, with alpha=0.18
z_max_144_PS = z_max_radio(alpha_high_144[ind_peaked], Power_144[ind_peaked], flux_lim_144_PS)
z_max_144[ind_peaked] = z_max_144_PS

# Check how many sources have redshifts higher than z_max --> tells us if limit is too conservative
wrong_z_144 = np.where((z_144 > z_max_144))
print('at 144GHz, # sources with too high z:',np.count_nonzero(wrong_z_144), "out of", (len(z_144)))

# Calculate corresponding V_max
V_max_radio = calc_volume(z_max_144) # Mpc^3
V_max_opt = calc_volume(z_max_opt)
# V_max_radio_PS = calc_volume(z_max_144_PS)

# Saving to a fits file
col1  = fits.Column(name='LoTSS_name', format = '34A', array = name[Power_144_ind])
col2  = fits.Column(name='RA', format = 'E', array = tbdata['RA'][ind_LF][Power_144_ind])
col3  = fits.Column(name='Dec', format = 'E', array = tbdata['DEC'][ind_LF][Power_144_ind])
col4  = fits.Column(name='alpha_low', format = 'E', array = alpha_low_144)
col5  = fits.Column(name='alpha_high', format = 'E', array = alpha_high_144)
col19  = fits.Column(name='e_alpha_low', format = 'E', array = e_alpha_low_144)
col20  = fits.Column(name='e_alpha_high', format = 'E', array = e_alpha_high_144)
col6 = fits.Column(name='LoTSS_flux', format = 'E', array = S_144[Power_144_ind])
col7 = fits.Column(name='e_LoTSS_flux', format = 'E', array = tbdata['e_LoTSS_flux'][ind_LF][Power_144_ind])
col8 = fits.Column(name='LoLSS_flux', format = 'E', array = tbdata['LoLSS_flux'][ind_LF][Power_144_ind])
col9 = fits.Column(name='e_LoLSS_flux', format = 'E', array = tbdata['e_LoLSS_flux'][ind_LF][Power_144_ind])
col10 = fits.Column(name='NVSS_flux', format = 'E', array = S_1400[Power_144_ind])
col11 = fits.Column(name='e_NVSS_flux', format = 'E', array = tbdata['e_NVSS_flux'][ind_LF][Power_144_ind])
col12 = fits.Column(name='z_reported', format = 'E', array = z_144)
# col13 = fits.Column(name='norm_high', format = 'E', array = norm_high[Power_144_ind])
col14 = fits.Column(name='z_max_144', format = 'F20', array = z_max_144)
col15 = fits.Column(name='V_max_radio', format = 'F20', array = V_max_radio)
col16 = fits.Column(name='Power_144', format = 'F20', array = Power_144)
# col17 = fits.Column(name='weights_144', format='F20', array = weights_144)
col18 = fits.Column(name='V_max_opt', format = 'F20', array = V_max_opt)

cols = fits.ColDefs([col1, col2, col3, col4, col19, col5,col20, col6, col7, col8, col9, col10, col11, col12, col14, col15, col16, col18])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/MPS_lum_func.fits', overwrite = True)