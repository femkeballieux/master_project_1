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
    z_range = np.arange(z_step,20,z_step)
    z_range_bigger = np.arange(z_step,150,z_step)

    for index, L_source in enumerate(tqdm(Power)):
        alpha_source = alpha[index]

        # Calculate the luminosity corresponding to the flux limit as function of redshift
        L = L_z(z_range, alpha_source, flux_lim)

        # Find the z-indices where luminosity < luminosity of source
        L_max_ind = np.where((L > L_source))

        try:
            #TODO: check of dit nog nodig is
            z_max[index] = np.min(z_range[L_max_ind])
        except:
            print('Trying with larger z_max')
            L_bigger = L_z(z_range_bigger, alpha_source, flux_lim)
            L_max_ind_bigger = np.where((L_bigger > L_source))
            try:
                z_max[index] = np.min(z_range_bigger[L_max_ind_bigger])
                print('Went right')
            except:
                print('Still nope')
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
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_redshift_SDSS.fits')
tbdata = hdulist[1].data
hdulist.close()

z = tbdata['z_best']

#We only us sources in our HF sample (LoLSS flux) that have redshift available
ind_HF = np.where((z>0.)&(tbdata['VLASS_flux']>0.))

name = tbdata['LoTSS_name'][ind_HF]
alpha_low = tbdata['alpha_low_VLASS'][ind_HF]
e_alpha_low = tbdata['e_alpha_low_VLASS'][ind_HF]
alpha_high = tbdata['alpha_high_VLASS'][ind_HF]
e_alpha_high = tbdata['e_alpha_high_VLASS'][ind_HF]
S_1400 = tbdata['NVSS_flux'][ind_HF]
# norm_high = tbdata['norm_high'][ind_HF]
S_144 = tbdata['LoTSS_flux'][ind_HF]
i_mag = tbdata['i'][ind_HF]
g_mag = tbdata['g'][ind_HF]
r_mag = tbdata['r'][ind_HF]
z = z[ind_HF]

# Compute the cosmic luminosity distance (in m)
dl = cosmo.luminosity_distance(z) # in Mpc

# Do the k-correction
k_corr = np.power((1 + z),-(1+alpha_high))

# Compute power at 1400MHz to check what happens without extrapolating
Power_1400 = (4*np.pi*pow(dl,2) * S_1400 * k_corr).to(u.m**2) * (10**-26)

# Remove sources with nan or 0 values, and redshifts > 3
Power_1400_ind = np.where((np.isfinite(Power_1400) & (Power_1400 != 0) & (z < 3.)\
                          & (i_mag > 0.) & (r_mag > 0.) & (g_mag > 0.)))
    
Power_1400 = Power_1400[Power_1400_ind].value
z_1400 = z[Power_1400_ind]
alpha_low_1400, alpha_high_1400 = alpha_low[Power_1400_ind], alpha_high[Power_1400_ind]
e_alpha_low_1400, e_alpha_high_1400 = e_alpha_low[Power_1400_ind], e_alpha_high[Power_1400_ind]
i_mag, g_mag, r_mag = i_mag[Power_1400_ind], g_mag[Power_1400_ind], r_mag[Power_1400_ind]

#This is the limit in janskys for NVSS
flux_lim_1400 = (3.24/1000)

#TODO: is this still up to date with out DR?
mag_lim_i = 21.3
#mag_lim_r = 17.77

# Calculate absolute i_band magnitude for all sources
g_i = g_mag - i_mag
g_r = g_mag - r_mag
i_MAG = i_MAG_z(z_1400, i_mag, g_i)
#r_MAG = i_MAG_z(z_1400, r_mag, g_r)

# Calculate maximum redshift at which a source can be detected
z_max_1400 = z_max_radio(alpha_high_1400, Power_1400, flux_lim_1400)
z_max_opt = z_max_opt(mag_lim_i, i_MAG, g_i)

# define PS sample
# ind_peaked = np.where((alpha_low_1400 > e_alpha_low_1400)\
#                 & (alpha_high_1400 < - e_alpha_high_1400))

#We do not handle the PS sources in a different way than the normal sources
# flux_lim_1400_PS = 0.013 #Jy, extrapolated from 95% complete LoLSS at 11mJy, with alpha=0.18
# z_max_1400_PS = z_max_radio(alpha_high_1400[ind_peaked], Power_1400[ind_peaked], flux_lim_1400_PS)
# z_max_1400[ind_peaked] = z_max_1400_PS

# Check how many sources have redshifts higher than z_max --> tells us if limit is too conservative
#TODO: NOPE, hier gaat iets goed mis. # sources with too high z: 31873 out of 31890
wrong_z_1400 = np.where((z_1400 > z_max_1400))
print('at 1400GHz, # sources with too high z:',np.count_nonzero(wrong_z_1400), "out of", (len(z_1400)))

# Calculate corresponding V_max
V_max_radio = calc_volume(z_max_1400) # Mpc^3
V_max_opt = calc_volume(z_max_opt)
# V_max_radio_PS = calc_volume(z_max_1400_PS)

# Saving to a fits file
col1  = fits.Column(name='LoTSS_name', format = '34A', array = name[Power_1400_ind])
col2  = fits.Column(name='RA', format = 'E', array = tbdata['RA'][ind_HF][Power_1400_ind])
col3  = fits.Column(name='Dec', format = 'E', array = tbdata['DEC'][ind_HF][Power_1400_ind])
col4  = fits.Column(name='alpha_low', format = 'E', array = alpha_low_1400)
col5  = fits.Column(name='alpha_high', format = 'E', array = alpha_high_1400)
col19  = fits.Column(name='e_alpha_low', format = 'E', array = e_alpha_low_1400)
col20  = fits.Column(name='e_alpha_high', format = 'E', array = e_alpha_high_1400)
col6 = fits.Column(name='LoTSS_flux', format = 'E', array = S_144[Power_1400_ind])
col7 = fits.Column(name='e_LoTSS_flux', format = 'E', array = tbdata['e_LoTSS_flux'][ind_HF][Power_1400_ind])
col8 = fits.Column(name='VLASS_flux', format = 'E', array = tbdata['VLASS_flux'][ind_HF][Power_1400_ind])
col9 = fits.Column(name='e_VLASS_flux', format = 'E', array = tbdata['e_VLASS_flux'][ind_HF][Power_1400_ind])
col10 = fits.Column(name='NVSS_flux', format = 'E', array = S_1400[Power_1400_ind])
col11 = fits.Column(name='e_NVSS_flux', format = 'E', array = tbdata['e_NVSS_flux'][ind_HF][Power_1400_ind])
col12 = fits.Column(name='z_reported', format = 'E', array = z_1400)
# col13 = fits.Column(name='norm_high', format = 'E', array = norm_high[Power_1400_ind])
col14 = fits.Column(name='z_max_1400', format = 'F20', array = z_max_1400)
col15 = fits.Column(name='V_max_radio', format = 'F20', array = V_max_radio)
col16 = fits.Column(name='Power_1400', format = 'F20', array = Power_1400)
# col17 = fits.Column(name='weights_1400', format='F20', array = weights_1400)
col18 = fits.Column(name='V_max_opt', format = 'F20', array = V_max_opt)

cols = fits.ColDefs([col1, col2, col3, col4, col19, col5,col20, col6, col7, col8, col9, col10, col11, col12, col14, col15, col16, col18])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/GPS_lum_func_1400.fits', overwrite = True)