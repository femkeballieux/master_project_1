# !/bin/python
# Plotting power distribution

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats

def poisson_errors(lum_func, err_lum_func, lum_func_counts):
    # Compute log poisson errors
    log_err_lum_func = np.zeros((2,np.shape(lum_func)[0]))
    
    for i, n in enumerate(lum_func_counts):
        if n < 5:
            # Use 84% confidence upper and lower limits based on Poisson statistics, tabulated in Gehrels (1986)
            m=n+1
            lambda_u = m*(1-1/(9*m)+1/(3*np.sqrt(m)))**3
            lambda_l = n*(1-1/(9*n)-1/(3*np.sqrt(n)))**3
            err_high = lum_func[i]*lambda_u/n - lum_func[i]
            err_low = lum_func[i] - lum_func[i]*lambda_l/n
            log_err_high = (1/np.log(10)) * (err_high / lum_func[i])
            log_err_low = (1/np.log(10)) * (err_low / lum_func[i])
            log_err_lum_func[:,i] = [log_err_low, log_err_high]
        else:
            max_val = lum_func[i] + err_lum_func[i]
            min_val = lum_func[i] - err_lum_func[i]
            log_min_val = np.log10(min_val)
            log_max_val = np.log10(max_val)
            log_err_high = log_max_val - np.log10(lum_func[i])
            log_err_low = -log_min_val + np.log10(lum_func[i])
            log_err_lum_func[:,i] = [log_err_low, log_err_high]
    return log_err_lum_func

plt.style.use('style.mplstyle')

hdulist = fits.open('data/lum_func_weighted_SDSS.fits')
tbdata = hdulist[1].data
hdulist.close()

name = tbdata['LoTSS_name']
alpha_low = tbdata['alpha_low']
alpha_high = tbdata['alpha_high']
S_1400 = tbdata['NVSS_flux']
norm_high = tbdata['norm_high']
S_144 = tbdata['LoTSS_flux']
z_max_144 = tbdata['z_max_144']
V_max_opt = tbdata['V_max_opt']
V_max_radio = tbdata['V_max_radio']
Power_144 = tbdata['Power_144']
z = tbdata['z_reported']
weights_144 = tbdata['weights_144']

# Compute 10log of power
log_Power_144 = np.log10(Power_144)


# define PS sample
ind_peaked = np.where((alpha_low >= 0.1) & (alpha_high <= 0.0))
V_max_radio_PS = V_max_radio[ind_peaked]


# LoLSS fractional area
str_to_deg2 = (180/np.pi)**2
frac_area_LoLSS = 740 / (4*np.pi * str_to_deg2)  #sr, LoLSS sky coverage (740 deg2)


# Calculate total volumes
V_max_144 = np.amin([V_max_radio,V_max_opt], axis = 0) * frac_area_LoLSS
V_max_144_PS = np.amin([V_max_radio_PS,V_max_opt[ind_peaked]], axis = 0) * frac_area_LoLSS
print(V_max_144)

# Compute luminosity function for 144MHz
# Define luminosity bins 
bin_edges_144 = np.arange(23., 29., 0.2)
n_144, bin_edges_144 = np.histogram(log_Power_144, bins=bin_edges_144)
print("# sources in each bin", n_144)
log_xerr_144 = np.diff(bin_edges_144) / 2


# Define luminosity bins PS sample
bin_edges_144_PS = np.arange(23., 29.4, 0.8)
n_144_PS, bin_edges_144_PS = np.histogram(log_Power_144[ind_peaked], bins=bin_edges_144_PS)
log_xerr_144_PS = np.diff(bin_edges_144_PS) / 2
print('# PS sources in each bin', n_144_PS)


# Find central bins & log x-error for plotting
central_bins_144  = np.array([(j+i)/2. for i, j in zip(bin_edges_144[:-1], bin_edges_144[1:])])
central_bins_144_PS = np.array([(j+i)/2. for i, j in zip(bin_edges_144_PS[:-1], bin_edges_144_PS[1:])])


# Compute the weighted luminosity function
# lum_func_144_w = (1/(np.diff(bin_edges_144))) * stats.binned_statistic(log_Power_144, weights_144/(V_max_144_w), 'sum', bins=bin_edges_144)[0]
# err_lum_func_144_w = np.sqrt( (1/(np.diff(bin_edges_144)))**2 * stats.binned_statistic(log_Power_144, (weights_144/V_max_144_w)**2, 'sum', bins=bin_edges_144)[0])
# # err_lum_func_weights = (1/(np.diff(bin_edges_144)))**2 * stats.binned_statistic(log_Power_144, (err_weights/V_max_144)**2, 'sum', bins=bin_edges_144)[0]
# log_lum_func_144_w = np.log10(lum_func_144_w)
# lum_func_counts = stats.binned_statistic(log_Power_144, weights_144/(V_max_144_w), 'count', bins=bin_edges_144)[0]
# log_err_lum_func_144_w = poisson_errors(lum_func_144_w, err_lum_func_144_w, lum_func_counts)
# print("log LF weighted master sample:", log_lum_func_144_w)


# Compute the SDSS luminosity function
lum_func_144_SDSS = (1/(np.diff(bin_edges_144))) * stats.binned_statistic(log_Power_144, 1/(V_max_144), 'sum', bins=bin_edges_144)[0]
err_lum_func_144_SDSS = np.sqrt( (1/(np.diff(bin_edges_144)))**2 * stats.binned_statistic(log_Power_144, (1/V_max_144)**2, 'sum', bins=bin_edges_144)[0])
err_lum_func_weights = (1/(np.diff(bin_edges_144)))**2 * stats.binned_statistic(log_Power_144, (1/V_max_144)**2, 'sum', bins=bin_edges_144)[0]
log_lum_func_144_SDSS = np.log10(lum_func_144_SDSS)
lum_func_counts_SDSS = stats.binned_statistic(log_Power_144, 1/(V_max_144), 'count', bins=bin_edges_144)[0]
log_err_lum_func_144_SDSS = poisson_errors(lum_func_144_SDSS, err_lum_func_144_SDSS, lum_func_counts_SDSS)
print("log LF SDSS master sample:", log_lum_func_144_SDSS)


# Compute the weighted luminosity function - PS sources
# lum_func_144_PS_w = (1/(np.diff(bin_edges_144_PS))) * stats.binned_statistic(log_Power_144[ind_peaked], weights_144[ind_peaked]/V_max_144_PS_w, 'sum', bins=bin_edges_144_PS)[0]
# err_lum_func_144_PS_w = np.sqrt( (1/(np.diff(bin_edges_144_PS))**2) * stats.binned_statistic(log_Power_144[ind_peaked], (weights_144[ind_peaked]/V_max_144_PS_w)**2, 'sum', bins=bin_edges_144_PS)[0])
# log_lum_func_144_PS_w = np.log10(lum_func_144_PS_w)
# lum_func_counts_PS = stats.binned_statistic(log_Power_144[ind_peaked], weights_144[ind_peaked]/(V_max_144_PS_w), 'count', bins=bin_edges_144_PS)[0]
# log_err_lum_func_144_PS_w = poisson_errors(lum_func_144_PS_w, err_lum_func_144_PS_w, lum_func_counts_PS)


# Compute the SDSS luminosity function - PS sources
lum_func_144_PS_SDSS = (1/(np.diff(bin_edges_144_PS))) * stats.binned_statistic(log_Power_144[ind_peaked], 1/V_max_144_PS, 'sum', bins=bin_edges_144_PS)[0]
err_lum_func_144_PS_SDSS = np.sqrt( (1/(np.diff(bin_edges_144_PS))**2) * stats.binned_statistic(log_Power_144[ind_peaked], (1/V_max_144_PS)**2, 'sum', bins=bin_edges_144_PS)[0])
log_lum_func_144_PS_SDSS = np.log10(lum_func_144_PS_SDSS)
lum_func_counts_PS_SDSS = stats.binned_statistic(log_Power_144[ind_peaked], 1/(V_max_144_PS), 'count', bins=bin_edges_144_PS)[0]
log_err_lum_func_144_PS_SDSS = poisson_errors(lum_func_144_PS_SDSS, err_lum_func_144_PS_SDSS, lum_func_counts_PS_SDSS)


# # LoTSS Luminosity function AGN values Sabater 2019
# log_L_central_bins_Sabater = np.arange(21.25,26.95,0.30)
# log_lum_func_Sabater = np.array([-3.75,-3.69,-3.73,-3.93,-4.03,-4.13,-4.38,-4.47,-4.65,-4.78,-5.01,-5.06,-5.12,-5.37,-5.55,-5.92,-6.07,-6.47,-6.55])
# log_lum_func_err_Sabater = np.array([[0.17,0.07,0.04,0.03,0.03,0.03,0.03,0.03,0.06,0.06,0.06,0.05,0.06,0.08,0.09,0.11,0.14,0.23,0.41],\
#                                     [0.12,0.06,0.04,0.03,0.03,0.03,0.03,0.03,0.05,0.05,0.05,0.05,0.05,0.06,0.07,0.09,0.10,0.15,0.20]])
# yerr_Sabater = np.log(10) * 10**(log_lum_func_Sabater) * log_lum_func_err_Sabater


# Heckman & Best lum function
P = 10**(np.linspace(19,np.max(log_Power_144)+2,1000))
P_0_14 = 10**(24.95)
alpha_extrap = -0.7
P_0 = P_0_14 * (144/1400)**alpha_extrap
a = 0.42
b = 1.66
rho_0 = 10**(-5.33)
rho = rho_0/( (P/P_0)**a + (P/P_0)**b)


# Plot luminosity function
fig = plt.figure(0)
gs = plt.GridSpec(1,1)
ax = plt.subplot(gs[0])


# ax.errorbar((central_bins_144), log_lum_func_144_w, xerr=log_xerr_144, yerr = log_err_lum_func_144_w, marker = 'o', linestyle='none', capthick = 2, label = 'Weighted Master Sample (z < 3)')
# ax.errorbar((central_bins_144_PS), log_lum_func_144_PS_w, xerr=log_xerr_144_PS, yerr = log_err_lum_func_144_PS_w, marker = 's', linestyle = 'none', capthick = 2, label = 'Weighted PS sample (z < 3)')
ax.errorbar((central_bins_144), log_lum_func_144_SDSS, xerr=log_xerr_144, yerr = log_err_lum_func_144_SDSS, marker = 'D', linestyle = 'none', capthick = 2, label = 'SDSS Master Sample (z < 3)')
ax.errorbar((central_bins_144_PS), log_lum_func_144_PS_SDSS, xerr=log_xerr_144_PS, yerr = log_err_lum_func_144_PS_SDSS, marker = 's', linestyle = 'none', capthick = 2, label = 'SDSS PS sample (z < 3)')
# ax.errorbar(log_L_central_bins_Sabater, log_lum_func_Sabater, xerr = 0.15, yerr = log_lum_func_err_Sabater, linestyle='none', marker = 's', markersize = 5, capthick = 2, label = 'LoTSS AGN - Sabater (2019) (150MHz)')
ax.plot(np.log10(P),np.log10(rho), label = 'Best et al. (2014) z < 0.3 AGN model', color = 'grey')


ax.set_xlim(22,np.max(log_Power_144)+1)
ax.set_ylim(-11.5,-3.5)
ax.tick_params(axis='both',which='both',top=True,right=True)

ax.set_xlabel(r'$\log_{10}$(L$_{144 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')
ax.set_ylabel(r'$\log_{10} ( \rho$ / Mpc$^{-3} \Delta \log_{10}$L$^{-1}$ )')
# ax.legend(loc = 'lower left')

# plt.show()
plt.savefig("plots/luminosity_function_144MHz_SDSS_k_corr.png")
plt.close()