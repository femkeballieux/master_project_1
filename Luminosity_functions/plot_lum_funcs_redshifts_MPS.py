"""
Code written by Martje Slob, edited by Femke Ballieux. 

whenever there is referred to index=5, I believe this plots the entire sample from 0<z<3.
We dont use it, so could be removed
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.gridspec as gridspec
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


def local_correction(z):
    h=cosmo.H(0).value/100
    gamma_s=1.66
    s0=3.76/h
    dl_limit=20.0/h
    dl=cosmo.luminosity_distance(z).value
    corr = np.zeros(np.shape(dl))
    if dl<=dl_limit:
        xi=(s0/dl)**gamma_s
        corr=1.0+3.0/(3.0-gamma_s)*xi
    else:
        corr=1.0

    return corr


def calc_volume(z_max):
    d = cosmo.comoving_distance(z_max).value
    corr = local_correction(z_max)
    return 4/3*np.pi*d**3 * corr # Mpc^3


def poisson_errors(lum_func, err_lum_func, lum_func_counts):
    # Compute log poisson errors
    log_err_lum_func = np.zeros((2,np.shape(lum_func)[0]))
    
    for i, n in enumerate(lum_func_counts):
        if n < 6:
            # Use 84% confidence upper and lower limits based on Poisson statistics, tabulated in Gehrels (1986)
            m=n+1
            lambda_u = m*(1-1/(9*m)+1/(3*np.sqrt(m)))**3
            lambda_l = n*(1-1/(9*n)-1/(3*np.sqrt(n)))**3
            err_low = lum_func[i]*lambda_u/n - lum_func[i]
            err_high = lum_func[i] - lum_func[i]*lambda_l/n
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

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.27)

plt.style.use('/net/vdesk/data2/bach1/ballieux/master_project_1/code/style.mplstyle')

# LoLSS fractional area
str_to_deg2 = (180/np.pi)**2
frac_area_LoLSS = 650 / (4*np.pi * str_to_deg2)  #sr, LoLSS sky coverage (650 deg2)

# Define different redshift bins, all bins have ~100+ sources
redshift_bins = [0.,0.1,0.5,1.,1.5,3.,3.]

idx1 = [0,0,1,1,2,2]
idx2 = [0,1,0,1,0,1]

#Cutoff where it becomes incmplete
idx_cutoff=[0,1,6,1,3]
idx_cutoff_PS = [0,0,1,0,0]

bin_size = [0.3,0.25,0.25,0.4,0.5]
bin_size_PS = [0.7,0.9,0.7,0.9,1.1]



'''Importing literature LFs'''
# Heckman & Best lum function
y_max = 28.
P = 10**(np.linspace(19,y_max+2,1000))
P_0_14 = 10**(24.95)
alpha_extrap = -0.7
P_0 = P_0_14 * (144/1400)**alpha_extrap
a = 0.42
b = 1.66
rho_0 = 10**(-5.33)
rho = rho_0/( (P/P_0)**a + (P/P_0)**b) 

# Pracy lum function
C_SF = 10**(-2.83)
P_0_SF = 10**(21.18)
P_0_SF = P_0_SF * (144/1400)**alpha_extrap
a_SF = 1.02
b_SF = 0.60  
rho_SF = C_SF*((P/P_0_SF)**(1-a_SF) * np.exp(-0.5*(np.log10(1+P/P_0_SF)/b_SF)**2))

# Williams 2018:
log_P_williams = np.arange(24.25,28.25,0.5)
z_05_1 = [-3.86,-4.71,-5.22,-5.71,-6.02,-6.29,-6.84,0.]
z_05_1_err = [[0.05,0.04,0.06,0.11,0.16,0.17,0.38,0.],[0.05,0.03,0.06,0.10,0.14,0.20,0.38,0.]]
z_1_15 = [0,-4.15,-5.23,-5.71,-6.09,-6.50,-6.87,-7.04]
z_1_15_err = [[0.,0.11,0.05,0.10,0.13,0.21,0.27,0.49],[0.,0.11,0.05,0.09,0.13,0.18,0.37,0.49]]
z_15_2 = [0,0,-5.04,-5.78,-6.21,-6.72,-7.42,-7.42]
z_15_2_err = [[0.,0.,0.10,0.10,0.14,0.26,0.70,0.69],[0.,0.,0.09,0.09,0.15,0.26,0.55,0.54]]

#Kondapally 2022
LF_Konda_05_1 = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/Konda_LF/LF_radio_excess_AGN_0.7_z_1.0_tab.csv', delimiter=',', skiprows = 1, usecols = (0,1,2,3))
LF_Konda_1_15 = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/Konda_LF/LF_radio_excess_AGN_1.3_z_1.7_tab.csv', delimiter=',', skiprows = 1, usecols = (0,1,2,3))
LF_Konda_15_2 = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/Konda_LF/LF_radio_excess_AGN_1.7_z_2.1_tab.csv', delimiter=',', skiprows = 1, usecols = (0,1,2,3))
# LF_Konda_2_25 = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/Konda_LF/LF_LERG_AGN_2_z_2.5_tab_quies.csv', delimiter=',', skiprows = 1, usecols = (0,1,2,3))

"""Literature Slob2022 
(Fuck my life met dit intypen) 
Take only those from SDSS selected sample
"""
log_L_SLob_01=[23.4,24.33] #Log P
log_L_SLob_01_err=[0.4,0.5] #Log P error
log_phi_slob_01=[-4.71,-5.58] #log_lumfunc
log_phi_slob_01_err=[[0.09,0.11],[0.12,0.15]] #Upper, lower

log_L_SLob_01_PS=[23.71] #Log P
log_L_SLob_01_err_PS=[0.5] #Log P error
log_phi_slob_01_PS=[-5.31] #log_lumfunc
log_phi_slob_01_err_PS=[[0.2],[0.38]] #Upper, lower

log_L_SLob_01_05=[23.98,24.38,24.78,25.18,25.58,25.98,26.38,26.99] #Log P
log_L_SLob_01_05_err=[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4] #Log P error
log_phi_slob_01_05=[-5.64,-5.9,-6.09,-6.09,-6.4,-6.78,-6.97,-7.62] #log_lumfunc
log_phi_slob_01_05_err=[[0.11,0.08,0.05,0.04,0.05,0.07,0.08,0.12],[0.15,0.09,0.06,0.04,0.05,0.08,0.10,0.17]] #Upper, lower

log_L_SLob_01_05_PS=[25.16,26.56] #Log P
log_L_SLob_01_05_err_PS=[0.7,0.7] #Log P error
log_phi_slob_01_05_PS=[-7.01,8.01] #log_lumfunc
log_phi_slob_01_05_err_PS=[[0.16,0.14],[0.26,0.21]] #Upper, lower

log_L_SLob_05_1=[25.34,25.74,26.14,26.54,26.94,27.34,27.92] #Log P
log_L_SLob_05_1_err=[0.2,0.2,0.2,0.2,0.2,0.2,0.4] #Log P error
log_phi_slob_05_1=[-6.34,-6.32,-6.39,-6.75,-7,-7.54,-7.8] #log_lumfunc
log_phi_slob_05_1_err=np.array([[0.05,0.04,0.03,0.04,0.07,0.09,0.1],[0.05,0.04,0.03,0.04,0.08,0.111,0.13]]) #Upper, lower

log_L_SLob_05_1_PS=[25.94,26.54,27.14,27.74] #Log P
log_L_SLob_05_1_err_PS=[0.3,0.3,0.3,0.3] #Log P error
log_phi_slob_05_1_PS=[-7.24,-7.92,-8.2,-8.38] #log_lumfunc
log_phi_slob_05_1_err_PS=np.array([[0.14,0.12,0.13,0.24],[0.21,0.17,0.19,0.42]]) #Upper, lower

log_L_SLob_1_15=[25.99,26.49,26.99,27.49,28.07] #Log P
log_L_SLob_1_15_err=[0.2,0.2,0.2,0.2,0.3] #Log P error
log_phi_slob_1_15=[-6.5,-6.45,-6.15,-6.82,-7.81] #log_lumfunc
log_phi_slob_1_15_err=[[0.12,0.10,0.22,0.12,0.15],[0.16,0.12,0.46,0.16,0.22]] #Upper, lower

log_L_SLob_1_15_PS=[26.81,27.61] #Log P
log_L_SLob_1_15_err_PS=[0.4,0.4] #Log P error
log_phi_slob_1_15_PS=[-7.7,-7.97] #log_lumfunc
log_phi_slob_1_15_err_PS=[[0.18,0.24],[0.31,0.42]] #Upper, lower

log_L_SLob_15_3=[26.51,27.11,27.71,28.44] #Log P
log_L_SLob_15_3_err=[0.3,0.3,0.3,0.4] #Log P error
log_phi_slob_15_3=[-7.1,-6.98,-7.56,-7.9] #log_lumfunc
log_phi_slob_15_3_err=[[0.1,0.09,0.09,0.13],[0.13,0.11,0.12,0.18]] #Upper, lower

log_L_SLob_15_3_PS=[27.1,28.1] #Log P
log_L_SLob_15_3_err_PS=[0.5,0.5] #Log P error
log_phi_slob_15_3_PS=[-7.89,-9.04] #log_lumfunc
log_phi_slob_15_3_err_PS=[[0.15,0.16],[0.24,0.25]] #Upper, lower


#Start the figure

fig, axs = plt.subplots(3, 2, figsize=(20,23), sharex = True, sharey = True)
plt.subplots_adjust(wspace=0, hspace=0)
ax6 = axs[2, 1]
ax6.set_visible(False)

fig2, ax1 = plt.subplots(figsize=(10,7.5))

for index, z_bin_min in enumerate(redshift_bins[:-2]):
    hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/MPS_lum_func.fits')
    tbdata = hdulist[1].data
    hdulist.close()

    alpha_low = tbdata['alpha_low']
    alpha_high = tbdata['alpha_high']
    e_alpha_low = tbdata['e_alpha_low']
    e_alpha_high = tbdata['e_alpha_high']
    z_max_144 = tbdata['z_max_144']
    V_max_opt = tbdata['V_max_opt']
    V_max_radio = tbdata['V_max_radio']
    Power_144 = tbdata['Power_144']
    z = tbdata['z_reported']

    # Calculate total volume
    V_max_144 = np.amin([V_max_radio,V_max_opt], axis = 0) * frac_area_LoLSS

    marker_face = 'k'
    marker_face_PS = 'crimson'
    label_name = 'SDSS'
    m_i = '21.3'
    
    # Calculate max volume of redshift range
    if index < 5:
        V_max_limit = calc_volume(redshift_bins[index+1]) * frac_area_LoLSS
    elif index == 5:
        V_max_limit = calc_volume(3.) * frac_area_LoLSS # for full z ranges, use z_max = 3
    
    # Define values for this redshift range
    if index == 5:
        #I believe this is only when plotting the entire thing function, so no redshift bins
        z_idx = np.where((z > 0.) & (z < 3.))
    elif index < len(redshift_bins) - 1:
        z_idx = np.where((z > z_bin_min) & (z < redshift_bins[index+1]) & (np.log10(Power_144) > 23.))
    else:
        break

    # Make array with correct size to find min of V_max 
    V_max_i = np.empty(np.shape(z_idx[0]))
    V_max_i.fill(V_max_limit)

    # Define values for this redshift range
    Power_144_z = Power_144[z_idx]
    z_max_144_z = z_max_144[z_idx]
    V_max_144_z = np.amin([V_max_144[z_idx],V_max_i], axis = 0)
    alpha_low_z = alpha_low[z_idx]
    alpha_high_z = alpha_high[z_idx]
    e_alpha_low_z = e_alpha_low[z_idx]
    e_alpha_high_z = e_alpha_high[z_idx]
    
    
    # Also take minimum z of range into account
    # 5 is redshift bin 0-3, 0 is redshift bin 0-0.1. 
    if index > 0 and index < 5: # only when z_min != 0
        V_min = calc_volume(z_bin_min) * frac_area_LoLSS
        V_i = V_max_144_z - V_min
        V_idx = np.where(V_i > 0)
        V_i = V_i[V_idx]
        alpha_low_z = alpha_low_z[V_idx]
        alpha_high_z = alpha_high_z[V_idx]
        e_alpha_low_z = e_alpha_low_z[V_idx]
        e_alpha_high_z = e_alpha_high_z[V_idx]
        Power_144_z = Power_144_z[V_idx]
    else:
        #0.1 is ongeveer 0
        V_min = calc_volume(0.01) * frac_area_LoLSS
        V_i = V_max_144_z - V_min

    # Compute 10log of power
    log_Power_144 = np.log10(Power_144_z)

    # define PS sample
    ind_peaked = np.where((alpha_low_z > e_alpha_low_z)\
                          & (alpha_high_z < - e_alpha_high_z))
    z_max_144_PS = z_max_144_z[ind_peaked]
    V_i_PS = V_i[ind_peaked]

    # Define luminosity luminosity bins for 144MHz, use SDSS limits
    bin_size_i = bin_size[index]
    if index == 0 or index == 5:
        low_limit_144 = 23 # don't use the lowest luminosities, incomplete there
    else:
        low_limit_144 = np.min(log_Power_144)
    factor = 1
    if index in [1,5]:
        factor = 2
    bin_edges_144 = np.hstack((np.arange(low_limit_144, np.max(log_Power_144)-factor*bin_size_i, bin_size_i),np.max(log_Power_144)))
    n_144, bin_edges_144 = np.histogram(log_Power_144, bins=bin_edges_144)
    
    # Define luminosity bins PS sample
    if index == 5:
        bin_edges_144_PS = np.arange(23.,29.4, bin_size_PS[index]) 
    else:
        bin_edges_144_PS = np.arange(np.min(log_Power_144[ind_peaked]), np.max(log_Power_144[ind_peaked])+bin_size_PS[index], bin_size_PS[index]) 
    n_144_PS, bin_edges_144_PS = np.histogram(log_Power_144[ind_peaked], bins = bin_edges_144_PS)
    

    # Find central bins for plotting
    central_bins_144  = np.array([(j+i)/2. for i, j in zip(bin_edges_144[:-1], bin_edges_144[1:])])
    central_bins_144_PS = np.array([(j+i)/2. for i, j in zip(bin_edges_144_PS[:-1], bin_edges_144_PS[1:])])

    # Compute the luminosity function
    lum_func_144 = (1/(np.diff(bin_edges_144))) * stats.binned_statistic(log_Power_144, 1/V_i, 'sum', bins=bin_edges_144)[0]
    err_lum_func_144 = np.sqrt( (1/np.diff(bin_edges_144))**2 * stats.binned_statistic(log_Power_144, (1/V_i)**2, 'sum', bins=bin_edges_144)[0])
    log_lum_func_144 = np.log10(lum_func_144)
    lum_func_counts_144 = stats.binned_statistic(log_Power_144, log_Power_144, 'count', bins=bin_edges_144)[0]
    log_err_lum_func_144 = poisson_errors(lum_func_144, err_lum_func_144, lum_func_counts_144)
    

    # Compute the PS luminosity function
    lum_func_144_PS = (1/np.diff(bin_edges_144_PS)) * stats.binned_statistic(log_Power_144[ind_peaked], 1/V_i_PS, 'sum', bins=bin_edges_144_PS)[0]
    log_lum_func_144_PS = np.log10(lum_func_144_PS)
    err_lum_func_144_PS = np.sqrt((1/np.diff(bin_edges_144_PS))**2 * stats.binned_statistic(log_Power_144[ind_peaked], (1/V_i_PS)**2, 'sum', bins=bin_edges_144_PS)[0])
    bin_count_PS = stats.binned_statistic(log_Power_144[ind_peaked], log_Power_144[ind_peaked], 'count', bins=bin_edges_144_PS)[0] # Number of sources per bin
    log_err_lum_func_144_PS = poisson_errors(lum_func_144_PS, err_lum_func_144_PS, bin_count_PS)

    # print("# sources in each bin", lum_func_counts_144)
    log_xerr_144 = np.diff(bin_edges_144) / 2
    # print("# PS sources in each bin", bin_count_PS)
    log_xerr_144_PS = np.diff(bin_edges_144_PS) / 2


    ax = axs[idx1[index], idx2[index]]

    # Plot luminosity function
    ax.errorbar(central_bins_144[idx_cutoff[index]:], log_lum_func_144[idx_cutoff[index]:],\
                xerr=log_xerr_144[idx_cutoff[index]:], yerr = log_err_lum_func_144[:,idx_cutoff[index]:], c = 'k', \
                marker = 'o', markersize = 10, linestyle ='none', mfc = marker_face, capthick = 2, label = 'LF Master Sample')#', $m_i\, <$'+ m_i)
    ax.errorbar(central_bins_144_PS[idx_cutoff_PS[index]:], log_lum_func_144_PS[idx_cutoff_PS[index]:],\
                xerr=log_xerr_144_PS[idx_cutoff_PS[index]:], yerr = log_err_lum_func_144_PS[:,idx_cutoff_PS[index]:], c = 'crimson',\
                marker = '^',  capthick = 2, markersize = 10, mfc = marker_face_PS, linestyle = 'none', label = 'MPS sample')#', $m_i\, <$'+ m_i)
    

    ax.set_xlim(22.5,29)
    ax.set_ylim((-11.1,-3.5))


    
    # if LF_name == 'SDSS':
    ax.plot(np.log10(P), np.log10(rho), linewidth = 2, linestyle = '--', color = 'k', zorder = -10, label = 'Heckman & Best (2014)\nz < 0.3')
    ax.annotate(str(len(z_idx[0]))+' SDSS sources, ' + str(len(z_max_144_PS)) + ' SDSS PS sources', xy = (23,-10.2), xycoords = 'data', fontsize = 15)

    ax.annotate(str(z_bin_min)+' < z <' + str(redshift_bins[index+1]), xy=(23,-9.8), xycoords = 'data', fontsize = 17)

    if idx1[index] == 0:                    
        if idx1[index] == 0 or idx1[index] == 1:
            ax.tick_params(left = True, labelleft = True, right = True, bottom = True, labelbottom = False, top = True)
        if idx2[index] == 0:
            ax.errorbar(log_L_SLob_01, log_phi_slob_01, xerr=log_L_SLob_01_err, \
                        yerr=log_phi_slob_01_err, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='#1f78b4')
            ax.errorbar(log_L_SLob_01_PS, log_phi_slob_01_PS, xerr=log_L_SLob_01_err_PS, \
                        yerr=log_phi_slob_01_err_PS, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='mediumorchid')
        elif idx2[index] == 1:
            ax.errorbar(log_L_SLob_01_05, log_phi_slob_01_05, xerr=log_L_SLob_01_05_err, \
                        yerr=log_phi_slob_01_05_err, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='#1f78b4')
            ax.errorbar(log_L_SLob_01_05_PS, log_phi_slob_01_05_PS, xerr=log_L_SLob_01_05_err_PS, \
                        yerr=log_phi_slob_01_05_err_PS, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='mediumorchid')            
    if idx1[index] == 1:
        if idx2[index] == 0:
            ax.set_ylabel(r'$\log_{10} ( \Phi$ / Mpc$^{-3}$ dex$^{-1}$ )')
            # ax.errorbar(log_P_williams, z_05_1, xerr=0.25, yerr = z_05_1_err, marker = 's', c = 'dodgerblue',\
            #             markersize = 10, mfc='none', capthick = 2, linestyle='none')
            # ax.errorbar(LF_Konda_05_1[:,0], LF_Konda_05_1[:,1], xerr=0.15, yerr = [LF_Konda_05_1[:,2],LF_Konda_05_1[:,3]], marker = 'D',\
            #           c = 'forestgreen', mfc='none', markersize = 10, capthick = 2, linestyle='none')
            ax.errorbar(log_L_SLob_05_1[1:], log_phi_slob_05_1[1:], xerr=log_L_SLob_05_1_err[1:], \
                        yerr=log_phi_slob_05_1_err[:,1:], marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='#1f78b4')
            ax.errorbar(log_L_SLob_05_1_PS, log_phi_slob_05_1_PS, xerr=log_L_SLob_05_1_err_PS, \
                        yerr=log_phi_slob_05_1_err_PS, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='mediumorchid')
                
                
        elif idx2[index] == 1:
            # ax.errorbar(log_P_williams, z_1_15, xerr=0.25, yerr = z_1_15_err, marker = 's', c = 'dodgerblue',\
            #             markersize = 10, mfc='none', capthick = 2, linestyle='none')
            # ax.errorbar(LF_Konda_1_15[:,0], LF_Konda_1_15[:,1], xerr=0.15, yerr = [LF_Konda_1_15[:,2],LF_Konda_1_15[:,3]], marker = 'D',\
            #           c = 'forestgreen', mfc='none', markersize = 10, capthick = 2, linestyle='none')
            ax.set_xlabel(r'$\log_{10}$(L$_{144 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')
            ax.errorbar(log_L_SLob_1_15, log_phi_slob_1_15, xerr=log_L_SLob_1_15_err, \
                        yerr=log_phi_slob_1_15_err, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='#1f78b4')
            ax.errorbar(log_L_SLob_1_15_PS, log_phi_slob_1_15_PS, xerr=log_L_SLob_1_15_err_PS, \
                        yerr=log_phi_slob_1_15_err_PS, marker="*", mfc='none', markersize = 10, \
                            capthick = 2, linestyle='none', color='mediumorchid')
    if idx1[index] == 2 and idx2[index] == 0:
        # ax.errorbar(log_P_williams, z_15_2, xerr=0.25, yerr = z_15_2_err, marker = 's',c = 'dodgerblue',\
        #             markersize = 10, capthick = 2, mfc='none', linestyle='none', label = 'Williams et al. (2018)\nSF + AGN')
        # ax.errorbar(LF_Konda_15_2[:,0], LF_Konda_15_2[:,1], xerr=0.15, yerr = [LF_Konda_15_2[:,2],LF_Konda_15_2[:,3]], marker = 'D',\
        #               c = 'forestgreen', markersize = 10, mfc='none', capthick = 2, linestyle='none', label = 'Kondapally et al. (2022)\nAGN')
        ax.set_xlabel(r'$\log_{10}$(L$_{144 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')
        ax.errorbar(log_L_SLob_15_3, log_phi_slob_15_3, xerr=log_L_SLob_15_3_err, \
                    yerr=log_phi_slob_15_3_err, marker="*", mfc='none', markersize = 10, \
                        capthick = 2, linestyle='none', color='#1f78b4', label='Slob et al. 2022 Master sample') 
        ax.errorbar(log_L_SLob_15_3_PS, log_phi_slob_15_3_PS, xerr=log_L_SLob_15_3_err_PS, \
                    yerr=log_phi_slob_15_3_err_PS, marker="*", mfc='none', markersize = 10, \
                        capthick = 2, linestyle='none', color='mediumorchid', label='Slob et al. 2022 PS sample')             
            
    if idx2[index] == 1:
        ax.tick_params(left = True, labelleft = False, right = True, bottom = True, labelbottom = True, top = True)

    if idx1[index] == 2 and idx2[index] == 0: 
        ax.legend(bbox_to_anchor = (1.05,0.8),frameon = True, handletextpad = 1.5, labelspacing = 0.6)

    # ax.tick_params(left = True, right = True, bottom = True, labelbottom = False, top = True)
    ax.tick_params(axis='both',which='both',top=True,right=True)


    # if index != 6:
    #     print('\n')
    #     print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; LF Sample} &")
    #     for i, LF in enumerate(log_lum_func_144[idx_cutoff[index]:]):
    #         print(np.round((central_bins_144[idx_cutoff[index]:])[i],2), "$\pm$",\
    #               np.round((log_xerr_144[idx_cutoff[index]:])[i],1), "&", (n_144[idx_cutoff[index]:])[i], \
    #                   "& $", np.round(LF,2), "^{+", np.round((log_err_lum_func_144[:,idx_cutoff[index]:])[1,i],2)\
    #                       , "}_{-", np.round((log_err_lum_func_144[:,idx_cutoff[index]:])[0,i],2), "}$&")

    #     print('\n')
    #     print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; MPS Sample} \ " + "\ ")
    #     for i, LF in enumerate(log_lum_func_144_PS[idx_cutoff_PS[index]:]):
    #         print(np.round((central_bins_144_PS[idx_cutoff_PS[index]:])[i],2), \
    #               "$\pm$", np.round((log_xerr_144_PS[idx_cutoff_PS[index]:])[i],1), \
    #                   "&", (n_144_PS[idx_cutoff_PS[index]:])[i], "& $", np.round(LF,2),\
    #                       "^{+", np.round((log_err_lum_func_144_PS[:,idx_cutoff_PS[index]:])[1,i],2),\
    #                           "}_{-", np.round((log_err_lum_func_144_PS[:,idx_cutoff_PS[index]:])[0,i],2), "}$ \ " + "\ ")



    # Plotting the residuals
    # Interpolate LF to get residuals from PS function
    markers = ["s", "o", "^", "D", "v"]
    colors = ['k', 'crimson', 'dodgerblue', 'darkorange', 'darkviolet']
    if index != 5:
        interp_LF = np.interp(central_bins_144_PS[idx_cutoff_PS[index]:], central_bins_144, log_lum_func_144)
        resid = interp_LF - log_lum_func_144_PS[idx_cutoff_PS[index]:]
        
        
        i_min_diff = []
        for x_PS in (central_bins_144_PS[idx_cutoff_PS[index]:]):
            i_min_diff.append(np.argmin(abs(x_PS - central_bins_144)))
            
        log_err_x = log_err_lum_func_144[:,i_min_diff]
        err_resid = np.sqrt(log_err_lum_func_144_PS[:,idx_cutoff_PS[index]:]**2 + log_err_x**2)
        # print("difference between PS LF & LF:", resid, '+/-', err_resid, "\n")
        # print("\n")
        
        err_resid_symm = np.sqrt( np.array(err_resid)[0,:] **2 + np.array(err_resid)[1,:] **2 ) # symmetric errors
        resid_weights = 1 / (np.array(err_resid_symm) ** 2) #weights 
        resid_mean = np.sum(resid * resid_weights) / np.sum(resid_weights) #weighted mean
        uncertainty_in_mean = np.sqrt(1/np.sum(resid_weights))

        print("{3}<z<{2} MPS {0:.4} \pm {1:.3}".format(resid_mean, uncertainty_in_mean,  str(redshift_bins[index+1]), str(redshift_bins[index])))
        
        # for i, r in enumerate(resid):
        #     print("& $", np.round(r,2), "^{+", np.round(err_resid[1,i],2), "}_{-", np.round(err_resid[0,i],2), "}$")

        # if LF_name == 'SDSS':
        mfc_color = colors[index]
        legend_lab = str(z_bin_min)+'< z <' + str(redshift_bins[index+1])
        ax1.errorbar(central_bins_144_PS[idx_cutoff_PS[index]:], resid, yerr = err_resid, linestyle = ':', marker = markers[index],\
                mfc = mfc_color, c = colors[index], markersize = 9, label = legend_lab)
        # else:
        #     print('\n')
        #     print("\multicolumn{3}{c}{$z < " + str(redshift_bins[index])+"$; Master Sample} \ " + "\ ")
        #     for i, LF in enumerate(log_lum_func_144):
        #         print(np.round(central_bins_144[i],2), "$\pm$", np.round(log_xerr_144[i],1), "&", lum_func_counts_144[i], "& $", np.round(LF,2),\
        #              "^{+", np.round(log_err_lum_func_144[1,i],2), "}_{-", np.round(log_err_lum_func_144[0,i],2), "}$")

        #     print('\n')
        #     print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; PS Sample} \ " + "\ ")

        #     for i, LF in enumerate(log_lum_func_144_PS):
        #         # for idx, r in enumerate(resid):
        #         r = resid[i]
        #         print(np.round(central_bins_144_PS[i],2), "$\pm$", np.round(log_xerr_144_PS[i],1), "&", bin_count_PS[i], "& $", np.round(LF,2),\
        #                 "^{+", np.round(log_err_lum_func_144_PS[1,i],2), "}_{-", np.round(log_err_lum_func_144_PS[0,i],2), "}$ & $", np.round(r,2),\
        #                 "^{+", np.round(err_resid[1,i],2), "}_{-", np.round(err_resid[0,i],2), "}$ \ " + "\ ")



        #     mfc_color = 'none'
        #     legend_lab = '_nolegend_'
        
        ax1.set_xlabel(r'$\log_{10}$(L$_{144 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')
        ax1.set_ylabel(r'$\Delta \log_{10}(\Phi)$')
        ax1.tick_params(axis='both',which='both',top=True,right=True)
        ax1.set_xlim(22.5,29)
        ax1.set_ylim(-0.1,2.5)
        ax1.legend()
    
fig.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/LF_lumfunc_revised.png')
fig2.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/LF_resid_lumfunc_revised.png')

