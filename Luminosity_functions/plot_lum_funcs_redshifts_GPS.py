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

# VLASS fractional area
str_to_deg2 = (180/np.pi)**2
frac_area_VLASS = 5634 / (4*np.pi * str_to_deg2)  #sr, LoLSS sky coverage (650 deg2)

# Define different redshift bins, all bins have ~100+ sources
redshift_bins = [0.,0.1,0.3,0.5,0.7,0.9,1.4,3.,3.]

idx1 = [0,0,1,1,2,2,3]
idx2 = [0,1,0,1,0,1,0]

#Cutoff where it becomes incomplete
#TODO: apply these as well to the residual plot
idx_cutoff=[0,0,0,1,1,2,3]
idx_cutoff_PS = [0,0,0,0,0,0,0] #last one 2

bin_size = [0.5,0.25,0.4,0.35,0.35, 0.3, 0.3]
bin_size_PS = [0.9,0.7,0.75,1,0.5, 0.6, 0.55] 


'''Importing literature LFs'''
# Heckman & Best 2014 lum function fit, 1400MHz
y_max = 27.
P = 10**(np.linspace(19,y_max+2,1000))
P_0_14 = 10**(24.95)
P_0 = P_0_14 #* (1400/1400)**alpha_extrap
a = 0.42
b = 1.66
rho_0 = 10**(-5.33)
rho = rho_0/( (P/P_0)**a + (P/P_0)**b) 

#Best & Heckman 2012 local lum function z<0.3
log_L_min=np.array([22,22.3,22.6,22.9,23.2,23.5,23.8,24.1,24.4,24.7,25,25.3,25.6,25.9,26.2])
log_L_max = np.array([22.3,22.6,22.9,23.2,23.5,23.8,24.1,24.4,24.7,25,25.3,25.6,25.9,26.2, 26.5])
log_L_central = (log_L_min + log_L_max)/2
all_radio = ([-3.09,-3.49,-3.87,-4.22, -4.56,-4.75,-4.9,-5.08,-5.25,-5.54,-5.82,-6.32,-6.58,-7.18,-7.78])
all_radio_err = np.array([[0.03,0.02,0.02,0.02,0.02,0.01,0.02,0.01,0.02,0.02,0.03,0.05,0.07,0.12,0.21],\
                         [0.03,0.03,0.02,0.02,0.02,0.01,0.02, 0.01,0.02, 0.02,0.03,0.06,0.08,0.17,0.43]])
starforming = [-3.11,-3.54,-4,-4.57,-5.13,-5.69,-6.35,-6.83,-7.43,0,0,0,0,0,0]
starforming_err = [[0.03,0.03,0.02,0.03, 0.03,0.04,0.06,0.11,0.23,0,0,0,0,0,0],\
    [0.04,0.03,0.03,0.03,0.03,0.05, 0.08, 0.15, 0.51]]
AGN_all = [-4.44,-4.45,-4.45,-4.48,-4.69,-4.8,-4.91,-5.09,-5.26,-5.54,-5.82,-6.32,-6.58,-7.18,-7.78]
AGN_all_err = [[0.15,0.05,0.04,0.02,0.02,0.01,0.02,0.01,0.02,0.02,0.03,0.05,0.07,0.12,0.21],[0.23,0.05,0.04,0.03,0.02,0.01,0.02,0.02,0.02,0.02,0.03,0.06,0.08,0.17,0.18]]

# Pracy lum function at 1400
#TODO: should we add these? 
C_SF = 10**(-2.83)
P_0_SF = 10**(21.18)
P_0_SF = P_0_SF #* (1400/1400)**alpha_extrap
a_SF = 1.02
b_SF = 0.60  
rho_SF = C_SF*((P/P_0_SF)**(1-a_SF) * np.exp(-0.5*(np.log10(1+P/P_0_SF)/b_SF)**2))

#Ceraj
Lum_AGN_01_05 = [21.82,22.62,23.42,24.22,25.02,25.82]
phi_AGN_01_05 = [-3.99,-4.09,-4.74,-5.08,-5.50,-5.68]
err_phi_AGN_01_05 = [0.15,0.04,0.06,0.10,0.22,0.30]

Lum_AGN_05_07 = [22.62,23.42,24.22,25.02]
phi_AGN_05_07 = [-4.07,-4.45,-4.88,-5.81]
err_phi_AGN_05_07 = [0.11,0.05,0.06,0.30]

Lum_AGN_07_09 = [22.96,23.76,24.56,25.36,26.16]
phi_AGN_07_09 = [-4.61,-4.93,-5.41,-6.26,-6.26]
err_phi_AGN_7_09 = [0.09,0.07,0.11,0.53,0.53]

Lum_AGN_07_09 = [22.96,23.76,24.56,25.36,26.16]
phi_AGN_07_09 = [-3.99,-4.35,-4.98,-5.66,-5.96]
err_phi_AGN_07_09 = [0.11,0.09,0.06,0.17,0.3]

Lum_AGN_09_11 = [23.21, 24.01, 24.81, 25.61, 26.41, 27.21]
phi_AGN_09_11 = [-4.25,-4.66,-5.35,-5.58,-6.06,-6.36]
err_phi_AGN_09_11 = [0.13,0.06,0.08,0.12,0.3,0.53]

Lum_AGN_11_14 = [23.42,24.22,25.02,25.82,26.62]
phi_AGN_11_14 = [-4.51,-4.78,-5.49,-5.91,-6.13]
err_phi_AGN_11_14 = [0.05,0.04,0.07,0.14,0.22]

Lum_AGN_17_21 = [23.85,24.65,25.45, 26.25]
phi_AGN_17_21 = [-4.67,-4.90,-5.62,-6.25]
err_phi_AGN_17_21 = [0.06,0.09,0.06, 0.16]

fig, axs = plt.subplots(4, 2, figsize=(20,28), sharex = True, sharey = True)
plt.subplots_adjust(wspace=0, hspace=0)
ax6 = axs[3, 1]
ax6.set_visible(False)

fig2, ax1 = plt.subplots(figsize=(10,7.5))

for index, z_bin_min in enumerate(redshift_bins[:-2]):

    hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/Luminosity_functions/GPS_lum_func_1400.fits')
    tbdata = hdulist[1].data
    hdulist.close()

    alpha_low = tbdata['alpha_low']
    alpha_high = tbdata['alpha_hiSh'] #for some reason there was a typo in the Vmax code. should be high
    e_alpha_low = tbdata['e_alpha_low']
    e_alpha_high = tbdata['e_alpha_high']
    z_max_1400 = tbdata['z_max_1400']
    V_max_opt = tbdata['V_max_opt']
    V_max_radio = tbdata['V_max_radio']
    Power_1400 = tbdata['Power_1400']
    z = tbdata['z_reported']

    # Calculate total volume
    V_max_1400 = np.amin([V_max_radio,V_max_opt], axis = 0) * frac_area_VLASS

    marker_face = 'k'
    marker_face_PS = 'crimson'
    label_name = 'SDSS'
    m_i = '21.3'
    
    # Calculate max volume of redshift range
    if index < 7:
        V_max_limit = calc_volume(redshift_bins[index+1]) * frac_area_VLASS
    elif index == 7:
        V_max_limit = calc_volume(3.) * frac_area_VLASS # for full z ranges, use z_max = 3
    
    # Define values for this redshift range
    if index == 7:
        #I believe this is only when plotting the entire thing function, so no redshift bins
        z_idx = np.where((z > 0.) & (z < 3.))
    elif index < len(redshift_bins) - 1:
        z_idx = np.where((z > z_bin_min) & (z < redshift_bins[index+1]) & (np.log10(Power_1400) > 23.))
    else:
        break

    # Make array with correct size to find min of V_max 
    V_max_i = np.empty(np.shape(z_idx[0]))
    V_max_i.fill(V_max_limit)

    # Define values for this redshift range
    Power_1400_z = Power_1400[z_idx]
    z_max_1400_z = z_max_1400[z_idx]
    V_max_1400_z = np.amin([V_max_1400[z_idx],V_max_i], axis = 0)
    alpha_low_z = alpha_low[z_idx]
    alpha_high_z = alpha_high[z_idx]
    e_alpha_low_z = e_alpha_low[z_idx]
    e_alpha_high_z = e_alpha_high[z_idx]
    
    
    # Also take minimum z of range into account
    # 5 is redshift bin 0-3, 0 is redshift bin 0-0.1. 
    if index > 0 and index < 7: # only when z_min != 0
        V_min = calc_volume(z_bin_min) * frac_area_VLASS
        V_i = V_max_1400_z - V_min
        V_idx = np.where(V_i > 0)
        V_i = V_i[V_idx]
        alpha_low_z = alpha_low_z[V_idx]
        alpha_high_z = alpha_high_z[V_idx]
        e_alpha_low_z = e_alpha_low_z[V_idx]
        e_alpha_high_z = e_alpha_high_z[V_idx]
        Power_1400_z = Power_1400_z[V_idx]
    else:
        #0.001 is ongeveer 0
        V_min = calc_volume(0.001) * frac_area_VLASS
        V_i = V_max_1400_z - V_min

    # Compute 10log of power
    log_Power_1400 = np.log10(Power_1400_z)

    # define PS sample
    ind_peaked = np.where((alpha_low_z > e_alpha_low_z)\
                          & (alpha_high_z < 0))
    # ind_peaked = np.where((alpha_low_z > e_alpha_low_z)\
    #                       & (alpha_high_z < 0))
    z_max_1400_PS = z_max_1400_z[ind_peaked]
    V_i_PS = V_i[ind_peaked]

    # Define luminosity luminosity bins for 1400MHz, use SDSS limits
    bin_size_i = bin_size[index]
    if index == 0 or index == 7:
        low_limit_1400 = 23 # don't use the lowest luminosities, incomplete there
    else:
        low_limit_1400 = np.min(log_Power_1400)
    factor = 1
    if index in [1,7]:
        factor = 2
    bin_edges_1400 = np.hstack((np.arange(low_limit_1400, np.max(log_Power_1400)-factor*bin_size_i, bin_size_i),np.max(log_Power_1400)))
    n_1400, bin_edges_1400 = np.histogram(log_Power_1400, bins=bin_edges_1400)
    
    # Define luminosity bins PS sample
    if index == 7:
        bin_edges_1400_PS = np.arange(23.,29.4, bin_size_PS[index]) 
    else:
        bin_edges_1400_PS = np.arange(np.min(log_Power_1400[ind_peaked]), np.max(log_Power_1400[ind_peaked])+bin_size_PS[index], bin_size_PS[index]) 
    n_1400_PS, bin_edges_1400_PS = np.histogram(log_Power_1400[ind_peaked], bins = bin_edges_1400_PS)
    

    # Find central bins for plotting
    central_bins_1400  = np.array([(j+i)/2. for i, j in zip(bin_edges_1400[:-1], bin_edges_1400[1:])])
    central_bins_1400_PS = np.array([(j+i)/2. for i, j in zip(bin_edges_1400_PS[:-1], bin_edges_1400_PS[1:])])

    # Compute the luminosity function
    lum_func_1400 = (1/(np.diff(bin_edges_1400))) * stats.binned_statistic(log_Power_1400, 1/V_i, 'sum', bins=bin_edges_1400)[0]
    err_lum_func_1400 = np.sqrt( (1/np.diff(bin_edges_1400))**2 * stats.binned_statistic(log_Power_1400, (1/V_i)**2, 'sum', bins=bin_edges_1400)[0])
    log_lum_func_1400 = np.log10(lum_func_1400)
    lum_func_counts_1400 = stats.binned_statistic(log_Power_1400, log_Power_1400, 'count', bins=bin_edges_1400)[0]
    log_err_lum_func_1400 = poisson_errors(lum_func_1400, err_lum_func_1400, lum_func_counts_1400)
    

    # Compute the PS luminosity function
    lum_func_1400_PS = (1/np.diff(bin_edges_1400_PS)) * stats.binned_statistic(log_Power_1400[ind_peaked], 1/V_i_PS, 'sum', bins=bin_edges_1400_PS)[0]
    log_lum_func_1400_PS = np.log10(lum_func_1400_PS)
    err_lum_func_1400_PS = np.sqrt((1/np.diff(bin_edges_1400_PS))**2 * stats.binned_statistic(log_Power_1400[ind_peaked], (1/V_i_PS)**2, 'sum', bins=bin_edges_1400_PS)[0])
    bin_count_PS = stats.binned_statistic(log_Power_1400[ind_peaked], log_Power_1400[ind_peaked], 'count', bins=bin_edges_1400_PS)[0] # Number of sources per bin
    log_err_lum_func_1400_PS = poisson_errors(lum_func_1400_PS, err_lum_func_1400_PS, bin_count_PS)

    # print("# sources in each bin", lum_func_counts_1400)
    log_xerr_1400 = np.diff(bin_edges_1400) / 2
    # print("# PS sources in each bin", bin_count_PS)
    log_xerr_1400_PS = np.diff(bin_edges_1400_PS) / 2


    ax = axs[idx1[index], idx2[index]]


    # Plot luminosity function
    ax.errorbar(central_bins_1400[idx_cutoff[index]:], log_lum_func_1400[idx_cutoff[index]:],\
                xerr=log_xerr_1400[idx_cutoff[index]:], yerr = log_err_lum_func_1400[:,idx_cutoff[index]:], c = 'k', \
                marker = 'o', markersize = 10, linestyle ='none', mfc = marker_face, capthick = 2, label = 'HF sample')#', $m_i\, <$'+ m_i)
    ax.errorbar(central_bins_1400_PS[idx_cutoff_PS[index]:], log_lum_func_1400_PS[idx_cutoff_PS[index]:],\
                xerr=log_xerr_1400_PS[idx_cutoff_PS[index]:], yerr = log_err_lum_func_1400_PS[:,idx_cutoff_PS[index]:], c = 'crimson',\
                marker = '^',  capthick = 2, markersize = 10, mfc = marker_face_PS, linestyle = 'none', label = 'GPS sample')#', $m_i\, <$'+ m_i)
    

    ax.set_xlim(22.5,30)
    ax.set_ylim((-12,-3.5))
    
    #Plot the Heckman&Best line
    #TODO: should we plot this at low redshift?
    ax.plot(np.log10(P), np.log10(rho), linewidth = 2, linestyle = '--', color = 'k', zorder = -10, label = 'Heckman & Best (2014)\nz ~ 0.1')    
    # ax.plot(np.log10(P), np.log10(rho_SF), linewidth = 2, linestyle = '--', color = 'k', zorder = -10, label = 'Test')    

    ax.annotate(str(len(z_idx[0]))+' SDSS sources, ' + str(len(z_max_1400_PS)) + ' SDSS PS sources', xy = (23,-10.2), xycoords = 'data', fontsize = 15)
    ax.annotate(str(z_bin_min)+' < z <' + str(redshift_bins[index+1]), xy=(23,-9.8), xycoords = 'data', fontsize = 17)

    if idx1[index] == 0:                    
        if idx1[index] == 0 or idx1[index] == 1: 
            ax.tick_params(left = True, labelleft = True, right = True, bottom = True, labelbottom = False, top = True)
            
        if idx2[index] == 0: #z<0.1
            ax.errorbar(log_L_central, all_radio, yerr=all_radio_err, marker = 'D', label='Best & Heckman (2012)\n [z=0-0.3]',\
                       c = 'forestgreen', mfc='none', markersize = 10, capthick = 2, linestyle='none')
        elif idx2[index] == 1: #0.1<z<0.3
            ax.errorbar(log_L_central, all_radio, yerr=all_radio_err, marker = 'D', label='Best & Heckman (2012)\n [z=0-0.3]',\
                       c = 'forestgreen', mfc='none', markersize = 10, capthick = 2, linestyle='none')


    if idx1[index] == 1:
        if idx2[index] == 0: #0.3<z<0.5
            ax.set_ylabel(r'$\log_{10} ( \Phi$ / Mpc$^{-3}$ dex$^{-1}$ )')
            ax.errorbar(Lum_AGN_01_05, phi_AGN_01_05, yerr=err_phi_AGN_01_05, marker = 's', label='AGN Ceraj et al. (2018)',\
                       c = 'deeppink', mfc='none', markersize = 10, capthick = 2, linestyle='none')

        elif idx2[index] == 1: #0.5<z<0.7
            ax.errorbar(Lum_AGN_05_07, phi_AGN_05_07, yerr=err_phi_AGN_05_07, marker = 's', label='AGN Ceraj et al. (2018)',\
                       c = 'deeppink', mfc='none', markersize = 10, capthick = 2, linestyle='none')

    if idx1[index] == 2:
        if idx2[index] == 0: #0.7<z<0.9
            ax.errorbar(Lum_AGN_07_09, phi_AGN_07_09, yerr=err_phi_AGN_07_09, marker = 's', label='AGN Ceraj et al. (2018)',\
                       c = 'deeppink', mfc='none', markersize = 10, capthick = 2, linestyle='none')
        elif idx2[index] == 1:   #0.9<z<1.4
                ax.errorbar(Lum_AGN_11_14, phi_AGN_11_14, yerr=err_phi_AGN_11_14, marker = 's', label='AGN Ceraj et al. (2018)',\
                           c = 'deeppink', mfc='none', markersize = 10, capthick = 2, linestyle='none')    

    if idx1[index] == 3: #1.4<z<3
        if idx2[index] == 0: #0.7<z<0.9
            ax.errorbar(0, 0, yerr=1, marker = 'D', label='Best & Heckman (2012)\nz < 0.3',\
                       c = 'forestgreen', mfc='none', markersize = 10, capthick = 2, linestyle='none')
            ax.errorbar(Lum_AGN_17_21, phi_AGN_17_21, yerr=err_phi_AGN_17_21, marker = 's', label='AGN Ceraj et al. (2018)',\
                       c = 'deeppink', mfc='none', markersize = 10, capthick = 2, linestyle='none')
            ax.set_xlabel(r'$\log_{10}$(L$_{1400 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')

        
    if idx2[index] == 1:
        ax.tick_params(left = True, labelleft = False, right = True, bottom = True, labelbottom = True, top = True)

    if idx1[index] == 3 and idx2[index] == 0: 
        ax.legend(bbox_to_anchor = (1.05,0.8), frameon = True, handletextpad = 1.5, labelspacing = 0.6)

    # ax.tick_params(left = True, right = True, bottom = True, labelbottom = False, top = True)
    ax.tick_params(axis='both',which='both',top=True,right=True)


    if index != 8:
        print('\n')
        print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; HF Sample} &")
        for i, LF in enumerate(log_lum_func_1400[idx_cutoff[index]:]):
            print(np.round((central_bins_1400[idx_cutoff[index]:])[i],2), "$\pm$",\
                  np.round((log_xerr_1400[idx_cutoff[index]:])[i],1), "&", (n_1400[idx_cutoff[index]:])[i], \
                      "& $", np.round(LF,2), "^{+", np.round((log_err_lum_func_1400[:,idx_cutoff[index]:])[1,i],2)\
                          , "}_{-", np.round((log_err_lum_func_1400[:,idx_cutoff[index]:])[0,i],2), "}$&")


        print('\n')
        print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; GPS Sample} \ " + "\ ")
        for i, LF in enumerate(log_lum_func_1400_PS[idx_cutoff_PS[index]:]):
            print(np.round((central_bins_1400_PS[idx_cutoff_PS[index]:])[i],2), \
                  "$\pm$", np.round((log_xerr_1400_PS[idx_cutoff_PS[index]:])[i],1), \
                      "&", (n_1400_PS[idx_cutoff_PS[index]:])[i], "& $", np.round(LF,2),\
                          "^{+", np.round((log_err_lum_func_1400_PS[:,idx_cutoff_PS[index]:])[1,i],2),\
                              "}_{-", np.round((log_err_lum_func_1400_PS[:,idx_cutoff_PS[index]:])[0,i],2), "}$ \ " + "\ ")

    # Plotting the residuals
    # Interpolate LF to get residuals from PS function
    markers = ["s", "o", "^", "D", "v", "*", "s"]
    colors = ['k', 'crimson', 'dodgerblue', 'darkorange', 'darkviolet', 'limegreen', 'deeppink']
    if index != 7:
        interp_LF = np.interp(central_bins_1400_PS[idx_cutoff_PS[index]:], central_bins_1400, log_lum_func_1400)
        resid = interp_LF - log_lum_func_1400_PS[idx_cutoff_PS[index]:]
        
        i_min_diff = []
        for x_PS in (central_bins_1400_PS[idx_cutoff_PS[index]:]):
            i_min_diff.append(np.argmin(abs(x_PS - central_bins_1400)))
            
        log_err_x = log_err_lum_func_1400[:,i_min_diff]
        err_resid = np.sqrt(log_err_lum_func_1400_PS[:,idx_cutoff_PS[index]:]**2 + log_err_x**2)
        
        # print("difference between GPS & HF:", resid, '+/-', err_resid, "\n")

        # for i, r in enumerate(resid):
        #     print("& $", np.round(r,2), "^{+", np.round(err_resid[1,i],2), "}_{-", np.round(err_resid[0,i],2), "}$")


        mfc_color = colors[index]
        legend_lab = str(z_bin_min)+'< z <' + str(redshift_bins[index+1])
        ax1.errorbar(central_bins_1400_PS[idx_cutoff_PS[index]:], resid, yerr = err_resid, linestyle = ':', marker = markers[index],\
                mfc = mfc_color, c = colors[index], markersize = 9, label = legend_lab)     # else:
        # print('\n')
        # print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; PS Sample} \ " + "\ ")
        # for i, LF in enumerate(log_lum_func_1400):
        #     print(np.round(central_bins_1400[i],2), "$\pm$", np.round(log_xerr_1400[i],1), "&", lum_func_counts_1400[i], "& $", np.round(LF,2),\
        #           "^{+", np.round(log_err_lum_func_1400[1,i],2), "}_{-", np.round(log_err_lum_func_1400[0,i],2), "}$")

        # print('\n')
        # print("\multicolumn{3}{c}{$"+str(z_bin_min)+' < z <' + str(redshift_bins[index+1])+"$; GPS Sample} \ " + "\ ")

        # for i, LF in enumerate(log_lum_func_1400_PS):
        #     # for idx, r in enumerate(resid):
        #     r = resid[i]
        #     print(np.round(central_bins_1400_PS[i],2), "$\pm$", np.round(log_xerr_1400_PS[i],1), "&", bin_count_PS[i], "& $", np.round(LF,2),\
        #             "^{+", np.round(log_err_lum_func_1400_PS[1,i],2), "}_{-", np.round(log_err_lum_func_1400_PS[0,i],2), "}$ & $", np.round(r,2),\
        #             "^{+", np.round(err_resid[1,i],2), "}_{-", np.round(err_resid[0,i],2), "}$ \ " + "\ ")



        #     mfc_color = 'none'
        #     legend_lab = '_nolegend_'
        
        ax1.set_xlabel(r'$\log_{10}$(L$_{1400 \, \mathrm{MHz}}$ / W Hz$^{-1}$)')
        ax1.set_ylabel(r'$\Delta \log_{10}(\Phi)$')
        ax1.tick_params(axis='both',which='both',top=True,right=True)
        ax1.set_xlim(22.5,29)
        ax1.set_ylim(-0.4,2.5)
        ax1.legend(bbox_to_anchor=(0.9, 0.53))
    
fig.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/HF_lumfunc_revised_alpha0.pdf')
fig2.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/HF_resid_lum_func_revised_alpha0.pdf')
