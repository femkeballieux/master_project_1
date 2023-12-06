#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A code for classifying sources as HERGs or LERGs. Method might be a bit convoluted but
it is the best I can do with the available data.

Start by crossmatching  in 2''the Best&Heckman, the DR12 and the DR8 catalogs to the master sample in topcat. Then run this code

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.27)

#read in the data
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_emission_lines.fits')
tbdata = hdulist[1].data
colnames = hdulist[1].columns
hdulist.close()


name=tbdata['LoTSS_name']
RA=tbdata['RA_1']
# print(colnames)
# print('hoi', len(RA[np.where(tbdata['H']==1)]))
# print('hoi', len(RA[np.where(tbdata['L']==1)]))
# print('hoi', len(RA[np.where((tbdata['L']==0) & (tbdata['H']==0))]))
#in which we store what the HERG/LERG classification is, and how it was defined
def L_z(z, alpha, flux):
    dl = cosmo.luminosity_distance(z).to(u.m).value
    k_corr = np.power((1 + z),-(1+alpha))
    lum = 4 * np.pi * dl**2 * flux * k_corr * 1e-26 # W Hz^-1 m^-2
    return lum

   
def L_optical(flux, z):
    """
    Gives luminosity of spectral line / luminosity of sun
    """
    dl_z = cosmo.luminosity_distance(z).to(u.cm) /u.cm # in cm  
    L = flux * 4 * np.pi * dl_z ** 2 # 1e-17 erg / s / cm^2 * cm^2 = 1e-17 erg/s
    L = L * 1e-17 # erg/s
    
    return L / 3.828e33

for PS in ['MPS', 'GPS']:
    classification = np.zeros_like(RA) #1 if LERG, 2 if HERG, 3 if unclassified
    classification_source = np.zeros_like(RA) #1 if from Best, 2 if from DR12 EI, 3 from DR12 EW, 4 from DR8 EW
    if PS == 'MPS':
        F = 'LF'
        sample='LoLSS'
    else:
        F = 'HF'
        sample = 'VLASS'
        print('')
        
    print('--------------------------------------------------------------------------')

    
    alpha_low = tbdata['alpha_low_'+sample]
    e_alpha_low = tbdata['e_alpha_low_'+sample]
    alpha_high = tbdata['alpha_high_'+sample]
    e_alpha_high = tbdata['e_alpha_high_'+sample]
    
    mask_F = (e_alpha_low>0.)
    # mask_PS = ((alpha_low > e_alpha_low) & (alpha_high < -e_alpha_high))
    mask_PS = ((alpha_low > e_alpha_low) & (alpha_high < 0))    
    #Masks for Bestman&Heck 2012 catalog
    Best_HERG_mask= (tbdata['H'] == 1)
    Best_LERG_mask= (tbdata['L'] == 1)
    not_HERG_mask= (tbdata['H']==0)
    not_LERG_mask= (tbdata['L']==0)
    
    # AGN_mask = (tbdata['A_2']==1)
    # starforming_mask = (tbdata['A_2']==0)
    # main_samp = (tbdata['M']==1)
    
    #Store the definition  
    classification[np.where(Best_LERG_mask)] = 1   #LERG
    classification[np.where(Best_HERG_mask)] = 2   #HERG
    classification[np.where(not_HERG_mask & not_LERG_mask)] = 3   #Unclassified
    classification_source[np.where(Best_LERG_mask)] = 1  
    classification_source[np.where(Best_HERG_mask)] = 1
    classification_source[np.where(not_HERG_mask & not_LERG_mask)] = 1
    
    A = len(name[np.where((classification==1) & mask_F & (classification_source==1))])
    B = len(name[np.where((classification==2) & mask_F & (classification_source==1))])
    C = len(name[np.where((classification==3) & mask_F & (classification_source==1))])
    print('Best&Heckman: there are {0} LERGS, {1} HERGS, and {2} unclassified in {3}'.\
          format(A, B, C, F))
        
    a = len(name[np.where((classification==1) & mask_PS & (classification_source==1))])
    b = len(name[np.where((classification==2) & mask_PS & (classification_source==1))])
    c = len(name[np.where((classification==3) & mask_PS & (classification_source==1))])
    
    print('Best&Heckman: there are {0} LERGS, {1} HERGS, and {2} unclassified in {3}'.\
          format(a, b, c, PS))
    print('')
    
    #mThese sources have not been classified in Heckman&Best
    not_yet_classified_mask = ((classification != 1) & (classification!=2))

                              
    #Read in the DR12 Portsmouth stuff
    H_alpha_flux_DR12 = tbdata['FLUX'][:,24]
    H_alpha_flux_err_DR12 = tbdata['FLUX_ERR'][:,24]
    
    H_beta_flux_DR12 = tbdata['FLUX'][:,15]
    H_beta_flux_err_DR12 = tbdata['FLUX_ERR'][:,15]
    
    N_II_flux_DR12 = tbdata['FLUX'][:,25] #6583.34A, from buttiglione
    N_II_flux_err_DR12 = tbdata['FLUX_ERR'][:,25]
    
    S_II_flux_DR12 = tbdata['FLUX'][:,26] # 6716.31A, there is another, but in the paper they use this one
    S_II_flux_err_DR12 = tbdata['FLUX_ERR'][:,26]
    
    O_III_flux_DR12 = tbdata['FLUX'][:,17] #  5006.77A
    O_III_flux_err_DR12 = tbdata['FLUX_ERR'][:,17]
    O_III_ew_DR12 = tbdata['EW'][:,17]
    O_III_ew_err_DR12 = tbdata['EW_ERR'][:,17]
    
    O_I_flux_DR12 = tbdata['FLUX'][:,22] # 6363.67A
    O_I_flux_err_DR12 = tbdata['FLUX_ERR'][:,22]
    
    #Calculate the excitation index
    EI = np.log10(O_III_flux_DR12 / H_beta_flux_DR12) - (1/3)*( np.log10(N_II_flux_DR12/H_alpha_flux_DR12)\
        +np.log10(S_II_flux_DR12 /H_alpha_flux_DR12) + np.log10(O_I_flux_DR12/H_alpha_flux_DR12))
    
    
    #define the necessary masks
    non_zero_mask = ((O_III_flux_DR12 > 0.) & (H_beta_flux_DR12 > 0.) & (H_alpha_flux_DR12>0.) & (N_II_flux_DR12>0.)&\
                     (O_I_flux_DR12>0.) & (S_II_flux_DR12>0.))
    some_zero_mask = ((O_III_flux_DR12 == 0.) | (H_beta_flux_DR12 == 0.) | (H_alpha_flux_DR12==0.) | (N_II_flux_DR12 == 0.)|\
                         (O_I_flux_DR12==0.) | (S_II_flux_DR12 ==0.))
        
    in_DR12_mask = (tbdata['Z_3']>0.)
    
    EI_LERG_mask = (EI<0.95) & (np.isfinite(EI))
    EI_HERG_mask = (EI>0.95) & (np.isfinite(EI))

    print('For those that do not have a classification yet:')
    
    
    print(len(classification[np.where(EI_LERG_mask & non_zero_mask &  (classification==2))]), 'LERG EI, HERG Best')
    print(len(classification[np.where(EI_HERG_mask & non_zero_mask &  (classification==1))]), 'HERG EI, LERG best')
    #If they have a detection in all lines and have not been classified yet then they might be HERG or LERG
    classification[np.where(EI_LERG_mask & non_zero_mask &  not_yet_classified_mask)] = 1 
    classification[np.where(EI_HERG_mask & non_zero_mask & not_yet_classified_mask)] = 2
    classification[np.where(in_DR12_mask & some_zero_mask &  not_yet_classified_mask)] = 3 

    classification_source[np.where(EI_LERG_mask & non_zero_mask & not_yet_classified_mask)] = 2 
    classification_source[np.where(EI_HERG_mask & non_zero_mask & not_yet_classified_mask)] = 2  
    classification_source[np.where(in_DR12_mask & some_zero_mask &  not_yet_classified_mask)] = 2
    
    D = len(name[(classification==1) & mask_F & (classification_source==2)])
    E = len(name[(classification==2) & mask_F & (classification_source==2)])
    G = len(name[(classification==3) & mask_F & (classification_source==2)])
    print('Portsmouth EI criterium: there are {0} LERGS, {1} HERGS, and {2} unclassified in {3}'.\
          format(D, E, G, F))
        
    d = len(name[(classification==1) & mask_PS & (classification_source==2)])
    e = len(name[(classification==2) & mask_PS & (classification_source==2)])
    g = len(name[(classification==3) & mask_PS & (classification_source==2)])
    print('Portsmouth EI criterium: there are {0} LERGS, {1} HERGS, and {2} unclassified in {3}'.\
          format(d, e, g, PS))
    print('')    
    print('For those that do not have a classification yet:')

    #if not yet classified, and oxygen line was detected, we can classify it 
    not_yet_classified_mask_2 = ((classification != 1) & (classification != 2))
    OIII_mask_dr12_HERG = (O_III_ew_DR12>5.) & (np.isfinite(O_III_ew_DR12))
    OIII_mask_dr12_LERG_unclassified = (O_III_ew_DR12<5.) & (np.isfinite(O_III_ew_DR12))
    
    # print(len(classification[np.where(OIII_mask_dr12_HERG & mask_PS &(classification==1) & (classification_source==2))]), 'PS: HERG OII, LERG otherwise')
    # print(len(classification[np.where(OIII_mask_dr12_HERG & mask_F &(classification==1) & (classification_source==2))]), 'F: HERG OII, LERG otherwise')
    
    classification[np.where(OIII_mask_dr12_HERG & not_yet_classified_mask_2)] = 2
    classification[np.where(in_DR12_mask & OIII_mask_dr12_LERG_unclassified & not_yet_classified_mask_2)] = 3  

    classification_source[np.where(OIII_mask_dr12_LERG_unclassified & not_yet_classified_mask_2)] = 3  
    classification_source[np.where(in_DR12_mask & OIII_mask_dr12_HERG & not_yet_classified_mask_2)] = 3

    H = len(name[(classification==2) & mask_F & (classification_source==3)])
    I = len(name[(classification==3) & mask_F & (classification_source==3)])
    print('Portsmouth OIII EW criterium: there are {0} HERGS an {1} unclassified or LERG in {2}'.\
          format(H, I, F))
        
    h = len(name[(classification==2) & mask_PS & (classification_source==3)])
    i = len(name[(classification==3) & mask_PS & (classification_source==3)])
    print('Portsmouth OIII EW criterium: there are {0} HERGS an {1} unclassified or LERG in {2}'.\
          format(h, i, PS))    

        
    print('')   
    # print('For those that do not have a classification yet:')
    # not_yet_classified_mask_3 = ((classification != 1) & (classification != 2))
    
    #Unfortunately, DR8 does not have all the emission lines, so only equivalent width
    #This is the data Heckman&Best uses
    O_III_flux_DR8 = tbdata['OIII_5007_FLUX'] #  5006.77A
    O_III_flux_err_DR8 = tbdata['OIII_5007_FLUX_ERR']
    O_III_ew_DR8 = tbdata['OIII_5007_EQW']
    O_III_ew_err_DR8 = tbdata['OIII_5007_EQW_ERR']
    

    L_14 = L_z(tbdata['z_2'], -0.8, tbdata['NVSS_flux']) 
    L_14_2 = L_z(tbdata['z_3'], -0.8, tbdata['NVSS_flux']) 
    L_14_3 = L_z(tbdata['z_4'], -0.8, tbdata['NVSS_flux']) 
    L_14 = np.where(np.isnan(L_14), L_14_2, L_14)
    L_14 = np.where(np.isnan(L_14), L_14_3, L_14)
    
    #L_OIII/L_sun, calculated either from DR8 or DR12.
    L_OIII = L_optical( O_III_flux_DR8 ,tbdata['Z_2'])
    L_OIII_2 = L_optical( O_III_flux_DR12 ,tbdata['Z_3'])    
    L_OIII = np.where(np.isnan(L_OIII),L_OIII_2, L_OIII)
    
    #Here we make the unclassified figure
    plt.figure(figsize=(10,8))
    plot_mask_HERG = ((classification == 2) & (mask_F))
    plot_mask_LERG = ((classification == 1) & mask_F)
    plot_mask_u = ((classification==3) & mask_F)
    plot_mask_HERG_PS = ((classification==2) & mask_PS) 
    plot_mask_LERG_PS = ((classification==1) & mask_PS) 
    plot_mask_u_PS = ((classification==3)&mask_PS) 
    
    if PS == 'GPS':
        plt.scatter(np.log10(L_14[plot_mask_HERG]), np.log10(L_OIII[plot_mask_HERG]) , alpha=0.5, zorder=1, color='red', s=4, label='HERG {}'.format(F)) 
        plt.scatter(np.log10(L_14[plot_mask_LERG]), np.log10(L_OIII[plot_mask_LERG]) , alpha=0.5, zorder=1, color='black', s=3, label='LERG {}'.format(F))
        # plt.scatter(L_14[plot_mask_u], L_OIII[plot_mask_u] , alpha=0.3, zorder=1, color='green', s=10, label='LERG {}'.format(F))
    
    if PS == 'MPS':
        plt.scatter(np.log10(L_14[plot_mask_HERG]), np.log10(L_OIII[plot_mask_HERG]) , alpha=0.7, zorder=1, color='red', s=4, label='HERG {}'.format(F)) 
        plt.scatter(np.log10(L_14[plot_mask_LERG]), np.log10(L_OIII[plot_mask_LERG]) , alpha=0.7, zorder=1, color='black', s=3, label='LERG {}'.format(F))
        # plt.scatter(L_14[plot_mask_u], L_OIII[plot_mask_u] , alpha=0.3, zorder=1, color='green', s=10, label='LERG {}'.format(F))
    
    plt.scatter(np.log10(L_14[plot_mask_HERG_PS]), np.log10(L_OIII[plot_mask_HERG_PS]) , alpha=1, zorder=20, color='red', s=40, marker='D', label='HERG {}'.format(PS)) 
    plt.scatter(np.log10(L_14[plot_mask_LERG_PS]), np.log10(L_OIII[plot_mask_LERG_PS]) , alpha=1, zorder=20, color='orange', s=40, marker='s', label='LERG {}'.format(PS))
    plt.scatter(np.log10(L_14[plot_mask_u_PS]), np.log10(L_OIII[plot_mask_u_PS]) , alpha=1, zorder=1, color='blue', marker='^', s=40, label='Unclassified {}'.format(PS))


    plt.ylim(5, 9)
    plt.xlim(22, 27)
    
    x_range = np.linspace(22, 27, 10000)
    plt.plot(x_range, (0.3009*(x_range)-0.72), color='black')
    plt.xlabel(r'$\log_{10}(L_{NVSS}/\rm{W\,Hz^{-1}})$', fontsize=20)
    plt.ylabel(r'$\log_{10}(L_{[OIII]}/L_{sun})$', fontsize=20)
    plt.legend()
    plt.savefig('unclassified_'+F+'.pdf')
    
    

    
    mask_belowline = (np.log10(L_OIII) < (0.3009 * np.log10(L_14) - 0.72))
    classification_source[np.where((L_14>0.) & (L_OIII>0.) & (mask_belowline) & (classification == 3))] = 4   
    classification[np.where((L_14>0.) & (L_OIII>0.) & (mask_belowline) & (classification == 3))] = 1

        
    X = len(name[  mask_F & (classification_source==4)])
    Y = len(name[np.where(mask_F & (L_14>0.) & (L_OIII>0.) & (~mask_belowline) & (classification == 3)) ])
    print('F: OIII lum: there are {0} LERGS, {1} above the line'.format(X, Y))
    x = len(name[mask_PS & (classification_source==4)])
    y = len(name[np.where(mask_PS & (L_14>0.) & (L_OIII>0.) & (~mask_belowline) & (classification == 3)) ])
    print('PS: OIII lum: there are {0} LERGS, {1} above the line'.format(x, y))
    print('')

    L = len(name[(classification==1) & mask_F])
    M = len(name[(classification==2) & mask_F])
    N = len(name[(classification==3) & mask_F])
    print('In total, there are {0} LERGS, {1} HERGS, and {2} unclassified or LERG in {3}'.\
          format(L, M, N, F))
        
    l = len(name[(classification==1) & mask_PS])
    m = len(name[(classification==2) & mask_PS])
    n = len(name[(classification==3) & mask_PS]) 
    print('In total, there are {0} LERGS, {1} HERGS, and {2} unclassified or LERG in {3}'.\
          format(l, m, n, PS))  
    
    print('')
    z_mask = ((tbdata['z_2']<0.3) & (tbdata['z_2']>0.)) | ((tbdata['Z_3']<0.3)& (tbdata['Z_3']>0.))# | ((tbdata['Z_4']<0.3)& (tbdata['Z_4']>0.)) 
    
    O = len(name[(classification==1) & mask_F & z_mask])
    P = len(name[(classification==2) & mask_F & z_mask])
    Q = len(name[(classification==3) & mask_F & z_mask])
    print('In total with z<0.3, there are {0} LERGS, {1} HERGS, and {2} unclassified or LERG in {3}'.\
          format(O, P, Q, F))
        
    o = len(name[(classification==1) & mask_PS & z_mask])
    p = len(name[(classification==2) & mask_PS & z_mask])
    q = len(name[(classification==3) & mask_PS & z_mask]) 
    print('In total with z<0.3, there are {0} LERGS, {1} HERGS, and {2} unclassified or LERG in {3}'.\
          format(o, p, q, PS))  
    print('')
    print('--------------------------------------------------------------------')

    

def ratio(x,y,z, x_err=0, y_err=0, z_err=0):
    R = x/(x+y+z)
    if x_err == 0:
        x_err = np.sqrt(x)
    if y_err == 0:  
        y_err = np.sqrt(x)
    if z_err==0:
        z_err = np.sqrt(x)    
    err_R = np.sqrt((x_err**2 * (y+z)**2 + x**2 * (y_err**2 + z_err**2))/(x+y+z)**4)
    return R, err_R

print(ratio( 6,10,2, y_err= 6.891, z_err=0.708, x_err=3.620))
      
      