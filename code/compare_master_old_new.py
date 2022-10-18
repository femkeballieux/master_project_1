"""
Written by Femke Ballieux
takes in the crossmatchd version of the old (Martjes) and new master sample.
Checks is all sources found in hers where also in mine, other way around"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.optimize as opt
import gpscssmodels


path_laptop = 'C:/Users/Femke/Documents/GitHub/master_project_1/data'
path_vdesk= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'

#read in the data
hdulist = fits.open(path_vdesk + '/crossmatch_old_new_master.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

my_name_array= np.array(tbdata['LoTSS_name_1'])
Martje_name_array = tbdata['LoTSS_name_2']

my_alpha_low = tbdata['alpha_low_1']
my_alpha_high = tbdata['alpha_high_1']
Martje_alpha_low = tbdata['alpha_low_2']
Martje_alpha_high = tbdata['alpha_high_2']

my_flux_LoLSS = tbdata['LoLSS_flux_1']
Martje_flux_LoLSS = tbdata['LoLSS_flux_2']
my_flux_LoTSS = tbdata['LoTSS_flux_1']
Martje_flux_LoTSS = tbdata['LoTSS_flux_2']
my_flux_NVSS = tbdata['NVSS_flux_1']
Martje_flux_NVSS = tbdata['NVSS_flux_2']
print(orig_cols)


my_counter = 0 #how many PS sources do I find
Martje_counter = 0 #How many of these does Martje find
nan_count_1 = 0 # how many do I find, that Martje does not, that can be explained by falling outside of footprint

print("These are in my sample, but not in Martje and are not Nan")
for i, name in enumerate(my_name_array):
    if (my_alpha_low[i] >= 0.1) & (my_alpha_high[i]<=0): #select when a source in my master sample is PS
        my_counter +=1
        if (Martje_alpha_low[i] >= 0.1) & (Martje_alpha_high[i]<=0):#Select when it is also PS in Martjes
            Martje_counter += 1
        else: #If it is not PS in martje
            if np.isnan(Martje_flux_LoLSS[i]): #Check if it is because of nan
                nan_count_1 +=1
            else: #In mine, not in martje, not because of nan
                print(i, name,"mine:", my_flux_LoLSS[i], ' Martje: ' , Martje_flux_LoLSS[i], 'ratio:', my_flux_LoLSS[i]/Martje_flux_LoLSS[i])
print('')
print ('In my sample there are ', my_counter ,' PS sources')
print(Martje_counter, ' of these sources are also in Martjes master sample' )
print(nan_count_1, ' of these remaining ', my_counter - Martje_counter, 'had nan flux in Martjes sample')
print('We thus need to explain', my_counter-Martje_counter-nan_count_1, 'Sources')


print("")
print("These are in Martjes sample but not in mine")
Martje_counter_2 = 0 #How many sources Martje finds
my_counter_2 = 0 #How many sources of these I find
nan_count_2 = 0 #How many of the difference can be explained by being nan = falling outside of footprint
for i, name in enumerate(Martje_name_array):

    if (Martje_alpha_low[i] >= 0.1) & (Martje_alpha_high[i]<=0): #select when a source in Martjes master sample is PS
        Martje_counter_2 +=1
        if (my_alpha_low[i] >= 0.1) & (my_alpha_high[i]<=0): #Check if also PS in my sample
            my_counter_2 += 1
        else: #If not PS in my sample, check if nan
            if np.isnan(my_flux_LoLSS[i]):
                nan_count_2 +=1
            else: #If not ps in my sample, but ps in her sample, but not nan, print flux values
                print(i, name, "mine: ", my_flux_LoLSS[i], ' Martje: ' , Martje_flux_LoLSS[i], 'ratio:', my_flux_LoLSS[i]/Martje_flux_LoLSS[i])

                flux_array=np.array([my_flux_LoTSS[i],my_flux_NVSS[i]] )
                freq_array=np.array([144., 1400.])
                plt.figure(figsize=(10,8))
                plt.scatter(freq_array, flux_array, label='LoTSS, NVSS', color='black')
                plt.scatter(54., Martje_flux_LoLSS[i], label='old LoLSS')
                plt.scatter(54., my_flux_LoLSS[i], label='new LoLSS')
                plt.xscale('log')
                plt.yscale('log')
                plt.title(my_name_array[i])
                plt.xlabel('freq in MHZ')
                plt.ylabel('flux in Jy')
                plt.legend()
            
                plt.savefig(path_vdesk + 'compare_sed_old_new/' + my_name_array[i] + '.pdf', bboxinches='tight')
                
print('')
print ('In Martjes sample there are ', Martje_counter_2 ,' PS sources')
print(my_counter_2, ' of these sources are also in my master sample' )
print(nan_count_2, 'of these remaining ', Martje_counter_2 - my_counter_2, 'had nan flux in my sample')
print("we thus need to explain", Martje_counter_2-my_counter_2-nan_count_2, 'sources')

 
