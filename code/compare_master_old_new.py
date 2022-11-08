"""
Written by Femke Ballieux
takes in the crossmatched version of the old (Martjes) and new master sample.
Checks is all sources found in hers where also in mine, other way around
So sorry for this mess of a code...
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#paths to run from 2 different machines
path_laptop = 'C:/Users/Femke/Documents/GitHub/master_project_1/data'
path_vdesk= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'

path = path_vdesk

#read in the data
hdulist = fits.open(path + '/crossmatch_old_new_master.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#get the right columns
my_name_array= np.array(tbdata['LoTSS_name_1'])
Martje_name_array = tbdata['LoTSS_name_2']

RA = tbdata['RA_1']
Dec= tbdata['Dec_1']

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

new_LoLSS_rms = tbdata['LoLSS_rms'] /1000 #here it is still mjy. In master_inbandLoLSS it is fixed
new_LoLSS_maj = tbdata['LoLSS_maj']
new_LoLSS_S_code = tbdata['LoLSS_S_code']
print(orig_cols)

"""
This for-loop runs over all my sources that are PS
"""
my_counter = 0 #how many PS sources do I find
Martje_counter = 0 #How many of these does Martje find
nan_count_1 = 0 # how many do I find, that Martje does not, that can be explained by falling outside of footprint

#Some empty list to store samples
ratio_in_mine = []
ratio_in_both = []
SNR_in_both = []
RA_in_mine=[]
Dec_in_mine=[]
SNR_in_mine = []
print("These are in my sample, but not in Martje and are not Nan")
for i, name in enumerate(my_name_array):
    if (my_alpha_low[i] >= 0.1) & (my_alpha_high[i]<=0): #select when a source in my master sample is PS
        my_counter +=1
        if (Martje_alpha_low[i] >= 0.1) & (Martje_alpha_high[i]<=0):#Select when it is also PS in Martjes
            Martje_counter += 1
            ratio_in_both.append(my_flux_LoLSS[i]/Martje_flux_LoLSS[i]) #These are PS in both
            SNR_in_both.append(my_flux_LoLSS[i]/new_LoLSS_rms[i])
        else: #If it is not PS in martje
            if np.isnan(Martje_flux_LoLSS[i]): #Check if it is because of nan
                nan_count_1 +=1
            else: #In mine, not in martje, not because of nan
                ratio_in_mine.append(my_flux_LoLSS[i]/Martje_flux_LoLSS[i])
                RA_in_mine.append(RA[i])
                Dec_in_mine.append(Dec[i])
                SNR_in_mine.append(my_flux_LoLSS[i]/new_LoLSS_rms[i])
                
print('')
print ('In my sample there are ', my_counter ,' PS sources')
print(Martje_counter, ' of these sources are also in Martjes master sample' )
print(nan_count_1, ' of these remaining ', my_counter - Martje_counter, 'had nan flux in Martjes sample')
print('We thus need to explain', my_counter-Martje_counter-nan_count_1, 'Sources')



"""
This figure plots a histogram of all the sources that are PS in my sample, 
not Martjes sample
"""
# plt.figure(figsize=(10,8))
# plt.hist(ratio_in_mine, bins=8)
# plt.title('PS in DR1, not in PDR')
# plt.xlabel('flux DR1/ flux PDR')
# plt.ylabel('Number of sources')
# plt.savefig(path+'/compare_old_new_LoLSS/in_mine_not_Martje.pdf', bboxinches='tight')
# print('mean:',np.mean(ratio_in_mine),'Std in mine:',np.std(ratio_in_mine))


"""
This for-loop runs over all Martjes sources that are PS
"""
print("")
print("These are in Martjes sample but not in mine")
Martje_counter_2 = 0 #How many sources Martje finds
my_counter_2 = 0 #How many sources of these I find
nan_count_2 = 0 #How many of the difference can be explained by being nan = falling outside of footprint
#some empty lists to store samples
ratio_in_Martje = []
RA_in_Martje=[]
Dec_in_Martje=[]
SNR_in_Martje = []

for i, name in enumerate(Martje_name_array):

    if (Martje_alpha_low[i] >= 0.1) & (Martje_alpha_high[i]<=0): #select when a source in Martjes master sample is PS
        Martje_counter_2 +=1
        if (my_alpha_low[i] >= 0.1) & (my_alpha_high[i]<=0): #Check if also PS in my sample
            my_counter_2 += 1
        else: #If not PS in my sample, check if nan
            if np.isnan(my_flux_LoLSS[i]):
                nan_count_2 +=1
            else: #If not ps in my sample, but ps in her sample, but not nan, print flux values
                #print(i, name, "mine: ", my_flux_LoLSS[i], ' Martje: ' , Martje_flux_LoLSS[i], 'ratio:', my_flux_LoLSS[i]/Martje_flux_LoLSS[i])
                ratio_in_Martje.append(my_flux_LoLSS[i]/Martje_flux_LoLSS[i])
                RA_in_Martje.append(RA[i])
                Dec_in_Martje.append(Dec[i])
                SNR_in_Martje.append(my_flux_LoLSS[i]/new_LoLSS_rms[i])
                
print('')
print ('In Martjes sample there are ', Martje_counter_2 ,' PS sources')
print(my_counter_2, ' of these sources are also in my master sample' )
print(nan_count_2, 'of these remaining ', Martje_counter_2 - my_counter_2, 'had nan flux in my sample')
print("we thus need to explain", Martje_counter_2-my_counter_2-nan_count_2, 'sources')


"""
This plot makes the positions of the discrepancy sources in the sky
"""
# plt.figure(0,figsize=(10,8))
# plt.title('RA, Dec', fontsize=14)
# plt.scatter(np.array(RA_in_mine), np.array(Dec_in_mine), label='Is PS using DR1, not using PDR')
# plt.scatter(np.array(RA_in_Martje), np.array(Dec_in_Martje), label='Is PS using PDR, not using DR1')
# plt.xlabel('RA', fontsize=14)
# plt.ylabel('Dec', fontsize=14)
# plt.legend()
# plt.savefig(path + '/compare_old_new_LoLSS/RA_Dec.pdf', bboxinches='tight')


"""
This plot makes a histogram for all sources that are PS in martjes, not in mine
"""
# plt.figure(figsize=(10,8))
# plt.hist(ratio_in_Martje, bins=17)
# plt.title('PS in PDR, non-PS in DR1')
# plt.xlabel('flux DR1/ flux PDR')
# plt.ylabel('Number of sources')
# plt.savefig(path+'/compare_old_new_LoLSS/in_Martje_not_mine.pdf', bboxinches='tight')
# print('mean:',np.mean(ratio_in_Martje),'Std In Martje:',np.std(ratio_in_Martje))
 

"""
This loop runs over all sources that are in both master samples and not nan. 
So more intrinsic to LoLSS, not only our PS sources
"""
ratio_isolated = []
my_flux_isolated = []
Martje_flux_isolated = []
SNR_isolated = []
for i, name in enumerate(my_name_array):
    if my_flux_LoLSS[i] >= 0.:
        if Martje_flux_LoLSS[i] >=0.:
            ratio_isolated.append(my_flux_LoLSS[i]/Martje_flux_LoLSS[i])
            my_flux_isolated.append(my_flux_LoLSS[i])
            Martje_flux_isolated.append(Martje_flux_LoLSS[i])
            SNR_isolated.append(my_flux_LoLSS[i]/new_LoLSS_rms[i])
            

ratio_isolated_array = np.array(ratio_isolated)           
print('isolated:', np.mean(ratio_isolated_array))
print('Std:',np.std(ratio_isolated_array))

"""
This plot shows a histogram for all isolated sources
"""
# plt.figure(figsize=(10,8))
# plt.hist(ratio_isolated_array, bins=70)
# plt.title('All isolated, unresolved sources present in both samples')
# plt.xlabel('flux DR1/ flux PDR')
# plt.ylabel('Number of sources')
# plt.savefig(path+'/compare_old_new_LoLSS/hist_isolated.pdf', bboxinches='tight')

"""
This plot shows the SNR-ratio plot for different populations of sources
"""
plt.figure(figsize=(10,8))
plt.scatter(SNR_in_mine, ratio_in_mine, zorder=1, label = 'Only PS in DR1')
plt.scatter(SNR_in_Martje, ratio_in_Martje, zorder=1, label='Only PS in PDR')
plt.scatter(SNR_isolated, ratio_isolated, zorder=0.1, alpha=0.5, label='all isolated sources')
plt.scatter(SNR_in_both, ratio_in_both, zorder=0.5, alpha=0.9, label='PS in both')
plt.title('All isolated sources with flux in PDR and DR1')
plt.xlabel('SNR DR1')
plt.xscale('log')
plt.ylabel('ratio flux DR1/PDR')
plt.hlines(1, xmin=0, xmax=5000)
plt.legend()
#plt.xlim(np.min(new_SNR),5000)
plt.savefig(path_vdesk+'/compare_old_new_LoLSS/DR1_SNR_ratio_Martje_vs_mine_both.pdf', bboxinches='tight')



"""
This plot shows all isolated, unresolved sources present in both samples
ratio vs flux DR1
"""
# plt.figure(figsize=(10,8))
# plt.scatter( ratio_isolated_array,my_flux_isolated,)
# plt.title('All isolated, unresolved sources present in both samples')
# plt.ylabel('flux DR1')
# plt.xlabel('ratio flux DR1/PDR')
# plt.vlines(1, ymin=0, ymax=18)
# plt.savefig(path+'/compare_old_new_LoLSS/DR1_flux_ratio.pdf', bboxinches='tight')


"""
This plot shows for all isolated, unresolved sources present in both samples
ratio vs flux PDR
"""
# plt.figure(figsize=(10,8))
# plt.scatter( ratio_isolated_array,Martje_flux_isolated,)
# plt.title('All isolated, unresolved sources present in both samples ')
# plt.ylabel('flux PDR')
# plt.xlabel('ratio flux DR1/PDR')
# plt.vlines(1, ymin=0, ymax=18)
# plt.savefig(path+'/compare_old_new_LoLSS/PDR_flux_ratio.pdf', bboxinches='tight')





