"""
Written by Femke Ballieux
takes in the crossmatched version of the preliminary and first LoLSS data release
checks for discrepancy
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#code can be run from eother laptop or pc
path_laptop = 'C:/Users/Femke/Documents/GitHub/master_project_1/data'
path_vdesk= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'

#read in the data
hdulist = fits.open(path_vdesk + '/compare_old_new_LoLSS/crossmatch_old_new_LoLSS.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()
print(orig_cols)

#read in the columns
new_name_array= np.array(tbdata['Source_name_1'])
new_flux_LoLSS = tbdata['Total_flux_LoLLS'] #all still in Mjy 
old_flux_LoLSS = tbdata['Total_flux_LoLSS']
new_Isl_rms = tbdata['Isl_rms_1']
old_Isl_rms = tbdata['Isl_rms_2']
S_code_new = tbdata['S_Code_1']
S_code_old = tbdata['S_Code_2']



#get the ratios of all the sources that are present in both samples
#SNR and index of the extremes as well
ratio_list = []
new_SNR = []
old_SNR = []
extreme_index=[]
for i, name in enumerate(new_name_array):
    if S_code_new[i] == 'S':
        ratio = new_flux_LoLSS[i]/old_flux_LoLSS[i]
        ratio_list.append(ratio)
        new_SNR.append(new_flux_LoLSS[i]/new_Isl_rms[i])
        old_SNR.append(old_flux_LoLSS[i]/old_Isl_rms[i])
        
        if ratio >= 3:
            extreme_index.append(i)
        elif ratio <= (1/3):
            extreme_index.append(i)

       
#Make a new table of just the extremes so we can check in topcat
col1 = fits.Column(name='LoTSS_name', format = '34A', array = new_name_array[np.array(extreme_index)])
col2 = fits.Column(name='RA', format = 'E', array= tbdata['RA_1'][np.array(extreme_index)])
col3 = fits.Column(name='Dec', format = 'E', array = tbdata['Dec_1'][np.array(extreme_index)])
col4 = fits.Column(name='old_flux_LoLSS', format = 'E', array = old_flux_LoLSS[np.array(extreme_index)])
col5 = fits.Column(name='new_flux_LoLSS', format = 'E', array = new_flux_LoLSS[np.array(extreme_index)])

cols = fits.ColDefs([col1, col2, col3,col4, col5])
tbhdu = fits.BinTableHDU.from_columns(cols)  
print("#----------------------------------------------------------#")
print('Saving to a fits file.')  

tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/compare_old_new_LoLSS/locations_extreme_outliers.fits', overwrite = True)


#Ratio histogram   
ratio_array= np.array(ratio_list)
plt.figure(figsize=(10,8))
plt.hist(ratio_array, bins=70)
plt.title('All LoLSS data present in both samples, only S')
plt.xlabel('flux DR1/ flux PDR')
plt.ylabel('Number of sources')
plt.savefig(path_vdesk+'/compare_old_new_LoLSS/hist_all_sources.pdf', bboxinches='tight')


plt.figure(figsize=(10,8))
plt.scatter( new_SNR, ratio_list, alpha=0.1)
plt.title('All of LoLSS, only S')
plt.xlabel('SNR DR1')
plt.xscale('log')
plt.ylabel('ratio flux DR1/PDR')
plt.hlines(1, xmin=0, xmax=5000)
plt.xlim(np.min(new_SNR),5000)
plt.savefig(path_vdesk+'/compare_old_new_LoLSS/DR1_SNR_ratio.pdf', bboxinches='tight')

plt.figure(figsize=(10,8))
plt.scatter( old_SNR, ratio_list, alpha=0.1)
plt.title('All of LoLSS, only S')
plt.xlabel('SNR PDR')
plt.xscale('log')
plt.ylabel('ratio flux DR1/PDR')
plt.hlines(1, xmin=0, xmax=1500)
plt.xlim(np.min(old_SNR),1200)
plt.savefig(path_vdesk+'/compare_old_new_LoLSS/PDR_SNR_ratio.pdf', bboxinches='tight')

print('max:',max(ratio_array), 'min:', min(ratio_array))
print('mean:', np.std(ratio_array))
print('std', np.mean(ratio_array))