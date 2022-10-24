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
hdulist = fits.open(path_laptop + '/compare_old_new_LoLSS/crossmatch_old_new_LoLSS.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

#read in the columns
new_name_array= np.array(tbdata['Source_name_1'])
new_flux_LoLSS = tbdata['Total_flux_LoLLS']
old_flux_LoLSS = tbdata['Total_flux_LoLSS']

#get the ratios of all the sources that are present in both samples
ratio_list = []
for i, name in enumerate(new_name_array):
    ratio = new_flux_LoLSS[i]/old_flux_LoLSS[i]
    ratio_list.append(ratio)

    
ratio_array= np.array(ratio_list)
plt.figure(figsize=(10,8))
plt.hist(ratio_array, bins=70)
plt.title('All LoLSS data present in both samples')
plt.xlabel('flux DR1/ flux PDR')
plt.ylabel('Number of sources')
plt.savefig(path_laptop+'/compare_old_new_LoLSS/hist_all_sources.pdf', bboxinches='tight')

print('max:',max(ratio_array), 'min:', min(ratio_array))
print('mean:', np.std(ratio_array))
print('std', np.mean(ratio_array))