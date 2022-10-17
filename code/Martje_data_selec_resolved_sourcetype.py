"""
Original code written by Joe Callingham and Martje Slob, adapted by Femke Ballieux
Takes in isolated sources from LoTTs, returns unresolved sources with single or multiple sourcetype
"""
#run on vdesk
import numpy as np
from astropy.io import fits
from astropy.table import Table

filepath = "/net/vdesk/data2/bach1/ballieux/master_project_1/data/isolated_test.fits"

# Read out fits files from data directory, place data and header in Astropy table
survey = fits.open(filepath)
data = Table(survey[1].data)
header = survey[1].header
survey.close()

# Get flux data + errors from data
peak_flux = data['Peak_flux']
E_peak_flux = data['E_Peak_flux']
total_flux = data['Total_flux']
E_total_flux = data['E_Total_flux']

# Calculate the SNR value for each source of the survey
SNR_1 = (total_flux/E_total_flux)

# Calculate R value for each source, defined in Shimwell (2022)
R = np.log(total_flux/peak_flux)

# Calculate 99.9% envelope for R-value, defined in Shmiwell (2022)
x = [0.65999588, 2.49058227, 96.57331031, 0.42495155]
R_99 = abs(x[3] + (x[0] + x[3])/(1+(SNR_1/x[2])**(x[1])))

# Find all sources with R outside 99.9% envelope, remove from data
extended_index = np.where((R >= R_99))

print(np.round(100*len(extended_index[0])/len(peak_flux),2), '% of sources (', len(extended_index[0]), ') are extended and removed')

data.remove_rows(extended_index)

# Also remove 'complex' sources
S_code = data['S_Code']

complex_source = np.where(S_code == 'C')

print(np.round(100*len(complex_source[0])/len(peak_flux),2), '% of sources (', len(complex_source[0]), ') are complex and removed')

#Now we have all good sources from LoTTS DR2
data.remove_rows(complex_source[0])
data.write('/net/vdesk/data2/bach1/ballieux/master_project_1/data/unresolved_isolated_S_source.fits', format='fits', overwrite = True)