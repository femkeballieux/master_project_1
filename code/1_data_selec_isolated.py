"""
Original code written by Joe Callingham and Martje Slob, adapted by Femke Ballieux
Takes in LoTTS fits file, returns a file with isolated sources
"""
#Run on vdesk
import numpy as np
from astropy.io import fits
from astropy.table import Table
import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm

filepath = "/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/LoTSS_DR2.fits"

# Read out LoTTS DR2 fits files from data directory, place data and header in Astropy table
survey = fits.open(filepath, memmap=True)
data = Table(survey[1].data)
survey.close()
print(data.colnames)
data.sort('RA')

#This is the peak flux, in later stages we use total flux 
peak_flux = data['Peak_flux']

# Get position data
sample = len(peak_flux) - 1
RA = np.array(data['RA'])
DEC = np.array(data['DEC'])

# Find sources within 47''
max_sep = 47 * u.arcsec # arcsec

unisolated_index = []
border = 2000

coords = SkyCoord(RA*u.deg, DEC*u.deg, frame = 'icrs')

#For_loop that checks if no 2 sources are within 47 arcsecs from eachother
start = time.time()
for index, i in enumerate(tqdm(RA)):
    # You can change sample size to test the code more quickly
    if index > sample:
        break

    # Take into account circular coordinates
    if index - border < 0:
        RA_i = np.concatenate((RA[index - border:], RA[:index + border]))
        DEC_i = np.concatenate((DEC[index - border:], DEC[:index + border]))
        peak_flux_i = np.concatenate((peak_flux[index - border:], peak_flux[:index + border]))
    elif index + border > len(RA):
        RA_i = np.concatenate((RA[index-border:], RA[:index+border-len(RA)]))
        DEC_i = np.concatenate((DEC[index-border:], DEC[:index+border-len(RA)]))
        peak_flux_i = np.concatenate((peak_flux[index-border:], peak_flux[:index+border-len(RA)]))
    else:
        RA_i = RA[index-border:index+border]
        DEC_i = DEC[index-border:index+border]
        peak_flux_i = peak_flux[index-border:index+border]
        # name_i = name[index-border:index+border]

    # Check if any source has an absolute distance <47'' in RA or DEC from a source
    coords_source = coords[index]
    coords_match = SkyCoord(RA_i*u.deg, DEC_i*u.deg, frame = 'icrs')
    d2d = coords_source.separation(coords_match)
    matchmsk = (d2d < max_sep)
    remove_self = (d2d != 0 * u.arcsec)
    isolate_mask = matchmsk * remove_self

    # remove sources from catalogue
    if np.any(isolate_mask):   
        flux_source = peak_flux[index]
        neighbour_flux = peak_flux_i[isolate_mask]
        for number, f in enumerate(neighbour_flux):
            if f/flux_source > 0.1:
                unisolated_index.append(index)
                break
    if index % 10000 == 0:
        end = time.time()
        print(index, end - start)
end = time.time()

print(np.round(100*len(unisolated_index)/sample,2), '% of sources (', len(unisolated_index), ') are unisolated and removed')
print('script takes', end - start, 'seconds to run')

#Remove the LoTTS data that overlaps, rewrite everything else to a new file
data.remove_rows(unisolated_index)
data.write('/net/vdesk/data2/bach1/ballieux/master_project_1/data/isolated_test.fits', format='fits', overwrite = True)
