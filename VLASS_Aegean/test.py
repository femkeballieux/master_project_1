import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.io import fits

path = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/vlass_catalog_10000.fits'



#Read in our mastersample for the coordinates
hdulist = fits.open(path)
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()

mask=np.where(tbdata['err_b']>=0.)
print(len(tbdata['err_b'][mask])/len(tbdata['err_b']))
print('you lose 20 percent approximately?')

