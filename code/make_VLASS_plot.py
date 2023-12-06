
import numpy as np
from astropy.io import fits
from astropy.table import Table
import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import matplotlib.pyplot as plt

print('hoi')
plt.style.use('/net/vdesk/data2/bach1/ballieux/master_project_1/code/style.mplstyle')
filepath = "/net/vdesk/data2/bach1/ballieux/master_project_1/data/official_VLASS_no_duplicates.fits"

survey = fits.open(filepath, memmap=True)
data = Table(survey[1].data)
survey.close()

NN=data['NN_dist']
maj=data['Maj']

plt.figure(figsize=(10, 6))

mask1 = np.where(NN > 0.)
mask2 = np.where((NN > 0.) & (maj > 3))
mask3 = np.where((NN > 0.) & (maj < 3))

bin_edges = np.logspace(np.log10(NN[mask1].min()), np.log10(NN[mask1].max()), num=101)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

plt.hist(NN[mask1], bins=bin_edges, histtype='step', label='All sources')
plt.hist(NN[mask2], bins=bin_edges, histtype='step', label='Major axis > 3\"')
plt.hist(NN[mask3], bins=bin_edges, histtype='step', label='Major axis < 3\"')

plt.ylabel('Number of sources')
plt.xlabel('Distance to nearest neighbour (\")')
plt.xscale('log')
#plt.vlines(2.9, 0, 62000, label='typical VLASS beam size', linestyle='dotted')
plt.vlines(47, 0, 67000, label='Isolation limit of 47\"', linestyle='dashed')
plt.legend(loc='upper left')
plt.xlim(2,1e3)
plt.ylim(0,65000)
plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/NN_dist.pdf')
plt.show()