# -*- coding: utf-8 -*-
"""
code written by Femke Ballieux, takes the crossmatch between master_sample and inband fluxes of LoLSS
and throws away all that is not in the master sample. Also renames the columns
run on computer
"""
import numpy as np
from astropy.io import fits
import gpscssmodels
import seds_plot_func_extreme
from tqdm import tqdm
import scipy.optimize as opt
import warnings
import time
warnings.filterwarnings("ignore")

#read in the data
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_master_inbandLoLSS.fits')
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
hdulist.close()# fit curvature models with rms cuts