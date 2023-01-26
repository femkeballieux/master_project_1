# !/bin/python
# Plotting sensitivity of surveys for detecting GPS sources

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import rc
from collections import OrderedDict
from astropy.modeling import models, fitting
import matplotlib.cm as cm
import gpscssmodels

#plt.style.use('style.mplstyle')
"""
for explanation on plot, see paper bu Martje, figure 1. Femkes adaptations:
    -I have added LoLSS DR1, by looking into the catalog and finding the lowest total fluxes reported
    -I have added VLASS (which is an estimate now, taken from https://arxiv.org/pdf/2102.11753.pdf). Need to still fix the 
    actual lowest observed flux. Ask Joe?
    -need to add the curve of my research once we know what we are going to be doing. 

"""
fig = plt.figure(0,figsize=(10, 10))
gs = plt.GridSpec(1,1)
ax = plt.subplot(gs[0])

#freq_survey_sing = np.array([74.,80.,85.5,150.,151,159.,160.,178.,330,365,408.,408.,843.,1400.,2700.,4850,4850.,20000.,144.,54.])
#sens_survey_sing = np.array([0.38,2.,6.,11.0e-3,0.2,10.,1.,0.7,0.02,0.4,2.,4.,0.006,2.e-3,2.,0.025,0.03,0.04,0.052e-3,13.0e-3]) * 1000
#survey_list = ['VLSSr','CCA','MSH','TGSS','7C','3C','CCA','4C','WENSS','TXS','MRC','PKS','SUMSS','NVSS','2Jy','87GB','PMN','AT20G','LoTSS', 'LoLSS']

survey_list = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/data/flux_limits_radio_surveys.csv', dtype = 'str', delimiter=',', unpack=True, skiprows = 1, usecols = 0)
freq_survey, flux_lim_survey = np.loadtxt('/net/vdesk/data2/bach1/ballieux/master_project_1/data/flux_limits_radio_surveys.csv', dtype = float, delimiter=',', unpack=True, skiprows = 1, usecols = [1,2])

flux_lim_survey *= 1000

freq_gleam = np.linspace(72,231,1000) 
sens_gleam = 20.189*freq_gleam**(-1.2411) * 1000

freq_gps = np.linspace(15,23000,4000)
# freq_gps2 = np.linspace(15,23000,1000)

flux_gps1 = gpscssmodels.singSSA(freq_gps,80.,2.4,100)
flux_gps2 = gpscssmodels.singSSA(freq_gps,300.,2.4,750)
flux_gps3 = gpscssmodels.singSSA(freq_gps,160.,2.4,190)
flux_gps4 = gpscssmodels.singSSA(freq_gps,40.,2.4,1000)

ax.plot(freq_survey,flux_lim_survey,marker='s',color='k',ls='none',markerfacecolor='none',markeredgewidth=2)
ax.plot(freq_gleam,sens_gleam)
ax.plot(freq_gps,flux_gps1,color='#377eb8',zorder=-1)
ax.plot(freq_gps,flux_gps2,ls='--',color='#ff7f00')
ax.plot(freq_gps,flux_gps3,ls='-.',color='#4daf4a')
# ax.plot(freq_gps,flux_gps4,color='darkgreen')

for i, txt in enumerate(survey_list):
    if txt in ['AT20G']:
        ax.annotate(txt, (freq_survey[i]-0.40*freq_survey[i],flux_lim_survey[i]+0.15*flux_lim_survey[i]), fontsize=12)
    elif txt in ['CCA']:
        ax.annotate(txt, (freq_survey[i]-0.30*freq_survey[i],flux_lim_survey[i]-0.27*flux_lim_survey[i]), fontsize=12)
    elif txt in ['WENSS','87GB','3C', 'MRC']:
        ax.annotate(txt, (freq_survey[i]+0.07*freq_survey[i],flux_lim_survey[i]-0.23*flux_lim_survey[i]), fontsize=12)
    elif txt in ['PMN']:
        ax.annotate(txt, (freq_survey[i]+0.07*freq_survey[i],flux_lim_survey[i]+0.23*flux_lim_survey[i]), fontsize=12)
    else:
        ax.annotate(txt, (freq_survey[i]+0.07*freq_survey[i],flux_lim_survey[i]), fontsize=12)

ax.annotate('GLEAM', (freq_gleam[120],sens_gleam[700]), rotation = -32)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([20,22000])
ax.set_ylim([0.02,1200])

ax.set_xlabel(r'Frequency (MHz)')
ax.set_ylabel(r'Limiting Flux Density (mJy)')
ax.tick_params(axis='both',which='both',top=True,right=True)
#plt.show()
plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/survey_limits_new.png')
