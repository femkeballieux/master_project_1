import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2


MPS = np.array([0.8084,1.26,1.277,1.372,1.719])
err_MPS = np.array([0.271,0.125,0.115,0.298,0.336])
GPS = np.array([0.9166,1.229,1.248,1.085,1.235])
err_GPS = np.array([0.18,0.0319,0.0286,0.185,0.103])
z_central = np.array([0.05,0.3,0.75,1.25,2.25])
z_err = (0.05,0.2,0.25,0.25,0.75)


def sigma(GPS,MPS, err_GPS, err_MPS):
    return np.sqrt((err_GPS/MPS)**2+(GPS * err_MPS / MPS**2)**2)


ratio = GPS/MPS
ratio_err=sigma(GPS,MPS, err_GPS, err_MPS)

print(ratio)
print(ratio_err)


plt.errorbar(z_central,MPS, yerr=err_MPS, xerr=z_err, fmt='', color='red', label='MPS', markersize='7', marker='o', linestyle='')
plt.errorbar(z_central,GPS, yerr=err_GPS, xerr=z_err, fmt='', label='GPS', markersize='7', marker='s', linestyle='')
plt.xlabel('z')
plt.ylabel(r'$\Delta \log_{10}(\Phi)$')


def y(x,a,b):
    return a * x + b


def linear(x, b):
    return 0* x+ b

popt_GPS, pcov_GPS = curve_fit(y, z_central, GPS, sigma = err_GPS)
popt_MPS, pcov_MPS = curve_fit(y, z_central, MPS, sigma = err_MPS)

print('Best fit for GPS is y ={0:.3} (+/- {2:.3})* x + {1:.3}(+/-{3:.3})'.format(popt_GPS[0],\
                popt_GPS[1], np.sqrt(np.diag(pcov_GPS))[0], np.sqrt(np.diag(pcov_GPS))[1]))
print('Best fit for MPS is y ={0:.3} (+/- {2:.3})* x + {1:.3}(+/-{3:.3})'.format(popt_MPS[0],\
                popt_MPS[1], np.sqrt(np.diag(pcov_MPS))[0], np.sqrt(np.diag(pcov_MPS))[1]))
    
    
chisquare_GPS = np.sum(((GPS - y(z_central, popt_GPS[0], popt_GPS[1]))/err_GPS)**2)
# p_GPS = 1 - chi2.cdf(chisquare_GPS, df=3)
print('chi^2 and normalized value for GPS', chisquare_GPS, chisquare_GPS/3)

chisquare_MPS = np.sum(((MPS - y(z_central, popt_MPS[0], popt_MPS[1]))/err_MPS)**2)
# p_MPS = 1 - chi2.cdf(chisquare_MPS, df=3)
print('chi^2 and normalized value for GPS', chisquare_MPS, chisquare_MPS/3)


popt_GPS_2, pcov_GPS_2 = curve_fit(linear, z_central, GPS, sigma = err_GPS)
chisquare_GPS_2 = np.sum((GPS - linear(z_central, popt_GPS_2[0]))**2/(err_GPS)**2)
p_GPS_2 = 1 - chi2.cdf(chisquare_GPS_2, df=4)

popt_MPS_2, pcov_MPS_2 = curve_fit(linear, z_central, MPS, sigma = err_MPS)
chisquare_MPS_2 = np.sum((MPS - linear(z_central, popt_MPS_2[0]))**2/(err_MPS)**2)
p_MPS_2 = 1 - chi2.cdf(chisquare_MPS_2, df=4)

print('Best fit for GPS is y ={0:.3}(+/-{1:.3})'.format(popt_GPS_2[0], np.sqrt(np.diag(pcov_GPS_2))[0]))
print('Best fit for MPS is y ={0:.3}(+/-{1:.3})'.format(popt_MPS_2[0], np.sqrt(np.diag(pcov_MPS_2))[0]))
print('chi^2 and normalized value for GPS', chisquare_GPS_2, chisquare_GPS/4)
print('chi^2 and normalized value for MPS', chisquare_MPS_2, chisquare_MPS/4)

z_range = np.linspace(0,3,1000)
plt.plot(z_range, y(z_range, popt_GPS[0], popt_GPS[1]), color='black', linestyle='-.', label='offset=az+b')
plt.plot(z_range, linear(z_range, popt_GPS_2[0]), color='black', linestyle=':', label='offset=b')
plt.plot(z_range, y(z_range, popt_MPS[0], popt_MPS[1]), color='red', linestyle='-.')
plt.plot(z_range, linear(z_range, popt_MPS_2[0]), color='red', linestyle=':')
plt.legend()
plt.savefig('median_offset.pdf')