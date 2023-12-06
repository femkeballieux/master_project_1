import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2


MPS = np.array([0.7452,1.232,1.235,1.37,1.708])
err_MPS = np.array([0.253,0.122,0.116,0.297,0.331])

GPS = np.array([0.7443,1.119,1.164,1.1,1.135])
err_GPS = np.array([0.156,0.028,0.0256,0.173,0.0942])
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
plt.plot(z_range, y(z_range, popt_GPS[0], popt_GPS[1]), color='black', linestyle='-.', label='Sloped fit')
plt.plot(z_range, linear(z_range, popt_GPS_2[0]), color='black', linestyle=':', label='Horizontal fit')
plt.plot(z_range, y(z_range, popt_MPS[0], popt_MPS[1]), color='red', linestyle='-.')
plt.plot(z_range, linear(z_range, popt_MPS_2[0]), color='red', linestyle=':')
plt.legend()
plt.savefig('median_offset.pdf')