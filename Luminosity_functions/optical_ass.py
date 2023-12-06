"""
Code written by Femke Ballieux, takes in the master sample and crossmatches it with optical components
"""
#run on laptop
import os
import astropy.units as u
# from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.io import fits
import time
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.pyplot as plt

"""
Here is how the entire optical association process works:
Making a combination of the LoTSS DR1 and DR2 optical catalogues, concatination in topcat
Added an empty column to DR1 with the same datatype as z_source for DR2, 
so that for both catalog the source can be written. I only use the redshift from this, 
but need the z_source to select on te different types of redshifts. Bit convoluted but this is the easy fix.
 
Crossmatch master sample to LoTSS optical ID's, in topcat as well. Exact matching to names 
You might need to change the string type such that they match .

In topcat also select the proper columns, and only rows with a redshift and save as the clean version
 as "master_sample_with_redshift"

Get SDSS i,g and r bands for our master sample with a redshift, in code. I chose radius = 2 arcsec, 
Did DR16 since this corresponds to the DR2 LoTSS optical catalog. Put this into a new table, easiest in topcat.
Then crossmatch this to the master sample with a radius of 2'', and after this do v_max.py code from Martje.

"""
def obtain_SDSS():
    plt.style.use('/net/vdesk/data2/bach1/ballieux/master_project_1/code/style.mplstyle')
    hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_sample_with_redshift.fits')
    tbdata = hdulist[1].data
    hdulist.close()
    
    length=len(tbdata['RA_1']) #How many sources
    coordinates=coords.SkyCoord(tbdata['RA_1'], tbdata['Dec_1'], frame='icrs', unit='deg') #All the coordinates in our file
    
    #Initialize arrays
    ra = []
    dec = []
    g = []
    err_g = []
    i = []
    err_i = []
    r = []
    err_r = []
    start=time.time() 
    batch_size = 40  # Set the number of coordinates to query in each batch
    
    for j in range(0, length, batch_size):
    
        if j+batch_size<=length:
            coords_batch = coordinates[j:j+batch_size]
        else:
            coords_batch = coordinates[j:length-1]
            
        # xid = SDSS.query_region(coords_batch, radius = 1 * u.arcsec, photoobj_fields=['ra', 'dec', 
        # 'psfMag_g', 'psfMagErr_g', 'psfMag_i', 'psfMagErr_i' ,'psfMag_r', 'psfMagErr_r'],
        #                          data_release=16)
        xid = SDSS.query_region(coords_batch, radius = 2.5 * u.arcsec, fields=['ra', 'dec', 
        'psfMag_g', 'psfMagErr_g', 'psfMag_i', 'psfMagErr_i' ,'psfMag_r', 'psfMagErr_r'], data_release=16)
    
        print(j)                          
                              
        # xid2 = SDSS.query_region(coords_batch, radius = 1 * u.arcsec,                           
        # specobj_fields = ['oiii_5008_flux', 'oiii_5008_flux_err', 'oiii_5008_eqw',
        # 'h_beta_flux', 'h_beta_flux_err','h_alpha_flux', 'h_alpha_flux_err',
        # 'nii_6584_flux', 'nii_6584_flux_err', 'sii_6717_flux', 'sii_6717_flux_err',
        # 'sii_6731_flux', 'sii_6731_flux_err', 'oi_6300_flux', 'oi_6300_flux_err'], data_release=17)
        # xid2 = SDSS.query_region(coords_batch, radius = 1 * u.arcsec,                          
        # specobj_fields = ['h_alpha_flux'], data_release=18)
        # print(xid2)
    
        try:
    
            if len(ra)==0: #For the first index, the new array is just this
                ra = np.array(xid['ra'].value)
                dec = np.array(xid['dec'].value)
                g = np.array(xid['psfMag_g'].value)
                err_g = np.array(xid['psfMagErr_g'].value)
                i = np.array(xid['psfMag_i'].value)
                err_i = np.array(xid['psfMagErr_i'].value)
                r = np.array(xid['psfMag_r'].value)
                err_r = np.array(xid['psfMagErr_r'].value)
            else: #For the next indeces, just concatenate the arrays
                ra = np.concatenate((ra, np.array(xid['ra'].value)))
                dec = np.concatenate((dec, np.array(xid['dec'].value)))
                g = np.concatenate((g, np.array(xid['psfMag_g'].value)))
                err_g = np.concatenate((err_g, np.array(xid['psfMagErr_g'].value)))
                i = np.concatenate((i, np.array(xid['psfMag_i'].value)))
                err_i = np.concatenate((err_i, np.array(xid['psfMagErr_i'].value)))
                r = np.concatenate((r, np.array(xid['psfMag_r'].value)))
                err_r = np.concatenate((err_r, np.array(xid['psfMagErr_r'].value)))
        except Exception as error:
            print("An exception occurred:", error) # An exception occurred: division by zero
    
    
    end=time.time()
    print("This process took {:.5} seconds".format(end-start))
    
    print(len(ra), len(dec), len(g), len(err_g))
    col2  = fits.Column(name='RA', format = 'E', array = np.array(ra))
    col3  = fits.Column(name='Dec', format = 'E', array = np.array(dec))
    col4  = fits.Column(name='g', format = 'E', array = np.array(g))
    col5  = fits.Column(name='err_g', format = 'E', array = np.array(err_g))
    col6  = fits.Column(name='i', format = 'E', array = np.array(i))
    col7  = fits.Column(name='err_i', format = 'E', array = np.array(err_i))
    col8  = fits.Column(name='r', format = 'E', array = np.array(i))
    col9  = fits.Column(name='err_r', format = 'E', array = np.array(err_i))
    
    
    
    cols = fits.ColDefs([ col2, col3, col4, col5, col6, col7, col8, col9])
    tbhdu = fits.BinTableHDU.from_columns(cols)  
    print("#----------------------------------------------------------#")
    print('Saving to a fits file.')  
    
    tbhdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/SDSS_sample.fits', overwrite = True)
    
hdulist = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_redshift.fits')
tbdata = hdulist[1].data
hdulist.close()
    #Make the redshift distribution for PS sources in 2 samples
z_best=tbdata['z_best']  
PS_MPS = np.where((tbdata['alpha_low_LoLSS'] > tbdata['e_alpha_low_LoLSS'])
                            & (tbdata['alpha_high_LoLSS'] < 0))    
PS_GPS = np.where((tbdata['alpha_low_VLASS'] > tbdata['e_alpha_low_VLASS'])
                            & (tbdata['alpha_high_VLASS'] <= 0))  


host = host_subplot(111)
twin = host.twinx()

# Plot the first histogram on the left y-axis
plt.hist(z_best[PS_MPS], bins=35, histtype='step', label='MPS', color='black')

plt.ylabel('Number of sources in MPS sample')


# Plot the second histogram on the right y-axis
twin.hist(z_best[PS_GPS], bins=35, histtype='step', label='GPS', color='red')

plt.xlabel('z')
twin.set_ylabel('Number of sources in GPS sample')

# Add legend for both histograms
plt.legend(['MPS', 'GPS'], loc='upper right', fontsize=14)

host.yaxis.get_label().set_color('red')

plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/redshift_dist.pdf')

print('MPS number of sources with redshift', len(z_best[PS_MPS]))
print('GPS number of sources with redshift', len(z_best[PS_GPS]))

# obtain_SDSS()



# def crossmatch_SDSS():
#     os.system('java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts \
#               tmatchn join1=always matcher=sky multimode=pairs nin=2 params=2 \
#         in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_redshift.fits values1="RA_1 DEC_1" \
#         in2=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/SDSS_sample.fits values2="ra dec" \
#         out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/master_redshift_SDSS.fits')
       
# crossmatch_SDSS()        