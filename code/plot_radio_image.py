# !/bin/python
# Plotting radio image of PSR B1508+55 and circular polarisation image as inset

import matplotlib
# mpl.style.use('classic')
import re
import aplpy
import matplotlib.pyplot as mpl
import numpy as np
from astropy.io import fits
from matplotlib import rc
from astropy import units as u
# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['serif']})

# position of B1508 in the obs

# ra_deg_src = 199.413307
# dec_deg_src = 41.262671

ra_deg_src = 199.411
dec_deg_src = 41.264

image_hdu = fits.open('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/B1315+415.fits')
image_hdu[0].data = image_hdu[0].data[:,:]*1000 #image from jy/beam to mjy per beam
image_hdu.writeto('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/B1315+415-mosaic_adjustedheader.fits',overwrite=True)
fig = mpl.figure(1,figsize=(10, 10))
f1 = aplpy.FITSFigure('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/B1315+415-mosaic_adjustedheader.fits',figure=fig)
f1.recenter(ra_deg_src,dec_deg_src,width=0.044,height=0.044)
f1.frame.set_linewidth(2)
# f1.show_colorscale(cmap='coolwarm',vmin=-0.5, vmax=1.3)
f1.show_colorscale(cmap='Blues',vmin=-0.4, vmax=1.3)
f1.add_colorbar()
f1.add_scalebar(30 * u.arcsec,lw=3)
f1.scalebar.set_corner('bottom left')
f1.scalebar.set_color('black')
f1.scalebar.set_label(r"30$''$")
f1.scalebar.set_font(size='20')
f1.tick_labels.set_xformat('hh:mm:ss.s')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.ticks.set_color('k')

f1.add_beam(fill=True,color='0.5',corner='top left')
# f1.beam(major=6./3600.,minor=6./3600.,angle=90,fill=True,color='0.5',corner='top left')
f1.beam.set_major(6./3600)
f1.beam.set_minor(6./3600)
f1.beam.set_angle(90)
f1.beam.set_color('0.2')
f1.beam.set_frame(True)

f1.axis_labels.set_xtext('Right Ascension (J2000)')
f1.axis_labels.set_ytext('Declination (J2000)')
f1.axis_labels.set_ypad(-1)
f1.tick_labels.set_font(size='16')
f1.axis_labels.set_font(size='20')
f1.ticks.set_linewidth(1)
f1.colorbar.set_axis_label_text(r'Flux density (mJy beam$^{-1}$)')
f1.colorbar.set_axis_label_font(size='20')
f1.colorbar.set_frame_linewidth(1)
f1.colorbar.set_font(size='14')
# f1.set_title('B1315+415', fontsize=20)
fig.savefig('/net/vdesk/data2/bach1/ballieux/master_project_1/plots/B1315+415_niceimage.pdf',bbox_inches='tight')
