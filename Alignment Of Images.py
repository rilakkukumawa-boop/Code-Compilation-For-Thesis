# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 14:49:04 2026

@author: lily zaidi
"""
from astropy.io import fits 
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib.colors import LogNorm
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord 

#files Ha and R 

file1 = fits.open(r'C:/Users/wafiy/T.W/Calibrated SII/MasterMedian_SII900s.fits')
file2 = fits.open(r'C:/Users/wafiy/T.W/Calibrated R/MasterMedian_R70s.fits')
img_h = file1[0].data 
img_h 

img_r = file2[0].data 
img_r 


plt.figure()
plt.imshow(img_h, origin='lower', cmap='magma', vmin=np.percentile(img_h,5), vmax=np.percentile(img_h,99))
plt.title("Ha")
plt.show()

plt.figure()
plt.imshow(img_r, origin='lower', cmap='gray', vmin=np.percentile(img_r,5), vmax=np.percentile(img_r,99))
plt.title("R")
plt.show()

#open given fits image and get data from file 
fileh = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/M33Ha_300s_WCS.fits')
img_ha = fileh[0].data 
img_ha 

#get the world co-ordinates from the header file 
wcs = WCS(fileh[0].header)

#using the cutout 2D function specify a centre position in wcs
#taken from centre on Ha image. 
position = SkyCoord(23.355575*u.deg, +30.61175972*u.deg, frame='icrs')

#size = (201,33) cutout 40%
#size = (895, 899) cutouts for Hanew and Rnew
size = (266,121)

cutout = Cutout2D(img_ha,position, size, wcs= wcs)

# --- save cutout as a FITS file ---
hdu = fits.PrimaryHDU(cutout.data, header=cutout.wcs.to_header())

outpath = r"C:/Users/wafiy/T.W/Aligned Photos/fourth_proper_cutout_px_Ha.fits"
hdu.writeto(outpath, overwrite=True)





