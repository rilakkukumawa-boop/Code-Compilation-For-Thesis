# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 21:07:41 2026

@author: lily zaidi
"""
#libraries 
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import matplotlib.pyplot as plt

# sii fits 
sii = fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_background_subtracted.fits')
#
plt.figure()
plt.figure(figsize=(8,8))
plt.imshow(sii, origin='lower', cmap='gray',
           vmin=np.percentile(sii,5), # contrast settings 
           vmax=np.percentile(sii,99))# contrast settings
plt.title('SII image')
plt.colorbar() 
plt.show()

#background model the old image

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground 
#Astropy Background Modelling 

sigma_clip = SigmaClip(sigma=2.5)
bkg_estimator = MedianBackground()
bkg = Background2D(sii, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
print("")

print(bkg.background_median)
print(bkg.background_rms_median)

plt.imshow(bkg.background, origin='lower', cmap='Greys_r',
           interpolation='nearest')

vmin = np.nanpercentile(sii, 0.5) #contrast settings
vmax = np.nanpercentile(sii, 99.8) #contrast settings 
plt.figure(figsize=(8,8))

fixed_sii = sii - bkg.background 

plt.imshow(fixed_sii, origin='lower', cmap='Greys_r', vmin =vmin, vmax=vmax)
plt.colorbar()
plt.show()

outpath = r'C:/Users/wafiy/Downloads/LilliesXD.fits'
header = fits.getheader(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_CS_WCS_cutout790px.fits')
fits.writeto(outpath,fixed_sii, header = header,overwrite=True)

