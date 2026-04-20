# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:09:33 2026

@author: lily zaidi
"""

#Ratio Analysis Matthewson and Clarke if [SII]/Ha > 0.4 = shock excited gas 
# if [SII]/Ha < 0.4 = photoionized region 

from astropy.io import fits #reading fits
import numpy as np #array
import matplotlib.pyplot as plt #plotting

SII_sub = fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_CS_WCS_cutout790px.fits')
Ha_sub =fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_CS_WCS_cutout790px.fits')

ratio= SII_sub / Ha_sub

plt.figure()
plt.imshow(ratio, origin='lower', cmap='plasma', vmin=0, vmax=1)
plt.colorbar(label='[SII]/Hα ratio')
plt.title('[SII]/Hα ratio map')
plt.show()

header = fits.getheader(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_CS_WCS_cutout790px.fits')
outpath= r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ratio_Maps_w_Haheader.fits'
fits.writeto(outpath, ratio ,header= header, overwrite=True)

