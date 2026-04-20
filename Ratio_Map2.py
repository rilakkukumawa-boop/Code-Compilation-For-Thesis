# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 19:03:52 2026

@author: lily zaidi
"""

#Ratio Analysis Matthewson and Clarke if [SII]/Ha > 0.4 = shock excited gas 
# if [SII]/Ha < 0.4 = photoionized region 

from astropy.io import fits #reading fits
import numpy as np #array
import matplotlib.pyplot as plt #plotting

#first read the data i.e 
SII_sub = fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_background_subtracted.fits')
Ha_sub =fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_background_subtracted.fits')

#ratio = SII_sub/Ha_sub

# only keep pixels where Hα is safely above zero
ha_sigma = 4.4 # replace with printed value
sii_sigma = 9.4 # replace with printed value



mask = (Ha_sub > 2*ha_sigma) & (SII_sub > 2*sii_sigma) #standard deviation of the backound box, and it keeps anything 
#above 2 standard deviations from the background std 

# compute ratio only for valid pixels
#so mask makes true or false array, true for if this condition is satisfied, and false if anything falls under
#3 sigma

#full-like(a,fill,dtype)
ratio = np.full_like(Ha_sub, np.nan, dtype=float)
#array of Ha subtracted is filled with np.nan
ratio[mask] = SII_sub[mask] / Ha_sub[mask]



plt.figure()
plt.imshow(ratio, origin='lower', cmap='plasma', vmin=0, vmax=1)
plt.colorbar(label='[SII]/Hα ratio')
plt.title('[SII]/Hα ratio map')
plt.show()

header = fits.getheader(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_CS_WCS_cutout790px.fits')
outpath= r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ratio_Maps_w_Haheader2NDCS.fits'
fits.writeto(outpath, ratio ,header= header, overwrite=True)




