# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:23:18 2026

@author: lily zaidi
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

Ha, ha_header = fits.getdata(
    r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/HA_CS_WCS_cutout790px.fits',
    header=True
)
Ha = Ha.astype(float)

SII = fits.getdata(r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_CS_WCS_cutout790px.fits').astype(float)

# full image
Ha_source = Ha
SII_source = SII

# background box
bx1, bx2 = 600,660
by1, by2 = 150,200


Ha_bg = Ha[by1:by2, bx1:bx2]
SII_bg = SII[by1:by2, bx1:bx2]

#Ha display
fig, ax = plt.subplots(figsize=(8,8))
ax.imshow(Ha, origin='lower', cmap='gray',
          vmin=np.percentile(Ha,5),
          vmax=np.percentile(Ha,99))

#overlay of cyan box on the Ha display 
import matplotlib.patches as patches
bg_rect = patches.Rectangle((bx1, by1), bx2-bx1, by2-by1,
                            edgecolor='cyan', facecolor='none', linewidth=2)
ax.add_patch(bg_rect)
plt.show()

#SII display
fig, ax = plt.subplots(figsize=(8,8))
ax.imshow(SII, origin='lower', cmap='gray', 
          vmin = np.percentile(SII, 5), 
          vmax = np.percentile(SII,99))
bg_rect2 = patches.Rectangle((bx1, by1), bx2-bx1, by2-by1,
                            edgecolor='cyan', facecolor='none', linewidth=2)
ax.add_patch(bg_rect2)
plt.show()




# check if background got subrtracted correctly 
Ha_image = Ha_source - np.mean(Ha_bg)
SII_image = SII_source - np.mean(SII_bg)

fig, ax = plt.subplots(1, 2, figsize=(12,6))

ax[0].imshow(Ha_image, origin='lower', cmap='gray',
             vmin=np.percentile(Ha_image, 5),
             vmax=np.percentile(Ha_image, 99))
ax[0].set_title("Hα")

ax[1].imshow(SII_image, origin='lower', cmap='gray',
             vmin=np.percentile(SII_image, 5),
             vmax=np.percentile(SII_image, 99))
ax[1].set_title("[SII]")

plt.show()

#-------------- TO CHECK BACKGROUND SUBTRACTION RAN PROPERLY --------------

print("")
print("Ha background mean of cutout=", np.mean(Ha_bg)) #avg background estimate
print("SII background mean cutout=", np.mean(SII_bg)) #avg background estimate 
print("")

print("Ha original mean value of image =", np.mean(Ha)) #avg value of sources + background !
print("Ha subtracted mean value of image=", np.mean(Ha_image)) #avg value of sources + background after background is subtracted 
print("")

print("SII original mean value of image =", np.mean(SII)) #avg value of sources + background !
print("SII subtracted mean value of image =", np.mean(SII_image)) #avg value of sources + background after background is subtracted 

# --- details of the background cutout ------ 
print("")
print("Ha background std=", np.std(Ha_bg)) #std of background 
print("SII background std=",np.std(SII_bg)) #std of background 

from astropy.io import fits

outpath1 = r"C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_background_subtracted.fits"
fits.writeto(outpath1, Ha_image, ha_header ,overwrite=True)


outpath2 = r"C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_background_subtracted.fits"
fits.writeto(outpath2, SII_image, ha_header ,overwrite=True)




