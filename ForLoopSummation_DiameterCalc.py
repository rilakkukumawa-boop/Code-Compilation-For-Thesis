# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 02:34:11 2026

@author: lily zaidi
"""
#code for summing pixel value over pixel value for each aperture 
#plus print out of its diameter
#printing column of co-ordinates, with ratio and diameter into a text file


import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions


region_file = r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/KAI_SA_BLOBS_FITTED.reg' #reg file
ha_fits = r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_background_subtracted.fits' #ha 
#sii_fits = r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/SII_background_subtracted.fits'#sii
sii_fits = r'C:/Users/wafiy/Downloads/LilliesXD.fits'
wcs_fits = r'C:/Users/wafiy/T.W/Pythonscripts_Tutorials/Continuum Subtractions/Ratio_Maps_26th/Ha_CS_WCS_cutout790px.fits'#wcs ha


regions = Regions.read(region_file, format='ds9')
print("")
print("Number of apertures in region file:", len(regions))
print("")

#images as floats
Ha = fits.getdata(ha_fits).astype(float)
SII = fits.getdata(sii_fits).astype(float)


# wcs of ha
header = fits.getheader(wcs_fits)
wcs = WCS(header)

#pixel scale
pc11 = header['PC1_1']
pc12 = header['PC1_2']
pc21 = header['PC2_1']
pc22 = header['PC2_2']

scale_x = np.sqrt(pc11**2 + pc12**2)
scale_y = np.sqrt(pc21**2 + pc22**2)


pix_scale_deg = (scale_x + scale_y) / 2.0
pix_scale_arcsec = pix_scale_deg * 3600.0
print("Pixel scale X_direction =", scale_x, "degree per pixel" )
print("Pixel scale Y_direction =", scale_y, "degree per pixel" )

print("")
print("Pixel Scale = ", pix_scale_deg, "degree per pixel")
print("converted                                  ↓")
print("Pixel scale =", pix_scale_arcsec, "arcsec per pixel")

print("")

for i, reg in enumerate(regions):
    pixel_reg = reg.to_pixel(wcs) #switches out sky coordinates to pixel coord
    
    # get diameter of circular aps (based on the region type)
    if hasattr(pixel_reg, 'radius'):
        diam_pix = pixel_reg.radius * 2 #takes the radius in pixels and multiplies it by 2 
    else:
        diam_pix = 0
        
    diam_arcsec = diam_pix * pix_scale_arcsec #multiplies the pixel diameter by 
    #pixel scale to give arcsec
    distance_pc = 820000  # <-- whatever your object's distance is
    diam_pc = diam_arcsec * distance_pc / 206265.0

    # mask is made to sum values 
    # changed 'method' to 'mode' here:
    mask = pixel_reg.to_mask(mode='center') 
    
    sii_values = mask.get_values(SII)
    ha_values = mask.get_values(Ha)
    
    # using np.nanmedian for each ratios because nansum a bit too high 
    sum_sii = np.nansum(sii_values)
    sum_ha = np.nansum(ha_values)
    
    # 4. Calculate Ratio
    if sum_ha != 0:
        ratio = sum_sii / sum_ha
    else:
        ratio = np.nan
        
    # getting the RA and Dec printed out again 
    # get RA/Dec from region (sky coords)
    
    if hasattr(reg, 'center'):
        ra_str = reg.center.ra.to_string(unit='hour', sep=':')
        dec_str = reg.center.dec.to_string(unit='deg', sep=':')
    else:
        ra, dec = np.nan, np.nan
        #getting the RA and Dec printed out again 
        # get RA/Dec from region (sky coords)
   
    if hasattr(reg, 'center'):
        ra_str = reg.center.ra.to_string(unit='hour', sep=':') #ra and dec hr mm second style
        dec_str = reg.center.dec.to_string(unit='deg', sep=':')
        # decimal degrees
        ra_deg = reg.center.ra.deg
        dec_deg = reg.center.dec.deg
    else:
        ra_str = 'nan'
        dec_str = 'nan'
        ra_deg = np.nan
        dec_deg = np.nan
        # getting the RA and Dec printed out again 
        # get RA/Dec from region (sky coords)
        #| RA: {ra_str} | Dec: {dec_str}

#convert distance to parsecs through 
# a = 206265 * linear diameter / distance 
#since arcsecs = linear size/ distance 
# converting arc secs to radians 360 * 3600/ 2 pi = 206265 

# linear diameter or size = (distance(pc) * arcsecs)/ 206265
# | RA(deg): {ra_deg:.6f} | Dec(deg): {dec_deg:.6f} | RA: {ra_str} | Dec: {dec_str} |
#    print(f'ID {i}|Ratio: {ratio:.3f} | Diam: {diam_arcsec:.2f}"| Diam: {diam_pc:.4f} pc')
        

    print(f'ID {i}| RA(deg): {ra_deg:.6f} | Dec(deg): {dec_deg:.6f} | RA: {ra_str} | Dec: {dec_str} |Ratio: {ratio:.3f} | Diam: {diam_pc:.4f} | Diam: {diam_pc:.4f} pc')
        
