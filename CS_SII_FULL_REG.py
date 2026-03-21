# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:55:34 2026

@author: lily zaidi
"""
#combining the regions that we are going to be using for a better scale factor 

# libraries + imports
from astropy.io import fits #read fits file 
import numpy as np #array 
import matplotlib.pyplot as plt #plotting 
from astrodendro import Dendrogram
from matplotlib.colors import LogNorm #normalize log
from astropy.io.fits import getdata #extract pixel value
from astropy.stats import sigma_clipped_stats #stats of data


# shows original SII image in full 
fileSII = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/M33SII_900s_WCS.fits')
img_SII = fileSII[0].data

plt.figure()
plt.imshow(img_SII, origin='lower', cmap='magma', vmin=np.percentile(img_SII,0.5), vmax=np.percentile(img_SII,99.8))
plt.title("SII")
plt.show()

# combining SII cutouts in one figure 
file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/sixth_proper_cutout_px_SII.fits')
file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fifth_proper_cutout_px_SII.fits')

file3 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/sixth_proper_cutout_px_R.fits')
file4 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fifth_proper_cutout_px_R.fits')

#sii hdu data
img1 = file1[0].data 
img1 
img2 = file2[0].data #can do hdul[0].header to open the header content 
img2 
#r hdu data 
img3 = file3[0].data
img3
img4 = file4[0].data
img4 

def combine_side_by_side(img1, img2, gap=20, pad_value=np.nan):
    h1, w1 = img1.shape
    h2, w2 = img2.shape

    new_h = max(h1, h2)
    new_w = w1 + gap + w2

    combined = np.full((new_h, new_w), pad_value) 

    combined[:h1, :w1] = img1
    combined[:h2, w1 + gap:w1 + gap + w2] = img2 #gap is between the image 
    #to prevent false detections, but we have to substitute np.nan for that gap
    #np.nan gap = not a number, gap

    return combined

img_s_combined = combine_side_by_side(img1, img2, gap=20)
plt.figure(figsize=(12,5))
plt.imshow(img_s_combined, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_s_combined, 0.5),
           vmax=np.nanpercentile(img_s_combined, 99.8))
plt.colorbar()
plt.title("Combined SII regional cutouts")
plt.show()

img_r_combined = combine_side_by_side(img3,img4, gap=20)
plt.figure(figsize=(12,5))
plt.imshow(img_s_combined, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_r_combined, 0.5),
           vmax=np.nanpercentile(img_r_combined, 99.8))
plt.colorbar()
plt.title("Combined R regional cutouts")
plt.show()




# get stats of the combined imgurrr :o 
mean_s,median_s,std_s = sigma_clipped_stats(img_s_combined, sigma = 3.0)
mean_r,median_r,std_r= sigma_clipped_stats(img_r_combined, sigma = 3.0)

print("mean_s: {},\nmedian_s: {}, \nstd_s: {}".format(mean_s,median_s,std_s))
print("")

# replacing the NaN in the gap with the median so that we don't have complications 
# with trying to subtract background or detect sources later
# DAOStarFinder will get error if it detects NaN
# DAO StarFinder will never detect the median as we are detecting stars levels above the 
# background median so we are good.

#empty pixel replacement with med
img_s_detect = np.nan_to_num(img_s_combined, nan=median_s)


from photutils.detection import DAOStarFinder

daofind = DAOStarFinder(fwhm=4, threshold=4*std_s)

sources_s = daofind(img_s_detect - median_s)

for col in sources_s.colnames:
    sources_s[col].info.format = '%.2f'

print(sources_s)
print("Total sources detected:", len(sources_s))

from photutils.aperture import CircularAperture, aperture_photometry

# make list of (x,y) from SII detections
positions = np.transpose((sources_s['xcentroid'], sources_s['ycentroid']))

# make apertures at those positions
apertures_s = CircularAperture(positions, r=3)

# do photometry using the apetures we take from sii
phot_s = aperture_photometry(img_s_detect, apertures_s) #using img s detect as detection space 

# we also have to do temporary detection for r, because the image combination on r, also needs a 
# detection of the same stars

img_r_detect = np.nan_to_num(img_r_combined, nan=median_r)

phot_r = aperture_photometry(img_r_detect, apertures_s)

plt.figure(figsize=(12,5))
plt.imshow(img_s_detect, origin='lower', cmap='magma')
plt.colorbar()
apertures_s.plot(color='red')
plt.show()

plt.figure(figsize=(12,5))
plt.imshow(img_r_detect, origin='lower', cmap='magma')
plt.colorbar()
apertures_s.plot(color='red')
plt.show()


# flux here is the direct measured estimate from the first 5 stars 
print('Flux SII')
print(phot_s['aperture_sum'][:5]) #photometry table 
print("")
print('Flux R')
print(phot_r['aperture_sum'][:5]) #photometry table 
# the aperture sum, finds the sum of the brightness within an individual aperture 
# for the first 5 stars 

print("")
print('Total sources detected for SII: {}'.format(len(phot_s))) # same sources = same output 
print('Total sources measured for R (from R list): {}'.format(len(phot_r))) # same sources = same output 


#sorted flux sum SII and R
print("")
sources_s.sort('flux') 
print(sources_s[:7]) #sorted the table by brightness 

print("")
from astropy.wcs import WCS
wcs_h = WCS(file1[0].header)
apertures_sky = apertures_s.to_sky(wcs_h)
print("Aperture_SII(RA/Dec)")
print("------------------")
print(apertures_sky)

print("")
wcs_r = WCS(file2[0].header)
print("Aperture_R_wcs(RA/Dec)")
print("------------------")
aperture_sky2 = apertures_s.to_sky(wcs_r)
print(aperture_sky2)
print("")

#------------------[REMOVING BACKGROUND LIGHT]---------------------------------
# we are removing background light within each aperture from the physical flux, 
# because we want only just the star brightness 

from photutils.aperture import CircularAnnulus
# We can use Circular Annulus to minus background light in each star 
#Circular Annulus has an empty hole in the middle so we minus the perimeter background light 
#within the aperture 
annulus = CircularAnnulus(positions, r_in=5, r_out=9)
#r in has to be bigger than the r aperture from CircularAperture



# photometry in the annulus
bkg_s = aperture_photometry(img_s_detect, annulus)
bkg_r = aperture_photometry(img_r_detect, annulus)

# we can get the areas of the annulus and the circular aperture
ap_area = apertures_s.area
ann_area = annulus.area 

# background per pixel
bkg_perpix_s = bkg_s['aperture_sum'] / ann_area
bkg_perpix_r = bkg_r['aperture_sum'] / ann_area

# background expected inside the aperture
bkg_in_ap_s = bkg_perpix_s * ap_area
bkg_in_ap_r = bkg_perpix_r * ap_area

# background-subtracted star fluxes
Fs = phot_s['aperture_sum'] - bkg_in_ap_s
Fr = phot_r['aperture_sum'] - bkg_in_ap_r

#=========[CORRECTED FLUX]================
print('Flux SII')
print(Fs[:5]) #photometry table 
print("")
print('Flux R')
print(Fr[:5]) #photometry table
print("")

# keep only good stars
Fs = np.array(Fs)
Fr = np.array(Fr)
good = np.isfinite(Fs) & np.isfinite(Fr) & (Fs > 0) & (Fr > 0) 
#^ this line removes NaN, inf, -inf  , keeps only positive fluxes of Fh 
#and Fr 

Fs_g = Fs[good] #array 
Fr_g = Fr[good] #array 


print("Number of good stars used:", len(Fs_g))

plt.figure()
plt.scatter(Fr_g,Fs_g , s= 10, color = 'black')
k, b = np.polyfit(Fr_g, Fs_g, 1)

x_line = np.linspace(Fr_g.min(), Fr_g.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color ='red', linewidth=2)

   # SII = k*R + b
print("k, b =", k, b)

plt.xlabel("R-band flux  ")
plt.ylabel("SII-band flux ")
plt.show()

# ----- load master images -----
file3 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/Failed Cutouts/cutout_800px_Rnew.fits')
img_R_band = file3[0].data


#misaligned need to fix either psf or alignment
file4 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/Failed Cutouts/cutout_800px_SIInew.fits')
img_SII_band = file4[0].data

# ----- safer subtraction (background removed) -----
img_s0 = img_SII_band - np.nanmedian(img_SII_band)
img_r0 = img_R_band - np.nanmedian(img_R_band)

new_SII_data = img_SII_band - (k* img_R_band)

vmin = np.nanpercentile(new_SII_data, 0.5) #contrast settings
vmax = np.nanpercentile(new_SII_data, 99.8) #contrast settings 
plt.figure()
plt.imshow(new_SII_data, origin="lower", cmap="magma", vmin=vmin, vmax=vmax)
plt.colorbar()
plt.show()







