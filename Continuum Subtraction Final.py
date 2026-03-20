# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 18:28:10 2026

@author: lily zaidi
"""
# libraries + imports
from astropy.io import fits #read fits file 
import numpy as np #array 
import matplotlib.pyplot as plt #plotting 
from astrodendro import Dendrogram
from matplotlib.colors import LogNorm #normalize log
from astropy.io.fits import getdata #extract pixel value
from astropy.stats import sigma_clipped_stats #stats of data

    
#------------------------------main objective----------------------------------
#main goal is to do a continuum subtraction 
#need to cut out a map from the fits image that matches both regions in Ha and R
#use detection on the the cut out image
#find the amount the image with the biggest amount of stars with no dust
#try subtract it for using an apeture
#find the flux line for the Rband and Ha 
#---------------------------------[ORIGINALS]----------------------------------

fileh = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/M33Ha_300s_WCS.fits')
img_ha = fileh[0].data

plt.figure()
plt.imshow(img_ha, origin='lower', cmap='magma', vmin=np.percentile(img_ha,0.5), vmax=np.percentile(img_ha,99.8))
plt.title("Ha")
plt.show()


# ------------------------ {FILES H & R}--------------------------------------

#file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/third_proper_cutout_px_Ha.fits')
#file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/third_proper_cutout_px_R.fits')

file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fourth_proper_cutout_px_Ha.fits')
file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fourth_proper_cutout_px_R.fits')


#these files are 2D-Grids
#file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_85px_Ha.fits')
#file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_85px_r.fits') #0.19
#open fits creating hdu lists python data units 

img_h = file1[0].data #to read header data/pixelvalues, generally you do hdul = fits.open()
img_h                 # then you do hdul[0].data to read
img_r = file2[0].data #can do hdul[0].header to open the header content 
img_r 


#--------------------[DISPLAY Ha & R] ----------------------------------------
plt.figure() #blank canvas 
plt.imshow(img_h, origin='lower', cmap='magma') #display the ha
plt.colorbar() #add colour bar for reference 
plt.show() # prompt plot to show

plt.figure()
plt.imshow(img_r, origin='lower', cmap='magma') #display the r
plt.colorbar() #add colour bar for reference 
plt.show()# prompt plot to show 
#note: orgin='lower' means co-ordinate [0,0] is @ bottom left 

#---------------------------- [STATS]  ---------------------------------------
mean_h,median_h,std_h = sigma_clipped_stats(img_h, sigma = 3.0) #stats
#we remove pixels that are 3 standard deviations away from the estimated background 
#mean pixel average, i.e gets clipped, outliers are removed 


print("mean_h: {},\nmedian_h: {}, \nstd_h: {}".format(mean_h,median_h,std_h))
print("")

mean_r,median_r,std_r = sigma_clipped_stats(img_r, sigma = 3.0)
print("mean_r: {},\nmedian_r: {}, \nstd_r: {}".format(mean_r,median_r,std_r))
#sigma_clipped_stats calculates the average background brightness of image 
#mean brightness , median brightness , and the std of the data 
print("")
#------------------------- [FINDING SOURCES]-----------------------------------
# to do photometry must identify a collection of sources 
# we use photoutils here and daostarfinder <- this automatically detects bright sources 
# in our image. 
from photutils.detection import DAOStarFinder
daofind =  DAOStarFinder(fwhm = 4.0, threshold = 3.0*std_h)
#builds a gaussian kernel and looks for stars that are 3 px wide. 
#finding stars that are about as big as 3 pixels 
#find stars that are above 4 standard deviations away from the background level brightness

# the DAOStarFinder has the DAOFind algorithm incorporated.
# it finds the local density maxima and peak ampiltude greater than the 
# threshold suggested in the image (note: we applied a condition to the threshold earlier)
# creates a 2D Gaussian Kernel, defined by fwhm, ratio, theta etc. 

# we can subtract the background of Ha, and we should have peaks that are n-sigma above the 
# average estimated background brightness
sources_h = daofind(img_h - median_h)
#creates a temporary image image_h <- background minus the median pixel value
#daofind you find sources on the new temporary image 


# we can yield a table with the information on the sources 
for col in sources_h.colnames:
    sources_h[col].info.format = '%.2f'

print(sources_h)    
print("")
# flux in here is different because its just a estimate 
# of the star as a gaussian blob 


# if we follow the same procedure in source detection up till now as a copy 
# for band-r, we will plausibly find that the sources in r and ha
# are detected differently, so to mitigate this problem we have to reuse the 
# same apeture and use the same list that we do in Ha 

from photutils.aperture import CircularAperture, aperture_photometry

# make list of (x,y) from Hα detections
positions = np.transpose((sources_h['xcentroid'], sources_h['ycentroid']))
#note: these are pixel co-ordinates, can convert these to sky co-ordinates and document
# them if prefered 

# make apertures at those positions
apertures_h = CircularAperture(positions, r=3)
#Cicular Aperture provides a measure of how bright everything is within the circle

# do photometry using the apetures we take from ha
phot_h = aperture_photometry(img_h, apertures_h) # this is raw flux, we dont use this 
#note: aperture_photometry function assumes the input data and background data have 
#been subtracted 

#if you background subtract, you assume the background is uniform which is not 
#true for ha emission so thats why ive left image as themselves 
phot_r = aperture_photometry(img_r, apertures_h) # this is raw flux, we dont use this

# flux here is the direct measured estimate from the first 5 stars 
print('Flux Ha')
print(phot_h['aperture_sum'][:5]) #photometry table 
print("")
print('Flux R')
print(phot_r['aperture_sum'][:5]) #photometry table 
# the aperture sum, finds the sum of the brightness within an individual aperture 
# for the first 5 stars 

print("")
print('Total sources detected for Ha: {}'.format(len(phot_h))) # same sources = same output 
print('Total sources measured for R (from Ha list): {}'.format(len(phot_r))) # same sources = same output 


#-------------[CORRECTED DISPLAY W/ APERTURES]--------------------------------
#note: only for visual representation, no subtraction is done here to the data 
plt.figure() #blank canvas for plot
plt.imshow(img_h, origin='lower', cmap='magma')
apertures_h.plot(color='red', lw=1.5)

plt.figure() #blank canvas for plot
plt.imshow(img_r, origin='lower', cmap='magma')
apertures_h.plot(color='red', lw=1.5 )



#-----------[FIND THE SKY CO-ORDINATES OF THE APERTURES]-----------------------
# 7 sources are detected in Ha, this is the same for R 
# here I will sort the sources out first 
print("")
sources_h.sort('flux') 
print(sources_h[:7]) #sorted the table by brightness 

print("")
from astropy.wcs import WCS
wcs_h = WCS(file1[0].header)
apertures_sky = apertures_h.to_sky(wcs_h)
print("Aperture_ha(RA/Dec)")
print("------------------")
print(apertures_sky)

print("")
wcs_r = WCS(file2[0].header)
print("Aperture_r_wcs(RA/Dec)")
print("------------------")
aperture_sky2 = apertures_h.to_sky(wcs_r)
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
bkg_h = aperture_photometry(img_h, annulus)
bkg_r = aperture_photometry(img_r, annulus)

# we can get the areas of the annulus and the circular aperture
ap_area = apertures_h.area
ann_area = annulus.area 

# background per pixel
bkg_perpix_h = bkg_h['aperture_sum'] / ann_area
bkg_perpix_r = bkg_r['aperture_sum'] / ann_area

# background expected inside the aperture
bkg_in_ap_h = bkg_perpix_h * ap_area
bkg_in_ap_r = bkg_perpix_r * ap_area

# background-subtracted star fluxes
Fh = phot_h['aperture_sum'] - bkg_in_ap_h
Fr = phot_r['aperture_sum'] - bkg_in_ap_r

#=========[CORRECTED FLUX]================
print('Flux Ha')
print(Fh[:5]) #photometry table 
print("")
print('Flux R')
print(Fr[:5]) #photometry table
print("")

# keep only good stars
Fh = np.array(Fh)
Fr = np.array(Fr)
good = np.isfinite(Fh) & np.isfinite(Fr) & (Fh > 0) & (Fr > 0) 
#^ this line removes NaN, inf, -inf  , keeps only positive fluxes of Fh 
#and Fr 

Fh_g = Fh[good] #array 
Fr_g = Fr[good] #array 


print("Number of good stars used:", len(Fh_g))

#Apply to your arrays (assuming Fs_g is your Ha flux here)
keep_condition = ~((Fr_g < 2500) & (Fh_g > 800)| (Fr_g < 25000) & (Fh_g >1800)|(Fr_g<1000) & (Fh_g>100))
Fr_g = Fr_g[keep_condition]
Fh_g = Fh_g[keep_condition]


plt.figure()
plt.scatter(Fr_g,Fh_g , s= 10, color = 'black')
k, b = np.polyfit(Fr_g, Fh_g, 1)

x_line = np.linspace(Fr_g.min(), Fr_g.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color ='red', linewidth=2)

   # Hα = k*R + b
print("k, b =", k, b)
plt.xlabel("R-band flux  ")
plt.ylabel("Hα-band flux ")
plt.show()

# ----- load master images -----
file3 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/Failed Cutouts/cutout_900px_Rnew.fits')
img_R_band = file3[0].data



file4 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/Failed Cutouts/cutout_900px_Hanew.fits')
img_H_band = file4[0].data

# ----- safer subtraction (background removed) -----
img_h0 = img_H_band - np.nanmedian(img_H_band)
img_r0 = img_R_band - np.nanmedian(img_R_band)

new_Ha_data = img_h0 - k* img_r0

vmin = np.nanpercentile(new_Ha_data, 0.5) #contrast settings
vmax = np.nanpercentile(new_Ha_data, 99.8) #contrast settings 
plt.figure()
plt.imshow(new_Ha_data, origin="lower", cmap="magma", vmin=vmin, vmax=vmax)
plt.colorbar()
plt.show()














