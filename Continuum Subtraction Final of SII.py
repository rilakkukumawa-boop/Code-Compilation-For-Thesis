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
#------------------------------------------------------------------------------
#---------------------------------[ORIGINALS]----------------------------------


fileSII = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/M33SII_900s_WCS.fits')
img_SII = fileSII[0].data

plt.figure()
plt.imshow(img_SII, origin='lower', cmap='magma', vmin=np.percentile(img_SII,0.5), vmax=np.percentile(img_SII,99.8))
plt.title("SII")
plt.show()




# ------------------------ {FILES H & R}---------------------------------------

#these files are 2D-Grids, 

#open fits creating hdu lists python data units 
file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/sixth_proper_cutout_px_SII.fits')
file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/sixth_proper_cutout_px_R.fits')
#file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fifth_proper_cutout_px_SII.fits')
#file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fifth_proper_cutout_px_R.fits')

img_s = file1[0].data #to read header data/pixelvalues, generally you do hdul = fits.open()
img_s                # then you do hdul[0].data to read
img_r = file2[0].data #can do hdul[0].header to open the header content 
img_r 


#--------------------[DISPLAY Ha & R] -----------------------------------------
plt.figure() #blank canvas 
plt.imshow(img_s, origin='lower', cmap='magma') #display the ha
plt.colorbar() #add colour bar for reference 
plt.show() # prompt plot to show 

plt.figure()
plt.imshow(img_r, origin='lower', cmap='plasma') #display the r
plt.colorbar() #add colour bar for reference 
plt.show()# prompt plot to show 
#note: orgin='lower' means co-ordinate [0,0] is @ bottom left 

#---------------------------- [STATS]  ----------------------------------------
mean_s,median_s,std_s = sigma_clipped_stats(img_s, sigma = 3.0) #stats
#we remove pixels that are 3 standard deviations away from the estimated background 
#mean pixel average, i.e gets clipped, outliers are removed 


print("mean_s: {},\nmedian_s: {}, \nstd_s: {}".format(mean_s,median_s,std_s))
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
daofind =  DAOStarFinder(fwhm = 4, threshold =4*std_s)
#builds a gaussian kernel and looks for stars that are 3 px wide. 
#finding stars that are about as big as 3 pixels 
#find stars that are above 5 standard deviations away from the background level brightness

# the DAOStarFinder has the DAOFind algorithm incorporated.
# it finds the local density maxima and peak ampiltude greater than the 
# threshold suggested in the image (note: we applied a condition to the threshold earlier)
# creates a 2D Gaussian Kernel, defined by fwhm, ratio, theta etc. 

# we can subtract the background of Ha, and we should have peaks that are n-sigma above the 
# average estimated background brightness
sources_s = daofind(img_s - median_s)
#creates a temporary image image_h <- background minus the median pixel value
#daofind you find sources on the new temporary image 


# we can yield a table with the information on the sources 
for col in sources_s.colnames:
    sources_s[col].info.format = '%.2f'

print(sources_s)    
print("")
# flux in here is different because its just a estimate 
# of the star as a gaussian blob 


# if we follow the same procedure in source detection up till now as a copy 
# for band-r, we will plausibly find that the sources in r and ha
# are detected differently, so to mitigate this problem we have to reuse the 
# same apeture and use the same list that we do in Ha 

from photutils.aperture import CircularAperture, aperture_photometry

# make list of (x,y) from Hα detections
positions = np.transpose((sources_s['xcentroid'], sources_s['ycentroid']))
#note: these are pixel co-ordinates, can convert these to sky co-ordinates and document
# them if prefered 

# make apertures at those positions
apertures_s = CircularAperture(positions, r=3)
#Cicular Aperture provides a measure of how bright everything is within the circle

# do photometry using the apetures we take from ha
phot_s = aperture_photometry(img_s, apertures_s) # this is raw flux, we dont use this 
#note: aperture_photometry function assumes the input data and background data have 
#been subtracted 

#if you background subtract, you assume the background is uniform which is not 
#true for ha emission so thats why ive left image as themselves 
phot_r = aperture_photometry(img_r, apertures_s) # this is raw flux, we dont use this

# flux here is the direct measured estimate from the first 5 stars 
print('Flux Ha')
print(phot_s['aperture_sum'][:5]) #photometry table 
print("")
print('Flux R')
print(phot_r['aperture_sum'][:5]) #photometry table 
# the aperture sum, finds the sum of the brightness within an individual aperture 
# for the first 5 stars 

print("")
print('Total sources detected for SII: {}'.format(len(phot_s))) # same sources = same output 
print('Total sources measured for R (from Ha list): {}'.format(len(phot_r))) # same sources = same output 


#-------------[CORRECTED DISPLAY W/ APERTURES]--------------------------------
#note: only for visual representation, no subtraction is done here to the data 
plt.figure() #blank canvas for plot
plt.imshow(img_s, origin='lower', cmap='magma')
apertures_s.plot(color='red', lw=1.5)

plt.figure() #blank canvas for plot
plt.imshow(img_r, origin='lower', cmap='magma')
apertures_s.plot(color='red', lw=1.5 )



#-----------[FIND THE SKY CO-ORDINATES OF THE APERTURES]-----------------------
# 7 sources are detected in Ha, this is the same for R 
# here I will sort the sources out first 
print("")
sources_s.sort('flux') 
print(sources_s[:7]) #sorted the table by brightness 

print("")
from astropy.wcs import WCS
wcs_h = WCS(file1[0].header)
apertures_sky = apertures_s.to_sky(wcs_h)
print("Aperture_ha(RA/Dec)")
print("------------------")
print(apertures_sky)

print("")
wcs_r = WCS(file2[0].header)
print("Aperture_r_wcs(RA/Dec)")
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
bkg_s = aperture_photometry(img_s, annulus)
bkg_r = aperture_photometry(img_r, annulus)

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

#---------------------[PLOTTING FOR SCALE FACTOR]----------------------------
# now we want to plot a flux of R vs Ha 
# so flux_Ha is the x-values 
# and flux_R is the y-values


#-----  R vs Ha ---------
# Keep only points where: (R < 1500) OR (SII > 400)
# This specifically targets that bottom-right "cluster" of outliers
#Fr_g = Fr_g[(Fr_g < 1500) | (Fs_g > 400)]
#Fs_g = Fs_g[(Fr_g < 1500) | (Fs_g > 400)]
#plt.figure()
#plt.scatter(Fs_g, Fr_g, s= 10, color = 'black')
#to find line of best fit for y = mx+c 
#m,c = np.polyfit(Fs_g, Fr_g,1)
#y = m*x +c 
#Fh_g = m*Fr_g + c
#xline = np.linspace(Fs_g.min(), Fs_g.max(), 200)
#yline =(m * xline) + c
#yline = best fit
#best_fit = m*x + c
#equation of line, i.e best fit
#plt.plot(xline, yline, color ='red', linewidth=2)

#plt.xlabel("SII-band flux  ")
#plt.ylabel("R-band flux ")#
#plt.title("Calibration: R vs SII (stars)")
#plt.show()
#print("m, c =", m,c)

# ----- SII vs R ? -----
# Keep points that are NOT in the 100-150 strip and not in the bottom-right corner
"""keep_condition = ~( ((Fr_g >= 500) & (Fr_g <= 2500) & (Fs_g >= 100) & (Fs_g <= 150)) |  # original strip
    ((Fr_g > 1200) & (Fs_g < 450))                                     # The 2 specific outliers
)

Fr_g = Fr_g[keep_condition]
Fs_g = Fs_g[keep_condition]

# then you can do plt fig, but keep in mind this runs for 85 px"""
#keep_condition = ~(((Fr_g <= 100) & (Fs_g <= 95)) )

#Fr_g = Fr_g[keep_condition]
#Fs_g = Fs_g[keep_condition]                 

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
file3 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_800px_Rnew.fits')
img_R_band = file3[0].data


#misaligned need to fix either psf or alignment
file4 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_800px_SIInew.fits')
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





#instructions 
#name apeture and save them 
#convert apertures to ra and dec 

#region in ha with stars and no diff emissions 
#find sources of stars in Ha 
#make sure the r is just right for the star 

# apetures that are found are summed
# list of apetures that i have in order 
# list of fluxes for each of the stars in R and Ha 
# for loop for R and Ha photometry extraction/ apeture photometry extraction 
# making a flux graph of R vs Ha
# y = mx + c <- after making flux graph of R vs Ha 
# y = mR + c <- R is the Rbandimage (can name this cal image )
# then you can do Ha = Ha -calimage














