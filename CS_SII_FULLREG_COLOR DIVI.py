# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:06:36 2026

@author: lily zaidi
"""

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
plt.imshow(img_SII, origin='lower', cmap='plasma', vmin=np.percentile(img_SII,0.5), vmax=np.percentile(img_SII,99.8))
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
plt.imshow(img_r_combined, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_r_combined, 0.5),
           vmax=np.nanpercentile(img_r_combined, 99.8))
plt.colorbar()
plt.title("Combined R regional cutouts")
plt.show()




# get stats of the combined imgurrr :o 
mean_s,median_s,std_s = sigma_clipped_stats(img_s_combined, sigma = 3.0)
mean_r,median_r,std_r= sigma_clipped_stats(img_r_combined, sigma = 3.0)
print("")
print("mean_s: {},\nmedian_s: {}, \nstd_s: {}".format(mean_s,median_s,std_s))
print("")

# replacing the NaN in the gap with the median so that we don't have complications 
# with trying to subtract background or detect sources later
# DAOStarFinder will get error if it detects NaN
# DAO StarFinder will never detect the median as we are detecting stars levels above the 
# background median so we are good.

#empty pixel replacement with med from combined SII image
img_s_detect = np.nan_to_num(img_s_combined, nan=median_s)


from photutils.detection import DAOStarFinder

daofind = DAOStarFinder(fwhm=4, threshold=4*std_s)

sources_s = daofind(img_s_detect - median_s)

for col in sources_s.colnames:
    sources_s[col].info.format = '%.2f' #column printed of detection 

print(sources_s)
print("Total sources of stars detected:", len(sources_s))

from photutils.aperture import CircularAperture, aperture_photometry

# make list of (x,y) from SII detections
positions = np.transpose((sources_s['xcentroid'], sources_s['ycentroid']))

# make apertures at those positions
apertures_s = CircularAperture(positions, r=3)

# do photometry using the apetures we take from sii
phot_s = aperture_photometry(img_s_detect, apertures_s) #using img s detect as detection space, 
#img_s_detect is just the temporary image (remember that this is the temporary image
#with a padded gap replaced from NaN => np.median)

# we also have to do temporary detection for r, because the image combination on r, also needs a 
# detection of the same stars

img_r_detect = np.nan_to_num(img_r_combined, nan=median_r) #replacing nan gap with 
#median_r

phot_r = aperture_photometry(img_r_detect, apertures_s) #detecting on temporary 
# image with same detections on s, these should be the same sources 
# this is raw flux, we dont use this 
#note: aperture_photometry function assumes the input data and background data have 
#been subtracted 

plt.figure(figsize=(12,5)) #scale of image
plt.imshow(img_s_detect, origin='lower', cmap='magma',vmin=np.nanpercentile(img_s_combined, 0.5),
vmax=np.nanpercentile(img_s_combined, 99.8)) #regions combined with 
#median in between
plt.colorbar() #color bar 
apertures_s.plot(color='red')
plt.title("SII apertures") # plotted circles aperture plot
plt.show()

plt.figure(figsize=(12,5)) #scale of image
plt.imshow(img_r_detect, origin='lower', cmap='magma',vmin=np.nanpercentile(img_r_combined, 0.5),
vmax=np.nanpercentile(img_r_combined, 99.8))
plt.colorbar()
apertures_s.plot(color='red')
plt.title("R apertures")
plt.show()


# flux here is the direct measured estimate from the first 5 stars 
print('Flux SII')
print(phot_s['aperture_sum'][:5]) #photometry table 
print("")
print('Flux R')
print(phot_r['aperture_sum'][:5]) #photometry table 
# the aperture sum, finds the sum of the brightness within an individual aperture 
# for the first 5 stars  
# this is raw flux, we dont use this 
#note: aperture_photometry function assumes the input data and background data have 
#been subtracted 

print("")
print('Total sources detected for SII: {}'.format(len(phot_s))) # same sources = same output 
print('Total sources measured for R (from R list): {}'.format(len(phot_r))) # same sources = same output 


#sorted flux sum SII and R
print("")
sources_s.sort('flux') 
print(sources_s[:7]) #sorted the table by brightness 

from astropy.wcs import WCS
gap = 20
w1 = img1.shape[1]   # width of left Ha cutout

x = positions[:, 0]
y = positions[:, 1]

#confirming if coordinates are the same for the matched cutouts
# ---------------- LEFT CUTOUT ----------------
left_mask = x < w1
pos_left = positions[left_mask]

if len(pos_left) > 0:
    apertures_left = CircularAperture(pos_left, r=3)
    wcs_left = WCS(file1[0].header)   # Ha left cutout WCS
    apertures_left_sky = apertures_left.to_sky(wcs_left)

    print("Aperture_SII_left (RA/Dec)")
    print("------------------------")
    print(apertures_left_sky)
    print("")

# ---------------- RIGHT CUTOUT Ha ----------------
right_mask = x >= (w1 + gap)
pos_right = positions[right_mask].copy()

# shift x positions back into the original right-cutout pixel system
pos_right[:, 0] = pos_right[:, 0] - (w1 + gap)

if len(pos_right) > 0:
    apertures_right = CircularAperture(pos_right, r=3)
    wcs_right = WCS(file2[0].header)   # Ha right cutout WCS
    apertures_right_sky = apertures_right.to_sky(wcs_right)

    print("Aperture_SII_right (RA/Dec)")
    print("-------------------------")
    print(apertures_right_sky)
    print("")
    
# ---------------- LEFT REGION (R WCS) ----------------
left_mask = x < w1
pos_left = positions[left_mask]

if len(pos_left) > 0:
    apertures_left_r = CircularAperture(pos_left, r=3)
    wcs_left_r = WCS(file3[0].header)   # R left cutout WCS
    apertures_left_r_sky = apertures_left_r.to_sky(wcs_left_r)

    print("Aperture_R_left (RA/Dec)")
    print("-----------------------")
    print(apertures_left_r_sky)
    print("")

# ---------------- RIGHT REGION (R WCS) ----------------
right_mask = x >= (w1 + gap)
pos_right = positions[right_mask].copy()

# shift x positions back into the original right-cutout pixel system
pos_right[:, 0] = pos_right[:, 0] - (w1 + gap)

if len(pos_right) > 0:
    apertures_right_r = CircularAperture(pos_right, r=3)
    wcs_right_r = WCS(file4[0].header)   # R right cutout WCS
    apertures_right_r_sky = apertures_right_r.to_sky(wcs_right_r)

    print("Aperture_R_right (RA/Dec)")
    print("------------------------")
    print(apertures_right_r_sky)
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
print(Fs[:5]) #photometry table = corrected flux
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

# also keep x positions for the same good stars
x_good = x[good]

# first scatter plot with left/right colours
left_good = x_good < w1
right_good = x_good >= (w1 + gap)

 #plt.scatter(Fr_g,Fs_g , s= 10, color = 'black') <- instead we can colour code these regions 
plt.figure()
plt.scatter(Fr_g[left_good], Fs_g[left_good], s=25, color='magenta', label='Left cutout')
plt.scatter(Fr_g[right_good], Fs_g[right_good], s=25, color='orange', label='Right cutout')

k, b = np.polyfit(Fr_g, Fs_g, 1)
x_line = np.linspace(Fr_g.min(), Fr_g.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color ='red', linewidth=1)

   # SII = k*R + b
print("k, b =", k, b)

plt.xlabel("R-band (pixel count) ")
plt.ylabel("SII-band (pixel count) ")
plt.title('Pixel Intensities of SII vs R')
plt.legend()
plt.grid()
plt.show()

#---------------- manual removal of outliers -----------------------------
remove_outliers = ~(
    ((Fr_g < 3000) & (Fs_g > 900)) |     # high SII for very small R
    ((Fr_g < 4000) & (Fs_g < 250))       # very low SII near origin
)

Fr_final = Fr_g[remove_outliers]
Fs_final = Fs_g[remove_outliers]
x_final = x_good[remove_outliers]


print("Number after outlier removal:", len(Fs_final))

# define left/right again after outlier removal
left_final = x_final < w1
right_final = x_final >= (w1 + gap)

plt.figure()
plt.scatter(Fr_final[left_final], Fs_final[left_final],
            s=25, color='magenta', label='Left cutout')

plt.scatter(Fr_final[right_final], Fs_final[right_final],
            s=25, color='orange', label='Right cutout')


k, b = np.polyfit(Fr_final, Fs_final, 1)

x_line = np.linspace(Fr_final.min(), Fr_final.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color='red', linewidth=1)

print("k, b =", k, b)
plt.xlabel("R-band (pixel count)")
plt.ylabel("SII-band (pixel count)")
plt.title('Pixel Intensities of SII vs R')
plt.legend()
plt.grid()
plt.show()


# ----- load master images -----
file5 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_800px_Rnew.fits')
img_R_band = file5[0].data


#misaligned need to fix either psf or alignment
file6 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_800px_SIInew.fits')
img_SII_band = file6[0].data

# ----- safer subtraction (background removed) -----
#img_s0 = img_SII_band - np.nanmedian(img_SII_band)
#img_r0 = img_R_band - np.nanmedian(img_R_band)

new_SII_data = img_SII_band - (k* img_R_band)

vmin = np.nanpercentile(new_SII_data, 0.5) #contrast settings
vmax = np.nanpercentile(new_SII_data, 99.8) #contrast settings 
plt.figure()
plt.imshow(new_SII_data, origin="lower", cmap="plasma", vmin=vmin, vmax=vmax)
plt.colorbar()  
plt.show()
#copyment of header
#header = file6[0].header.copy()
#header transfer and new fits made: 
#fits.writeto("SII_subtracted.fits", new_SII_data, header=header, overwrite=True)
