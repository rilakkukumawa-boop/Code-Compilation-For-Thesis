# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 19:15:25 2026

@author: lily zaidi
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:02:43 2026

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

#shows full Ha science image:

fileh = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/M33Ha_300s_WCS.fits')
img_ha = fileh[0].data

plt.figure()
plt.imshow(img_ha, origin='lower', cmap='magma',
           vmin=np.percentile(img_ha,0.5),
           vmax=np.percentile(img_ha,99.8))
plt.title("Ha")
plt.show()

# combining Ha cutouts in one figure 

file1 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/third_proper_cutout_px_Ha.fits')
file2 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fourth_proper_cutout_px_Ha.fits')

file3 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/third_proper_cutout_px_R.fits')
file4 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/fourth_proper_cutout_px_R.fits')

#ha hdu
img1 = file1[0].data 
img2 = file2[0].data 

#r hdu 
img3 = file3[0].data
img4 = file4[0].data

def combine_side_by_side(img1, img2, gap=20, pad_value=np.nan):
    h1, w1 = img1.shape
    h2, w2 = img2.shape

    new_h = max(h1, h2)
    new_w = w1 + gap + w2

    combined = np.full((new_h, new_w), pad_value) 

    combined[:h1, :w1] = img1
    combined[:h2, w1 + gap:w1 + gap + w2] = img2 #gap is between the image 
    #to prevent false detections, but we have to substitute np.nan for that gap
    
    return combined

img_ha_combined = combine_side_by_side(img1, img2, gap=20)
plt.figure(figsize=(12,5))
plt.imshow(img_ha_combined, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_ha_combined, 0.5),
           vmax=np.nanpercentile(img_ha_combined, 99.8))
plt.colorbar()
plt.title("Combined HA regional cutouts")
plt.show()

img_r_combined = combine_side_by_side(img3, img4, gap=20)
plt.figure(figsize=(12,5))
plt.imshow(img_r_combined, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_r_combined, 0.5),
           vmax=np.nanpercentile(img_r_combined, 99.8))
plt.colorbar()
plt.title("Combined R regional cutouts")
plt.show()

# get stats of the combined image
mean_h, median_h, std_h = sigma_clipped_stats(img_ha_combined, sigma=3.0)
mean_r, median_r, std_r = sigma_clipped_stats(img_r_combined, sigma=3.0)

print("")
print("mean_h: {},\nmedian_h: {}, \nstd_h: {}".format(mean_h, median_h, std_h))
print("")

#empty pixel replacement with median
img_h_detect = np.nan_to_num(img_ha_combined, nan=median_h)

#source detection
from photutils.detection import DAOStarFinder

daofind = DAOStarFinder(fwhm=4, threshold=4*std_h)
sources_h = daofind(img_h_detect - median_h)

for col in sources_h.colnames:
    sources_h[col].info.format = '%.2f'

print(sources_h)
print("Total sources detected:", len(sources_h))

from photutils.aperture import CircularAperture, aperture_photometry

# make list of (x,y) from ha detections
positions = np.transpose((sources_h['xcentroid'], sources_h['ycentroid']))

# make apertures at those positions
apertures_h = CircularAperture(positions, r=3)

# do photometry using the apertures we take from ha
phot_h = aperture_photometry(img_h_detect, apertures_h)

# temporary R detection image with nans converted to median
img_r_detect = np.nan_to_num(img_r_combined, nan=median_r)
phot_r = aperture_photometry(img_r_detect, apertures_h)

plt.figure(figsize=(12,5))
plt.imshow(img_h_detect, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_ha_combined, 0.5),
           vmax=np.nanpercentile(img_ha_combined, 99.8))
plt.colorbar()
apertures_h.plot(color='red')
plt.title("Ha apertures")
plt.show()

plt.figure(figsize=(12,5))
plt.imshow(img_r_detect, origin='lower', cmap='magma',
           vmin=np.nanpercentile(img_r_combined, 0.5),
           vmax=np.nanpercentile(img_r_combined, 99.8))
plt.colorbar()
apertures_h.plot(color='red')
plt.title("R apertures")
plt.show()

# flux here is the direct measured estimate from the first 5 stars 
print('Flux Ha')
print(phot_h['aperture_sum'][:5])
print("")
print('Flux R')
print(phot_r['aperture_sum'][:5])

print("")
print('Total sources detected for Ha: {}'.format(len(phot_h)))
print('Total sources measured for R (from R list): {}'.format(len(phot_r)))

#sorted flux sum ha
print("")
sources_h.sort('flux') 
print(sources_h[:7])

print("")
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

    print("Aperture_Ha_left (RA/Dec)")
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

    print("Aperture_Ha_right (RA/Dec)")
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
# we are removing background light within each aperture from the physical flux

from photutils.aperture import CircularAnnulus

# Circular annulus to estimate local background
annulus = CircularAnnulus(positions, r_in=5, r_out=9)

# photometry in the annulus
bkg_h = aperture_photometry(img_h_detect, annulus)
bkg_r = aperture_photometry(img_r_detect, annulus)

# areas of annulus and circular aperture
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
print(Fh[:5])
print("")
print('Flux R')
print(Fr[:5])
print("")

# keep only good stars
Fh = np.array(Fh)
Fr = np.array(Fr)
good = np.isfinite(Fh) & np.isfinite(Fr) & (Fh > 0) & (Fr > 0)

Fh_g = Fh[good]
Fr_g = Fr[good]

print("Number of good stars used:", len(Fh_g))

# also keep x positions for the same good stars
x_good = x[good]

# first scatter plot with left/right colours
left_good = x_good < w1
right_good = x_good >= (w1 + gap)

plt.figure()
plt.scatter(Fr_g[left_good], Fh_g[left_good], s=25, color='blue', label='Left cutout')
plt.scatter(Fr_g[right_good], Fh_g[right_good], s=25, color='green', label='Right cutout')

k, b = np.polyfit(Fr_g, Fh_g, 1)
x_line = np.linspace(Fr_g.min(), Fr_g.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color='red', linewidth=1)

print("k, b =", k, b)
plt.title('Pixel Intensities of Ha vs R')
plt.xlabel("R-band (pixel count)")
plt.ylabel("Hα-band (pixel count)")
plt.grid()
plt.legend()
plt.show()

# ------------------ manual removal of outliers ----------------- 
remove_outliers = ~(((Fr_g < 3000) & (Fh_g > 800)) |
                    ((Fr_g < 4000) & (Fh_g < 500)))

Fr_final = Fr_g[remove_outliers]
Fh_final = Fh_g[remove_outliers]
x_final = x_good[remove_outliers]

print("Number after outlier removal:", len(Fh_final))

# define left/right again after outlier removal
left_final = x_final < w1
right_final = x_final >= (w1 + gap)

plt.figure()
plt.scatter(Fr_final[left_final], Fh_final[left_final],
            s=25, color='blue', label='Left cutout')

plt.scatter(Fr_final[right_final], Fh_final[right_final],
            s=25, color='green', label='Right cutout')

k, b = np.polyfit(Fr_final, Fh_final, 1)

x_line = np.linspace(Fr_final.min(), Fr_final.max(), 200)
y_line = k * x_line + b
plt.plot(x_line, y_line, color='red', linewidth=1)

print("k, b =", k, b)
plt.title('Pixel Intensities of Ha vs R')
plt.xlabel("R-band (pixel count)")
plt.ylabel("Ha-band (pixel count)")
plt.grid()
plt.legend()
plt.show()

# ----- load master images -----
file5 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_900px_Rnew.fits')
img_R_band = file5[0].data

file6 = fits.open(r'C:/Users/wafiy/T.W/Aligned Photos/cutout_900px_Hanew.fits')
img_H_band = file6[0].data

# ----- subtraction -----
new_Ha_data = img_H_band - k * img_R_band

vmin = np.nanpercentile(new_Ha_data, 0.5)
vmax = np.nanpercentile(new_Ha_data, 99.8)

plt.figure()
plt.imshow(new_Ha_data, origin="lower", cmap="magma", vmin=vmin, vmax=vmax)
plt.colorbar()
plt.title("Continuum-subtracted Hα image")
plt.show()

# optional save
# header = file6[0].header.copy()
# fits.writeto("Ha_subtracted.fits", new_Ha_data, header=header, overwrite=True)