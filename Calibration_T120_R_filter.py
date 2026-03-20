# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:34:35 2025

@author: lily zaidi
"""#neils MasterMedian_R.fits 11-11-25
from astropy import units as u
import numpy as np
from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits

bias1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0001.fits')
bias2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0002.fits')
bias3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0003.fits')
bias4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0004.fits')
bias5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0005.fits')
bias6 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0006.fits')
bias7 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0007.fits')
bias8 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0008.fits')
bias9 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0009.fits')
bias10 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120_Cal/BIAS-0010.fits')

bias1 = CCDData(bias1, unit='adu')
bias2 = CCDData(bias2, unit='adu')
bias3 = CCDData(bias3, unit='adu')
bias4 = CCDData(bias4, unit='adu')
bias5 = CCDData(bias5, unit='adu')
bias6 = CCDData(bias6, unit='adu')
bias7 = CCDData(bias7, unit='adu')
bias8 = CCDData(bias8, unit='adu')
bias9 = CCDData(bias9, unit='adu')
bias10 = CCDData(bias10, unit='adu')

MasterBias = ccdproc.combine([bias1, bias2, bias3, bias4, bias5,\
                              bias6, bias7, bias8, bias9, bias10],
'MasterBias.fit', 'median')


flat1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat R-0001R.fits')
flat2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat R-0001R.fits')
flat3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat R-0001R.fits')
flat4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat R-0001R.fits')
flat5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat R-0001R.fits')

flat1 = CCDData(flat1, unit ='adu')
flat2 = CCDData(flat2, unit ='adu')
flat3 = CCDData(flat3, unit ='adu')
flat4 = CCDData(flat4, unit ='adu')
flat5 = CCDData(flat5, unit ='adu')

MasterFlat = ccdproc.combine([flat1, flat2, flat3, flat4, flat5],\
'MasterFlatRcopy.fit', 'median')
    
dark1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Dark-0001.fits')
dark2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Dark-0002.fits')
dark3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Dark-0003.fits')
dark4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Dark-0004.fits')
dark5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Dark-0005.fits')

dark1 = CCDData(dark1, unit = 'adu')
dark2 = CCDData(dark2, unit = 'adu')
dark3 = CCDData(dark3, unit = 'adu')
dark4 = CCDData(dark4, unit = 'adu')
dark5 = CCDData(dark5, unit = 'adu')

MasterDark = ccdproc.combine([dark1, dark2, dark3, dark4, dark5],\
'MasterDark60scopy.fit', 'median')
    
#---------------------------------------------(MasterBias)
input_filename1 = 'MasterBias.fit'
output_filename1 = 'MasterBias_bin2.fit'

def rebin_array(arr, factor):
    """Rebin a 2D array by an integer factor."""
    shape = (arr.shape[0] // factor, factor, arr.shape[1] // factor, factor)
    return arr[:shape[0]*factor, :shape[2]*factor].reshape(shape).mean(-1).mean(1) #* factor**2  # you may need to remove the * factor**2 depending on how the data are binned at OHP
# Read FITS file
with fits.open(input_filename1) as hdul:
    
    hdu = hdul[0]
    data = hdu.data
    header = hdu.header.copy()
    
# Rebin by factor of 2
factor = 2
rebinned_data = rebin_array(data, factor)

# Update header for new binning
if 'CRPIX1' in header:
    header['CRPIX1'] = (header['CRPIX1'] - 0.5) / factor + 0.5

if 'CRPIX2' in header:
    header['CRPIX2'] = (header['CRPIX2'] - 0.5) / factor + 0.5

if 'CDELT1' in header:
    header['CDELT1'] *= factor

if 'CDELT2' in header:
    header['CDELT2'] *= factor

if 'NAXIS1' in header:
    header['NAXIS1'] = rebinned_data.shape[1]

if 'NAXIS2' in header:
    header['NAXIS2'] = rebinned_data.shape[0]

# Save rebinned image
fits.writeto(output_filename1, rebinned_data, header, overwrite=True)


sci1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/R_Filter_M33/M33-0001R.fits')
sci2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/R_Filter_M33/M33-0002R.fits')
sci3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/R_Filter_M33/M33-0003R.fits')
sci4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/R_Filter_M33/M33-0004R.fits')
sci5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/R_Filter_M33/M33-0005R.fits')

sci1 = CCDData(sci1, unit = 'adu')
sci2 = CCDData(sci2, unit = 'adu')
sci3 = CCDData(sci3, unit = 'adu')
sci4 = CCDData(sci4, unit = 'adu')
sci5 = CCDData(sci5, unit = 'adu')

MasterBias_bin2 = CCDData.read('MasterBias_bin2.fit', unit='adu')

BiasSubDark = ccdproc.subtract_bias(MasterDark, MasterBias_bin2)
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias_bin2)
BiasSubSci1 = ccdproc.subtract_bias(sci1, MasterBias_bin2)
BiasSubSci2 = ccdproc.subtract_bias(sci2, MasterBias_bin2)
BiasSubSci3 = ccdproc.subtract_bias(sci3, MasterBias_bin2)
BiasSubSci4 = ccdproc.subtract_bias(sci4, MasterBias_bin2)
BiasSubSci5 = ccdproc.subtract_bias(sci5, MasterBias_bin2)

DarkSubFlat = ccdproc.subtract_dark(BiasSubFlat, BiasSubDark,
                                    dark_exposure = (60 * u.second),\
                                        data_exposure = (1* u.second),
                                        scale = True)

DarkSubSci1 = ccdproc.subtract_dark(BiasSubSci1, BiasSubDark,
                                    dark_exposure = (60 * u.second),\
                                        data_exposure = (120 * u.second),
                                        scale = True)

DarkSubSci2 = ccdproc.subtract_dark(BiasSubSci2, BiasSubDark,
                                    dark_exposure = (60* u.second),\
                                        data_exposure = (120 * u.second),
                                        scale = True)
    

DarkSubSci3 = ccdproc.subtract_dark(BiasSubSci3, BiasSubDark,
                                    dark_exposure = (60* u.second),\
                                        data_exposure = (120 * u.second),
                                        scale = True)
                                    

DarkSubSci4 = ccdproc.subtract_dark(BiasSubSci4, BiasSubDark,
                                    dark_exposure = (60* u.second),\
                                        data_exposure = (120 * u.second),
                                        scale = True)
    

DarkSubSci5 = ccdproc.subtract_dark(BiasSubSci5, BiasSubDark,
                                    dark_exposure = (60* u.second),\
                                        data_exposure = (120 * u.second),
                                        scale = True)

FinalSci1 = ccdproc.flat_correct(DarkSubSci1, DarkSubFlat)
FinalSci2 = ccdproc.flat_correct(DarkSubSci2, DarkSubFlat)
FinalSci3 = ccdproc.flat_correct(DarkSubSci3, DarkSubFlat)
FinalSci4 = ccdproc.flat_correct(DarkSubSci4, DarkSubFlat)
FinalSci5 = ccdproc.flat_correct(DarkSubSci5, DarkSubFlat)


FinalSci1 = np.asarray(FinalSci1)
FinalSci1 = (FinalSci1/2.0)
FinalSci1 = CCDData(FinalSci1, unit = 'adu')

FinalSci2 = np.asarray(FinalSci2)
FinalSci2 = (FinalSci2/2.0)
FinalSci2 = CCDData(FinalSci2, unit = 'adu')

FinalSci3 = np.asarray(FinalSci3)
FinalSci3 = (FinalSci3/2.0)
FinalSci3 = CCDData(FinalSci3, unit = 'adu')

FinalSci4 = np.asarray(FinalSci4)
FinalSci4 = (FinalSci4/2.0)
FinalSci4 = CCDData(FinalSci4, unit = 'adu')

FinalSci5 = np.asarray(FinalSci5)
FinalSci5 = (FinalSci5/2.0)
FinalSci5 = CCDData(FinalSci5, unit = 'adu')

MasterSci = ccdproc.combine([FinalSci1, FinalSci2, FinalSci3, FinalSci4,FinalSci5],\
'MasterMedian_R.fits',\
'median')

MasterScience = ccdproc.cosmicray_lacosmic(MasterSci, sigclip=3.0)
MasterScience.write('MasterScienceR10th.fits')

