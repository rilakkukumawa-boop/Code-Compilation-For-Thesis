# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:12:44 2025

@author: lily zaidi
"""

from astropy import units as u
import numpy as np
from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits


flat1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Ha Flat-0002Ha.fits')
flat2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Ha Flat-0003Ha.fits')
flat3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Ha Flat-0004Ha.fits')


flat1 = CCDData(flat1, unit ='adu')
flat2 = CCDData(flat2, unit ='adu')
flat3 = CCDData(flat3, unit ='adu')

MasterFlat = ccdproc.combine([flat1, flat2, flat3],\
'MasterFlatHacopy.fit', 'median')

sci1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Ha_Filter_M33/M33-0001Ha.fits')
sci2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Ha_Filter_M33/M33-0002Ha.fits')
sci3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Ha_Filter_M33/M33-0003Ha.fits')


sci1 = CCDData(sci1, unit = 'adu')
sci2 = CCDData(sci2, unit = 'adu')
sci3 = CCDData(sci3, unit = 'adu')
#
   
MasterDark = CCDData.read('MasterDark60scopy.fit', unit='adu')
MasterBias_bin2 =CCDData.read('MasterBias_bin2.fit', unit='adu')


BiasSubDark = ccdproc.subtract_bias(MasterDark, MasterBias_bin2)
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias_bin2)
BiasSubSci1 = ccdproc.subtract_bias(sci1, MasterBias_bin2)
BiasSubSci2 = ccdproc.subtract_bias(sci2, MasterBias_bin2)
BiasSubSci3 = ccdproc.subtract_bias(sci3, MasterBias_bin2)

DarkSubFlat = ccdproc.subtract_dark(BiasSubFlat, BiasSubDark,
                                    dark_exposure = (290* u.second),\
                                        data_exposure = (150* u.second),
                                        scale = True)

DarkSubSci1 = ccdproc.subtract_dark(BiasSubSci1, BiasSubDark,
                                    dark_exposure = (290 * u.second),\
                                        data_exposure = (300 * u.second),
                                        scale = True)

DarkSubSci2 = ccdproc.subtract_dark(BiasSubSci2, BiasSubDark,
                                    dark_exposure = (290* u.second),\
                                        data_exposure = (300 * u.second),
                                        scale = True)
    

DarkSubSci3 = ccdproc.subtract_dark(BiasSubSci3, BiasSubDark,
                                    dark_exposure = (290* u.second),\
                                        data_exposure = (300 * u.second),
                                        scale = True)
                                    


FinalSci1 = ccdproc.flat_correct(DarkSubSci1, DarkSubFlat)
FinalSci2 = ccdproc.flat_correct(DarkSubSci2, DarkSubFlat)
FinalSci3 = ccdproc.flat_correct(DarkSubSci3, DarkSubFlat)

FinalSci1 = np.asarray(FinalSci1)
FinalSci1 = (FinalSci1/2.0)
FinalSci1 = CCDData(FinalSci1, unit = 'adu')

FinalSci2 = np.asarray(FinalSci2)
FinalSci2 = (FinalSci2/2.0)
FinalSci2 = CCDData(FinalSci2, unit = 'adu')

FinalSci3 = np.asarray(FinalSci3)
FinalSci3 = (FinalSci3/2.0)
FinalSci3 = CCDData(FinalSci3, unit = 'adu')

MasterSci = ccdproc.combine([FinalSci1, FinalSci2, FinalSci3],\
'MasterMedianHa.fits',\
'median')

MasterScience = ccdproc.cosmicray_lacosmic(MasterSci, sigclip=3.0)
MasterScience.write('MasterScienceHa10th.fits')
