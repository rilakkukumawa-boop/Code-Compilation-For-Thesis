# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:14:04 2025

@author: lily zaidi
"""

from astropy import units as u
import numpy as np
from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits

flat1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat V-0001V.fits')
flat2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat V-0002V.fits')
flat3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat V-0003V.fits')
flat4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat V-0004V.fits')
flat5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/Calibration Files/T120 Cal/Flat V-0005V.fits')

flat1 = CCDData(flat1, unit ='adu')
flat2 = CCDData(flat2, unit ='adu')
flat3 = CCDData(flat3, unit ='adu')
flat4 = CCDData(flat4, unit ='adu')
flat5 = CCDData(flat5, unit ='adu')

MasterFlat = ccdproc.combine([flat1, flat2, flat3, flat4, flat5],\
'MasterFlatVcopy.fit', 'median')
    
sci1 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/V_Filter_M33/M33-0001Vr.fits')
sci2 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/V_Filter_M33/M33-0002Vr.fits')
sci3 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/V_Filter_M33/M33-0003Vr.fits')
sci4 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/V_Filter_M33/M33-0004Vr.fits')
sci5 = fits.getdata('C:/Users/wafiy/Desktop/CopyLILYUSB/V_Filter_M33/M33-0005Vr.fits')

sci1 = CCDData(sci1, unit = 'adu')
sci2 = CCDData(sci2, unit = 'adu')
sci3 = CCDData(sci3, unit = 'adu')
sci4 = CCDData(sci4, unit = 'adu')
sci5 = CCDData(sci5, unit = 'adu')

MasterDark = CCDData.read('MasterDark60scopy.fit', unit='adu')
MasterBias_bin2 =CCDData.read('MasterBias_bin2.fit', unit='adu')

BiasSubDark = ccdproc.subtract_bias(MasterDark, MasterBias_bin2)
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias_bin2)
BiasSubSci1 = ccdproc.subtract_bias(sci1, MasterBias_bin2)
BiasSubSci2 = ccdproc.subtract_bias(sci2, MasterBias_bin2)
BiasSubSci3 = ccdproc.subtract_bias(sci3, MasterBias_bin2)
BiasSubSci4 = ccdproc.subtract_bias(sci4, MasterBias_bin2)
BiasSubSci5 = ccdproc.subtract_bias(sci5, MasterBias_bin2)    

DarkSubFlat = ccdproc.subtract_dark(BiasSubFlat, BiasSubDark,
                                    dark_exposure = (60 * u.second),\
                                        data_exposure = (0.59999999999999998* u.second),
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

MasterSci = ccdproc.combine([FinalSci1, FinalSci2, FinalSci3, FinalSci4, FinalSci5],\
'MasterMedianV.fits',\
'median')

MasterScience = ccdproc.cosmicray_lacosmic(MasterSci, sigclip=3.0)
MasterScience.write('MasterScienceV10th.fits')

    

                                    
   
