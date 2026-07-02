#!/usr/bin/env python3
#@file catalogMaker.py
#@brief Creates catalogs readable by phosim from WIYN_ODI files
#
#@brief  Created by:
#@author Caleb T. Remocaldo (Purdue)

from astropy.io import fits
import sys, getpass, glob

print("Building catalogs")

filterType = ['u', 'g', 'r', 'i', 'z']
user = getpass.getuser()

path = sys.argv[1] #Passed from ../validation_4F_pipeline
wiynFiles = glob.glob(path+"/*.fits.fz")

for j in range(0, len(wiynFiles)):
        fitsFile = fits.open(wiynFiles[j], memmap=True)

        #Phosim catalogs associate each filter with a number
        for i in range(0, 5):
                if filterType[i] in fitsFile[0].header['PHOTFILT']:
                        filt = str(i)
                        break
                
        fileName = "/Users/" + user  + "/phosim_core/validation/validation_4F_" + "{num:03}".format(num=j) + "_catalog"

        #Converting RA and DEC to decimal degrees
        RA = fitsFile[0].header['RA'].split(':')
        rightAscension = float(RA[0]) + (float(RA[1])/60) + (float(RA[2])/3600)

        DEC = fitsFile[0].header['DEC'].split(':')
        declination = float(DEC[0]) + (float(DEC[1])/60) + (float(DEC[2])/3600)

        altitude = 90 - fitsFile[0].header['ZD']

        #Writing values in
        with open(fileName, 'w+') as f:
                f.write("rightascension\t" + str(rightAscension) + "\n" +
                        "declination\t" + str(declination) + "\n" +
                        "sunalt\t" + str(fitsFile[0].header['SUN__ALT']) + "\n" +
                        "dist2moon\t" + str(fitsFile[0].header['MOON_D']) + "\n" +
                        "moonphase\t" + str(fitsFile[0].header['MOONPHSE']) + "\n"+
                        "moonalt\t" + str(fitsFile[0].header['MOON_ALT']) + "\n" +
                        "altitude\t" + str(altitude) + "\n" +
                        "azimuth\t" + fitsFile[0].header['AZMAP'] + "\n" +
                        "mjd\t" + str(fitsFile[0].header['MJD-MID']) + "\n" +
                        "nsnap\t1" + "\n" +
                        "vistime\t60" + "\n" +
                        "filter\t" + filt + "\n" +
                        "obshistid\t" + "{num:03}".format(num=j)  + "\n" +
                        "object\t0\t" + str(rightAscension) + "\t" + str(declination) + " 35 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none")

        fitsFile.close()

print("Catalogs complete")
