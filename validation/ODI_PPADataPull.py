#!/usr/bin/env python3
#@file ODI_PPADataPull.py
#@brief: Pulls image data from ODI PPA one image at a time, then writes relevant data to catalogs and realData.txt
#@author: Caleb Remocaldo (Purdue)

from subprocess import call
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import glob, getpass

def main():
    user = getpass.getuser()
    path = '/users/' + user + '/phosim_core/validation/temp/'
    chipMedians = []
    filterType = ['u', 'g', 'i', 'r', 'z']
    
    with open('fileList.txt') as f:
        j = 0
        for link in f:
            # Download image
            call('wget ' + link[:-1] + ' -P ' + path, shell=True)

            # Open image
            fitsFile = fits.open(path + link[76:-1])
            
            #Phosim catalogs associate each filter with a number
            for i in range(5):
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

            j += 1

            print('Calculating statistics...\n')
            with open('/users/' + getpass.getuser() + '/phosim_core/validation/realData.txt', 'a') as f:
                for x in range(1, 31):
                    mean, median, std = sigma_clipped_stats(fitsFile[x].data)
                    chipMedians.append(median)

                #intensity = ((np.mean(chipMedians)*100.0/EXPTIME)>1.0)*EXPTIME #Takes the maximum value of calculated intensity and 1
                intensity = np.mean(chipMedians)*100.00

                f.write(str(fitsFile[0].header['MJD-OBS']) + '\t' + str(intensity) + '\t' + fitsFile[0].header['PHOTFILT'] + '\n')

            fitsFile.close()
            
            call('rm -f ' + path + '*fits.fz', shell=True)
            
main()
