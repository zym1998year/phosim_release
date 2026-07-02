#!/usr/bin/env python3
#@file realDataMaker.py
#@brief: Create .txt file of real data for faster validation
#@author: Caleb Remocaldo (Purdue)

from astropy.io import fits
import sys, os.path, getpass, glob, numpy as np

chipMedians = []
path = sys.argv[1] #Passed from ../validation_4F_pipeline or command line

while not (os.path.exists(path)):
    path = str(raw_input('Path not found, please enter a different one: '))

print('Writing observation data to /validation/realData.txt')

realFiles = glob.glob(path+"*.fits.fz")
dest = '/users/' + getpass.getuser() + '/phosim_core/validation/realData.txt'
    
with open(dest, 'w+') as f:
    for i in range(0, len(realFiles)):
        fitsFile = fits.open(realFiles[i], memmap=True)

        #Tried to use a generator instead of list to improve speed, but most of the time
        #is spent calculating the mean, which is done regardless of using list of generator
        for j in range(1, 31):
            median = np.nanmedian(fitsFile[j].data)
            chipMedians.append(median)

        #intensity = ((np.mean(chipMedians)*100.0/EXPTIME)>1.0)*EXPTIME #Takes the maximum value of calculated intensity and 1
        intensity = np.mean(chipMedians)*100.00

        f.write(str(fitsFile[0].header['MJD-OBS']) + '\t' + str(intensity) + '\t' + fitsFile[0].header['PHOTFILT'] + '\n')
        fitsFile.close()

        print( str(i+1) + ' of ' + str(len(realFiles)) + ' completed\r')
