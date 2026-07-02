#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:10:12 2022

@author: anirban
"""

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import sys,os,shutil


#Read ra dec loaction of supernovae 
ra_nov = str(sys.argv[2])
dec_nov = str(sys.argv[3])
mag_nov = str(sys.argv[4])
req_exp = str(sys.argv[5])

#Accept the input image location
img_loc  = str(sys.argv[1])


#Read the input image
f=fits.open(img_loc)
data = np.array(f[0].data)
w = WCS(f[0].header)
hdr = f[0].header
if('CRPIX1' not in hdr ):
    print ('Warning ! Incomplete WCS. Significant error is likely to occur')
    
    
if('TOTEXPT' in hdr):
    exp = (f[0].header)['TOTEXPT']
else:
    print ('Could not find TOTEXPT keyword in header. Using default exposure of 60 sec')
    exp = 60

if('MAGZP' in hdr):
    zp = (f[0].header)['MAGZP']
else:
    print ('Could not find MAGZP keyword in header. Using default ZP=20')
    zp=20
f.close()


#Find the central ra dec
(y,x) =np.shape(data)
central_wcs = w.wcs_pix2world([[y/2, x/2]], 0)[0]
print ('RA DEC of central pixel is:',central_wcs)
#Conver nans to 0
data[np.isnan(data)] = 0
data[data<0] = 0

#Find background and gain
mean, median, sigma = sigma_clipped_stats(data, cenfunc= np.median)
bkg = median
gain = (sigma**2)/mean


#Subtract bkg 
data=data- bkg
data[data<0] = 0

#Make the filename  and copy file
#Copy is needed to to get the header 
dst_file = os. getcwd() +'/data/images/temp.fits'
fits.writeto(dst_file, data, hdr, overwrite= True)


counts = np.sum(data)
mag=-2.5*np.log10(counts/gain/exp) + zp
print ('Magnitude brightness is :',mag)


#Create run file 
run_file = os. getcwd()+'/examples/temp.txt'
f=open(run_file, 'w+')
f.write('rightascension '+str(central_wcs[0])[0:10] +' \n')
f.write('declination '+str(central_wcs[1])[0:10] +' \n')
f.write('altitude 89 '+'\n')
f.write('vistime '+str(req_exp) +'\n')
f.write('nsnap 1 '+'\n')
f.write('moonalt -90.0 '+'\n')
f.write('obshistid 1001 '+'\n')
f.write('object 0 '+str(central_wcs[0])[0:10]+ ' '+str(central_wcs[1])[0:10]+ ' '+str(mag)[0:6]+' ../sky/sed_flat.txt 0 0 0 0 0 0 temp.fits -1 0 none none '+'\n')
f.write('object 1 '+str(ra_nov)+' ' +str(dec_nov)+' '+str(mag_nov)+' ../sky/sed_flat.txt 0 0 0 0 0 0 star none none '+'\n')
 
f.close()