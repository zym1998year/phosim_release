
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 18:44:15 2022

@author: anirban
"""
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import sys,os,shutil
zpArr=['ZP', 'PHOTZP', 'MAGZP']
expArr=['EXP', 'EXPTIME', 'EXPOSURE', 'EXPT', 'TOTEXPT' ]
hubbleHdr = ['PHOTFLAM', 'PHOTPLAM']
argLen = len(sys.argv)
#print (len(sys.argv))
#Accept the input image location
img_loc  = str(sys.argv[1])
req_exp = float(sys.argv[2])
if(argLen == 3):
    mode = 0
elif(argLen == 6):
    mode=1
else:
    mode= None
    
print (mode)
#Test cases
#img_loc = '/home/anirban/m35_40min_red.fits'
#req_exp = 60
#mode = 0
#zp = 0
#End test cases 

#Read the input image
f=fits.open(img_loc)
data = np.array(f[0].data)
w = WCS(f[0].header)
hdr = f[0].header
if('CRPIX1' not in hdr ):
    print ('Warning ! Incomplete WCS. Significant error is likely to occur')
f.close()

###########################################################################################
if(mode == 0):
    #Search for Exposure time in image. If not available ask user
    exp=''
    for key in expArr:
        if(key not in hdr):
            continue
        if('sec' in hdr.comments[key]):
            exp = float(hdr[key])
        elif('min' in hdr.comments[key] or 'minute' in hdr.comments[key]):
            exp = float(hdr[key])*60
        elif('hour' in hdr.comments[key] ):
            exp = float(hdr[key])*3660
        else:
            val = str(input("Could not determine if exposure time is in seconds or minutes or hours. Do you wish to use seconds as default [y/n]?"))
            if(val =='y' or val=='Y'):
                exp = float(hdr[key])
            elif(val =='n' or val=='N'):
                exp = float(input("Enter exposure time(in seconds): "))
            else:
                print ('Invalid Input')
                sys.exit()
            
    if (exp == ''):   
        exp = float(input("Could not find exposure time of input image. Enter exposure time(in seconds): "))
    
    #Search for Zero Point in image. If not available ask user
    hubbleFlag = 0 #This flag turns to 1 for hubble image
    zp=''
    if ('PHOTFLAM' in hdr and 'PHOTPLAM' in hdr):
        hubbleFlag = 1
        print ('Hubble Image Detected')
        zp = -2.5*np.log10(hdr['PHOTFLAM']) - 5*np.log10(hdr['PHOTPLAM']) - 2.408
    
    #If not hubble image look for other keys
    if(zp == ''):
        for key in zpArr:
            if(key not in hdr):
                continue
            if(float(hdr[key])>40 or float(hdr[key])<0 ):
                continue
            zp = float(hdr[key])
    
       
    if (zp == ''): 
        valKey = str(input("Could not find Zero Point of input image. Press a to input Aperture or z to enter ZP: "))
        if(valKey == 'A' or valKey == 'a'):
            val = str(input("Enter Aperture(in cm) : "))
            dia = float(val)
            area = np.pi* (dia/2)**2
            h_cgs = 6.62e-27
            zp = -2.5 *np.log10((h_cgs*500)/ (100*area))-48.6
        elif(valKey == 'z' or valKey == 'Z'):
            val = str(input("Enter ZP : "))
            zp = float(val)
        else:
            print ('Invalid Input')
            sys.exit()

elif(mode == 1) :
    exp = float(sys.argv[3])
    valKey = str(sys.argv[4])
    val = str(sys.argv[5])
    if(valKey == 'A' or valKey == 'a'):
        dia = float(val)
        area = np.pi* (dia/2)**2
        h_cgs = 6.62e-27
        zp = -2.5 *np.log10((h_cgs*500)/ (100*area))-48.6
    elif(valKey == 'z' or valKey == 'Z'):
        zp = float(val)
    else:
        print ('Invalid Input. Use a or z. Check Documentation.')
        sys.exit()
else:
    print ('Invalid numer of arguments. Check Documentation. ')
    sys.exit()
####################################################################################    
        
print (zp , exp) 
#Find the central ra dec
(y,x) =np.shape(data)
central_wcs = w.wcs_pix2world([[y/2, x/2]], 0)[0]
print ('RA DEC of central pixel is:',central_wcs)

print(os. getcwd())

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
dst_file = os. getcwd()+'/data/images/temp.fits'
fits.writeto(dst_file, data, hdr, overwrite= True)
print(os. getcwd())

counts = np.sum(data)
mag=-2.5*np.log10(counts/gain/exp) + zp


#Create run file 
run_file = os. getcwd()+'/examples/temp.txt'
f=open(run_file, 'w+')
f.write('rightascension '+str(central_wcs[0])[0:10] +' \n')
f.write('declination '+str(central_wcs[1])[0:10] +' \n')
f.write('vistime '+str(req_exp) +'\n')
f.write('obshistid 1001 '+'\n')
f.write('object 0 '+str(central_wcs[0])[0:10]+ ' '+str(central_wcs[1])[0:10]+ ' '+str(mag)[0:6]+' ../sky/sed_flat.txt 0 0 0 0 0 0 temp.fits -1 0 none none '+'\n')
 
f.close()



