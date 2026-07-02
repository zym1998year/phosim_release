#!/usr/bin/env python
"""
  @ @package phosim
  @ @file material_generate
  @ @brief script to genrate optical material file using Sellmeier1 formula
  @
  @ @brief Created by:
  @ @author Glenn Sembroski (Purdue)
  @
  @
  @warning This code is not fully validated
        and not ready for full release.  Please
        treat results with caution.
  @
  @ Usage Example:
  @ python material_generate.py <glass name> <glass database file>
  @
  @ This code will return: ubk7.txt
"""
import sys, os, math


#Impliments calculation of index of refraction using sellmeier 1 formula
#see: https://neurophysics.ucsd.edu/Manuals/Zemax/ZemaxManual.pdf page 590
#lambda in um.
def Sellmeier1(lmbda, kparam, lparam ) :
    lambda2=lmbda * lmbda
    n2 = 1.0 + ( ( kparam[0] * lambda2) / ( lambda2 - lparam[0] ) ) + ( ( kparam[1] * lambda2) / ( lambda2 - lparam[1] ) ) + ( ( kparam[2] * lambda2) / ( lambda2 - lparam[2] ) )
    n = math.sqrt(n2)
    return n


#######################
#
#Start of Main code#
#
#######################
#We are going to create a phosim compatable optical 'glass.txt' type file for
#the glass material name given in the first argument (typically this for a lens
#material, but not always).
#Output file name will be <glass name>.txt
#format of <glass name>.txt output file:

#  wavelentgh(um)   Index of refraction  Imaginary index of refraction(opacity)

#Each row is for a particular wavelength. Wavelengths go from lambda_mim (um)
#to lambda_max(um). The step size in lambda has been arbitrarliy been set to
#~5 nm (.005 um). lambdaMin and lambdaMax will be found in the
#<glass data base file>,

#The index of refraction is cauluated using the Sellmeier 1 formula (See ref
#above), which uses 6 coefficients. These coeficients will be found in the
#<glass data base file>.

#The name of the material type is the first argument.
material = sys.argv[1]
#Name of the glass data base file is in the second argument
glassDataBaseFile = sys.argv[2]

#Find line in the glass data base file for the specified glass
for line in open(glassDataBaseFile).readlines() :
    #Find glass
    if material == line.split()[0] :
        #print ( line )
        lambda_min = float( line.split()[1] )
        lambda_max = float( line.split()[2] )
        lambda_step = 0.005          # in um
        lambdaNum = int ( ( lambda_max - lambda_min) / lambda_step ) 

        print('Lambda Min, LambdaMax, LambdaStep (um), lambdaNum:' +
              str(lambda_min) + ', '+ str(lambda_max) + ', ' + str(lambda_step) +
              ', ' + str(lambdaNum) )

        #Get the Sellmeier 1 coeficients for this material
        kparam = []
        lparam = []
        kparam.append( float( line.split()[3] ) )
        lparam.append( float( line.split()[4] ) )
        kparam.append( float( line.split()[5] ) )
        lparam.append( float( line.split()[6] ) )
        kparam.append( float( line.split()[7] ) )
        lparam.append( float( line.split()[8] ) )

        #Create and open the file for output
        outputFile=material + '.txt'
 
        #We will replace any existing file of this name
        try:
            os.remove(outputFile)
        except OSError:
            pass

        out = open(outputFile, 'a' )
        lmbda = (lambda_min - lambda_step)

        for i in range (1,lambdaNum) :
            lmbda = lambda_min + i * lambda_step
            n = Sellmeier1( lmbda, kparam , lparam)

            #Note imag part of index of refraction always set to 0.0. (mainly
            #because I don't know how to generate it)
            out.write('%-4f  %4f %s \n ' % ( lmbda, n, str(0.0) ) )

        out.close()
        print ('File: ' + outputFile + ' created.' )
        sys.exit()


#Failed to find the glass in the glass DataBase file
print('Failed to find ' + material + ' in ' + glassDataBaseFile )
sys.exit()


    

