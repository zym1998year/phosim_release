#!/usr/bin/env python

"""
  @package phosim
  @file generatePlateScale.py
  @brief python script to read in an optics_*.txt file, determine the plate
  @scale (and f number), 
   
  @brief Created by:
  @author Glenn Sembroski (Purdue)

Usage: python generatePlateScale optics_*.txt
Aumption is that this is run from the phosim home directory.
Ex:
$> python data/tools/generatePlateScale.py data/lsst/optics_0.txt

This program uses the paraxial formalism to determine the plate scale of an
optical system defined by an optics_*.txt file. It also determines the
actual position of the focal plane and the f number. The main source for this
can be found in: Introduction to Optics, 3rd Edition, Pedrotti,F et. al.,
17 August 2019, ISBN: 9781108552493, Ch 18, pg 396-418.

The paraxial method has several assumptions and approximations.
1: Optical system is axially symetric.
2. Initial meduium is air.= (not space, fix later)
3. Ignore central holes in all optics.
4. Paraxial assumtion: All light angles are close to 0 deg ( << 10 deg)
5. Paraxial assumtion: Any non sprerical suface corrections can be ignored.
6. Paraxial assumption: Small angle approximation everywhere.
"""

import sys, os, subprocess, shutil
import numpy as np


def paraxialTranslationMatrix(disz) :
    m = np.array( [ [1.0, disz] , [ 0.0, 1.0] ] )
    return m

def paraxialReflectanceMatrix( R ) :
    if ( R == 0.0 ) :
        #Flat mirror
        m = np.array( [ [1.0, 0.0] , [0, 1.0] ] )
    else :
        m = np.array( [ [1.0, 0.0] , [2.0/R, 1.0] ] )
    return m

def paraxialRefractionMatrix( R, n0, n1 ) :
    indexRatio = n0/n1
    #print("Index Ratio: " + str(indexRatio) )
    
    if ( R == 0.0 ) :    #Phosim indicates flat mirror (radius of curavtuire
                         #is infinite) with R == 0.0
        #flat lens/filter surface
        m = np.array( [ [1.0, 0.0] , [ 0, indexRatio] ] )
    else :
        m = np.array( [ [1.0, 0.0] ,
                                 [( indexRatio - 1.0 ) / R, indexRatio] ] )
    return m

def indexAtLambda(fileName, lambda_um, opticsPath ) :
    #The following is a  messy subroutine.
    #Problem is getting the index file itself. May be in:
    #1: local direcotry as M (if specified with *.txt format)
    #2: data/material/optical/ as M ((if specified with *.txt format and not
    #found in local dir)
    #3: data/matrial/optical/ as M.txt
    #4: as 'air' (we generate index from function)
    #5: as 'vacuum' n=1.

    #Get the index of refraction from the file(except for air and vacuum)
    n = 0.
    if (fileName == "vacuum" ) :
        n=1.0
        
    elif (fileName == "air") :
       #See:https://refractiveindex.info/tmp/database/data/other/mixed%20gases/air/nk/Ciddor.html
       x_nm = lambda_um*1000.
       n=1+0.05792105/(238.0185-x_nm**-2)+0.00167917/(57.362-x_nm**-2)
    else :
        #We look for the index or refraction file in the local directory
        filePath = opticsPath + fileName
        #print(filePath)
        if ( not  os.path.exists(filePath) ):
            #Try in the matrials/optical directory
            filePath = "data/material/optical/" + fileName
            #print(filePath)
            if ( not  os.path.exists(filePath) ):
                #Try with adding a .txt file extention
                filePath=filePath + ".txt"
                #print(filePath)
                if ( not  os.path.exists(filePath) ):
                    #Couldn't fine it . quit.
                    print("Fatal: Medium file for : " + fileName
                          + "from optics_0.txt file not found" )
                    sys.exit()
        #Found the file. Search for the index and interpolate
        firstTableIndex = True
        print(f"Attemping to read material file: {filePath}")
        for line in open(filePath,'r', encoding='utf-8-sig').readlines(): #To resolve some BOM issues with the optical materials file
            l=line.split()
            if ( (float( l[0] ) > lambda_um and firstTableIndex ) or
                 ( float( l[0] ) == lambda_um ) ):
                    n = float( l[1] )
                    break
            elif (float( l[0] ) > lambda_um ) :
                #We found the intervel , interoplate within it.
                lambdaNew = float( l[0] )
                indexNew = float( l[1] )
                #Interplolate
                f = ( lambda_um - lambdaOld ) / ( lambdaNew - lambdaOld )
                n = f * ( indexNew - indexOld ) + indexOld
                break
            else :
                firstTableIndex = False
                lambdaOld = float( l[0] )
                indexOld = float( l[1] )
                n = float( l[1] )   #This for when requested index is
                                    #beyond table
    return n

def sysMatrixGen( name,  curv, direction, disz, n0, n1, alphaDeg)  :
    #Iterate over the surfaces listed in the optics_*.txt file up to the last
    #one before the "det"

    sysMatrix = np.array( [ [0,0], [0,0] ] )

    #Create test vector to test for angle. Save max angle
    maxAlphaRad = 0.
    alphaRad = alphaDeg * np.pi / 180.
    vector = np.array( [ 0.0, alphaRad ] ) #We start on axis.

    #print ("len(name): " + str(len(name)) )
    for i in range( len( name ) - 1 ):
        #Now we generate the translation matrix to this surface and the 
        #matrix for the surface itself. If this is the first surface, then 
        #there is no translation. We start at the origen at the first
        #surface.
        #print(" name: " + name[ i ] + "\n curv: " + str( curv[ i ] )
        #      + "\n disz: " + str( disz[ i ] ) + "\n n0: "
        #      + str( n0[ i ] ) + "\n n1: " + str( n1[ i ] )
        #      + "\n alphaDeg: " + str( alphaDeg ) )
        if ( i != 0 ) : #No translation matrix for first surface
            translationMatrix = paraxialTranslationMatrix(disz[ i ])
            #print("translation Matrix:")
            #print(translationMatrix)
            #print('  Det: ' + str( np.linalg.det( translationMatrix ) ) )
            #Add it to the system matrix
            sysMatrix = np.matmul(translationMatrix, sysMatrix)
            #print("sysMatrix")
            #print(sysMatrix)
            #print('  Det: ' + str( np.linalg.det(sysMatrix) ) )
            result_vector = np.matmul( sysMatrix, vector )
            #print("Result vector:")
            #print( "[" + str(result_vector[0] )  + "  " +
            #       str(result_vector[ 1 ] * 180. / np.pi) + "]" )

            
            
        #Create the surface matrix
        #Now we have to determine if this surface is concave or convex.
        #Paraxial convention is concave urface has R negative and convex
        #surface ahs R positive
        #Phosim convention is that a positive curv is concave when direction
        #is negative(downwards) and convex when curv is negative
        #Thus
        #If phosim  direction is negative (direction == -1) and phosim
        #curv is positive this is a concave surface. Then Paraxial R == curv
        #If photon direction is positive and phosim cuvature is + this is a
        #paraxial concave surface. Paraxial  R = - curv

        #Determine sign of R
        R = -direction[i] * curv[i]
        #print("\n R: " + str( R ) )
        
        if (name[i] == "mirror") :
           surfaceMatrix = paraxialReflectanceMatrix( R )
        elif ( ( name[i] == "lens" ) or ( name[i] == "filter" ) ) :
            surfaceMatrix = paraxialRefractionMatrix( R, n0[i], n1[i] )

        #print("surfaceMatrix:")
        #print(surfaceMatrix)
        #print('  Det: ' + str(np.linalg.det(surfaceMatrix) ) )
        if ( i == 0 ) :
            sysMatrix = surfaceMatrix
        else :
            sysMatrix = np.matmul( surfaceMatrix, sysMatrix )

        #We now have the system matrix to this surface.
        #print( "sysMatrix:" )
        #print( sysMatrix )
        #print('  Det: ' + str(np.linalg.det(sysMatrix) ) )


        #Get the exit angle with our test vector as input.
        result_vector = np.matmul( sysMatrix, vector )
        #print("Result vector:")
        #print( "[" + str(result_vector[0] )  + "  " +
        #           str(result_vector[ 1 ] * 180. / np.pi) + "]" )

        #Save exit angle if max.
        exitAlphaRad = result_vector[ 1 ]
        if ( abs( result_vector[ 1 ] )  > maxAlphaRad ) :
             maxAlphaRad = abs( result_vector[ 1 ]  )

        #end of for loop

    maxAlphaDeg =   maxAlphaRad * 180. / np.pi
    return maxAlphaDeg, sysMatrix
                                    
                                    
#######################
#
#Start of Main code#
#
#######################
#Check usage and input file existance.

if len(sys.argv) != 2:
    print("Usage: python GeneratePlateScale <optics_*.txt>")
    sys.exit(1)

# Initial input file from user
optics_file = sys.argv[1]
l = optics_file.split('/')
opticsPath = ""
for i in range( len (l)-1 ) :
    opticsPath = opticsPath + l[i] + "/"


 
if not os.path.isfile(optics_file):
    print(f"Error: File '{input_file}' not found.")
    sys.exit(1)
    
firstLine = True
firstSurface = True
sysMatrix = np.array( [ [0,0], [0,0] ] )

#Go though the optics_*.txt file surface by surface, collecting the needed
#values.
name = []
curv = []
disz = []
direction = []
aperture = []
n0 = []
n1 = []

direct =-1.0   #Inciates photon Z direction through the optics. Starts -
#(to match phosim convention)

for line in open(optics_file).readlines():
    #print(line)
    l=line.split()
    firstWord = l[0]
    
    if( firstWord[0]  == "#") :
        continue

    if(firstLine) :
        #The first line in the optics_*.txt file is not a surface but holds
        #various useful values for Phosim optimization. This include the
        #approximate plate scale (which may be a dummy value). Save this line.
        #We will replace with an identical line but with an updated platescale
        #which will be determined.
        properties=line   #save first line  for later restoration with correct
        #platescale
      
        print(optics_file + " properties:" )
        print("     Configuration  Name: " + l[1] )
        print("              lambda min: " + l[2] + " (um)" )
        print("              lambda max: " + l[3] + " (um)" )
        print("          Central lambda: " + l[4] + " (um)" )
        centerLambda= float(l[4])
        print(" Approximate Plate Scale: " + l[5] + " (um)/deg" )
        initialPlateScale = float(l[5])
        print("        X Focalplane Dir: " + l[6] )
        print("        Y Focalplane Dir: " + l[7] )
        firstLine = False
        n = indexAtLambda( "air", centerLambda, opticsPath )  #Fix later, Might be vacuum

    else :
        
        #Now fill up the surace specification arrays from the optics_*.txt
        #file.
        #(Note  Test at each surface the angle of the photon of the test
        #ray to insure it doesen't exced some maximum value (where the
        #paraxil approximations may start to be violated.)

        if ( l[1] == "mirror"  ) or (l[1] == "lens") or (l[1] == "filter") or ( l[1] == "det"  ) :
            #Good  surface, Save optical values
            name.append( l[1] )
            #print( "len(name): " + str(len(name)) )

            #Curvature of surface in mm. Note PHOSIM has curvature
            #set to 0 if it is a flat surface(curvature=infinity)
            #We will follow this convention
            #The PHOSIM  convention for curvature sign is concave in the + z
            #direction is +,  and - for convex.
            #Curvature for Paraxial is + for concave when Phosim photon
            # direction is -. For Phosim direction + concave needs - curvature
            #Thus we need to keep track of photon direction through the optics
            #and change PHOSIM mirror curvature sign as needed.

            curv.append( float( l[2] ) )

            direction.append(direct)
            #Photon direction changes sign after a mirror
            if (l[1] == "mirror" ) :
                direct = direct * -1.0   #For next surface
                     
            #Distance from previous surface.
           
            disz.append( float( l[3] ) )  
                     #distance from previous surface in mm. Sign may indicate
                     #ray direction unless this is 0.0. Thus the use of the
                     #direction list

            aperture.append(float ( l[4] )  )
                     #outer diameter of this surface in mm
            #Get the index of refraction that before and after this surface at
            #the centerLambda
            n0.append( n )        #previous index of refraction
            
            #Note here that the plate scale if for the central wavelength
            n = indexAtLambda( l[23] , centerLambda, opticsPath)  #Next index of refracrtion
            n1.append( n )

            if ( l[1] == "det") :  #Quit when we reach the detector
               break
#End of optics_*.txt read 'for' loop
        
#Generate system matrix to last surface before "det"
#Check that our off axis ray is always less then 1 deg at all
#surfaces. Reduce inital angle (at first surface) until it meets this
#criteria (paraxial assumption)

#alphaDeg = .01   #First guess
alphaDeg = .02   #First guess
maxAlphaDeg, sysMatrix = sysMatrixGen(name, curv, direction, disz, n0, n1,
                                      alphaDeg)
#returns maximum alpha at any surface when we start with alphaDeg at the
#first surface  and it also returns the System Matix to the surface before
#the detector

#Test that this was a good starting alphDeg. This will be true if
#maxAlphaDeg at all surfaces is less than 1 deg. If not we can scale
#alphaDeg and try again.
if (maxAlphaDeg >1.0 ) :
     alphaDeg = alphaDeg * ( 1.0 / maxAlphaDeg )
     #I know this recalculating the sysMatrix is inefficent but we do it
     #to find the maxAlpha to make sure we remain within the paraxil
     #assumptions.
     maxAlphaDeg,sysMatrix = sysMatrixGen(name, curv, direction, disz, n0,
                                          n1, alphaDeg)
#print("For AlphaDeg: " + str(alphaDeg) + " maxAlphaDeg is: " +
#       str(maxAlphaDeg) )

        
#1.Determine the distance to the actual focal plane which may be
#  different from the detector placement. TBD
#1.Complete the system Matrix to the detector and determine the
#  platescale
#3.Use the aperture of the first surface to come up with a '1/f'
#  number for the system.TBD

#Determine distance to focal plane from last surface.

#Print the system matrix at last surface
#print(sysMatrix)
#print('  Det: ' + str(np.linalg.det(sysMatrix) ) )
#Distance to focal plane occurs when distance from last surface to the
#focal plan causes the first element of sysMatrix with the translation
#to the focal plane causes A (first element of sysMatrix to be 0.
#i.e. A + disz(focalplane) = 0 Thus disz(focalplane)=-A

detIndex = len( name ) - 1
#The q is the distance from the output plane to the focal plane
q  = direction[ detIndex ]* sysMatrix[1,0] / sysMatrix[1,0]   #-(A/C)

#print("Best focal plane at " + str(q) +
#      "mm from last surface")
#print("Detector at: " + str( disz[ detIndex ] ) + " mm")

#determine f number
#Focal equivalrent focal length is -1/C for sysMatrix = [ [A,B], [C, D] see
#reference top of page
#print( "sysMatrix[1,0]; " + str(sysMatrix[1,0])  )
       
equivFocalLength_mm = direction[detIndex] * 1./sysMatrix[1,0]
                                                          #-(1/C)
# f unmber is equivalent focal length/diameter of aperture
#print("Aperture: " + str( aperture[0] ) ) 
fNumber = equivFocalLength_mm/ (2 * aperture[0] )
#print( "Equivalent focal length: " + str(equivFocalLength_mm) + " mm")
#print( "F number: f/" + str(fNumber))

#Translation to detector surface.
#print("At Detector:")

translationMatrix = paraxialTranslationMatrix( disz[ detIndex ] )
#print("translationMatrtx:")
#print(translationMatrix)

sysMatrix = np.matmul(translationMatrix, sysMatrix)
#print(sysMatrix)
#print('  Det: ' + str(np.linalg.det(sysMatrix) ) )
#We are now ready to determine the plate scale.
#We assume here that the system is axialy symetric and that a photon
#enering the system on the axis with 0 inclination will arive at 0,0
#when it reaches the detector.
#We then can see where a photon on axis but with an inclination of
#alphaDeg 
alphaRad = ( alphaDeg / 180. ) * np.pi

#Create input vector
vector = np.array([0.0, alphaRad])
#vector = np.array([1000.0, alphaRad])
        
#find ouput vector
result_vector = np.dot(sysMatrix, vector)
#print( "[" + str(result_vector[0] )  + "  " +
#                   str(result_vector[ 1 ] * 180. / np.pi) + "]" )

#Height this photon lands on detector in um is thus
y_um=result_vector[0]*1000.

#determine platescale
plateScale = y_um/alphaDeg
#and thats all we need
print("Old plate scale: " + str(initialPlateScale) + "um/deg" )
print("Found plate scale: " + str(plateScale) + "um/deg" )

#Thick
#f1=1./(((1.46232-1.0002726)/1.0002726)*(1./100. + 1./100. +(1.0002726/1.46232-1.)*10./( 100.*100.)))
#Thin
#f1=1./(((1.46232-1.0002726)/1.0002726)*(1./100. + 1./100.))

#print("f1= " + str(f1) )
#print("f2= " + str(-f1) )
print("done")
