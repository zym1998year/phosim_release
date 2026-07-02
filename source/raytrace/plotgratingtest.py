#!/usr/bin/env python3

import sys
import os
import math
import matplotlib.pyplot as plt
import numpy as np

######################################################################
# Read in the file "gratingtest.dat" from local directory and plot the
# test result distrubutions.
# Lines that start with "#" should be ignored.
# Test 1: As function of anglebetween plot the acos of dot product of
# input vector and reconstructed output vector which shoul be identicle.
# Each line in the test data file  with element == 1.1  anglebetween is
# the 11th item on each line.
# Test 2.1 :
# probability intensity at central vavelength, with anglein ~= 0
# and angleout. Maybe skipped if Cumulative prob tablke was read form a file.
# Test 2.2 :
# Cumulative probability at central vavelength, with anglein ~= 0
# and angleout.
######################################################################

def GenerateDataFloatArray( dataAsStringList, useLimits = False, \
                            lowLimit = 0, highLimit = 1) :
    #Convert the line list to a numpy array of strings
    #then to a numpy array of floats
    dataAsStringArray = np.array(dataAsStringList)
    dataAsFloatArray = dataAsStringArray.astype(np.double)
    if ( not useLimits) :
        print('# in array: ' + str( len( dataAsStringList) ) )
        return dataAsFloatArray
    #set limits on the final array. Delete elements out of rangs
    dataAsFloatArrayAbove = np.delete(dataAsFloatArray, \
                            np.where(dataAsFloatArray < lowLimit) )
    dataAsFloatArrayInRange = np.delete(dataAsFloatArrayAbove, \
                              np.where(dataAsFloatArrayAbove > highLimit) )

    print('# in full array: ' + str( len( dataAsStringList) ) + \
             ' # in range: ' + str(dataAsFloatArrayInRange.size) )
    return dataAsFloatArrayInRange
################

def GenerateArrayPairInFirstRange( dataFirstArray, dataSecondArray, \
                            lowLimit, highLimit) :
    #Remove elements in first array that are not in range.
    dataFirstAbove = np.delete(dataFirstArray, \
                          np.where(dataFirstArray < lowLimit) )
    dataFirstInRange = np.delete(dataFirstAbove, \
                          np.where(dataFirstAbove > highLimit) )

    #Now remove the same elements in the Second arry so that the array pairs
    #are still the same
    dataSecondAbove = np.delete(dataSecondArray, \
                          np.where(dataFirstArray < lowLimit) )
    dataSecondInRange = np.delete(dataSecondAbove, \
                          np.where(dataFirstAbove > highLimit) )
    print('# in full array: ' + str( dataFirstArray.size)  + \
             ' # in range: ' + str(dataFirstInRange.size) )
    return dataFirstInRange, dataSecondInRange

#######################
def GenerateHiLoOfArray(AngleOutDegStr, width) :
    AngOutInRangeFull =  GenerateDataFloatArray(AngleOutDegStr)
    y, binEdges = np.histogram(AngOutInRangeFull, bins=100000)
    peakIndex=findPeakIndexSlidingWindow(y,windowSize)
    #peakIndex = np.argmax(y)
    bincenter = 0.5 * (binEdges[peakIndex] + binEdges[peakIndex-1])
    print('Bin Center at: ' + str(bincenter) + 'deg\n')
    peakLo=bincenter - width
    peakHi=bincenter + width
    return peakLo, peakHi
#######################################################################    
def findPeakIndexSlidingWindow(array,windowSize) :
    #Using a sligin window of width "windowSize" bins find peak index
    #(center of window)
    #This will be brute force for now.
    n = len(array)

    # Initialize result
    max_sum = 0
    max_index=0
    # Consider all blocks
    # starting with i.
    for i in range(n - windowSize + 1):
        current_sum = 0
        for j in range(windowSize):
            current_sum = current_sum + array[i + j]
 
        # Update result if required.
        if current_sum > max_sum :
            max_sum = max(current_sum, max_sum)
            #rounding down here. Best if windowSize is odd
            max_index=i + math.floor(windowSize/2)
    return max_index
# #########################



#############################
# Start of Main
#############################

#get the "gratingtest.dat" file path and read it in
testResultsPath = "./gratingtest.dat"
#Make sure this file exists and if so read it in.
if os.path.isfile(testResultsPath)  :
    testFile = open(testResultsPath,"r")
    testData = testFile.readlines()
    testFile.close()
else :
    print("Fatal-- ", testResultsPath, "file fails to exist")
    sys.exit()

# Set up the plotting
# Extract all the distibutions
#Test1:
dot = []           #Dot product of reconstructed output vector with expected
                   #output vector
#Test2:
probIntensity = [] #Probability intensity at central vavelength, anglen ~=0
                   #and along angleout.
intensityAngleOut = []
cumProb = []       # Cumulative probability, same as above
cumProbAngleOut = []

#Test3: For middle angleIN index to grating
minIndexAngleOutDeg = [] #Min wavelength (uses index, ie along table axis)
midIndexAngleOutDeg = [] #Middle wavelwngth Index
midPlus1IndexAngleOutDeg = [] # mid wavelength index +1
midPlus2IndexAngleOutDeg = [] # mid wavelength index +2
maxIndexAngleOutDeg = [] # maximum wavelength index in lookup table

#For various wavength stepped from 500 nm, for middle angle in index
AngleOutDeg41 = []
AngleOutDeg42 = []
AngleOutDeg43 = []
AngleOutDeg44 = []

#For wavelength of 500.75 nm , angleIn stepped from value at mid index.
AngleOutDeg50a = []
AngleOutDeg50b = []
AngleOutDeg51 = []
AngleOutDeg52 = []
AngleOutDeg53 = []
AngleOutDeg54 = []

Interference50a = []
Envelope50a = []
Interference50b = []
Envelope50b = []

AngleOutXDeg71 = []
AngleOutXDeg72 = []

WavelengthNM71 = []
WavelengthNM72 = []

windowSize = 5;   #best if odd

for line in testData :
    #ignore comments
    lineItems = line.split(' ')
    if not  lineItems[0].startswith('#'):     #Skips comments.
        if lineItems[0] == "1.1" :
            dot.append(lineItems[11])
 #Test2:
        elif lineItems[0] == "2.1" :
            intensityAngleOut.append(lineItems[3])
            probIntensity.append(lineItems[4])
            Wavelength21 = lineItems[1]
            angleInDeg21 = lineItems[2]
        elif lineItems[0] == "2.2" :
            cumProbAngleOut.append(lineItems[3])
            cumProb.append(lineItems[4])
            Wavelength22 = lineItems[1]
            angleInDeg22 = lineItems[2]
            
 #Test3:
        elif lineItems[0] == "3.1" :
            minIndexAngleOutDeg.append(lineItems[5])
            Wavelength31 = lineItems[2]
            angleInDeg31 = lineItems[3]
            if ( abs( float(angleInDeg31) ) < 1.e-5 ) :
                 angleInDeg31 = "0.0"
        elif lineItems[0] == "3.2" :
            midIndexAngleOutDeg.append(lineItems[5])
            Wavelength32=lineItems[2]
            angleInDeg32 = lineItems[3]
            if ( abs( float(angleInDeg32) ) < 1.e-5 ) :
                 angleInDeg32 = "0.0"
        elif lineItems[0] == "3.3" :
            midPlus1IndexAngleOutDeg.append(lineItems[5])
            Wavelength33=lineItems[2]
            angleInDeg33 = lineItems[3]
            if ( abs( float(angleInDeg33) ) < 1.e-5 ) :
                 angleInDeg33 = "0.0"
        elif lineItems[0] == "3.4" :
            midPlus2IndexAngleOutDeg.append(lineItems[5])
            Wavelength34=lineItems[2]
            angleInDeg34 = lineItems[3]
            if ( abs( float(angleInDeg34) ) < 1.e-5 ) :
                 angleInDeg34 = "0.0"
        elif lineItems[0] == "3.5" :
            maxIndexAngleOutDeg.append(lineItems[5])
            Wavelength35=lineItems[2]
            angleInDeg35 = lineItems[3]
            if ( abs( float(angleInDeg35) ) < 1.e-5 ) :
                 angleInDeg35 = "0.0"

 #Test4
        elif lineItems[0] == "4.1" :
            AngleOutDeg41.append(lineItems[4])
            Wavelength41 = lineItems[1]
            angleInDeg41 = lineItems[2]
            if ( abs( float(angleInDeg41) ) < 1.e-5 ) :
                 angleInDeg41 = "0.0"
        elif lineItems[0] == "4.2" :
            AngleOutDeg42.append(lineItems[4]) 
            Wavelength42 = lineItems[1] 
            angleInDeg42 = lineItems[2]
            if ( abs( float(angleInDeg42) ) < 1.e-5 ) :
                 angleInDeg42 = "0.0"
        elif lineItems[0] == "4.3" :
            AngleOutDeg43.append(lineItems[4]) 
            Wavelength43 = lineItems[1] 
            angleInDeg43 = lineItems[2]
            if ( abs( float(angleInDeg43) ) < 1.e-5 ) :
                 angleInDeg43 = "0.0"
        elif lineItems[0] == "4.4" :
            AngleOutDeg44.append(lineItems[4]) 
            Wavelength44 = lineItems[1] 
            angleInDeg44 = lineItems[2]
            if ( abs( float(angleInDeg44) ) < 1.e-5 ) :
                 angleInDeg44 = "0.0"
            
 #Test5
        elif lineItems[0] == "5.0a" :
            AngleOutDeg50a.append( lineItems[3] )
            Interference50a.append( lineItems[4] )
            Envelope50a.append( lineItems[5] )
            Wavelength50a=lineItems[1]
            angleInDeg50a = lineItems[2]
            if ( abs( float(angleInDeg50a) ) < 1.e-5 ) :
                 angleInDeg50a = "0.0"
        elif lineItems[0] == "5.0b" :
            AngleOutDeg50b.append( lineItems[3] )
            Interference50b.append( lineItems[4] )
            Envelope50b.append( lineItems[5] )
            Wavelength50b=lineItems[1]
            angleInDeg50b = lineItems[2]
            if ( abs( float(angleInDeg50b) ) < 1.e-5 ) :
                 angleInDeg50b = "0.0"
        elif lineItems[0] == "5.1" :
            AngleOutDeg51.append(lineItems[4]) 
            Wavelength51=lineItems[1]
            angleInDeg51 = lineItems[2]
            if ( abs( float(angleInDeg51) ) < 1.e-5 ) :
                 angleInDeg51 = "0.0"
        elif lineItems[0] == "5.2" :
            AngleOutDeg52.append(lineItems[4]) 
            Wavelength52=lineItems[1] 
            angleInDeg52 = lineItems[2]
            if ( abs( float(angleInDeg52) ) < 1.e-5 ) :
                 angleInDeg52 = "0.0"
        elif lineItems[0] == "5.3" :
            AngleOutDeg53.append(lineItems[4]) 
            Wavelength53=lineItems[1] 
            angleInDeg53 = lineItems[2]
            if ( abs( float(angleInDeg53) ) < 1.e-5 ) :
                 angleInDeg53 = "0.0"
        elif lineItems[0] == "5.4" :
            AngleOutDeg54.append(lineItems[4]) 
            Wavelength54=lineItems[1] 
            angleInDeg54 = lineItems[2]
            if ( abs( float(angleInDeg54) ) < 1.e-5 ) :
                 angleInDeg54 = "0.0"
 #Test 7
        elif lineItems[0] == "7.1" :
            AngleOutXDeg71.append(lineItems[4]) 
            WavelengthNM71.append(lineItems[1])
            angleInXDeg71 = lineItems[2]
            if ( abs( float(angleInXDeg71) ) < 1.e-5 ) :
                 angleInXDeg71 = "0.0"
        elif lineItems[0] == "7.2" :
            AngleOutXDeg72.append(lineItems[4]) 
            WavelengthNM72.append(lineItems[1])
            angleInXDeg72 = lineItems[2]
            if ( abs( float(angleInXDeg72) ) < 1.e-5 ) :
                 angleInXDeg72 = "0.0"
            
        
#Now we plot stuff using matplotlib and numpy

nbins=10000
width = 0.1
#Test1 and Test2.
fig1_2, ( (test1_1, test2_1a), (test2_2a, test2_1b), (test2_2b, test2_1c) ) = plt.subplots(3, 2, figsize=(8, 11.0))
fig1_2.suptitle('test_gradiant Results: Test1 and Test2')

print('Test11')
dotndarray = GenerateDataFloatArray(dot)
test1_1.hist(dotndarray, bins=50)
test1_1.set_title('Output Vector Generation Check')
test1_1.set_xlabel('Difference of Initail Normal to Final')
test1_1.set_ylabel('Number')

if ( len(probIntensity) != 0 ) :
    print('Test21a')
    intnsAngOut = GenerateDataFloatArray( intensityAngleOut )
    probIntns   = GenerateDataFloatArray( probIntensity )
    test2_1a.plot(intnsAngOut,probIntns)
    title21a = "Probability Intensity Full\n" + 'Log Scale ' + Wavelength21 + ' nm, ' + angleInDeg21 + ' deg'
    test2_1a.set_title(title21a )
    test2_1a.set_yscale("log")    #This line is difference from test21b
    test2_1a.set_xlabel('AngleOut (Deg)')
    test2_1a.set_ylabel('Relative Intensity')

    print('Test21b')
    test2_1b.plot(intnsAngOut,probIntns)
    title21b = "Probability Intensity, Full\n" + Wavelength21 + ' nm, ' + angleInDeg21 + ' deg'
    test2_1b.set_title(title21b )
    test2_1b.set_xlabel('AngleOut (Deg)')
    test2_1b.set_ylabel('Relative Intensity')

    print('Test21c')
    peakIndex=findPeakIndexSlidingWindow(probIntns,windowSize)
    #peakIndex=np.argmax(probIntns)
    peakAngle21c=intnsAngOut[peakIndex]

    peakAngleLo21c = peakAngle21c - width
    peakAngleHi21c = peakAngle21c + width

    print('Test 2.1 Probability Max at angleOut: ' + str(peakAngle21c) + ' deg')
    intnsAngOutInRange, probIntnsInRange = GenerateArrayPairInFirstRange( \
        intnsAngOut, probIntns, peakAngleLo21c, peakAngleHi21c)
    test2_1c.plot(intnsAngOutInRange,probIntnsInRange)
    title21c = "Probability Intensity Mode 2\n " + Wavelength21 + ' nm, ' + angleInDeg21 + ' deg'
    test2_1c.set_title(title21c)
    test2_1c.set_xlabel('AngleOut (Deg)')
    test2_1c.set_ylabel('Relative Intensity')

    print('Test22a')
    cmAngOut =  GenerateDataFloatArray( cumProbAngleOut )
    cmPrb =  GenerateDataFloatArray( cumProb)
    test2_2a.plot( cmAngOut,cmPrb)
    title22a = "Cumulative probability\n" +  Wavelength22 + ' nm, ' + angleInDeg22 + ' deg'
    test2_2a.set_title(title22a)
    test2_2a.set_xlabel('AngleOut (Deg)')
    test2_2a.set_ylabel('Probability')

    print('Test22b')
    cmAngOutInRange, cmPrbInRange = GenerateArrayPairInFirstRange( \
        cmAngOut, cmPrb, peakAngle21c, peakAngleHi21c)
    test2_2b.plot( cmAngOutInRange,cmPrbInRange)
    title22b = "Cumulative Probability Mode 2\n" +  Wavelength22 + ' nm, ' + angleInDeg22 + ' deg'
    test2_2b.set_title(title22b)
    test2_2b.set_xlabel('AngleOut (Deg)')
    test2_2b.set_ylabel('Probability')
else :
    print("No Test2. Table read from file")
    
print("Save Fig1_2")
sys.stdout.flush()
fig1_2.subplots_adjust( hspace=.4 )
fig1_2.savefig('tests1_2.dat.png')

#Test 3:

width = 0.05
#width = 10.0
fig3, ( (test3_1, test3_2), (test3_3, test3_4), (test3_5, test3_6) ) = \
      plt.subplots(3, 2,  figsize=(8, 11.0))
fig3.suptitle('test3: tests grating.findAngleOutRadFromTables method \n MidAngleIn(Index), Various Min,Mid,Max Wavelengths(index)\n Mode2')
fig3.subplots_adjust( hspace=.4 )

print('Test31;')
#Set Limits
peakLo, peakHi = GenerateHiLoOfArray(minIndexAngleOutDeg, width)
minAngOutInRange =  GenerateDataFloatArray(minIndexAngleOutDeg, True, \
                                           peakLo,peakHi)
y31, binEdges = np.histogram(minAngOutInRange, bins=nbins)
bincenters31 = 0.5 * (binEdges[1:] + binEdges[:-1])
test3_1.plot(bincenters31, y31, '-')
title31 = angleInDeg31+ ' deg(index) ' + Wavelength31 + 'nm(Index)'
test3_1.set_title(title31)
test3_1.set_xlabel('AngleOut (Deg)')
test3_1.set_ylabel('Number')

print('test32;')
peakLo, peakHi = GenerateHiLoOfArray(midIndexAngleOutDeg, width)
midAngOutInRange =  GenerateDataFloatArray(midIndexAngleOutDeg, True, \
                                           peakLo, peakHi)
y32, binEdges = np.histogram(midAngOutInRange, bins=nbins)
bincenters32 = 0.5 * (binEdges[1:] + binEdges[:-1])
test3_2.plot(bincenters32, y32, '-')
title32 = angleInDeg32+ ' deg(index) ' + Wavelength32 + 'nm(Index)'
test3_2.set_title(title32)
test3_2.set_xlabel('AngleOut (Deg)')
test3_2.set_ylabel('Number')

print('test33;')
peakLo, peakHi = GenerateHiLoOfArray(midPlus1IndexAngleOutDeg, width)
midP1AngOutInRange = GenerateDataFloatArray(midPlus1IndexAngleOutDeg, True,\
                                           peakLo, peakHi)
y33, binEdges = np.histogram(midP1AngOutInRange, bins=nbins)
bincenters33 = 0.5 * (binEdges[1:] + binEdges[:-1])
test3_3.plot(bincenters33, y33, '-')
title33 = angleInDeg33+ ' deg(index) ' + Wavelength33 + 'nm(Index)'
test3_3.set_title(title33)
test3_3.set_xlabel('AngleOut (Deg)')
test3_3.set_ylabel('Number')

print('test34;')
peakLo, peakHi = GenerateHiLoOfArray(midPlus2IndexAngleOutDeg, width)
midP2AngOutInRange =  GenerateDataFloatArray(midPlus2IndexAngleOutDeg, True, \
                                           peakLo, peakHi)

#test3_4.hist(midP2AngOutInRange, bins=1000)
y34, binEdges = np.histogram(midP2AngOutInRange, bins=nbins)
bincenters34 = 0.5 * (binEdges[1:] + binEdges[:-1])
test3_4.plot(bincenters34, y34, '-')
title34 = angleInDeg34+ ' deg(index) ' + Wavelength34 + 'nm(Index)'
test3_4.set_title(title34)
test3_4.set_xlabel('AngleOut (Deg)')
test3_4.set_ylabel('Number')

print('test35;')
peakLo, peakHi = GenerateHiLoOfArray(maxIndexAngleOutDeg, width)
maxAngOutInRange = GenerateDataFloatArray(maxIndexAngleOutDeg, True, \
                                           peakLo, peakHi)
y, binEdges = np.histogram(maxAngOutInRange, bins=nbins)
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
test3_5.plot(bincenters, y, '-')
title35 = angleInDeg35+ ' deg(index) ' + Wavelength35 + 'nm(Index)'
test3_5.set_title(title35)
test3_5.set_xlabel('AngleOut (Deg)')
test3_5.set_ylabel('Number')

print('test36;')
test3_6.plot(bincenters32, y32, '-')
test3_6.plot(bincenters33, y33, '-')
test3_6.plot(bincenters34, y34, '-')
legend32 = Wavelength32 + 'nm'
legend33 = Wavelength33 + 'nm'
legend34 = Wavelength34 + 'nm'
test3_6.legend( [legend32, legend33, legend34], \
                loc="upper left")
title36 =  angleInDeg31+ ' deg(index),Mid Wavelengths(Index)'
test3_6.set_title(title36)
test3_6.set_xlabel('AngleOut (Deg)')
test3_6.set_ylabel('Number')

print("Save Fig3")
sys.stdout.flush()
fig3.savefig('tests3.dat.png')


#Test4:
width = 0.05
fig4, ( (test4_1a, test4_1b), (test4_2, test4_3), (test4_4, test4_5) ) =  \
                                    plt.subplots(3, 2,  figsize=(8, 11.0))
title4 = 'Test4: Tests grating.calculateAngleOutRadFromTable method\n ' +  'AngleIn: ' + angleInDeg41 + 'deg,  Various Wavelengths, Mode 2'
fig4.suptitle(title4)
fig4.subplots_adjust( hspace=.4 )

print('test41;')
AngOutDeg41a =  GenerateDataFloatArray(AngleOutDeg41, \
                                         True, 0.0, 30.0)
y41a, binEdges = np.histogram(AngOutDeg41a, bins=nbins)
bincenters41a = 0.5 * (binEdges[1:] + binEdges[:-1])
test4_1a.plot(bincenters41a, y41a, '-')
title41a = 'Lambda: '+ Wavelength41 +' nm'
test4_1a.set_title(title41a)
test4_1a.set_xlabel('AngleOut (Deg)')
test4_1a.set_ylabel('Number')
#test4_1a.set_yscale("log")


peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg41, width)
AngOutDeg41b =  GenerateDataFloatArray(AngleOutDeg41, \
                                         True, peakLo, peakHi)
y41b, binEdges = np.histogram(AngOutDeg41b, bins=nbins)
bincenters41b = 0.5 * (binEdges[1:] + binEdges[:-1])
test4_1b.plot(bincenters41b, y41b, '-')
title41b = 'Lambda: '+ Wavelength41 +' nm'
test4_1b.set_title(title41b)
test4_1b.set_xlabel('AngleOut (Deg)')
test4_1b.set_ylabel('Number')


#Test42
peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg42, width)
AngOutDeg42 =  GenerateDataFloatArray(AngleOutDeg42, \
                                         True, peakLo, peakHi)
y42, binEdges = np.histogram(AngOutDeg42, bins=nbins)
bincenters42 = 0.5 * (binEdges[1:] + binEdges[:-1])
test4_2.plot(bincenters42, y42, '-')
title42 = 'Lambda: '+ Wavelength42 +' nm'
test4_2.set_title(title42)
test4_2.set_xlabel('AngleOut (Deg)')
test4_2.set_ylabel('Number')

print('test43;')
peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg43, width)
AngOutDeg43 =  GenerateDataFloatArray(AngleOutDeg43, \
                                         True, peakLo, peakHi)
y43, binEdges = np.histogram(AngOutDeg43, bins=nbins)
bincenters43 = 0.5 * (binEdges[1:] + binEdges[:-1])
test4_3.plot(bincenters43, y43, '-')
title43 = 'Lambda: '+ Wavelength43 +' nm'
test4_3.set_title(title43)
test4_3.set_xlabel('AngleOut (Deg)')
test4_3.set_ylabel('Number')

print('test44;')
peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg44, width)
AngOutDeg44 =  GenerateDataFloatArray(AngleOutDeg44, \
                                         True, peakLo, peakHi)
y44, binEdges = np.histogram(AngOutDeg44, bins=nbins)
bincenters44 = 0.5 * (binEdges[1:] + binEdges[:-1])
test4_4.plot(bincenters44, y44, '-')
title44 = 'Lambda: '+ Wavelength44 +' nm'
test4_4.set_title(title44)
test4_4.set_xlabel('AngleOut (Deg)')
test4_4.set_ylabel('Number')

print('test45;')
test4_5.plot(bincenters41b, y41b, '-')
test4_5.plot(bincenters42, y42, '-')
test4_5.plot(bincenters43, y43, '-')
test4_5.plot(bincenters44, y44, '-')
legend41b = Wavelength41 +  'nm'
legend42 = Wavelength42 +  'nm'
legend43 = Wavelength43 +  'nm'
legend44 = Wavelength44 +  'nm'
test4_5.legend( [legend41b, legend42, legend43, legend44], \
                loc="upper left")
title45 = Wavelength41 + ', ' + Wavelength42 + ', ' +  Wavelength43 + ', ' \
          +  Wavelength44 + ' (nm)'
test4_5.set_title(title45)
test4_5.set_xlabel('AngleOut (Deg)')
test4_5.set_ylabel('Number')

print("Save Fig4")
sys.stdout.flush()
fig4.savefig('tests4.dat.png')

#Test5:
width = 0.05
fig5a, ( (test50aI, test50bI) , (test50aE, test50bE ), (test50a, test50b ))  = \
                                    plt.subplots(3, 2,  figsize=(8, 11.0))
title50 = 'Test5.0: Tests calculateSingleSlitEnvelope() and \n'+ 'calculateInterferenceIntensity()\n' + 'at ' + Wavelength50a + ' nm.'
fig5a.suptitle(title50)
fig5a.subplots_adjust( hspace=.4 )

print('test5_0a;')

AngOutDeg50a = GenerateDataFloatArray( AngleOutDeg50a , True, -2.0, 30.0)
Intrf50a   = GenerateDataFloatArray( Interference50a )
test50aI.plot(AngOutDeg50a,Intrf50a)
title50aI = "Interference term\n" + Wavelength50a + ' nm, ' + angleInDeg50a + ' deg'
test50aI.set_title(title50aI )
test50aI.set_xlabel('AngleOut (Deg)')
test50aI.set_ylabel('Relative Interference')

Env50a = GenerateDataFloatArray( Envelope50a , True, -2.0, 30.0)
test50aE.plot(AngOutDeg50a,Env50a)
title50aE = "Envelope term\n" + Wavelength50a + ' nm, ' + angleInDeg50a + ' deg'
test50aE.set_title(title50aE )
test50aE.set_xlabel('AngleOut (Deg)')
test50aE.set_ylabel('Relative Envelope')

Intns50a = np.multiply( Env50a, Intrf50a)
test50a.plot(AngOutDeg50a, Intns50a)
title50a = "E*I= Intensity\n" + Wavelength50a + ' nm, ' + angleInDeg50a + ' deg'
test50a.set_title(title50a )
test50a.set_xlabel('AngleOut (Deg)')
test50a.set_ylabel('Relative Intensity')

print('test5_0b;')

AngOutDeg50b = GenerateDataFloatArray( AngleOutDeg50b , True, -2.0, 30.0)
Intrf50b   = GenerateDataFloatArray( Interference50b )
test50bI.plot(AngOutDeg50b,Intrf50b)
title50bI = "Interference term\n" + Wavelength50b + ' nm, ' + angleInDeg50b + ' deg'
test50bI.set_title(title50bI )
test50bI.set_xlabel('AngleOut (Deg)')
test50bI.set_ylabel('Relative Interference')

Env50b = GenerateDataFloatArray( Envelope50b , True, -2.0, 30.0)
test50bE.plot(AngOutDeg50b,Env50b)
title50bE = "Envelope term\n" + Wavelength50b + ' nm, ' + angleInDeg50b + ' deg'
test50bE.set_title(title50bE )
test50bE.set_xlabel('AngleOut (Deg)')
test50bE.set_ylabel('Relative Envelope')

Intns50b = np.multiply( Env50b, Intrf50b)
test50b.plot(AngOutDeg50b, Intns50b)
title50b = "E*I= Intensity\n" + Wavelength50b + ' nm, ' + angleInDeg50b + ' deg'
test50b.set_title(title50b )
test50b.set_xlabel('AngleOut (Deg)')
test50b.set_ylabel('Relative Intensity')

print("Save Fig5a")
sys.stdout.flush()
fig5a.savefig('tests5a.dat.png')



fig5b, ( (test5_1a, test5_1b), (test5_2a, test5_2b), (test5_3a, test5_3b), (test5_4a, test5_4b) ) = \
                                    plt.subplots(4, 2,  figsize=(8, 11.0))
title5 = 'Test5b: Tests grating.calculateAngleOutRadFromTable method \n' + 'Various AngleIn and Lambda at ' + Wavelength51 + ' nm'
fig5b.suptitle(title5)
fig5b.subplots_adjust( hspace=.4 )


print('test51a;')
AngOutMidDeg51a =  GenerateDataFloatArray(AngleOutDeg51, True, -2.0, 30.0)
#AngOutMidDeg51a =  GenerateDataFloatArray(AngleOutDeg51, True, 13.2, 13.3)
y, binEdges = np.histogram(AngOutMidDeg51a, bins=nbins)
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_1a.plot(bincenters, y, '-')
title51a = angleInDeg51+ ' deg(index) ' + Wavelength51, 'nm'
test5_1a.set_title(title51a)
test5_1a.set_xlabel('AngleOut (Deg)')
test5_1a.set_ylabel('Number')

peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg51, width)
AngOutMidDeg51b =  GenerateDataFloatArray(AngleOutDeg51, \
                                       True, peakLo, peakHi)
y51b, binEdges = np.histogram(AngOutMidDeg51b, bins=nbins)
bincenters51b = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_1b.plot(bincenters51b, y51b, '-')
title51b = angleInDeg51+ ' deg(index) ' + Wavelength51, 'nm'
test5_1b.set_title(title51b)
test5_1b.set_xlabel('AngleOut (Deg)')
test5_1b.set_ylabel('Number')

print('test52;')
AngOutMidDeg52a =  GenerateDataFloatArray(AngleOutDeg52, True, -2.0, 30.0)
#AngOutMidDeg52a =  GenerateDataFloatArray(AngleOutDeg52, True, 13.2, 13.3)
y, binEdges = np.histogram(AngOutMidDeg52a, bins=nbins)
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_2a.plot(bincenters, y, '-')
title52a = angleInDeg52+ ' deg(index) ' + Wavelength52, 'nm'
test5_2a.set_title(title52a)
test5_2a.set_xlabel('AngleOut (Deg)')
test5_2a.set_ylabel('Number')

peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg52, width)
AngOutMidDeg52b =  GenerateDataFloatArray(AngleOutDeg52, \
                                       True, peakLo, peakHi)
y52b, binEdges = np.histogram(AngOutMidDeg52b, bins=nbins)
bincenters52b = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_2b.plot(bincenters52b, y52b, '-')
title52b =  angleInDeg52+ ' deg(index) ' + Wavelength52, 'nm'
test5_2b.set_title(title52b)
test5_2b.set_xlabel('AngleOut (Deg)')
test5_2b.set_ylabel('Number')

print('test53;')
AngOutMidDeg53a =  GenerateDataFloatArray(AngleOutDeg53, True, -2.0, 30.0)
#AngOutMidDeg53a =  GenerateDataFloatArray(AngleOutDeg53, True, 13.2, 13.3)
y, binEdges = np.histogram(AngOutMidDeg53a, bins=nbins)
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_3a.plot(bincenters, y, '-')
title53a = angleInDeg53+ ' deg(index) ' + Wavelength53, 'nm'
test5_3a.set_title(title53a)
test5_3a.set_xlabel('AngleOut (Deg)')
test5_3a.set_ylabel('Number')

peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg53, width)
AngOutMidDeg53b =  GenerateDataFloatArray(AngleOutDeg53, \
                                       True, peakLo, peakHi)
y53b, binEdges = np.histogram(AngOutMidDeg53b, bins=nbins)
bincenters53b = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_3b.plot(bincenters53b, y53b, '-')
title53b =  angleInDeg53+ ' deg(index) ' + Wavelength53, 'nm'
test5_3b.set_title(title53b)
test5_3b.set_xlabel('AngleOut (Deg)')
test5_3b.set_ylabel('Number')

print('test54;')
AngOutMidDeg54a =  GenerateDataFloatArray(AngleOutDeg54, True, -2.0, 30.0)
#AngOutMidDeg54a =  GenerateDataFloatArray(AngleOutDeg54, True, 13.2, 13.3)
y, binEdges = np.histogram(AngOutMidDeg54a, bins=nbins)
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_4a.plot(bincenters, y, '-')
title54a = angleInDeg54+ ' deg(index) ' + Wavelength54, 'nm'
test5_4a.set_title(title54a)
test5_4a.set_xlabel('AngleOut (Deg)')
test5_4a.set_ylabel('Number')

peakLo, peakHi = GenerateHiLoOfArray(AngleOutDeg54, width)
AngOutMidDeg54b =  GenerateDataFloatArray(AngleOutDeg54, \
                                       True, peakLo, peakHi)
y54b, binEdges = np.histogram(AngOutMidDeg54b, bins=nbins)
bincenters54b = 0.5 * (binEdges[1:] + binEdges[:-1])
test5_4b.plot(bincenters54b, y54b, '-')
title54b =  angleInDeg54+ ' deg(index) ' + Wavelength54, 'nm'
test5_4b.set_title(title54b)
test5_4b.set_xlabel('AngleOut (Deg)')
test5_4b.set_ylabel('Number')


print("Save Fig5b")
sys.stdout.flush()
fig5b.savefig('tests5b.dat.png')

#SummryTest3457:
fig3457, ( (test3Sum, test4Sum), (test5Sum, test7Sum) ) = \
                                    plt.subplots(2, 2,  figsize=(8, 11.0))
title3457 = 'Summary test3,4,5,7 Summary of tests'
fig3457.suptitle(title3457)
fig3457.subplots_adjust( hspace=.4 )

test3Sum.plot(bincenters32, y32, '-')
test3Sum.plot(bincenters33, y33, '-')
test3Sum.plot(bincenters34, y34, '-')
legend32 = Wavelength32 + 'nm'
legend33 = Wavelength33 + 'nm'
legend34 = Wavelength34 + 'nm'
test3Sum.legend( [legend32, legend33, legend34], \
                loc="upper left")
title3Sum =  angleInDeg31+ ' deg(index),Mid Wavelengths(Index)'
test3Sum.set_title(title3Sum)
test3Sum.set_xlabel('AngleOut(Deg)')
test3Sum.set_ylabel('Number')

test4Sum.plot(bincenters41b, y41b, '-')
test4Sum.plot(bincenters42, y42, '-')
test4Sum.plot(bincenters43, y43, '-')
test4Sum.plot(bincenters44, y44, '-')
legend41b = Wavelength41 +  'nm'
legend42 = Wavelength42 +  'nm'
legend43 = Wavelength43 +  'nm'
legend44 = Wavelength44 +  'nm'
test4Sum.legend( [legend41b, legend42, legend43, legend44], \
                loc="upper left")
title4Sum = Wavelength41 + ', ' + Wavelength42 + ', ' + Wavelength43 + ', ' \
          +  Wavelength44 + ' (nm)'
test4Sum.set_title(title4Sum)
test4Sum.set_xlabel('AngleOut (Deg)')
test4Sum.set_ylabel('Number')

test5Sum.plot(bincenters51b, y51b, '-')
test5Sum.plot(bincenters52b, y52b, '-')
test5Sum.plot(bincenters53b, y53b, '-')
test5Sum.plot(bincenters54b, y54b, '-')
legend51b =  "%.2g" % float(angleInDeg51) +  'deg'
legend52b =  "%.2g" % float(angleInDeg52) +  'deg'
legend53b =  "%.2g" % float(angleInDeg53) +  'deg'
legend54b =  "%.2g" % float(angleInDeg54) +  'deg'
test5Sum.legend( [ legend51b, legend52b, legend53b, legend54b ], \
                loc="upper center")
title5Sum = Wavelength51 + 'nm, (' + "%.2g" % float(angleInDeg51) + '. ' + \
          "%.2g" % float(angleInDeg52) + ', ' + \
          "%.2g" % float(angleInDeg53) + ', ' + \
          "%.2g" % float(angleInDeg54) + ' deg)'
test5Sum.set_title(title5Sum)
test5Sum.set_xlabel('AngleOut (Deg)')
test5Sum.set_ylabel('Number')
print("Save Fig3457")
sys.stdout.flush()
fig3457.savefig('tests3457.dat.png')


#Test7:
nbins = 100000
fig7, ( (test7_1a, test7_2a), (test7_1b, test7_2b), (test7_1c, test7_2d) ) =  \
                                    plt.subplots(3, 2,  figsize=(8, 11.0))
title7 = 'Test7: Tests grating.Diffract method at angleIn of  ' + angleInXDeg71 + 'deg,\n ' + ' Flat spectrum and stepped wavelengths'
fig7.suptitle(title7)
fig7.subplots_adjust( hspace=.7 )

print('test71;')
AngOutXDeg71 =  GenerateDataFloatArray(AngleOutXDeg71)
y71a, binEdges = np.histogram(AngOutXDeg71, bins=nbins)
bincenters71a = 0.5 * (binEdges[1:] + binEdges[:-1])
test7_1a.plot(bincenters71a, y71a, '-')
title71a = 'Flat Spectrum: Diffract Test\n' + 'AngleIN: '+ angleInXDeg71 +' deg'
test7_1a.set_title(title71a)
test7_1a.set_xlabel('AngleOut (Deg)')
test7_1a.set_ylabel('Number')

print('test71b;')
Wavelngth71 = GenerateDataFloatArray(WavelengthNM71)
y71b, binEdges = np.histogram(Wavelngth71, bins=nbins)
bincenters71b = 0.5 * (binEdges[1:] + binEdges[:-1])
test7_1b.plot(bincenters71b, y71b, '-')
title71b = 'Flat Spectrum: Diffract Test\n' + 'AngleIN: '+ angleInXDeg71 +' deg'
test7_1b.set_title(title71b)
test7_1b.set_xlabel('Wavelength (nm)')
test7_1b.set_ylabel('Number')


print('test72;')
AngOutXDeg72a =  GenerateDataFloatArray(AngleOutXDeg72, True, -2.0, 30.0)
y72a, binEdges = np.histogram(AngOutXDeg72a, bins=nbins)
bincenters72a = 0.5 * (binEdges[1:] + binEdges[:-1])
test7_2a.plot(bincenters72a, y72a, '-')
title72a = 'Lines with  5 nm Wavelength seperation\n' + 'AngleIN: '+ angleInXDeg71 +' deg'
test7_2a.set_title(title72a)
test7_2a.set_xlabel('AngleOut (Deg)')
test7_2a.set_ylabel('Number')

AngOutXDeg72b =  GenerateDataFloatArray(AngleOutXDeg72, True, 15, 18)
y72b, binEdges = np.histogram(AngOutXDeg72b, bins=nbins)
bincenters72b = 0.5 * (binEdges[1:] + binEdges[:-1])
test7_2b.plot(bincenters72b, y72b, '-')
title72b = 'Lines with  5 nm Wavelength seperation\n' + 'Mode 2, AngleIN: '+ angleInXDeg71 +' deg'
test7_2b.set_title(title72b)
test7_2b.set_xlabel('AngleOut (Deg)')
test7_2b.set_ylabel('Number')

print("Save Fig7")
sys.stdout.flush()
fig7.savefig('tests7.dat.png')

