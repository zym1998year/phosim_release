#!/usr/bin/env python

"""
  @package phosim
  @file ZMXtoPhosim.py
  @brief python script to convert a ZEMAX file to PHOSIM optics_*.txt(sort of)
 
  @brief Created by:
  @author En-Hsin Peng (Purdue)
  @Updated Glenn Sembrosk(GHS)i 2023-02-14 (GHS)
  @Updated Glenn Sembroski(GHS) 2025-02-25 (GHS)
       (make use of input.lis file optional)
  @Updated Cory Hoo 2025-05-23
        (converts .zmx file from unicode to ascii and generates simple
        focalplanelayout.txt and segmentation file) (suggested use:
        copy python script to new telescope directly and run in the new
        directory) 


  @warning This code is not fully validated
        and not ready for full release.  Please
        treat results with great caution.

Usage: python ZMXtoPhosim.py zemax.zmx input.lis(optional)

In general the ZEMAX.ZMX files as used by the ZEMAX programs will be in
unicode(but see update 2025-05-23 above).
This python script expects an input zemax.zmx file to be in ascii.
Thus you need to convert the ZEMAX.ZMX unicode file to a zemax.zmx ascii
file. In linux that looks like:
    $> cat ZEMAX.ZMX | iconv -f utf-16 > zemacs.zmx
    Note: we could probably add this to this script,(but see update
    2025-05-23 above).

    NOTE:(GHS)  My and En-Hsin's interpretation of how to translate the
    zemac.zmx file is basically an educated guess. As far as I can tell
    there is no non-proprietary manual that describes the various options
    in the zemacx.zmx file. I (and En-Hsin before me) are making our
    best guess.
    NOTE:(GHS,CH) Making requirement of input.lis file optional requires
    determining the information of what was in the input.lis file from the
    input ZEMAX.zmx file. Esentially, if no input.lis file is specified on
    the command line a new input.lis file is generated from the zemax.zmx
    file itself. That new input.lis file is then used to generate the phosim
    optics_*.txt file.

    input.lis is an optional, user provided file or it can be a file created
    by this prograsm. It contains a line for each surface to be included in
    the phosim optics.*.txt file with the following format:
    input.lis:
       Column 1: Zemax surface number
       Column 2: Surface name (M1, L1, etc) 
       Column 3: Surface type (mirror, lens, det)
                                 (in future? also filter and  grating) 
       Column 4: Coating file
       Column 5: Medium file
       Column 6: Material file

       The input.lis  file can be provided by the user.  To create this file
       the user needs to look through the zemax.zmx file looking for the
       'SURF' entries. The user then picks out those surfaces that are of
       interest to include in the optics_*.txt file and includes a line for
       each such surface in the input.lis file as described above. The
       values in the input.lis file will be used in the optics_*.txt file
       that this program generates unless they are found to be superseded in
       the zemax.zmx file.

       If an input.lis file is supplied on the command line:
           All supplied values (Columns 2-6) as specified in the file will
           be used in creating the optics_*.txt files.

       If this program creates an input.lis file from values found in the
       zemax.zmx file, then the various columns will be:
       Column2:  M1, L1, L1E, Det etc. as determined
       Column3:  mirror, lens or det as determined
       Column4:  none  (always)
       Column5:  air or for lend first surface, the GLAS command material
       Column6:  none

       If the user finds the values in the generated input.lis file need
       to be edited( most likely sittuation), ZMXtoPHOSIM.py should be rerun
       with the edited input.lis file specified on the command line. 

       Example input.lis file (From LSST) e.g.,
       23   M1    mirror  m1_protAl_Ideal.txt  air  none
       26   M2    mirror  m2_protAl_Ideal.txt  air  none
       29   M3    mirror  m3_protAl_Ideal.txt  air  none
       30   none  none    none                 air  none
       31   L1    lens    lenses.txt           silica_dispersion.txt  none
       32   L1E   lens    lenses.txt           air  none
       33   L2    lens    lenses.txt           silica_dispersion.txt  none
       34   L2E   lens    lenses.txt           air  none
       35   F     filter  filter_x.txt         silica_dispersion.txt  none
       36   FE    filter  none                 air  none
       37   L3    lens    lenses.txt           silica_dispersion.txt  none
       38   L3E   lens    lenses.txt           air  none
       39   D     det     detectorar.txt       air  none

   Note:1 This script assumes that all the values in the zemax.zmx file
        which are measurments are in mm  unless a 'UNIT' option is specified
        in the zemax.zmx file. In that case, the values found in zemax.zmx
        file will be converted to mm as necessary . PHOSIM requires units to
        be in mm.
   Note:2 This program will overwrite any existing optics_*.txt in the
        current folder. 

"""

import sys, os, subprocess, shutil

class Surface(object):            #Values for a line in the optics_*.txt file
    def __init__(self,curv,disz,z,rout,rin,coni,an,med):
        self.curv=curv   #Spherical curvature of the surface 0.0 defined as
                         #flat in PHOSIM,otherwise as normal!
        self.disz=disz   #Offset to the next surface
        self.z=z         #Absolute z of next surface (I think)
        self.rout=rout   #Outer radius of surface
        self.rin=rin     #Inner radius of surface (hole)
        self.coni=coni   #Conical value for the surface(1.0 is parabola)
        self.an=an       #polynomial constants starting at 2nd degree
        self.medium=med  #Medium. Type of material folloing the surface

#This function Finds the next surface (surf) as specified in the input.lis
#file.
#Collect all the values for this surface including its relative offset
#from the previous input.lis specified surface.
def findSurface(zmx,surf,surf0,flt,med,zprev=0.0):
    
    curv=0.0        #Curvature of surface in mm. Note PHOSIM wants curvature
                    #set to 0 if it is a flat surface(curvature=infinity)
    rout=0.0        #Outer radius
    rin=0.0         #Inner radius
    coni=0.0        #Conic constant (==1 for parabolic)
    an=[0.0 for i in range(17)] #polynomial constants starting at 2nd degree

    medium = med   #Default medium type


    ################################################
    #The following is mainly about filling in the absolute z values of the
    #surfaces.

    disz=[0.0 for i in range(100)]  #relative offsets between surfaces.
    pzup=[[-1.0,-1.0,-1.0] for i in range(100)]  #Not sure what this is(GHS)
    readCurv=False

    #read DISZ, update THIC, CRVT
    #Restart at the beginning of the file. This is not particularly
    #efficient but files are small so this will take little time.
    #Note: I think this fails for flt>0. i.e. more than one versions of an
    #instrument
    for line in open(zmx).readlines():  
        if 'SURF ' in line:                #Beginning of the next surface
            s=int(float(line.split()[1]))  #Pickup zemacs surface number
        elif 'DISZ' in line:
            #collect offset to the next surface (following this one), 
            if line.split()[1]!='INFINITY': #Skip first surface,defaults
                                            #to 0.0
                disz[s]=float(line.split()[1]) #Save offset to next surface,
                                               #indexing by surface number
        elif 'THIC' in line:       #No idea what THIC is(GHS)
            l=line.split()
            if float(l[2])-1==flt:
                disz[int(float(l[1]))]=float(l[3])
        elif 'CRVT' in line:   #This must be an inverse radius of  curvature
            l=line.split()     #This is part of the multi telescope version
                               #analysis. I have no idea what its doing(GHS).
            if float(l[2])-1==flt and float(l[1])==surf:
                curv=-(1/float(l[3]))
                readCurv=True

        elif 'PZUP' in line:   #Not sure what this is.(GHS)
            l=line.split()
            pzup[s][0]=float(l[1])
            pzup[s][1]=float(l[2])
            pzup[s][2]=float(l[3])

    #PZUP
    #Nor sure how this is used. Anyway, pzup was defined as originally
    #loaded with all "-1.0" So if no PZUP were found above the following
    #is skipped and disz is left untouched.
    for i in range(surf):#Not sure whats this is about. (GHS)
        if pzup[i][0]>0:
            disz[i]=disz[int(pzup[i][0])]*pzup[i][1]+pzup[i][2]


    #Now we are going to use the disz to determine the absolute z location
    #of this surface.
    z=0.0
    if surf>surf0:  #Surf0 is first surface listed in the input.lis file.
                    #Going to let its z location be 0.
        for i in range(surf): #Note i's  this range is from 0 to surf-1. Thus
                              # at the end of this loop, z is the absolute
                              #position of the surface(With + z going up.
            if i==surf0:      #Make z relative to surf[0]
                z0=z          #Grab z of surf0 to use as a reference
            z-=disz[i]   #For PHOSIM invert direction. The sky is at +z.
                 
        z-=z0            #Remove base
    #So we end up with z as the distance between surf0 and surf
    #######################################    

    #We will us this value of z to determine the offset between this surface
    #and the previous input.lis specified surface as needed for the PHOSIM
    #optics_*.txt file.

#Read spherical parameters and size of this surface (and any inner hole in
#them)
    found=False
    #Reread the zemax.zmx file from the beginning. Not the most efficient but
    #we really don't care.
    for line in open(zmx).readlines():
        if 'SURF '+str(surf) in line:   #look for surface 'surf' in the
                                        #zemacs.zmx file.
            found=True
            continue
        if found:                       #We are in the surf block we want.
            if 'CURV' in line and readCurv==False: 
                if float(line.split()[1])!=0:      #Pick up spherical radius of
                    curv=-1/float(line.split()[1]) #curvature. Zemax CURV
                                                   #uses the negative inverse
                                                   #of the PHOSIM convention.
  
            #In the following we determine the inner and outer radius of the
            #circular surface. My guess is that the DIAM value is the active
            #outer radius of the surface. CLAP (CLear APerture?)or FLAP(?)
            #option values are the the actual radii of the surface (Inner and
            #Outer) and we only need the inner (may not be the active inner
            #radius but its what we have). The OBSC (OBSCuration?) seems also
            #to specify an inner radius and it may actually be the active
            #inner radius. DIAM (DIAMeter? maybe not since the value is for a
            #radius) options gives the outer active radius limit of the
            #surface. We always seem to have a DIAM option for a surface and
            # a CLAP,FLAP or OBCS option if we have an inner hole. Thus, for
            #the PHOSIM optics_*.txt file we want to use the DIAM radius for
            #the outper radius and the CLAP, FLAP or OBSC values for the
            #inner radius (if we need one).

            elif 'DIAM' in line:
                rout=float(line.split()[1])   #Actual outer radius of surface
                                              #Used if no CLAP of FLAP option
                rin=0.0

            elif 'CLAP' in line or 'FLAP' in  line : 
                rin=float(line.split()[1])    #Inner radius (of the hole) in
                                              #the surface

            elif 'OBSC' in line: 
 
                rin=float(line.split()[2])    #Inner radius (of the hole) in
                                              #the surface

            #Conical value of the surface (EX: for parabola == 1.0)
            elif 'CONI' in line:
                coni=float(line.split()[1])
            
            #polynomial constants starting at 2nd degree of surface
            elif 'PARM' in line:  #parm 1: a2, parm 2: a4
                lstr=line.split()
                an[int(float(lstr[1]))*2-1]=float(lstr[2])/1e3
            elif 'GLAS' in line:  #Get the glass type
                lstr=line.split()
                #If this is a MIRROR should the medium type be set to 'air' or
                #'none' or left to the input.lis file value ? Lets do that!
                if lstr[1] != 'MIRROR' :
                    medium = lstr[1]+'.txt'
            elif 'SURF' in line or 'CONF' in line:
                    break   #Quit when we finish reading this  surfaces
                            #options.

    #Load up the surface structure
    surface=Surface(curv,(z-zprev),z,rout,rin,coni,an,medium)
    return surface


def printOptics(output,surface,name,typ,coating,medium,material,flt):
    out=open(output,'a')
    #Setup default size of detector.
    if typ=='det' and surface.rout==0:
        surface.rout=400.0/unitsConv
#    print(' rout,rin:'+str(surface.rout) + ' ' + str(surface.rin))
    out.write('%-7s %7s %8.1f %10.4f %10.4f %6.1f %9.6f '  % (name, typ,
               surface.curv*unitsConv, surface.disz*unitsConv,
               surface.rout*unitsConv, surface.rin*unitsConv, surface.coni))
    for i in range(15):
        out.write('%13.6e ' % (-surface.an[i+2]))
    if coating!='none' and typ=='filter':
        coating='filter_'+str(flt)+'.txt'
    out.write('%s %s %s\n' % (coating,medium,material))
    out.close()


def printInputList(inputList, ZMXSurface, name, typ, coating, medium, material):
    out=open(inputList, 'a' )
    out.write('%-7s %7s  %7s %12s %12s %12s \n'  % (str(ZMXSurface), name,
              typ, coating, medium, material ))
    out.close()

def is_utf16(filename):
    with open(filename, 'rb') as f:
        header = f.read(2)
        return header in [b'\xff\xfe', b'\xfe\xff']


#######################
#
#Start of Main code#
#
#######################

if len(sys.argv) != 2:
    print("Usage: python ZMXtoPhosim.py <zmxfile>")
    sys.exit(1)

# Initial input file from user
input_file = sys.argv[1]
if not os.path.isfile(input_file):
    print(f"Error: File '{input_file}' not found.")
    sys.exit(1)

# Check encoding and convert if needed
if is_utf16(input_file):
    base, ext = os.path.splitext(input_file)
    ascii_file = f"{base}_ascii{ext}"
    print(f"Converting {input_file} from UTF-16 to ASCII format...")
    try:
        with open(ascii_file, "w") as outfile:
            subprocess.run(
                ["iconv", "-f", "utf-16", "-t", "ascii//TRANSLIT", input_file],
                stdout=outfile,
                check=True
            )
        print(f"Converted file saved as: {ascii_file}")
        zmxFile = ascii_file
    except subprocess.CalledProcessError as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)
else:
    print(f"{input_file} is already ASCII.")
    # Truncate at the first '.' and construct the new filename
    first_dot_index = input_file.find('.')
    if first_dot_index != -1:
        base = input_file[:first_dot_index]
    else:
        base = os.path.splitext(input_file)[0]
    ascii_file = f"{base}_ASCII.txt"
    shutil.copy(input_file, ascii_file)
    print(f"Copied original ASCII file to: {ascii_file}")
    zmxFile = ascii_file

"""It is possible that the zemax.zmx file will have specifications for a
   number of versions of the instrument. This may correspond to the different
   optics_*.txt files (different filters for example).
"""
#Count the number of instrument versions. (All the surface commands for a
#particular version end with the MNUM option in the zemac.zmx file.

fltNum=0.0
print('count number of MNUM instances in zemax file\n')
with open(zmxFile,"r") as fp:
    lines=fp.read()    # The .read() reads the whole file into a single string
    fltNum=lines.count("MNUM")    #Counts number of telescope versions
print(' Number of telescope versions found (MNUM instances): ' + str(fltNum))

#Find if we have a :"UNIT" option specified in the zemax.zmx file
# and specfiy the conversion factor to mm
units="mm"
unitConv=1.0
for line in open(zmxFile).readlines():
    if 'UNIT' in line:
        print('UNIT line: ' + line)
        units=line.split()[1]
        print('units: ' + units)
        if "IN" in units:
            print('Set conversion to mm from IN')
            unitsConv = 25.4
            break
        elif units == "MM":
            unitsConv=1.0   #mm
            break
        else :
            print('ZMXtoPHOSIM.py: Unknown units type specified in input .zmx file.\n')
            exit(1)
            
    else:
        #Assume units in zemax.zmx file in mm.
        unitsConv=1.0

print('Measurment units conversion factor: '+str(unitsConv)+'\n')

###################################################################
#Check to see if input.lis file was specfied in command line. If not
#determine the information directly from the zemax.zmx file (at least as
# well as we can). 
###################################################################
inputLisExists=False
if ( ( len(sys.argv) - 1) > 1) :
    inputLisExists=True
    inputList=sys.argv[2]

#Setup arrays for the input.lis file

zmxSurface=[]  #ZEMAX file number id for each surface
name=[]        #surface names
typ=[]         #Surface types
coating=[]     #Coating file names
medium=[]      #Medium file names
material=[]    #material file names

#If input.lis file was specfied, read in the input.lis values. Only do this
#once for all telescope versions.
if (inputLisExists) :
    for line in open(inputList).readlines():
        l=line.split()
        zmxSurface.append(int(float(l[0])))
        name.append(l[1])
        typ.append(l[2])
        coating.append(l[3])
        medium.append(l[4])
        material.append(l[5])
 
    ####################################################################
    #If no input.lis file is specified on the command line, generate one. Its
    #determine by values from  the zemax.zmx file. Again we only do this once.
    #We use the first telescope version to get this info. We thus assume here
    #that all telescope versions will have the same surfaces of interest.
    # (May need to fix this later).

    """ iterate through the zemax file starting at the beginning. We are
      searching for a line that begins with 'SURF'.  We quit if we find
      a line that begins with 'MNUM'
      Surface search algorithm: As we process the zemax.zmx file line by
      line, use the finding of a 'SURF' command (start of a new surface) to
      indicate we need to first process the previous surface to see if it
      should be added to the input.lis file we are generating.
      We only add 4 types of surfaces to the input .lis file:

      1: Previous surface was found to be a mirror: argument to 'GLAS'
      command = 'MIRROR'.
       goodMirror==True
      
      2: Previous surface was found to be front surface of a lens: argument
      to GLAS command is not MIRROR , but instead glass type(medium)  of lens.
      goodFrontLens == True
      
      3: Previous surface was found to be the back surface of a lens. First
      surface following a lens front surface and no GLAS command
      goodBackLens == True

      4: Previous surface is at the focal plane and is thus a detector
      surface. This will be last surface specified before MNUM command.
    """ 
else :
    goodMirror = False
    goodFrontLens = False
    goodBackLens = False
    mirrorIndex=0
    lensIndex=0
    surfID=0
    surfName='M'
    surfTyp = 'mirror'
    surfCoating = 'none'
    surfMaterial = 'air'
    surfMedium = 'none'
    
    for line in open(zmxFile).readlines() :  
#        print(' line: ' + line + ' \n')
        if 'SURF' in line:
            surfID = int(float(line.split()[1]))  #Pickup zemacs surface number
#            print('SURF in line,goodMirror,goodFrontLens,goodBackLens: ' +
#                        str(goodMirror) + ' ' + str(goodFrontLens) + ' ' +
#                        str(goodBackLens)) 
            #Beginning of the next surface
            #thus also end of previous
            #If previous surface was good keep its info.
            #ignore surface 0 (always seems to be at infinity)
            if (goodMirror or goodFrontLens ) and surfID > 0 :
                zmxSurface.append(surfID-1)
                name.append(surfName)
#                print(' surfName: ' + surfName )
                typ.append(surfTyp)
                coating.append(surfCoating)
                medium.append(surfMedium)
                material.append(surfMaterial)

                #setup for start of next surface
            elif goodBackLens and surfID > 0 :
                zmxSurface.append(surfID-1)
                name.append(surfName)
                typ.append(surfTyp )
                coating.append(surfCoating)
                medium.append(surfMedium)
                material.append(surfMaterial)

            #setup for start of next surface
            if goodMirror :
#                print('2 goodMirror,goodFrontLens,goodBackLens: ' +
#                    str(goodMirror) + ' ' + str(goodFrontLens) + ' ' +
#                    str(goodBackLens)) 
                goodMirror = False
                goodFrontLens = False

            goodBackLens = False

            if goodFrontLens :
                goodBackLens = True #Flag we are expecting the back
                                    #surface of a lens Could be another
                                    #Front surface)
                surfName = 'L' + str(lensIndex) + 'E'
                surfTyp = 'lens'
#                print('3 goodMirror,goodFrontLens,goodBackLens: ' +
#                    str(goodMirror) + ' ' + str(goodFrontLens) + ' ' +
#                    str(goodBackLens)) 
                goodFrontLens = False
                #Pickup zemacs surface ID number
                surfID=int(float(line.split()[1]))
                                    #defaults for next suface (if it is used)
                surfCoating  = 'lenses.txt'
                surfMedium   = 'air'
                surfMaterial = 'none'
               
        elif 'MNUM' in line:
            #End of first intrument version found
            #Assume this surface was the final focal plane.
            #place detector here
            zmxSurface.append(surfID)  #surfID is tricky here
            name.append( 'Det' )
            typ.append( 'det'  )
            coating.append( 'detector_ar' )
            medium.append( 'air' )
            material.append( 'none')
            #we are done.
            break
            #Process the surface
        elif 'GLAS' in line:
            #This always exist for the mirror and  front lens surfaces
            surfMedium = ( line.split()[1])
            if surfMedium == 'MIRROR' : 
#                print(' surfMaterial:' + surfMaterial)      
                mirrorIndex += 1
                surfName = 'M' + str(mirrorIndex)
                surfTyp = 'mirror'
                surfCoating = 'protected_al2o3'
                surfMedium = 'air '
                surfMaterial = 'none'
                goodMirror = True
                goodFrontLens = False
                goodBackLens  = False
            else :
                #Its a the front surface of a lens (maybe a filter,
                #treat as front surface of a lens)
#                print(' surfMedium:' + surfMedium)      
                lensIndex += 1
                surfName = 'L' + str(lensIndex)
                surfTyp = 'lens'
                surfCoating = 'lenses.txt'
                surfMaterial = 'none'
                goodMirror = False
                goodFrontLens = True
                goodBackLens  = False
#                print('4 goodMirror,goodFrontLens,goodBackLens: ' +
#                    str(goodMirror) + ' ' + str(goodFrontLens) + ' ' +
#                    str(goodBackLens)) 

               #Ignore all other commands.

    #All input.lis data is in arrays. Save to a new input.lis file. Any
    #existing input.lis file get overwritten.
    inputList = 'input.lis'
    try:
        os.remove(inputList)
    except OSError:
        pass

    numSurfaces= len(name)
    print('Number of Surfaces: ' + str(numSurfaces))

    #Save the newly generated input.lis file to disk (locally) replacing any
    #existing input.lis file.
    for i in range(len(name)):
        printInputList(inputList, zmxSurface[i], name[i], typ[i], coating[i],
                       medium[i], material[i])

    #And we are done generating and saveing the input.lis file. 
      

      
#Iterate over the various telescope versions. (optics_*.txt files)
#where * goes from 0 to fltNum-1 
for flt in range(fltNum):
    output='optics_'+str(flt)+'.txt'
    try:
        os.remove(output)
    except OSError:
        pass

    zprev=0.0       #Init the previous absolute  z  value.

    #Iterate over the surfaces listed in the input.lis file
    for i in range(len(name)):
        #Find the next surface that is listed in the input.lis file. We
        #have to keep track of the cumulative offset in z (position) from
        #the last listed surface. PHOSIM wants that value in the
        #optics_*.txt file
        
        #Fill in what we can for this surface from the zmemax.zmx file.
        surface=findSurface(zmxFile, zmxSurface[i], zmxSurface[0], flt, medium[i], zprev)
                            #Note: surface.medium now holds medium type.

        #Save previous absolute z of this surface
        zprev=surface.z

        #Write this surface out to the optics_*.txt file for this version
        #of the telescope
        printOptics(output,surface,name[i],typ[i],coating[i],surface.medium,
                    material[i],flt)
    print('New '+ output +' file created, but it still needs an OPTICS ' +
              'line at the begining to be added before it can be used!')
    print('Segmentation and focalplane files may need to be constructed')
    
