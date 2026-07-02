///
/// @package phosim
/// @file GalaxyComplex.h
/// @brief Generates phosim sersicComplex galaxy parameters
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
/// @author 
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef GALAXYCOMPLEX_H
#define GALAXYCOMPLEX_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>    //std::min_element, std::max_element, std:lower_bound
#include <sstream>

#include "../ancillary/random.h"
const double kPI = 3.1415926535897932384626433832795028841971693993751;
const int kNumValues = 10000;   //Number of Mabsolute,Mapparent and z 
                                //divisions of range.

class GalacticComplexObject
{
  // **********************************************************************
  //Structure to hold sersicComplex values for phosimcatgen GalaxyComplex class
  // **********************************************************************
 public:
  std::string ID;
  double RADeg;
  double DecDeg;
  double M;                                         //magnitude
  std::string SEDFileName;
  double z;                                         //Red Shift

  //sersicComplex parameters
  double a;                      //size of axis 1 in arcseconds
  double b;                      //size of axis 2 in arcseconds
  double c;                      //size of axis 3 in arcseconds
  double alpha;                  //polar angle in Deg.
  double beta;                   //position angle in Deg.
  double gamma;                  //rotation angle in Deg.
  int    index;                  //sersic index
  double clumpFrac;              //fraction of light in clumps
  double clump;                  //number of clumps
  double clumpWidth;             //gaussian clump size in arcseconds
  double spiralFrac;             //fraction of light in spiral
  double alphaSpiral;            //winding angle of spiral in degrees
  double bar;                    //spiral bar size in arcseconds
  double spiralWidth;            //spiral gaussian width in arcseconds
  double phi0;                   //position angle of spiral in deg.
};
// *****************************************************************

class GalaxyComplex
// ***************************************************************
// This class uses IDL code from John Peterson (2018-02-15), converted to C++
// by Glenn Sembroski. (Not sure of definitions or souces of any of this) 
// It produces a random galaxy, both  the disk and bulge magnitudes, 
// the red shift z
// and it fills the sersicComplex class for each (disk and bulge) in order for 
// the calling program to generate a galaxy object(2 objects, disk and bulge
//  actually) for a phosim catalog.
// ****************************************************************
{
 public:
  GalaxyComplex();
  void Init(); 
  bool GenerateGalaxyForCells(GalacticComplexObject& Bulge, 
                              GalacticComplexObject& Disk, double MinimumM, 
                              double MaximumM,  
                              std::vector < double >& galaxyUniformRandom, 
                              Uint64 j);
  double GetGalaxyDensityDeg2(){return fGalaxyDensityDeg2;};
  std::string  GenerateSersicComplexObjectString(
											   GalacticComplexObject& Object);
  void WriteObject(GalacticComplexObject& Object, std::ofstream* ofs);

private:
  double fGalaxyDensityDeg2;
  std::vector<double >::iterator fVectorIteratorUpper;

 public:
  std::vector < double > fMabsoluteValues;
  std::vector < double > fCumulativeMabsoluteProb;
  std::vector < double > fZValues;
  std::vector < double > fCumulativeZProb;
};
// *********************************************************************

#endif
