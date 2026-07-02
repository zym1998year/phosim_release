// @brief CGParms
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef CGPARMS_H
#define CGPARMS_H
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <algorithm>    // std::min_element, std::max_element
#include <sstream>
#include <dirent.h>     // UIsed to check if catalog dirrctory exists.

#include "sphericalcell.h"
#include "../ancillary/readtext.h"
#include "../ancillary/random.h"

const double PAL__DPI  = 3.1415926535897932384626433832795028841971693993751;

// **************************************************************************
class  CGParams
{
public:
  CGParams(int argc, char* argv[]);
  void SetDefaults();
  void PrintOptionValues();
  bool makeStars(){return fMakeStars;};
  bool makeStarGrid(){return fMakeStarGrid;};
  bool makeGalaxies(){return fMakeGalaxies;};

  std::string fObservationID;
  uint32_t fRandomSeed;
  double   fRADegOrigin;       // Telescope pointing
  double   fDecDegOrigin;      
  double   fCatalogFOVRADeg;   // Catalog pointing. Mmay be different from  
  double   fCatalogFOVDecDeg;  // telescope pointing.

  std::string fSEDCatalogCollectionPath;


  //command arguments;
  double fStarsMinimumMagnitude;
  double fStarsMaximumMagnitude;
  double fStarsDiameterFOVDeg;

  double fStarGridSpacingDeg;
  double fStarGridMagnitude;
  double fStarGridWidthDeg;

  double fGalaxiesMinimumMagnitude;
  double fGalaxiesMaximumMagnitude;
  double fGalaxiesDiameterFOVDeg;



private:
  bool fMakeStars;
  bool fMakeStarGrid;
  bool fMakeGalaxies;

public:
  int           fCatalogObjectID;
  int           fCatalogObjectCount;
  double        fStarDensityPerDeg2;
  double        fGalaxyDensityPerDeg2;
  std::string   fStarGridSEDFileName;
  std::string   fCatFileName;
  std::ofstream fCatFile;
};
// ***********************************************************************
 
#endif
