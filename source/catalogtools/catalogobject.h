// @brief CatalogObject
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

#ifndef CATALOGOBJECT_H
#define CATALOGOBJECT_H

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
#include <dirent.h>     // UIsed to check if catalog dir exist.

#include "cgparams.h"
#include "sphericalcell.h"
#include "../ancillary/readtext.h"
#include "../ancillary/random.h"

//Constants for pal* and era* functions (and also used as needed otherwise)
const double ERFA_D2PI = 6.283185307179586476925287;
const int kPrintFreq = 10000;

class CatalogObject
{
 public:
  //CatalogObject(CGParams* pCGPar , Random* pRandom);
  CatalogObject(CGParams* pCGPar );
  void   GenerateSingleRandomRaDec();
  void   GenerateGridRaDec(double RARad,  double DecRad, 
                           double FOVDeg, double StepSizeDeg);
  double GetRADeg(){return fRARad*180.0/PAL__DPI;};
  double GetDecDeg(){return fDecRad*180.0/PAL__DPI;};
  void   SetMagnitude(double mag ){fMagnitude=mag;return;};
  double GetMagnitude(){return fMagnitude;};
  std::string GetSEDFileName(){return fCatalogSEDFileName;};
  void   GetSphericalCellsInFOV(double FOVDeg); 
  int    ClosestElement(std::vector<int>* pVec, double value);
  double ClosestElement(std::vector<double>* pVec, double value);
  bool   fDebugPrint;

    double angularSep(double ra1, double dec1, double ra2, double dec2);

  // ***********************************************************************
  // Utility function definitions for Sofa routines (see refereneces below)
  // ***********************************************************************
  double eraAnp(double a);
  double eraPdp(double a[3], double b[3]);
  double eraPm(double p[3]);
  void   eraPxp(double a[3], double b[3], double axb[3]);
  void   eraS2c(double theta, double phi, double c[3]);
  double eraSepp(double a[3], double b[3]);
  double eraSeps(double al, double ap, double bl, double bp);
  void   palDtp2s ( double xi, double eta, double raz, double decz,
                    double *ra, double *dec );
  void   palDh2e ( double az, double el, double phi, double *ha, double *dec);

  // ***********************************************************************

 protected:
  CGParams* pCGParams;
  Random fRandom;
  std::vector < double > fRightAscension;
  std::vector < double > fDeclination;

  std::vector < SphericalCell > fCells;
  int fNumThetaCells;
  int fNumPhiCells;
  //Object vaiables
  double fRARad;
  double fDecRad;
  double fMagnitude;
  double fZenithMaxRange;
  std::string fCatalogSEDFileName;
 
  std::vector < double > lsstChipCenterXYDeg;

}; 
// ***************************************************************************

#endif
