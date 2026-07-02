// @brief GalaxyObject
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

#ifndef GALAXYOBJECT_H
#define GALAXYOBJECT_H

#include "catalogobject.h"
#include "galaxycomplex.h"
const double kHubbleConstant = 70.5; //Units km/sec/MPC
const double kKMperMPC = 30856775812799586000.0;
const double kSecPerYear = 365.25 * 24. * 60. * 60.;
const double kAgeUniverseYears  = 
                           (1.0/kHubbleConstant) * (kKMperMPC / kSecPerYear); 
const double kMadauHighMin = 3.0;
const double kGalacticMetallicityMean  = .1;
const double kGalacticMetallicitySigma = .1;
const double kZsolar = 0.0196;
const double kBulgeOffsetLogTMean  = .12;
const double kBulgeOffsetLogTSigma = .15;

class GalaxyObject :public CatalogObject
// **********************************************************
// Note this class derived from CatalogObject and thus has 
// CatalogObject's methods
// **************************************************
{
public:
  //GalaxyObject(CGParams* pCGPar, std::ofstream& objFile, Random* pRandom);
  GalaxyObject(CGParams* pCGPar, std::ofstream& objFile);
  void GenerateComplexGalaxiesForCells();
  void GalaxySetSEDFileName()
              {fCatalogSEDFileName=pCGParams->fStarGridSEDFileName;return;};
  void GenerateGalaxySEDFileNames(std::vector <double>& galaxyRandomUniform,
                                  std::vector <double>& galaxyRandomNormal,
                                  GalacticComplexObject& disk,
                                  GalacticComplexObject& bulge);
  double GetRandomRedShiftFromSRF(double slopeIndex,double zGalaxy,double zMax,
                                  double u1);
  void BuildGalaxySEDName(std::string Catalog,  double DiskLogT,
                          double BulgeLogT, double Metallicity, 
                          std::string& BulgeSEDName, std::string& DiskSEDName);
private:
  std::ofstream* fOfs;
  std::vector < double > fPopStarKZ;
  std::vector < double > fPopStarKLogT;
  std::string ConstructSED(std::string catalogDirectory, 
                           std::string catalogSpecifier, int closestMeta,
                           double closestLogT);
};
// *****************************************************************

#endif
