// @brief StarObject
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

#ifndef STAROBJECT_H
#define STAROBJECT_H


#include "catalogobject.h"

const double kStellarMetallicitySigma = 0.19;
const double kStellarMetallicityMean = -0.09;
const double kIMFAlpha = -2.3;         //For range .5Msun < M
const double kStellarMassMinimum = .5; //Other ranges are possible,check the 
                                       //references listed in 
                                       //GenerateStarSEDFileName()
const double kSolarTeff = 5778;
const double kLoggLogTeffLoSlope = 16.6113;  //For Teff <  5814(K)
const double kLoggLogTeffLoIntercept = -8.6113;
const double kLoggLogTeffHiSlope = -0.7299;  //For Teff >= 5814(K)
const double kLoggLogTeffHiIntercept = 4.646;
const double kTeffBreak = 5814.;

class StarObject :public CatalogObject
// **********************************************************
// Note this class derived from CatalogObject and thus has 
// CatalogObject's methods
// **************************************************
{
 public:
  //StarObject(CGParams* pCGPar, std::ofstream& objFile, Random* prandom);
  StarObject(CGParams* pCGPar, std::ofstream& objFile);
  void GenerateGridStars();
  void GenerateRaDecFromStarXYMFile(std::string fileName);
  void GenerateGalacticStars(); 
  std::string  GenerateStarSEDFileName(Random& random);
  void   StarsSetSEDFileName()
              {fCatalogSEDFileName=pCGParams->fStarGridSEDFileName;return;};
  void   StarGridSetSEDFileName()
              {fCatalogSEDFileName=pCGParams->fStarGridSEDFileName;return;};
private:
  bool BuildStarSEDName(std::string Catalog, double Metallicity, double Teff, 
                  double Logg, std::string& SEDName);

  std::ofstream* fOfs;

  // Catalog Star SED acceptable values (sorted!!)
  std::vector< int > fCK04Metallicity;
  std::vector< int > fCK04Teff;
  std::vector< int > fCK04Logg;
  std::vector< int > fK93Metallicity;
  std::vector< int > fK93Teff;
  std::vector< int > fK93Logg;

};
// *****************************************************************


#endif
