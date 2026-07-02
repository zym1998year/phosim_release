// @brief CGParams
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

#include "cgparams.h"

 
CGParams::CGParams(int argc, char* argv[])
{
  // *****************************************************************
  // This code has little to no error checking. Assumed this program is run 
  // by another program (phosim.py usually)  which knows how things should
  // be setup
  // *****************************************************************

  SetDefaults();

  // This constructer parses the parametersa form the input (obs_*.pars) 
  // file comming in though stdin.

  std::string cmdArguments;

  // Read parameters from stdin.
  readText pars(std::cin);

  for (size_t t(0); t < pars.getSize(); t++) {
	std::string line(pars[t]);
	readText::get(line, "obshistid",   fObservationID);
	readText::get(line, "obsseed",     fRandomSeed);
	readText::get(line, "pointingra",  fRADegOrigin);
	readText::get(line, "pointingdec", fDecDegOrigin);
	readText::get(line, "seddir",      fSEDCatalogCollectionPath);

  }

  // have to run through twice in case pointing appears at the end
  for (size_t t(0); t < pars.getSize(); t++) {
	std::string line(pars[t]);

	std::istringstream iss(line);
	std::string keyName;
	iss >> keyName;

	if ( keyName == "stars") {
	  fMakeStars = true;
	  iss >> fStarsMinimumMagnitude >> fStarsMaximumMagnitude 
		  >> fStarsDiameterFOVDeg;
	  //We may have center of FOV different from telescop pointing.
      //This will be in last 2 arguments
      if( iss.good()) {
        iss >> fCatalogFOVRADeg >> fCatalogFOVDecDeg;
      }
      else {
        fCatalogFOVRADeg  = fRADegOrigin;
        fCatalogFOVDecDeg = fDecDegOrigin;
      }
    }

  	if ( keyName == "stargrid") {
	  fMakeStarGrid = true;
	  iss >> fStarGridSpacingDeg >> fStarGridWidthDeg >> fStarGridMagnitude;
	  //We may have center of FOV different from telescop pointing.
      //This will be in last 2 arguments
      if( iss.good()) {
        iss >> fCatalogFOVRADeg >> fCatalogFOVDecDeg;
      }
      else {
        fCatalogFOVRADeg  = fRADegOrigin;
        fCatalogFOVDecDeg = fDecDegOrigin;
      }
        
	}

  	if ( keyName == "galaxies") {
	  fMakeGalaxies = true;
	  iss >> fGalaxiesMinimumMagnitude >> fGalaxiesMaximumMagnitude 
		  >> fGalaxiesDiameterFOVDeg;
	  //We may have center of FOV different from telescop pointing.
      //This will be in last 2 arguments
      if( iss.good()) {
        iss >> fCatalogFOVRADeg >> fCatalogFOVDecDeg;
      }
      else {
        fCatalogFOVRADeg  = fRADegOrigin;
        fCatalogFOVDecDeg = fDecDegOrigin;
      }
    }
  }

}
// **********************************************************************

void CGParams::PrintOptionValues() {
    
    std::cout << "Observation ID:                                             " << fObservationID << std::endl;
    std::cout << "Catalog ra (deg):                                           " << fCatalogFOVRADeg  << std::endl;
    std::cout << "Catalog dec (deg):                                          " << fCatalogFOVDecDeg << std::endl;
    

// For debug reason print out galactic coords
  double latDeg;
  double longDeg;
  SphericalCell tempCell;
  tempCell.GetGalacticFromRaDec(fCatalogFOVRADeg, fCatalogFOVDecDeg, longDeg,
                                       latDeg);
  std::cout << "Galactic latitude (deg):                                    " << latDeg << std::endl;
  std::cout << "Galactic longitude (deg):                                   " << longDeg << std::endl;
   
  std::cout << "Random number generator seed:                               " << static_cast<long>(fRandomSeed) << std::endl;
  std::cout << "SED path:                                                   " << fSEDCatalogCollectionPath << std::endl;
  if (fMakeStars) {
      std::cout << "Make stars:                                                 True"      << std::endl;
      std::cout << "Min stars magnitude:                                        " << fStarsMinimumMagnitude 
                << std::endl;
      std::cout << "Max stars magnitude:                                        " << fStarsMaximumMagnitude 
                << std::endl;
      std::cout << "Stars FOV diameter (deg):                                   " << fStarsDiameterFOVDeg 
                << std::endl;
      std::cout << "Stars FOV Area (sq. deg):                                   " << PAL__DPI * 
                                  pow(fStarsDiameterFOVDeg/2.,2) << std::endl;
  } else {
      std::cout << "Make stars:                                                 False" << std::endl;
  }

  if (fMakeStarGrid) {
      std::cout<< "Make grid of stars:                                         True" << std::endl;
      std::cout<< "Star grid spacing (deg):                                    " << fStarGridSpacingDeg 
               << std::endl;
      std::cout<< "Star grid magnitude (deg):                                  " << fStarGridMagnitude 
               << std::endl;
      std::cout<< "Star grid FOV diameter (deg):                               " << fStarGridWidthDeg 
               << std::endl;
  } else {
      std::cout<< "Make star grid:                                             False"<< std::endl;
  }

  if (fMakeGalaxies) {
      std::cout<< "Make galaxies:                                              True" << std::endl;
      std::cout<< "Min galaxies magnitude:                                     " << fGalaxiesMinimumMagnitude 
               << std::endl;
      std::cout<< "Max galaxies magnitude:                                     " << fGalaxiesMaximumMagnitude 
               << std::endl;
      std::cout<< "Galaxies FOV diameter (deg):                                " << fGalaxiesDiameterFOVDeg 
               << std::endl;
      std::cout<< "Galaxies FOV Area (sq. deg):                                " << PAL__DPI * 
                              pow(fGalaxiesDiameterFOVDeg/2.,2) << std::endl;
  } else {
      std::cout<< "Make galaxies:                                              False" << std::endl;
  }
  std::cout << "---------------------------------------------------------------------"
      "-------------------------------" << std::endl;
  return;
}
// **********************************************************************

void CGParams::SetDefaults()
{
  fMakeStars    = false;
  fMakeStarGrid = false;
  fMakeGalaxies = false;
  fStarGridSEDFileName = "../sky/sed_flat.txt";

  return;
}
// ************************************************************************

