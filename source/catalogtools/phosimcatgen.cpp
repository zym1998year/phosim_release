///
/// @package phosim
/// @file phosimcatgen.cpp
/// @brief Generates phosim Object catalog and files with 
///        random or grid of stars and/or galaxies.
///
/// @brief Created by: 
/// @author Glenn Sembroski (Purdue)
/// @brief Modified by:
/// @author 
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

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
#include "catalogobject.h"
#include "starobject.h"
#include "galaxyobject.h"

#include "../ancillary/readtext.h"
#include "../ancillary/random.h"

int main(int argc, char* argv[]) {

  std::cout << "-----------------------------------------------------------"
	           "-----------------------------------------" << std::endl;
  std::cout << "Catalog Generator" << std::endl;
  std::cout << "-----------------------------------------------------------"
	           "-----------------------------------------" << std::endl;

  // First handle the input options.
  CGParams* pCGParams = new CGParams(argc,argv);
  pCGParams->PrintOptionValues();

  //Random*  pRandom = new Random();

  // Need to set random seed  before we use random numbers
  //pRandom->setSeed32( pCGParams->fRandomSeed); 
  
  // **********************************************************************
  // Initalize object IDs to start at zero. While this may duplicate object 
  // ids in other catalogs thats ok.
  // **********************************************************************
  pCGParams->fCatalogObjectID=0;    
  pCGParams->fCatalogObjectCount=0;    


  // *************************************************************************
  // Now we start to produce objects for the catalog.  We may have a combined 
  // catalog of stars, galaixies and a stargrid as chosen by user.
  // *********************************************************************
  // Set up the catalog file name and open the catalog
  if (pCGParams->makeStars() || pCGParams->makeStarGrid() || 
	  pCGParams->makeGalaxies()) {
	

	pCGParams->fCatFileName = "catgen_" + pCGParams->fObservationID +  ".cat";
	pCGParams->fCatFile.open(pCGParams->fCatFileName.c_str(), std::ios::out);
  }
  else{
	std::cout<<" No catalog request found. No catalog made"<< std::endl;
	pCGParams->fCatFile.close();
	exit(1);
  }
  // *************************************************************************
  // First are random stars.All use default SED.
  // *************************************************************************
  if ( pCGParams->makeStars() ) {
	StarObject starObj(pCGParams,pCGParams->fCatFile);
    // *********************
    // Generate stars with galactic latitude density distribution and 
    // magnitude. No longitude distribution as of yet.
    // This method also used the sphericalCell alogrithum for reproducibility.
    // *********************
    starObj.GenerateGalacticStars();
}

  // ***************************
  // Grid of stars. All have same magnitiude and defualt SED
  // ***************************
  if( pCGParams->makeStarGrid() ) {
	StarObject starGridObj(pCGParams,pCGParams->fCatFile);
	starGridObj.GenerateGridStars();
  }

  // ************************************************************************
  // Random Complex Galaxies (in ra/dec cells)
  // ************************************************************************
  if ( pCGParams->makeGalaxies() ) {
	GalaxyObject galaxyObj(pCGParams,pCGParams->fCatFile);
	galaxyObj.GenerateComplexGalaxiesForCells();
  }

  // Display number of objects created (stars + grid + galaxies) in file.
  std::cout << "----------------------------------------------------------------------------------------------------"
            <<std::endl;
  std::cout<< "Total number of simulated objects:                          "<< pCGParams->fCatalogObjectCount 
           << std::endl; 
   return 0;
}			  

