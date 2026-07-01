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

#include "starobject.h"

StarObject::StarObject(CGParams* pCGPar, std::ofstream& objFile): 
                         CatalogObject(pCGPar)
  //StarObject::StarObject(CGParams* pCGPar, std::ofstream& objFile, 
  //Random* prandom):
  //                            CatalogObject(pCGPar,prandom)  
{
  //fZenithMaxRange used to generate  ra/dec for random stars.
  fZenithMaxRange = ( 1. - cos( ( pCGParams->fStarsDiameterFOVDeg / 2. )
                                * PAL__DPI / 180.0 ) ) / 2.0;
  fOfs = &objFile;

  // ***********************************************
  // Load CastillKurucz2004 allowed parameter value vectors;
  // ***********************************************
  // fCK04Metalicity AND fCK04Logg values loaded with vector assign method.
  // Only way I see to do it pre c++-11
  int ck04Marray[] = { -25, -20, -15, -10, -5, 0, 2, 5 };
  fCK04Metallicity.clear();
  for (int i = 0; i< 8; i++) {
      fCK04Metallicity.push_back(ck04Marray[i]);
  }

  fCK04Teff.clear();
  for ( int i = 3500; i <= 13000; i = i + 250 ) {
      fCK04Teff.push_back(i);
  }
  for ( int i = 14000; i <= 50000; i = i + 1000 ) {
      fCK04Teff.push_back(i);
  }

  fCK04Logg.clear();
  int ck04Logg[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50}; 
  for (int i = 0; i < 11; i++) {
      fCK04Logg.push_back(ck04Logg[i]);
  }

  // ***********************************************
  // Load Kurucz1993 allowed parameter value vectors;
  // ***********************************************
  int F93M[] = { -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, -3, -2, -1,
                 0, 1, 2, 3, 5, 10 }; 
  fK93Metallicity.clear();
  for( int i = 1; i < 19; i++) {
      fK93Metallicity.push_back(F93M[i]);
  }

  fK93Teff.clear();
  for ( int i = 3500; i <= 10000; i = i + 250 ) {
      fK93Teff.push_back(i);
  }
  for ( int i = 10500; i <= 13000; i = i + 500 ) {
      fK93Teff.push_back(i);
  }

  for ( int i = 14000; i <= 34000; i = i + 1000 ) {
    fK93Teff.push_back(i);
  }
  for ( int i = 35000; i <= 50000; i = i + 2500 ) {
      fK93Teff.push_back(i);
  }

  int F93Logg[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};  
  fK93Logg.clear();
  for (int i = 0; i < 11; i++ ) {
      fK93Logg.push_back(F93Logg[i]);
  }
}

void StarObject::GenerateGridStars()
{
  // *******************************************************
  // Using a single default (for now) SED for all grid stars
  // *******************************************************
  StarGridSetSEDFileName();              //CatalogObject Method

  // ***************************************************************
  // Get Ra/Dec of where the catalog is centered. Maynot be the same as 
  // telescope pointing.
  // ***************************************************************
  double raRad = pCGParams->fCatalogFOVRADeg * PAL__DPI / 180.0;
  double decRad  = pCGParams->fCatalogFOVDecDeg * PAL__DPI / 180.0;


  // **************************************************************
  //Generate the grid locations, centered on raRad,decRad.
  // Places ra/dec values for grid in rightAscension and declination
  // vectors.
  // **************************************************************
  GenerateGridRaDec(raRad, decRad, pCGParams->fStarGridWidthDeg, 
					       pCGParams->fStarGridSpacingDeg);   
                                                 //CatalogObject method
  SetMagnitude(pCGParams->fStarGridMagnitude);
  
  int numRaDecPoint=fRightAscension.size();

  // ************************************
  // Write star object to catalog file
  // ************************************
  for ( int i = 0; i < numRaDecPoint; i++ ) {
    pCGParams->fCatalogObjectID++;
    *fOfs<< "object " << pCGParams->fCatalogObjectID << " " << std::fixed 
       << std::setprecision(8) << fRightAscension.at(i) * 180.0 / PAL__DPI 
       << " " << fDeclination.at(i) * 180.0 / PAL__DPI 
       << "  " << std::setprecision(6) << GetMagnitude()
       << " "  <<  GetSEDFileName() 
       << " 0.0 0.0 0.0 0.0 0.0 0.0 star none none" 
       << std::endl;
  }
  pCGParams->fCatalogObjectCount += numRaDecPoint;
  std::cout<< "Number of stellar grid objects added to catalog:            " 
           << numRaDecPoint << std::endl;
  return;
}
// ****************************************************************

void StarObject::GenerateGalacticStars()
// ****************************************************************
// Generate stars randomly (but repeatably) with a density that is
// galactic latitude dependent (may include longitude dependency later)
// Do this by dividing celestial sphere into ra and dec cells whithin which 
// density is aproximatly constant. Use density and area of cell to get 
// mean number in cell, fluxtuate it with approprate poison distribution and 
// create and wirte object stars to catalog file.
// ***************************************************************
{
  // *********************************************************************
  // First test that the CK04 (CastilliKurucz2004) and K93 (Kurucz1993) 
  // catalogs exist. If not we will default to the sed_flat SED file is used
  // *********************************************************************
  int  numStarsInRange=0;
  bool hasCK04Catalog = false;
  bool hasK93Catalog  = false;

  StarsSetSEDFileName();//This just sets default (no catalog) sed_flat.txt 

  if ( !pCGParams->fSEDCatalogCollectionPath.empty() ) {
    std::string CK04catalogPath = 
      pCGParams->fSEDCatalogCollectionPath +  "/CastelliKurucz2004";
    if (std::filesystem::is_directory(CK04catalogPath)) {
      /* Directory exists. */
      hasCK04Catalog = true;
      std::string K93catalogPath = pCGParams->fSEDCatalogCollectionPath +
        "/Kurucz1993";
      if (std::filesystem::is_directory(K93catalogPath)) {
        hasK93Catalog = true;
      }
    }
  
    if (!hasCK04Catalog || !hasK93Catalog ) {
      std::cout << "****************************************************"
                <<std::endl;
      std::cout << "* WARNING! Required SED stellar catalogs missing from " 
                << pCGParams->fSEDCatalogCollectionPath << std::endl;
      std::cout << "* WARNING! Defaulting to all SED's as sed_flat.txt" 
                << std::endl;
      std::cout << "* WARNING! See: https://bitbucket.org/phosim/"
                   "phosim_release/wiki/Using%20PhoSim" << std::endl;

      std::cout << "****************************************************"
                <<std::endl;
    }
  } // End of check that catalog collection Path exist

  // ****************************************************
  // Find all the cells in our FOV (and then some for safty's sake)
  // Places cells in vector fCells
  // ****************************************************
  GetSphericalCellsInFOV(pCGParams->fStarsDiameterFOVDeg);

  // *****************************************************
  // Iterate through the cells
  int numCells = fCells.size();
  for (int i = 0; i < numCells; i++) {
    // We will now need to iterate though the Magnitude range since density is 
    // a function of magnitude. Use deltaM=1 steps
    // Since our densitiy is for magnitude bins with integer low edge 
    // and high edge, handle partial lowest and highest bins.

    int lowMEdge=int(pCGParams->fStarsMinimumMagnitude);  //Round DOWN!

    Random random;
    // Determine random seed for this cell. Using a different but cell 
    // index dependent seed for each cell insures we can reproduce same 
    // stars independent of origin of FOV
    // Set the seed in this cell random number generator
    //Note this way every cell will always have the same stars in it.
    int cellSeed = fCells.at(i).fPhiIndex +
                                      fCells.at(i).fThetaIndex * fNumPhiCells; 
    random.setSeed32(cellSeed);
    //std::cout << "i, cellSeed: " << i << " " << cellSeed <<std::endl;

    while ( lowMEdge < pCGParams->fStarsMaximumMagnitude ) {
      // ***************************************************
      // Mean number of stars in this Magnitude bin for this cell
      // Using low edge of bin here. 
      // Possible lowest and hists bins are fractional. Handle low edge first
      // We only reset m here if lowMEdge is not the minimum magnitude in this
      // bin. Otherwise we use lowMEdge for m.
      // ****************************************************
      double mBinFraction = 1;
      double m = (double) lowMEdge;
      if (double(lowMEdge) < pCGParams->fStarsMinimumMagnitude ) {
        m=pCGParams->fStarsMinimumMagnitude;
        mBinFraction = (lowMEdge+1)-pCGParams->fStarsMinimumMagnitude;
      }
      else  if (lowMEdge+1. >pCGParams->fStarsMaximumMagnitude ) {
        //Partial top Magnitude cell. Scale things as needed
        mBinFraction = (pCGParams->fStarsMaximumMagnitude-lowMEdge);
        //use defulat low edge (above)
      }

      double density = fCells.at(i).CalculateStarDensity(m);
      double meanNumStars = density * fCells.at(i).fAreaDeg2 * mBinFraction;

      //Possion fluxtuate this number to get correct statistical number
      int numStars = random.poisson(meanNumStars);
      
      //Now generate these stars and place in the catalog
      if (numStars >0 ) {
        for (int j = 0; j < numStars;  j++ ) {
          double RADeg=0.;
          double decDeg=0.;
          //fCells.at(i).GenerateRandomRADec(RADeg, decDeg);
          fCells.at(i).GenerateCellRandomRADec(RADeg, decDeg, random);
          
          // At this point we need to check that this particular star is within 
          // the exact FOV. (not the modified). This saves time and space.
          double angDistDeg = eraSeps( 
                                 RADeg * (PAL__DPI/180.), 
   
                              decDeg * (PAL__DPI/180.), 
                                 pCGParams->fCatalogFOVRADeg *(PAL__DPI/180.), 
                                 pCGParams->fCatalogFOVDecDeg * 
                                        (PAL__DPI/180.) ) * (180./PAL__DPI);

          // *************************************************************
          // In order for the random number stream to remain consistent for
          // this cell we need to pre calculate the random nnumber we would 
          // normaly use after the  following test
          double u = random.uniform();
          if ( angDistDeg <= pCGParams->fStarsDiameterFOVDeg/2.) { 
            //Do an even distribution within magnitude bin for M.
            //double magnitude= m * fCells.at(i).fRandom.uniform() * 
            //Handle special cases of first and last bin by setting m and 
            // mbinFraction as appropriate above.
            double magnitude;
            magnitude = m + (u * mBinFraction);

            //Debug print
            //double longDeg;
            //double latDeg;
            //fCells.at(i).GetGalacticFromRaDec(RADeg, decDeg, longDeg, latDeg);
            //std::cout << RADeg << " " << decDeg << " " << magnitude << " " 
            //          << mBinFraction << " " << latDeg << " " << longDeg 
            //          << std::endl;
            
            // *******************
            // SED Generation
            // If catalogs don't exist use flat SED
            // *******************

            std::string starSEDFileName;
            if ( hasCK04Catalog && hasK93Catalog) {
              starSEDFileName = GenerateStarSEDFileName(random);
            }
            else {
              // No catalogs, use sed_flat.txt
              starSEDFileName =  GetSEDFileName();
            }

            // Write star to catalog 
            pCGParams->fCatalogObjectID++;
            *fOfs<< "object " << pCGParams->fCatalogObjectID << " " 
                 << std::fixed << std::setprecision(8) << RADeg << " " 
                 << decDeg 
                 << "  " << std::setprecision(6) << magnitude << " "  
                 <<  starSEDFileName 
                 << " 0.0 0.0 0.0 0.0 0.0 0.0 star none none" << std::endl;
            pCGParams->fCatalogObjectCount++;        
          
            numStarsInRange++;
          }
        }
      }
      //Bump to next magnitude bin.
      lowMEdge++;
    }
  }
  std::cout<< "Number of stellar objects added to catalog:                 " << numStarsInRange << std::endl;
  return;
}  
// **************************************************************************

std::string  StarObject::GenerateStarSEDFileName(Random& random)
// *************************************************************************
// Get a SED for this star. If SED libraries exits ( we use
// CastelliKurucz2004 (call "CK04" here))  and KURUCZ1993 ( called "k93" here)
// generate from  distributions
// [Fe/H], Teff and Logg for this star. If synthetic the catalogs 
// (CastelliKurucz2004 and KURUCZ1993) are unavailable default to the use 
// of the flat SED. 
// ************************************************************************
// CASTELLI AND KURUCZ 2004
// The atlas contains about 4300 stellar atmosphere models for a wide range of 
//  metallicities, effective temperatures and gravities. These LTE models 
// with no convective overshooting computed by Fiorella Castelli, have 
// improved opacities and abundances upon previously used by Kurucz (1990). 
// The main improvements from previous opacity distribution functions listed 
// in: Castelli & Kurucz 2003 (IAU Symposium 210, Modelling of Stellar 
// Atmospheres, Uppsala, Sweden, eds. N.E. Piskunov, W.W. Weiss. and 
// D.F. Gray, 2003, ASP-S210)
// Found at: 
// http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html
// *************************
// KURUCZ 1993
// This  ATLAS contains about 7600 stellar atmosphere models covering a 
// wide range of metallicities, effective temperatures and gravities. The 
// original atlas was created on August 22, 1993 by Dr. Kurucz 
// Found at:
// http://www.stsci.edu/hst/observatory/crds/astronomical_catalogs.
// html#kurucz1993
// *************************************************************************
// This will be a very simple "toy" spectra generation model. No claim to any 
// great accuracy is made here. 
// First a metallicity is picked from a Gaussian distribution.  Then a 
// star mass if picked from a simple IMF model. This Stellar Mass is used to 
// determine a Teff. Then using this Teff a Logg is generated. 
// Mostly I've tried to stay with distributions determined for the local Solar 
// neighborhood. Thus a lack of high metalicities as would be expected in the 
// bar of the core. I've also limited the mass to a minimum of .5 solar 
// masses for Teff and thus Logg.
// Various improvements to this model can be made in a striaght forward mannor 
// for all 3 parameters.. This may happen if requested.(ie. If color color 
// plots look funny I can add some more randomness to smear things out, go to 
// lower Teff etc.
// *************************************************************************
// References:
// Metallicity: [Fe/H]: 
// "S4N: A Spectroscopic Survey of Stars in the Solar Neighborhood"
// Carlos Allende Prieto ,et al, Astron.Astrophys. 420 (2004) 183-205,
//  arXiv:astro-ph/0403108v1,  https://arxiv.org/abs/astro-ph/0403108v1
//
// IMF:
//"On the variation of the initial mass function". Kroupa, Pavel (2001). 
// MNRAS. 322 (2): 231–246. arXiv:astro-ph/0009005. 
// Bibcode:2001MNRAS.322..231K.,  doi:10.1046/j.1365-8711.2001.04022.x
//
// Logg:
// "The atmospheric parameters and spectral interpolator for the MILES stars",
// Ph. Prugniel, I. Vauglin, and M. Koleva,A&A 531, A165 (2011), 
// https://www.aanda.org/articles/aa/pdf/2011/07/aa16769-11.pdf
// from plot on page 7.
// *************************************************************************
// If we are unable to genrate an acceptable file name we will make another 
// attemp. This routine thus always return an existing SED file name.
// *************************************************************************
// random comes as an argument so that we are deterministic(can reporduce 
// exactly if same seed used).
// *******************************************************************88
{  
  bool hasCK04File = false;
  bool hasK93File  = false;
  double metallicity = 0.;
  double Teff = 0.;
  double logg = 0.;
  std::string SEDName;

  while ( !hasCK04File && !hasK93File) {
    // ********************
    // Pick a metalicity from a Gaussian(normal) distribution (solar 
    // neighborhood)
    metallicity = random.normal()*kStellarMetallicitySigma + 
                                                    kStellarMetallicityMean;

    // Pick a stellar mass (in units of Solar mass) from simple IMF.
    // n ~ M**(kIMFAlpha) for M > .5Msun
    // Note we may want to put an upper limit cut here (stellarMass<10Msun?)
    double r = random.uniform();
    double stellarMass = pow( (1.-r), (1./(1.+kIMFAlpha))) * 
                                                         kStellarMassMinimum;

    // Use solar mass and temp to calibrate how we generate Teff from 
    // stellarMass.  Using: Teff~M**(3/4).
    Teff = kSolarTeff*pow(stellarMass,(3.0/4.0) );

    // Logg is function of Teff. Approximate with 2 lines in Logg vs 
    // log10(Teff/1000)
    // plot. One for Teff>5814 and another for Teff < 5844.  We could improve 
    // this somewhat by producing a logg and then spread it out a liltle. See 
    // plot in Prugniel on page 7. For now we won't.

    if (Teff > kTeffBreak ) {
      logg = log10(Teff/1000.)* kLoggLogTeffHiSlope + kLoggLogTeffHiIntercept;
      //Derived directly from plot (very crudely)
    }
    else {
      logg = log10(Teff/1000)* kLoggLogTeffLoSlope + kLoggLogTeffLoIntercept;
      //Derived directly from plot (very crudely)
    }

    // ***************************************************
    // We now have a [Fe/H], Teff, and Logg for a star. These values are 
    // encoded within the SED file name. We will make up a file name here 
    // for the SED file. SED files for the CK04 catalog all start with 
    // "ck". All SED files in the k93 catalog start with "k". 
    // Note that that the CK04 SEDs are newer and improved over the same SEDs 
    // in the k93 catalog. However the SEDs in the CK04 catalog are a limited 
    // subset of the SEDs in the k93 catlog. So our algorithm is to look first
    // for a SED in the CK04 catalog and if its not there look for it in the 
    // k93 catalog.
    // ***************************************************
    // Generate the CK04 base file name:
    // *****************************
  
    // First try the CK04 catalog. It has priority.
    // Its SEDs are salight improvment over K93' SDEs.
    // If that didn't work check the K93 catalog
    hasCK04File  = BuildStarSEDName("CK04",metallicity, Teff, logg, SEDName);
    if( !hasCK04File ) {
      hasK93File = BuildStarSEDName("K93", metallicity, Teff, logg, SEDName);
    }

  }
  //std::cout <<metallicity << " " << Teff << " " << logg << std::endl;
  return  SEDName;
 
}
// ***********************************************************************

bool StarObject::BuildStarSEDName(std::string Catalog, double Metallicity, 
                            double Teff, double Logg, std::string& SEDName)
// *********************************************************************
// Generate, if we can, the SED file name (includes path) to use as the SED
// for the star with given values of the parameters.
// Only CK04 (CastelliKurucz2004) and F93 (Kurucz1993) catalogs are allowed 
// for now.
// If the parameters are out of range for the selected catalog ( specified in 
// string variable "Catalog") return false.
// **********************************************************************
// File names in CK04 catalog have the form:
// "ck"metallicity*10"_"Teff".fits._g"logg""_"Teff".gz"
// or
// "k"metallicity*10"_"Teff".fits._g"logg""_"Teff".gz"

// Example: ckm15_5250.fits_g10_5250.gz
// This has [Fe/H] = -1.5, Teff = 5250, Logg =1.0  
// Yes the Teff comes in twice!
// ***********************************************************************
{
  //Set up vectors
  std:: vector < int>* pMetallicity;
  std:: vector < int>* pTeff;
  std:: vector < int>* pLogg;
  std::string catalogSpecifier;
  std::string catalogDirectory;
  std::string ck04Str = "CK04";
  std::string k93Str = "K93";
  if ( Catalog != ck04Str && Catalog != k93Str ) {
    std::cout << "Specified Catalog: " << Catalog << " not:" << ck04Str 
              << " and not: "<< k93Str  << std::endl;
    exit(1);
  }
  // *********************************************************************
  // Load up pointers to allowed values. Vectors setup in constructor of 
  // this class.
  // ********************************************************************* 
  if (Catalog == ck04Str ) {
    pMetallicity = &fCK04Metallicity; 
    pTeff        = &fCK04Teff;
    pLogg        = &fCK04Logg;
    catalogSpecifier  = "ck";
    catalogDirectory  = "CastelliKurucz2004";
  }
  if (Catalog  == k93Str ) {
    pMetallicity = &fK93Metallicity; 
    pTeff        = &fK93Teff;
    pLogg        = &fK93Logg;
    catalogSpecifier = "k";
    catalogDirectory  = "Kurucz1993";
  }
  // Test if within Catalog Range:
  if (Metallicity >= ( (double)pMetallicity->front() / 10.0 ) && 
     Metallicity <= ( (double)pMetallicity->back() / 10.0   ) &&
     Teff >= (double)pTeff->front() && 
     Teff <= (double)pTeff->back()  &&
     Logg >= ( (double)pLogg->front() / 10.)  && 
     Logg <= ( (double)pLogg->back() / 10. )  )   {
    // Paramters are within the SED catalog limits
    // Build file name and see if it exists
    
    // Metallicity:
    // The value  in the SED file name is encoded as the multiplicity * 10.0
    // Only the ordered values in vector pointed to by pMetalicity 
    // are allowed.
    // (Logg is handled simularly)
    int closestMeta = ClosestElement( pMetallicity, Metallicity*10.0);
    std::string metallicitySign = "p";
    if (closestMeta < 0 ) {
      closestMeta = -closestMeta;
      metallicitySign = "m";
    }
   
    //Now Teff:
    // The allowed Teff values are listed in sorted vector pTeff
    // Teff are always positive
    int closestTeff = ClosestElement( pTeff, Teff);
    
    // Logg:
    // Logg always positive.
    int closestLogg = ClosestElement( pLogg, Logg*10.0);
    
    // Build the SED file name
    std::ostringstream SEDss;
    SEDss.str("");
    SEDss << catalogDirectory << "/"
          << catalogSpecifier << metallicitySign
          << std::setw(2) << std::setfill('0') << closestMeta 
          << "_" << closestTeff 
          << ".fits_g" << std::setw(2) << std::setfill('0') << closestLogg 
          << "_" << closestTeff <<".gz";
    // Now we have the file name. Check that such a file exists in the 
    // this catalog.(not all cominations of closestMeta, closestTeff and 
    // closestLogg have exist in CK04 or K93.
    SEDName = "./" + SEDss.str();  //This is what raytrace wants.
    
    std::string fullSEDName = pCGParams->fSEDCatalogCollectionPath + "/" + 
      SEDss.str();
    std::ifstream infile(fullSEDName.c_str());
    if ( infile.good()) {
      // The SED  File exists. We are done.
      infile.close();
      return true;
    }       
  } 
  //This set of parameterts does not have an SED in the CK04 catalog
  return false;  
}
// *********************************************************************

