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

#include "galaxyobject.h"

// ************************************************************************

GalaxyObject::GalaxyObject(CGParams* pCGPar, std::ofstream& objFile) :
  CatalogObject(pCGPar)  
  //GalaxyObject::GalaxyObject(CGParams* pCGPar, std::ofstream& objFile, 
  //                           Random* pRandom):
  //  CatalogObject(pCGPar, pRandom)  
{

  //fZenithMaxRange used to generate  ra/dec
  fZenithMaxRange = ( 1. - cos( ( pCGParams->fGalaxiesDiameterFOVDeg / 2. )
                                * PAL__DPI / 180.0 ) ) / 2.0;
  fOfs = &objFile;

  // ***********************************************
  // Load PopStarK catalog allowed parameter value vectors;
  // ***********************************************
  // fPopStarKZ (Metallicity) and fPopStarKLogT (Age) values loaded with 
  // vector assign method. Only way I see to do it pre c++-11
  fPopStarKZ.clear();
  double popStarKZarray[] = {.0001, .0004, .004, .008, .02, .05 }; 
  for (int i = 0; i< 6; i++) {
    fPopStarKZ.push_back(popStarKZarray[i]);
  } 

  fPopStarKLogT.clear();
  double popStarKLogT[] = { 5.00, 5.48, 5.70,
    5.85, 6.00, 6.10, 6.18, 6.24, 6.30, 6.35, 6.40, 6.44, 6.48, 6.51, 6.54,
    6.57, 6.60, 6.63, 6.65, 6.68, 6.70, 6.72, 6.74, 6.76, 6.78, 6.81, 6.85,
    6.86, 6.88, 6.89, 6.90, 6.92, 6.93, 6.94, 6.95, 6.97, 6.98, 6.99, 7.00,
    7.04, 7.08, 7.11, 7.15, 7.18, 7.20, 7.23, 7.26, 7.28, 7.30, 7.34, 7.38,
    7.41, 7.45, 7.48, 7.51, 7.53, 7.56, 7.58, 7.60, 7.62, 7.64, 7.66, 7.68,
    7.70, 7.74, 7.78, 7.81, 7.85, 7.87, 7.90, 7.93, 7.95, 7.98, 8.00, 8.30,
    8.48, 8.60, 8.70, 8.78, 8.85, 8.90, 8.95, 9.00, 9.18, 9.30, 9.40, 9.48,
    9.54, 9.60, 9.65, 9.70, 9.74, 9.78, 9.81, 9.85, 9.90, 9.95, 10.00, 
    10.04, 10.08, 10.11, 10.12, 10.13, 10.14, 10.15, 10.18};
  for (int i = 0; i < 106; i++) {
    fPopStarKLogT.push_back(popStarKLogT[i]);
  }
}
// **************************************************************************

void GalaxyObject::GenerateComplexGalaxiesForCells()
// ****************************************************************
// Generate sersicComplex galaxies randomly (but repeatably) using idl code 
// from John Peterson converted to C++
// **************************************************************
// Do this by dividing celestial sphere into ra and dec cells  Use density and
// area of cell to get  mean number of galaxies in a cell, fluctuate it with 
// appropriate poison distribution and create and write object disk and bulge 
// of galaxies to catalog file.
// ***************************************************************
{
  GalaxySetSEDFileName();//This just flat for now until we find a better way. 

  // *********************************************************************
  // First test that the CK04 (CastilliKurucz2004) and K93 (Kurucz1993) 
  // catalogs exist. If not we will default to the sed_flat SED file is used
  // *********************************************************************

  bool hasPopStarKCatalog = false;

  if ( !pCGParams->fSEDCatalogCollectionPath.empty() ) {
    std::string PopStarKcatalogPath = 
      pCGParams->fSEDCatalogCollectionPath +  "/POPSTARKroupa";
    DIR* dir = opendir(PopStarKcatalogPath.c_str());
    if (dir) {
      /* Directory exists. */
      hasPopStarKCatalog = true;
      closedir(dir);
    }
    else{
      std::cout << "****************************************************"
                <<std::endl;
      std::cout << "* WARNING! Required SED Galaxy catalog for POPSTARKroupa "
                   " missing from " << pCGParams->fSEDCatalogCollectionPath 
                << std::endl;
      std::cout << "* WARNING! Defaulting to all SED's as sed_flat.txt" 
                << std::endl;
      std::cout << "* WARNING! See: https://bitbucket.org/phosim/"
                   "phosim_release/wiki/Using%20PhoSim" << std::endl;
      std::cout << "****************************************************"
                <<std::endl;
        }
  } // End of check that catalog collection Path exist


  // ****************************************************
  // Find all the cells in our FOV (and then some for safety's sake)
  // Places cells in vector fCells
  // ****************************************************
  GetSphericalCellsInFOV(pCGParams->fGalaxiesDiameterFOVDeg);

  // Define a GalaxyComplex object which sets up all the cumulative 
  // probability arrays for galaxy density (and other important things). 

  GalaxyComplex complexGalaxy;
  double densityDeg2 =  complexGalaxy.GetGalaxyDensityDeg2();

  //std::cout << "densityDeg2:" << densityDeg2 <<std::endl;


  GalacticComplexObject disk;
  GalacticComplexObject bulge;

  int numGalaxiesInRange=0;

  // *****************************************************
  // Iterate through the ra/dec cells that were found in the FOV
  int numCells = fCells.size();
  for (int i = 0; i < numCells; i++) {
    // Generate a unique but repeatable random number sequence generator for
    // all the galaxies in this cell
    Random random;
    
    // Determine random seed for this cell. In order to avoid correlations 
    // withing the random number generator always generate a unique set of 
    // random number that MAY be used. Thus we generate a unique set for each 
    // galaxy even if we don;t use them. This keeps the random number 
    // generator honest and repeatable. 
    // I know this is very inefficient but the random generator is not random 
    // enough (free of correlations) to do it any other way that I can think of.
    // This will  insure we can reproduce the same galaxies independent of 
    // (tracking) origin and diameter of the FOV.
    // Set the seed in this cell's  random number generator
    // Note: You find fNumPhiCells in CatalogObject.
    Uint64 cellSeed =   fCells.at(i).fPhiIndex 
      + fCells.at(i).fThetaIndex * fNumPhiCells;
    random.setSeed64(cellSeed);
    
    // Determine for this cell the number of galaxies over the complete 
    // Absolute Magnitude range. We will cut on the Apparent Magnitude later.

    double meanNumGalaxies = densityDeg2 * fCells.at(i).fAreaDeg2;
    //std::cout << i << " " << fCells.at(i).fAreaDeg2 << std::endl;

    //Possion fluctuate this number to get correct statistical number
    int numGalaxies = random.poisson(meanNumGalaxies);
 
    // Now generate these galaxies and if they are within the specified 
    // magnitude range, we write to the disk and bulge objects to the catalog 
    // file
    if (numGalaxies >0 ) {
      double RADeg=0.;
      double decDeg=0.;
      std::vector < double > galaxyUniformRandom;
      std::vector < double > galaxyNormalRandom;
      std::vector < double > galaxySEDRandomUniform;
      std::vector < double > galaxySEDRandomNormal;
      for (int j = 0; j < numGalaxies;  j++ ) {
        fCells.at(i).GenerateCellRandomRADec(RADeg, decDeg, random);
           
        // At this point we need to check that this particular Galaxy is 
        // within the exact FOV. (not the modified). 
        // This saves time and space.
        // double angDistDeg = eraSeps( 
        //                           RADeg * (PAL__DPI/180.), 
        //                           decDeg * (PAL__DPI/180.), 
        //                           pCGParams->fCatalogFOVRADeg *(PAL__DPI/180.), 
        //                           pCGParams->fCatalogFOVDecDeg *(PAL__DPI/180.) ) *
        //                                                          (180./PAL__DPI);
        double angDistDeg = angularSep(RADeg * (PAL__DPI/180.), 
                                  decDeg * (PAL__DPI/180.), 
                                  pCGParams->fCatalogFOVRADeg *(PAL__DPI/180.), 
                                        pCGParams->fCatalogFOVDecDeg *(PAL__DPI/180.) ) *
            (180./PAL__DPI);
        // printf("%lf %lf\n",angDistDeg,angDistDeg2);

        // We need 17 random values for each galaxy we generate in 
        // GenerateGalaxyForCells(and don't generate).
        galaxyUniformRandom.clear();
        for (int i=0; i<17; i++) {
          galaxyUniformRandom.push_back(random.uniform());
        }
        // Also, we need 3 uniform random values and 2 normal distributed 
        // random values for each galaxy's SEDs that will 
        // be generated (or are not generated). Doing it this way insures 
        // reproducability.
        galaxySEDRandomUniform.clear();
        for (int i=0; i<3; i++) {
          galaxySEDRandomUniform.push_back(random.uniform());
        }
        galaxySEDRandomNormal.clear();
        for (int i=0; i<2; i++) {
          galaxySEDRandomNormal.push_back(random.normal() );
        }

        if ( angDistDeg <= pCGParams->fGalaxiesDiameterFOVDeg/2.) { 


          //Now try to generate the galaxy objects (2 of them, a disk and 
          // a bulge). Reject if the galaxy Mapparent magnitude is not in
          // the range. GenerateGalaxyFor Cells does this.
          // ************************

          bulge.RADeg  = RADeg;   
          bulge.DecDeg = decDeg;
          disk.RADeg   = RADeg;
          disk.DecDeg  = decDeg;
   
          bool weHaveAGalaxy = 
            complexGalaxy.GenerateGalaxyForCells( 
                                      bulge, disk,
                                      pCGParams->fGalaxiesMinimumMagnitude,
                                      pCGParams->fGalaxiesMaximumMagnitude,
                                      galaxyUniformRandom, j);
          if( weHaveAGalaxy) {
            pCGParams->fCatalogObjectID++;   //Bump up the object ID
 
            // *************************************************************
            // Get the SEDs for this galaxy. Using the pre-generated random
            // numbers needed to do this. This pre-generation is a little 
            // inefficient but keeps the random numbers reproducible 
            // independent of the users choice of FOV and magnitude limits.
            // **************************************************************
            if ( hasPopStarKCatalog) {
              GenerateGalaxySEDFileNames(galaxySEDRandomUniform, 
                                         galaxySEDRandomNormal, bulge, disk);
            }
            else {
              // No catalogs, use sed_flat.txt
               bulge.SEDFileName = GetSEDFileName();
               disk.SEDFileName  = GetSEDFileName();
            }

            // ********************
            //Fill in the rest 
            std::ostringstream osBulgeID;
            osBulgeID << pCGParams->fCatalogObjectID << ".1";
            bulge.ID=osBulgeID.str();
            complexGalaxy.WriteObject( bulge, fOfs);

            std::ostringstream osDiskID;
            osDiskID << pCGParams->fCatalogObjectID << ".2";
            disk.ID = osDiskID.str();
            complexGalaxy.WriteObject( disk,  fOfs);
              
            numGalaxiesInRange++;
          }   //end of magnitude test
        }  //end of ra/dec in FOV test
	  } //end of galaxies in cell loop
    }
  }   //end of cell loop

  //count the objects we add to catalog.
  //Each galaxy requires 2 objects(disk and bulge)
  std::cout<< "Number of galaxies added to catalog:                        " << numGalaxiesInRange << std::endl;
  
  pCGParams->fCatalogObjectCount += 2*numGalaxiesInRange;
  return;
}
// ***********************************************************************

void GalaxyObject::GenerateGalaxySEDFileNames(
                                std::vector <double>& galaxySEDRandomUniform,
                                std::vector <double>& galaxySEDRandomNormal,
                                GalacticComplexObject& disk,
                                GalacticComplexObject& bulge)
// **********************************************************************
// If the galaxy SED library exists generate a disk and bulge SED file name.
// We use the PopStar-Kroupa( called PopStarK here) catalog. If the synthetic 
// catalog PopStar-Kroupa is unavailable default to the use of the flat SED. 
// ************************************************************************
// PopStar-Kroupa found at :
// http://svo2.cab.inta-csic.es/theory/newov2/index.php?models=pop_kroupa
// *********
// POPSTAR with Kroupa IMF
// PopStar Evolutionary synthesis models. Using IMF from Kroupa (2002). This 
// grid of Single Stellar Populations covers a wide range in both, age and 
// metallicity. The models use the most recent evolutionary tracks together 
// with the use of new NLTE atmosphere models.
// Kroupa P. 2002, ASPC in Modes of Star Formation and the Origin of Field 
// Populations, 285, 86
// Molla et al 2009, MNRAS, 398,451
// *************************************************************************
// This will be a very simple "toy" spectra generation model. No claim to any 
// great accuracy is made here. We randomly pick an age following the  Madau 
// SFR distribution (within redShift limits). This makes a little sense but 
// is obviously just a toy model. Metallicities for the disk 
// and bulge of a galaxy are chosen from a Gaussian centered at <[Fe/H]> = .1
// We then build valid PopStarK SED file names from these values.
// ******************************************************************
// We can also very very crudely account for disk and bulge differences in 
// ages and metallicity. 
// From: "Ages of Galaxy Bulges and Disks from
// Optical and Near-Infrared Colours", Peletier,R.,F., Balcells,M., 1996,
// arXiv:astro-ph/9509123v2   
// "For the sample, we find that the average difference in log Z will 
// be ~0.12 +/-   0.16 (1 sigma), and ~0.11 +/- 0.15 in log(t)."
// With either: "bulge is at most 30% older than its disk", or  disk should 
// be younger than its bulge. 
// Because our SEDS have very crude steps in Metallicity (probably won't see
// a difference) and fairly fine steps in age we chose to vary the LogT age
// for the bulge with it being the disk LogT age plus a normally distributed 
// offset with a mean .11 and sigma of  .15. 
// **************************************************************************
{
  // ***********************************************************************
  // First we have to generate this galaxies age. We will assume that the z
  // of the absolute effective birth of this galaxy follows the Madau Star 
  // Formation Rate (SFR)  curve of the universe. Get Madau plot from:
  // "Cosmic Star Formation History",Piero Madau, Mark Dickinson, 
  // arXiv:1403.0007,pg 48,49.
  // The text says we can divide this curve into 3 regions:
  // "cosmic SFH:a rising phase, scaling as SRF(z) ~(1 + z)**-2.9  
  // at 3 < z < 8, slowing and peaking at some  point probably between z = 2 
  // and 1.5, when the Universe was ~ 3.5 Gyr old, followed by a gradual 
  // decline to the present day, roughly as SFR(z) ∝ (1 + z)**2.7."
  // Using these formulas we can randomly pick a z(birth) that follows the 
  // Madau SFR distribution by using the standard inversion method for the 
  // use of random number to generate a particular distribution.
  // ***********************************************************************
  // Find region
  double zGalaxy = bulge.z;
  double zT0;
  double u1;
  double u2;
  
  double SFR3Area = (  pow((1.0 + 8.0),(1.0 - 2.9)) - 
                       pow((1.0 + kMadauHighMin),(1.0 - 2.9)) ) /(1.0-2.9);
  double SFR1FinalValue = pow((1.0 + 1.5 ),(+2.7));
  double SFR3InitialValue = pow((1.0 + kMadauHighMin ),(-2.9));

  if (zGalaxy >= kMadauHighMin) {
    // We are in the third region.  The zT0 (birth) of this galaxy which we
    // get from the SFR curve should be greater than zGalaxy. We will thus 
    // need the area of zGalaxy < zT0 <8.   SRF ~ (1+z)**(-2.9)

    u1 = galaxySEDRandomUniform.at(0);
    zT0 = GetRandomRedShiftFromSRF(-2.9, zGalaxy, 8.0, u1);
  }
  else if(zGalaxy >= 1.5) {
    //Middle SRF region is approximated as flat. 1.5<z<kMadauHighMin
    // But we start at zGalaxy in middle region of SRF
    // Get its area. Make it match start of region 3 (z=kMadauHighMin)
    // We need relative areas.
    double SFR2Area =  SFR3InitialValue * (kMadauHighMin - zGalaxy);
    
    //thow random to see which area we are in.
    u1 = galaxySEDRandomUniform.at(0);
    u2 = galaxySEDRandomUniform.at(1);
    if ( u1  < (SFR2Area / (SFR2Area + SFR3Area)) ) {
      //In middle region.
      zT0 =zGalaxy + u2 *(kMadauHighMin-zGalaxy);
    }
    else {
      //; In region 3
      zT0 = GetRandomRedShiftFromSRF(-2.9, kMadauHighMin, 8.0, u2);
    }
  }
  else {
    // ***********************
    // zGalaxy is in the first region. Now it gets tricky to determine which 
    // region the zT0 is in.  We will need the relative areas of all 3 regions
    // (with the first going from zGalaxy to 1.5). To do that they should 
    // match at z=1.5 and z=kMadauHighMin and thus should scale relative to 
    // one of the regions (lets say 1)
    double SFR1Area = ( pow((1.0 + 1.5),(1.0 + 2.7)) - 
                        pow((1.0 + zGalaxy),(1.0 + 2.7)) ) / (1+2.7);
    double SFR2Area = SFR1FinalValue * (kMadauHighMin - 1.5);

    // Scale SFR3Area:
    double SFR3AreaScaled = SFR3Area *(SFR1FinalValue / SFR3InitialValue);

    //Now pick where we are.
    u1 = galaxySEDRandomUniform.at(0);
    u2 = galaxySEDRandomUniform.at(1);
    if (u1 >  (SFR1Area + SFR2Area) / (SFR1Area + SFR2Area + SFR3AreaScaled) ) {
      //In region 3
      zT0 = GetRandomRedShiftFromSRF(-2.9, kMadauHighMin, 8.0, u2);    
    }
    else if (u1 > (SFR1Area / (SFR1Area + SFR2Area + SFR3AreaScaled) ) ) {
      //in region 2
      zT0 =1.5 + u2 *(kMadauHighMin-1.5);        
    }
    else {
      //In region 1 (but start at zGalaxy)
      zT0 = GetRandomRedShiftFromSRF(2.7, zGalaxy, 1.5, u2);    
    }
  }
  // Now convert this to an age.
  // First: time since big bang  at zGalaxy
  double tzGalaxy = kAgeUniverseYears/(1.0 + zGalaxy);
  
  // and time since big bang of the galaxy at  zT0;
  double tzT0 =kAgeUniverseYears/(1.0 + zT0); 
  
  // and finally age of galaxy at time zGalaxy:
  double ageGalaxy = tzGalaxy - tzT0;
  
  // We need log10 of the age
  // Max value of zT0 should be 8, Min value of zGalaxy = 0.
  // So max value of ageGalaxy should be 12.32E9 years or logt = 10.09
  // Min value of age would be 0. We need to limit it to 100,000 log10T=5.0
  // Acceptable values in PopStarK catalog are 5.0 to 10.18
  // We need to limit it to 100,000 log10T=5.0
  double diskLogT =  log10(ageGalaxy);
  if (diskLogT < 5.0) {
    diskLogT = 5.0;
  }

  // Now create the Bulge LogT which is norm distributed diskLogT + mean=.12
  // +/- .15
  double diskBulgeOffsetLogT = 
    (galaxySEDRandomNormal.at(0) * kBulgeOffsetLogTSigma) + 
    kBulgeOffsetLogTMean;
  double  bulgeLogT = diskLogT +  diskBulgeOffsetLogT;
  if (bulgeLogT < 5.0) {
    bulgeLogT = 5.0;
  }
  if (bulgeLogT > 10.18) {
    bulgeLogT = 10.18;
  }

  // ******************************************************************
  // Now generate a metallicity.  
  // 
  // Allowed PopStarK range of metallicity Zgalaxy is .0001 <= Zgalaxy <= 0.05;
  // The [Fe/H] value is related to the metallicity Z:
  // [Fe/H] = log10(Z/0.0196):"New solar metallicity measurements",Sunny 
  // Vagnozzi, arXiv:1703.10834 
  // If we approximate(very very crudely) the the distribution of metallicity
  // values as Gaussian with a log(Z/Zsolar) = .1 => Z=.02467
  // with a width of ~ Sigma(log(Z/Zsolar))  = .1
  // our max acceptable would be Z=.05 or log(Z/Zsolar)=.4 Ie 4 sigma
  // Min acceptable would be Z=.0001 or log(Z/Zsolar)= -2.29 ~20 sigma
  // So we will have to keep an eye on the range of the Z we generate.
  // This distribution is from(see pg 11): 
  // "The ages and metallicities of galaxies in the local universe",
  // Gallazzi, Anna, et al.,Mon. Not. R. Astron. Soc. 362, 41–58 (2005),
  // arXiv:astro-ph/0506539v1, pg 11
  // This is again a very crude model. Just really to get something thats not 
  // completly bonkers.
  // *******************************************************************
  // Note the random value here is "normally" distributed pre calculated for 
  // reproducability).
  double metallicity = galaxySEDRandomNormal.at(1) * kGalacticMetallicitySigma
                                    + kGalacticMetallicityMean;
  // Convert to Z 
  // Note: Convention. metallicity is an upper case "Z"  the lower case "z" 
  // is redShift.
  double Z = pow(10, metallicity)*kZsolar;

  // Debug
  //std::cout << logT << " " << Z << " " << zT0 << std::endl;

  // ******************************************************************
  // We are now ready to compose the galaxy SED.
  // ***************************************************
  // We now have an Age and a Metallicity (logT,Z) for this galaxy (disk). 
  // These values are to be encoded within the closets PopStarK SED file name.
  // SED files for the PopStarK catalog that we will use all start with 
  // "kroupa_spneb_kro_0.15_100_". The 0.15_100 is the mass range used in the 
  // IMF.
  // ***************************************************
  BuildGalaxySEDName("PopStarK", diskLogT, bulgeLogT, Z, bulge.SEDFileName, 
                     disk.SEDFileName);

  return;
}

// *********************************************************************

double GalaxyObject::GetRandomRedShiftFromSRF(double slopeIndex, 
                                              double zGalaxy,
                                              double zMax, double u1)
// ********************************************************************
// Pick a z between zGalaxy and zMax from a function that goes like
// (1+z)**(slopeIndex) using random number u1.
// ********************************************************************
{
  // Integral of (1+z)**(slopeIndex) is proportional to  (1_z)**(1+slopeIndex)
  // Evaluate at zGalaxy and ZMax
  double SFRInitialAreaValue = pow((1.0 +zGalaxy ),(1.0 + slopeIndex));
  double SFRFinalAreaValue = pow((1.0 +zMax ),(1.0 + slopeIndex));
  double SFRArea = SFRFinalAreaValue -  SFRInitialAreaValue;
                                     
  // Find z where area fraction from ZGalaxy to z is u1:
  double SRF = SFRInitialAreaValue + ( u1 * SFRArea);
  double zT0 = pow(SRF,(1/(1.0+slopeIndex)) ) - 1.0;
  return zT0;
}

// ***********************************************************************
void GalaxyObject::BuildGalaxySEDName(std::string Catalog, double DiskLogT,
                                      double BulgeLogT,
                                      double Metallicity, 
                                      std::string& BulgeSEDName,
                                      std::string& DiskSEDName)
// *********************************************************************
// Generate the SED file name (includes path) to use as the SED
// for the star with given values of the parameters.
// Only PopStarK (PopStar Kroupa) catalog implemented for now.
// If the parameters are out of range for the selected catalog, bring them to
// the range. Thus we should always get a valid SED. This is checked(SED 
// file exists).
// **********************************************************************
// File names in PopStarK catalog have the form:
// kroupa_spneb_kro_0.15_100_z"Metallicity"_t"LogT"_total.dat.txt.gz
// Notes:
// 1. We use the IMF made with the range of star masses of .15 to 100 solar m.
// 2. We use the total (star + nebular) SED.
// Example: kroupa_spneb_kro_0.15_100_z0500_t9.00_total.dat.txt.gz
// This has Z = .05 and logT = 9.00
// ***********************************************************************
{
  //Set up vectors
  std:: vector < double >* pMetallicity = &fPopStarKZ;//Just init to something.
  std:: vector < double >* pLogT = &fPopStarKLogT;
  std::string catalogSpecifier;
  std::string catalogDirectory;
  std::string popStarKStr = "PopStarK";
  if ( Catalog != popStarKStr ) {
    std::cout << "Fatal--Specified Catalog: " << Catalog << " not:" 
              << popStarKStr << std::endl;
    exit(1);
  }
  // *********************************************************************
  // Load up pointers to allowed values. Vectors setup in constructor of 
  // this class. (Do it this inefficient way for now in case we add other 
  // catalogs later)
  // Note we have previously checked that a legal catalog has been chosen
  // ********************************************************************* 
  if (Catalog == popStarKStr ) {
    pMetallicity = &fPopStarKZ; 
    pLogT        = &fPopStarKLogT;
    catalogSpecifier  = "kroupa_spneb_kro_0.15_100";
    catalogDirectory  = "POPSTARKroupa";
  }
  else{
    // Should never get here.
    std::cout << " Fatal -- GalaxyObject::BuildGalaxySEDName : Illegal "
                 " catalog selected: " << Catalog << " Should never get here!" 
              << std::endl;
    exit(1);
  }

  // *********************************************************************
  // Build file name and see if it exists.
  // *************************************
  // Metallicity:
  // The value  in the SED file name is encoded as the Metallicity * 10000
  // Only the ordered values in vector pointed to by metallicity 
  // are allowed. ClosestElement will set any out range value to nearest 
  // valid value (and not tell you about it!).
  int closestMeta = (int)(ClosestElement( pMetallicity, Metallicity) *10000);
   
  // LogT:
  // The allowed LogT  values are listed in sorted vector pLogT
  // ******************************************************************
  // We can also very very crudely account for disk and bulge differences in 
  // ages and metallicity. From: "Ages of Galaxy Bulges and Disks from
  // Optical and Near-Infrared Colours", Peletier,R.,F., Balcells,M., 1996,
  // arXiv:astro-ph/9509123v2   
  // "For the sample, we find that the average difference in log Z will 
  // be ~0.12 +/-   0.16 (1 sigma), and ~0.11 +/- 0.15 in log(t)."
  // With either: "bulge is at most 30% older than its disk", or  disk should 
  // be younger than its bulge. 
  // Because our SEDS have very crude steps in Metallicity (probably we 
  // wouldn't see a difference) and we have fairly fine steps in age. We '
  // will choose an older age for the bulge from  a normal distribution , 
  // mean .11 and sigma .15.
  // ******************************************************************
 
  double closestDiskLogT  = ClosestElement( pLogT, DiskLogT);
  double closestBulgeLogT = ClosestElement( pLogT, BulgeLogT);

  //Debug
  //std:: cout << closestMeta      << "\t" << Metallicity << "\t" 
  //           << closestDiskLogT  << "\t" << DiskLogT    << "\t"
  //           << closestBulgeLogT << "\t" << BulgeLogT   << std::endl;

  // Build the SED file names
  // Disk first
  DiskSEDName  = ConstructSED(catalogDirectory, catalogSpecifier, closestMeta,
                              closestDiskLogT);
  BulgeSEDName = ConstructSED(catalogDirectory, catalogSpecifier, closestMeta,
                              closestBulgeLogT);
  return;  
}
// **************************************************************************

std::string GalaxyObject::ConstructSED(std::string catalogDirectory, 
                                       std::string catalogSpecifier, 
                                       int closestMeta, double closestLogT)
{
  std::ostringstream SEDss;
  SEDss.str("");
  SEDss << catalogDirectory << "/"
        << catalogSpecifier << "_z" << std::setw(4) << std::setfill('0') 
        << closestMeta << "_t" << std::fixed << std::setprecision(2) 
        << closestLogT << "_total.dat.txt.gz";
  // Now we have the file name. Check that such a file exists in the 
  // this catalog.
  std::string SEDName = "./" + SEDss.str();  //This is what raytrace wants.
    
  std::string fullSEDName = pCGParams->fSEDCatalogCollectionPath + "/" + 
    SEDss.str();
  std::ifstream infile(fullSEDName.c_str());
  if ( !infile.good()) {
    // A SED file of this name does not exist
    // Should never get here.
    std::cout<< "Fatal-A non-existant SED file name was created: " 
             << fullSEDName << std::endl;
    exit(1);
  }
  infile.close();
  return SEDName;
}
