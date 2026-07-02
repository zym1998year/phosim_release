///
/// @package phosim
/// @file GalaxyComplex.cpp
/// @brief Generates phosim sersicComplex galaxy parameters randomly (but 
///   repeatably) using IDL code provided by John Peterson converted to C++
// **************************************************************
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
/// Original algorithm from IDL code provided by John Peterson.

#include "galaxycomplex.h"


// *********************************************************************
GalaxyComplex::GalaxyComplex()
{
  Init();
}

void GalaxyComplex::Init()
{
  // MagnitueValues and zValues are lookup tables for more efficiency.
  fMabsoluteValues.resize( kNumValues, 0.0);
  fCumulativeMabsoluteProb.resize( kNumValues, 0.0);

  std:: vector < double > densityMabsolute;
  densityMabsolute.resize( kNumValues, 0.0);

  std:: vector < double > MabsoluteProb;
  MabsoluteProb.resize( kNumValues, 0.0);


  double mmin = -17.0;   //This is an lowest absolute magnitude we are 
                         //interested in
  double densityMabsoluteTotal = 0;
  for(int i = 0; i < kNumValues; i++) {
    // Create the Absolute Magnitude array. This has as M values starting at 
    // M=-25.  Has steps of .00001 M. Last value is M~=-17. Milky way is -20.8
    // Mabsolute. A look up table, sort of.
	fMabsoluteValues.at(i) = -25.0 + ((double)i / (double) kNumValues ) *
	                                                          (mmin + 25.0);
	
	// determine densityMabsolute.  It is distributed as a function of Mabsolute
    double ppower =  0.4 * ( -20.3 - fMabsoluteValues.at(i));
	densityMabsolute.at(i)= 0.4 * log(10.0) * 0.005 * 
			 pow(10., (ppower * ( -1.2 + 1.0))  ) * exp(- pow(10.,  ppower )) ;
    //Total density over Mabsolute range.
    densityMabsoluteTotal += densityMabsolute.at(i); 
  }

  // Cumulative Mabsolute Probability distribution
  // **********************************************************
  // Normalize densityMabsolute. Makes it into a probability. Create 
  // cumulative densityMabsolute prob. This will be used to find a 
  // galaxy with a random Mabsolute. 
  // **********************************************************
  fCumulativeMabsoluteProb.at(0) = MabsoluteProb.at(0);  //Initial value
  for( int i = 1; i < kNumValues; i++) { 
	MabsoluteProb.at(i) = densityMabsolute.at(i)/densityMabsoluteTotal;
    fCumulativeMabsoluteProb.at(i) = fCumulativeMabsoluteProb.at(i-1) + 
                                                             MabsoluteProb.at(i);
  }
  // ************************************************************

  //Redshift Distribution
  fZValues.resize(kNumValues,0.0);
  fCumulativeZProb.resize(kNumValues,0.0);

  std::vector < double > zProb;
  zProb.resize(kNumValues,0.0);


  double zc = 1.0;  //Some type of z scale?

  double totalZProb = 0.0;
  for (int i=0; i<kNumValues; i++) {
	//z=5 is max z we are interested in.
	fZValues.at(i) = ((double)i/(double)kNumValues)*5.0;
	zProb.at(i) = pow( fZValues.at(i), 2) * exp(-fZValues.at(i) / zc); 
	totalZProb += zProb.at(i);
  }

  //Cumulative Redshift Probability Distribution
  fCumulativeZProb.at(0) = zProb.at(0);
  for (int i=1; i< kNumValues ; i ++ ) {
	zProb.at(i) = zProb.at(i) /	totalZProb;
    fCumulativeZProb.at(i) = fCumulativeZProb.at(i-1) + zProb.at(i);
  }

  // Density of galaxies
  // ***********************************************************************
  // This version is for the spherical cell algorithm (cells may be of 
  // differing areas ( in deg**2). So leave inclusion of areas to later)
  // To get a number in a cell we will apply area in deg2 later. So now add 
  // conversions to get denGalaxyDeg2 the correct units(the 
  // "pow(60., 2)" deg**2 to arcmnin**2 term in denGalaxyDeg2).
  // Don't forget John's fudge factor which he says we will address later.
  // Note: We use a double for the densityMabsolute here. We calculate 
  // complete number later (we will multiply by an area in deg2 and fluctuate 
  // appropriately (Poisson).).
  // ***********************************************************************
  // Breaking up density calculation into terms for computation clarity.
  // But where the breaks should be physics wise (ie what terms mean what) 
  // is unknown to me

  double term1 = ( mmin + 25.0 ) / kNumValues;
  double term2 = 4. * kPI * pow( (3e5 / 70.0), 3) ;
  double term3 = ( ( 4. * kPI) / pow( (kPI/180.), 2) )  * 60.0 * 60.0;
  double denGalaxyDeg2 = densityMabsoluteTotal * term1 * (term2 / term3) * 
                                      pow(60., 2);
  
  // I think  this is densityMabsolute of galaxies for the complete possible 
  // Mabsolute range!
  fGalaxyDensityDeg2 = ( denGalaxyDeg2 * (totalZProb/kNumValues) * 5.0 );  

  return;
}
// ***************************************************************************

std::string GalaxyComplex::GenerateSersicComplexObjectString(
												 GalacticComplexObject& Object)
{
  // **********************************************************************
  // Generate string that goes into phosim object catalog for the 
  // sersicComplex part
  // **********************************************************************
  std::ostringstream os;
  os << "sersicComplex" 
	 << std::fixed << std::setw(9) << std::setprecision(4)
	 << Object.a << " " << Object.b << " " << Object.c << " "
	 << Object.alpha << " " << Object.beta << " " << Object.gamma << " " 
     << Object.index << " " << Object.clumpFrac << " " << Object.clump << " " 
	 << Object.clumpWidth << " " << Object.spiralFrac << " " 
	 << Object.alphaSpiral << " " << Object.bar << " "
	 << Object.spiralWidth << " " << Object.phi0 << " ";

  return os.str();
}
// ************************************************************************

void GalaxyComplex::WriteObject( GalacticComplexObject& Object, 
								 std::ofstream* ofs)
// ************************************************************************
// Write the object line to the catalog file
// ************************************************************************
{
  // *********************************
  // Get the sersic part of the line
  // *********************************
  std::string sersicCmplxStr=GenerateSersicComplexObjectString( Object);

   *ofs<< "object " << Object.ID << " " << std::fixed 
       << std::setprecision(8) << Object.RADeg  << " " <<Object.DecDeg
       << " " << std::setprecision(2) << Object.M << " " 
       << Object.SEDFileName << " " << std::setprecision(3) << Object.z 
       << " 0.0 0.0 0.0 0.0 0.0 " << sersicCmplxStr << "none none " 
       << std::endl;
   return;
}
// *************************************************************************

bool GalaxyComplex::GenerateGalaxyForCells(GalacticComplexObject& Bulge, 
                                           GalacticComplexObject& Disk,
                                           double MinimumM, double MaximumM,
                                 std::vector < double >& galaxyUniformRandom,
                                           Uint64 galaxyIndex)
// ****************************************************************************
// This method generates a Mapparent for a galaxy and if within range of 
// Mminimum -Mmaximum it creates the Magnitudes, sizes and sersicComplex  
// parameters for a random galaxy disk and bulge following simple 
// distributions of such.
// ****************************************************************************
// galaxyIndex used for debug prints only
// ****************************************************************************
// Note: a 17 value random sequence is used that was generated in the calling 
// program. Note the 17 values was hardwired to make sure we don't screw this 
// up. 
// and put into the  galaxyUniformRandom vector
// ****************************************************************************
// GalacticComplexObject defined in GalaxyObject class
// ****************************************************************************
// Original algorithm from IDL code provided by John Peterson.
{
  // ***************************************************************
  // To determine the Mapparent of a galaxy with this method we need at the 
  // minimum Mabsolute and the distance.
  // **************************************************************
  // For distance first pick a z from our z distribution.
  // Find index to z in the z arrays. Find index where the 
  // cumulative probability exceeds the random selection  
  
  int indexRandom = 0;
  // From:http://www-numi.fnal.gov/offline_software/srt_public_context/
  //  WebDocs/Companion/cxx_crib/increment.html
  // "If ++ follows the variable, e.g. counter++, the value returned is the 
  // value in counter before it has been incremented."

  double u1 = galaxyUniformRandom.at(indexRandom++); 

  // To find the bin in a Cumulative probability distribution vector
  // (naturally sorted with first bin >= 0 (and maybe a few more as 0) and with
  // the last bin < 1.0) which holds the random 
  // probability u1 we will use the C++ standard library algorithm 
  // "std::upper_bound". upper_bound does a binary search:
  // upper_bound(ForwardIterator first, ForwardIterator last, const T& val);
  // "Returns an iterator pointing to the first element in the range 
  // (first,last) which compares greater than val" ie. its the first bin who's 
  // value is   > val. This works for all values of u1 even u1 == 0.0.
  // Note: We start range at bin 1 not 0. Why? Think about it and you will see.

  fVectorIteratorUpper = std::upper_bound(fCumulativeZProb.begin()+1, 
                                          fCumulativeZProb.end(), u1);
  //Find bin index for this value of u1. Back up one bin
  int g1 = fVectorIteratorUpper - fCumulativeZProb.begin() - 1; 
  
  // Thus:
  double redShift = fZValues.at(g1);     //redShift in range: 0 -> 5.0
  double distance=3e5/(70.0)*redShift;  //Units?


  // Now pick a random Mabsolute magnitude from its distribution
  double u2 = galaxyUniformRandom.at(indexRandom++);

  fVectorIteratorUpper = std::upper_bound(fCumulativeMabsoluteProb.begin()+1, 
                                          fCumulativeMabsoluteProb.end(), u2);
  //back up 1 bin
  int g2 = fVectorIteratorUpper - fCumulativeMabsoluteProb.begin() - 1;

  //Get absolute magnitude of galaxy
  double MAbsolute = fMabsoluteValues.at(g2);

  // Correct for distance. Generate total Apparent magnitude
  double MApparent = MAbsolute + (5.0 * log10( distance * 1e6/10.0 ) ); 

  // Debug
  //std::cout << redShift << "   " << distance << "  " << u1 << " " << u2 
  //          << " " << MApparent << " " << MAbsolute << " " << Bulge.RADeg 
  //          << " " << Bulge.DecDeg << " " << galaxyIndex << std::endl; 
 
  //Now is the soonest we can check that this MApparent is within the range
  // we specified

  if (MApparent >MaximumM || MApparent < MinimumM) {
    return false;
  }
  else {

    // OK its a good one! Make the bulge and disk objects!
    // Determine total flux and then convert to disk and bulge magnitude
    // Now the angular size of the galaxy
    double ssize =  pow( (pow( 10.0, (-0.4*(MAbsolute+21.0)))  ), .333) *0.01 /
      ( distance * (kPI /(180.*3600.0) ) );  //Units? (arcsec?)
    
    
    //Divide up MApparent between disk and bulge randomly
    double bulgefrac =  galaxyUniformRandom.at(indexRandom++);
    
    Bulge.M = -2.5 * log10(    bulgefrac   * pow(10., (-0.4*MApparent)) );
    Disk.M  = -2.5 * log10( (1.-bulgefrac) * pow(10., (-0.4*MApparent)) );
    Bulge.z = redShift;
    Disk.z  = redShift;
    
    //debug
    //std::cout << MApparent <<" " << Bulge.M <<" " << Disk.M << std::endl;
    
    // Now fill in the sersic values for sersicComplex.
    Bulge.a=ssize*(0.8+0.2*galaxyUniformRandom.at(indexRandom++))*0.75;
    Bulge.b=ssize*(0.8+0.2*galaxyUniformRandom.at(indexRandom++))*0.75;
    Bulge.c=ssize*(0.8+0.2*galaxyUniformRandom.at(indexRandom++))*0.75;
                                               // Polar angle deg
    Bulge.alpha = galaxyUniformRandom.at(indexRandom++) * 360.;
                                               //position angle deg
    Bulge.beta  = ( acos( 2.0 * galaxyUniformRandom.at(indexRandom++) - 1.0) 
                           - (kPI/2.0) ) * 180. / kPI; 
                                               // Rotation angle deg
    Bulge.gamma = galaxyUniformRandom.at(indexRandom++) * 360.;
    Bulge.index = 4;
    Bulge.clumpFrac   = galaxyUniformRandom.at(indexRandom++) * 0.5; 
    Bulge.clump       = 300;
    Bulge.clumpWidth  = ssize*0.1;
    Bulge.spiralFrac  = 0;
    Bulge.alphaSpiral = 0;
    Bulge.bar         = 0;
    Bulge.spiralWidth = 0;
    Bulge.phi0        = 0;
                                                           //axis 1 in arcsec
    Disk.a=ssize*(0.8+0.2*galaxyUniformRandom.at(indexRandom++))*1.25; 
                                                           //axis 2 in arcsec
    Disk.b=ssize*(0.8+0.2*galaxyUniformRandom.at(indexRandom++))*1.25; 
                                                           //axis 3 in arcsec
    Disk.c=ssize*(0.2*galaxyUniformRandom.at(indexRandom++))*1.25;
    Disk.alpha = Bulge.alpha;        //polar angle in deg, same as bulge
    Disk.beta  = Bulge.beta;         //position angle deg, same as bulge
    Disk.gamma = Bulge.gamma;        //Rotation angle deg, same as bulge
    Disk.index = 1;                                       //sersic index
                                                   //fraction in clumps
    Disk.clumpFrac   = galaxyUniformRandom.at(indexRandom++) * 0.5;
    Disk.clump       = 300;                        //number of clumps
    Disk.clumpWidth  = ssize*0.1;                  //gaussian clump size arcsec
                                           //fraction in spiral
    Disk.spiralFrac  = galaxyUniformRandom.at(indexRandom++);
                                           //winding angle of spiral in deg
    Disk.alphaSpiral = galaxyUniformRandom.at(indexRandom++) * 60.0 - 30.0;
    Disk.bar         = (ssize*0.5);        //spiral bar size in arcsec
    Disk.spiralWidth = (ssize*0.1);        //spiral Gaussian width in arcsec
                                           //position angle of spiral in deg.
    Disk.phi0        = galaxyUniformRandom.at(indexRandom++) * 360.0 ;
  
    //Debug
    // Bulge.M = 18.;
    //Disk.M = 16.0;
    //Bulge.a = 4.5;
    //Bulge.b = 4.5;
    //Bulge.c = 4.5;
    //Bulge.alpha = 0.0;
    //Bulge.beta =  90.0;

    //Disk.a = 7.8;
    //Disk.b = 7.8;
    //Disk.c = 0.6;
    //Disk.alpha = 0.0;
    //Disk.beta =  90.0;

    
    return true; //return total Magnitude of galaxy. Used to select M range
  }
}
// ************************************************************************

