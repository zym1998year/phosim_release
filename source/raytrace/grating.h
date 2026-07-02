/// @brief Grating Class
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#ifndef GRATING_H
#define GRATING_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <utility>

#include "constants.h"   //Defines DEGREE (degree to radians)
                         // This is a Phosim file in soucre/raytrace
#include "../ancillary/random.h"

// ****************
// Grating definition: Defaults.(mostly used in testing)
// See for intensity tutorial:
// ht7.09tp://hyperphysics.phy-astr.gsu.edu/hbase/phyopt/gratint.html
// ****************
double const kBlazeAngleDeg = 7.09;// Angle between face of groove and 
                                   // plane of grating
//int    const kNumSlits = 30;    
int    const kNumSlits = 600;      // Number of illuminated grooves in the 
                                   // grating/object 
double const kDeltaNm  = 3333.33;  // Distance between grooves  (nm). Inverse 
                                   // is the groove density

int    const kTableVersion = 1;
// *****************
// Cumulative Probablity table indexes parameters
// See Grating constructor methods for wavelength range and number of steps
// specifications.
// ***********************************************
// These values have been optimized for table length, time of construction,
// and table accuracy. Change at your own peril
double const kWavelengthStepSizeNm = .25;
double const kAngleInStepSizeDeg = 180.0/10000.0;
double const kAngleOutStepSizeDeg = 180.0/30000.0;
double const kAngleOutProbabilityStep = .00005; //Min probability steps in 
// Set value used to test for possible table failure
double const kMaxAllowedAngleOutDiffDeg = 5.0;  //Aprox sigma of angleout 
                                                //peaks
// ********************

// ********************
// Defaults: 
// Using values different that these (which is perfectly acceptable and in 
// most case required) will change table length and time of construction.
double const kAngleInMinDeg = -2.0;
double const kAngleInMaxDeg = 2.0;
double const kAngleOutMinDeg = -45.0;
double const kAngleOutMaxDeg = +45.0;
// *************************************************************************

class Grating 
// **********************************************************************
// Class to determine photon diffraction from a diffraction grating
// **********************************************************************
{
 public:
  // Use defualt values for grating parameters
  Grating(std::string instrDir, Random* pRandm, double Minwavelength, 
          double Maxwavelength);
  Grating(std::string instrDir, Random* pRandm, double Minwavelength, 
          double Maxwavelength,   int TestToEnable);

  // Use values from "surface." structure. Was read from "optics_*.txt" file.
  // Placed in a vector to avoid circular Grating/Surface definitions
  Grating(std::string instrDir, Random* pRandm, 
          std::vector < double > &  gratingParams ); 
 
  ~Grating();
  
  void setupGrating(Random* pRandm, double Minwavelength, double Maxwavelength);
  void init(std::string InstrDir);
  
  void diffract(double vxIn, double vyIn, double vzIn,
                double& vxOut, double& vyOut, double& vzOut,
                double wavelengthNm);
  
  
  int    testToEnable;

  //structs used for reading and writing to cumulative lookup table file
  //The grating definition
  struct gratingDefinition {
    double fAngleBlazeRad;
    double fDeltaNm;
    int    fNumSlits;
  public:
    bool cmp(const gratingDefinition& g) {
      if ( fAngleBlazeRad == g.fAngleBlazeRad &&
           fDeltaNm       == g.fDeltaNm &&
           fNumSlits      == g.fNumSlits ) {
        return true;
      }
      else {
        return false;
      }
    }
    void print() {
        std::cout << "    Grating definition:" << std::endl;
        std::cout << "      Blaze Angle(Deg): " << fAngleBlazeRad / DEGREE 
                  << std::endl;
        std::cout << "      Slit spacing(nm): " << fDeltaNm << std::endl;
        std::cout << "      Number slits illuminated: " << fNumSlits 
                  << std::endl;
    }
  };


  gratingDefinition fGDef;      // Grating definition
  gratingDefinition fGTestDef;  // Grating definition as read from file

  //The table structure definition
  struct tableDefinition {
    int    fTableVersion;
    double fWavelengthMinNm ;
    double fWavelengthMaxNm ;
    double fWavelengthStepSizeNm;
    int    fNumWavelengthSteps;

    double fAngleInMinRad;
    double fAngleInMaxRad;
    double fAngleInStepSizeRad;
    int    fNumAngleInSteps;

    double fAngleOutMinRad;
    double fAngleOutMaxRad;
    double fAngleOutStepSizeRad;
    int    fNumAngleOutSteps;
    double fAngleOutProbabilityStep;

  public:
    bool cmp(const tableDefinition& t) {
      if (fTableVersion            == t.fTableVersion &&
          fWavelengthMinNm         == t.fWavelengthMinNm &&
          fWavelengthMaxNm         == t.fWavelengthMaxNm &&
          fWavelengthStepSizeNm    == t.fWavelengthStepSizeNm &&
          fNumWavelengthSteps      == t.fNumWavelengthSteps &&
          fAngleInMinRad           == t.fAngleInMinRad &&
          fAngleInMaxRad           == t.fAngleInMaxRad &&
          fAngleInStepSizeRad      == t.fAngleInStepSizeRad &&
          fNumAngleInSteps         == t.fNumAngleInSteps &&
          fAngleOutMinRad          == t.fAngleOutMinRad &&
          fAngleOutMaxRad          == t.fAngleOutMaxRad &&
          fAngleOutStepSizeRad     == t.fAngleOutStepSizeRad &&
          fNumAngleOutSteps        == t.fNumAngleOutSteps && 
          fAngleOutProbabilityStep == t.fAngleOutProbabilityStep) {
        return true;
      }
      else {
        return false;
      }
    }
    void print() {
        std::cout << "    Grating Table V." << fTableVersion
 << " Limits:" 
                  << std::endl;
        std::cout << "      " <<  fWavelengthMinNm << " < Wavelength(nm) < " 
                  << fWavelengthMaxNm << std::endl;
        std::cout << "      " <<  fAngleInMinRad / DEGREE
                  << " < AngleIn(deg) < " << fAngleInMaxRad / DEGREE 
                  << std::endl;
        std::cout << "      " <<  fAngleOutMinRad / DEGREE
                  << " < AngleOut(deg) < " << fAngleOutMaxRad / DEGREE 
                  << std::endl;
    }
  };
  tableDefinition fTDef;     // Table defintion
  tableDefinition fTTestDef; // Table definition as read from file
    
  //Table management methods:
  std::string fGratingTableFileName = "/grating.tbl";
  void generateCumulativeProbTableParameters();
  void defineCumulativeProbTable();
  void makeCumulativeProbTable();
  bool loadCumulativeProbTable(std::string InstrDir);
  void saveCumulativeProbTableToFile(std::string InstrDir);

  //Methods(Those that start with "_" are for original Yu Dong algorithm)
  double calculateFunctionTransmissionGrating(double AngleInRad, 
                                              double AngleOutRad,
                                              double WavelengthNm);
  //A transmission intensity sub functions
  double calculateInterferenceIntensity(double AngleInRad, 
                                        double AngleOutRad,
                                        double WavelengthNm);
  double calculateSingleSlitEnvelope(double AngleInRad, 
                                     double AngleOutRad,
                                     double WavelengthNm);

  double calculateSpecificAngleOutRad( double WavelengthNM, double AngleInRad, 
                                       double Probability );

  // ********************************************************
  // New algorithm: For efficent (time wise not space waise) processing we
  // will make a 3d table (wavelength,angle in, angle out) 
  // Define table and function to build the 3d table   
  // ********************************************************
  Random* pRandom;
  std::vector < std::vector < std::vector < std::pair <float,float> > > >
                                                      fCumulativeProbTable;
  double interpolate(double x1, double y1, double x2, double y2, double x);
  float interpolate(float x1, float y1, float x2, float y2, float x);
  double findAngleOutRadFromTables( int iWavelength, int iAngleIn, 
                                               double probability);
  double calculateAngleOutRadFromTable(double wavelengthNm, double angleInRad,
                                       double probability);
  double fCosBlazeAngle;

 public:
  void build3DVector(double angleXOutRad, double angleYOutRad, double& vxOut,
                     double& vyOut, double& vzOut);
  int iTableFailCount1;
  int iTableFailCount2;
  int iTableFailCountAnd;
  int iTableFailCountOr;


// **********************************************************
};
#endif
