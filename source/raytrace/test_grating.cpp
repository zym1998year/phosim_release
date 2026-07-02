///
/// @package phosim
/// @file test_grating.cpp
/// @brief Debug tests for grating class.
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

#include "grating.h"
#include "ancillary/random.h"
#include <ctime>
#include <string>
#include <fstream>
#include <iostream>
#include <string>

double const kPi = 3.1415926536;

void lookupPhotToTestFile( std::string TestID, int photonID,
                              std::ofstream& TestFile, 
                              Random& random, int WavelengthIndex, 
                              int AngleInIndex, Grating& grating)
{
  double probability = random.uniform();
  double angleOutDeg = grating.findAngleOutRadFromTables( WavelengthIndex,
                                                          AngleInIndex, 
                                                          probability) 
                                                                 / DEGREE;
  double WavelengthNm = grating.fTDef.fWavelengthMinNm + 
                             WavelengthIndex * grating.fTDef.fWavelengthStepSizeNm;
  double angleInDeg = (grating.fTDef.fAngleInMinRad + 
                       AngleInIndex * grating.fTDef.fAngleInStepSizeRad ) / DEGREE;
  if ( abs(angleInDeg) < 1e-5 ) {
    angleInDeg = 0.0;
  }

  TestFile << TestID << " " << photonID << " " << WavelengthNm << " " 
           << angleInDeg << " " << probability << " " << angleOutDeg 
           << std::endl;
  return;
 }

void usage()
{
  std::cout << " test_grating usage: " << std::endl;
  std::cout << "    $> test_grating < test ID >  "<< std::endl;
  std::cout << " where <test ID > : " << std::endl;
  std::cout << "  '1'  : Test grating.build3Dvector method." << std::endl;
  std::cout << "  '2'  : Test grating.MakeCumulativeProbabilityTable method." 
            << std::endl;
  std::cout << "  '3'  : Test grating.findAngleOutRadFromTables method at" 
            << std::endl;
  std::cout << "  '4'  : Test grating.calculateAngleOutRadFromTable method at" 
            << std::endl;
  std::cout << "              various fixed wavelengths" << std::endl;
  std::cout << "  '5'  : Test grating.calculateAngleOutRadFromTable at" 
            << std::endl;
  std::cout << "              various fixed angleIn" << std::endl;
  std::cout << "  '6'  : Test timing ( duration) of  grating.diffract" 
            << std::endl;
  std::cout << "  '7'  : End to end test of diffract()" << std::endl;
  std::cout << "  'All' or 'all' : Execute all the above tests" << std::endl;
  std::cout << " Ex: './test_grating 3' or './test_grating all' " 
            << std::endl;
  std::cout << " Results of test are appended to file: ./gratingtest.dat" 
            << std::endl;
  std::cout << " Results can be plotted by running plotgratingtest.py" 
            << std::endl;
  return;
}


// **************************************************************************
int main( int argc, char** argv)
{
  // If no arguments print usage
  if ( argc == 1 ) {
    usage();
    return 0;
  }
  
  //Get the first argument which is the test to execute.
  int testToEnable;
  if ( argc >= 2 ) {
    // See if we are to do all tests
    std::string argument = argv[1];
    if ( argument== "All" ||  argument == "all" ) {
      testToEnable = 0; 
    }
    else {   
      testToEnable = std::stoi(argv[1]);
    }
  }
  if (testToEnable < 0 || testToEnable > 7 ) {
    std::cout << "test_grating: Illegal argument:" << argv[1] << std:: endl;
    usage();
    return 0;
  }

  //Setup for tests
  int nPhotons = 100;
  
  double vxIn;
  double vyIn;
  double vzIn;
  double vxOut;
  double vyOut;
  double vzOut;
  
  Random random;

  // Mostly scopio values
  //  double testWavelengthNM=  757.0;
  double testWavelengthNM=  503.0;
  //Landon's test values:
  //double minWavelengthNM = 691.0;
  //double maxWavelengthNM = 818.0;
  double minWavelengthNM = 450.0;
  double maxWavelengthNM = 800.0;
 
  // Define Output Test file:
  std::string instrDir = "./";
  std::ofstream testFile;
  testFile.open ("gratingtest.dat", std::ios_base::out);
  if (!testFile) {
      std::cout << "Failed to open new gratingtest.dat output file!" 
                <<std::endl;
      return 1;
    }
  std::cout<< "Existing gratingtest.dat file opened successfully." 
           << std::endl;

  // Generate Grating object. Don't generate the cumulative probability table
  // if this is test 1 (Thats handled in the Grating class).
  std::clock_t start;
  double duration;
  start = std::clock();


  Grating grating(instrDir, &random, minWavelengthNM, maxWavelengthNM, 
                  testToEnable);

  if (testToEnable != 1) {
    std::cout << "Grating Cumulative Table Acquisition Complete" 
            << std::endl;
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Table acquisition duration: "<< duration/60 << " min" 
             << std::endl;
  }

  // Test 1
  // test build3Dvector code
  if ( testToEnable == 0 || testToEnable == 1 ) {
    // Tests grating.build3DVector method
    // Tests that the generation of the output photon vector in the local
    // coord 
    // system of the grating works when the input vector's angleXOutRad is 
    // different from the original angleXInRad because of the gratings 
    // diffration.
    std::cout << "Test1." << std::endl;
    nPhotons = 100;   
    
    // Above we have instantiated the grating class so we can use the build3DVector
    // method.( But we dont generate the table if this is only test 1).
    
    //Put out a heading for the output results table
    testFile << "#Test 1:" << std::endl;
    testFile << "#theta:phi:vxI:vyI:vzI:angleXOut:angleYOut:vxO:vyO:vzO:"
                "angleBetween" << std::endl; 
    double sumDifferenceDeg = 0;
    for (int i = 0; i < nPhotons; i++) {
      double thetaRad = random.uniform()* kPi/2.;
      double phiRad = random.uniform()* 2*kPi;

      // Make a "reversed" photon vector:
      double vxRIn = cos(phiRad) * sin(thetaRad);
      double vyRIn = sin(phiRad) * sin(thetaRad);
      double vzRIn = cos(thetaRad);   //Positive GOING PHOTON
      //reverse its direction to negative going photon(which is what we want)
      vxIn = -vxRIn;
      vyIn = -vyRIn;
      vzIn = -vzRIn;

      //We want the out going angles relative to negative going photon direction 
      // I.e - z axis which is reason for the - signs(use of VxRIn etc is just for 
      // ease of use here).
      double angleXOutRad = -asin( vxRIn / sqrt(vxRIn*vxRIn + vzRIn*vzRIn) );
      double angleYOutRad = -asin( vyRIn / sqrt(vyRIn*vyRIn + vzRIn*vzRIn) );
 
      grating.build3DVector(angleXOutRad, angleYOutRad, vxOut, vyOut, vzOut);
      
      //acos of dot product (angle beteen):
      double cosAngleBetween = (vxIn*vxOut + vyIn*vyOut + vzIn*vzOut);
      if (cosAngleBetween > 1.0) {
        cosAngleBetween = 1.0;
      }
      double angleBetweenDeg =  acos( cosAngleBetween ) * 180 / kPi;
      
      testFile << "1.1 " <<thetaRad*180/kPi << ' ' << phiRad*180/kPi << ' ' 
                << vxIn 
                << ' ' << vyIn << ' ' << vzIn << ' ' << angleXOutRad*180/kPi 
                << " " <<  angleYOutRad*180/kPi << " " << vxOut << " " 
                << vyOut 
                << ' ' << vzOut << ' ' << angleBetweenDeg << std::endl;
      sumDifferenceDeg += angleBetweenDeg;
    }
    testFile << "#SumDifferenceDeg: " << sumDifferenceDeg 
             << "#for nPhotons: " << nPhotons << std::endl;
    testFile << "#Test 1 Complete." << std::endl;
    std::cout << "Test 1 Complete!" <<std::endl;
  }
  
  // Test 2
  // Tests grating.MakeCumulativeProbabilityTable method.
  if ( testToEnable == 0 || testToEnable == 2 ) {
    // Do this by making a text file of the intensity distribution along 
    // angleOut for the approximate central wavelength in the Probability 
    // table and angleIn ~= 0.0. This is a "line" in the table. We will 
    // then write to this file 2 sets of data. 
    // The first (flagged by first element == 2.1 (1) ) will be the 
    // differntial probabilty intensty (element 5) for the closest central 
    // wavelength (2), anglein ~ 0.0 (3)  and along the 'fine ' angleout (4) 
    // values. This distribution will be generated on the fly by this 
    // program. 
    // Following this set will be he second set (flaged by first 
    // element == 2 (1) ) will be for the same wavelength (2), anglein (3)) 
    // and variable angle out (4) of the cumulative proability (5).
    // This cumulative probability distribution for this line will be 
    // retrieved from the table.
 
    testFile << "#Test 2:" << std::endl;

    int angleOutSize = grating.fTDef.fNumAngleOutSteps+1;
    std::vector < float > angleOutTest2Intensity(angleOutSize,0.0);
    std::vector < float > angleOutTest2Deg(angleOutSize,0.0);

    // Set up the patrcular test2 table indices
    int test2WavelengthIndex;
    int test2AngleInIndex;
    float angleInDeg;
    // Find angleIn Index (near angleIn ~=  0 deg) and wavelength (near central 
    // wavelength) if none specified.
    // testWavelengthNM = -1 if flag that none was specified.
    // in that case use the middle wavelingth index
    if( testWavelengthNM < 0 ) {
      test2WavelengthIndex = grating.fTDef.fNumWavelengthSteps/2;
    }
    else {
      test2WavelengthIndex  = (testWavelengthNM - grating.fTDef.fWavelengthMinNm) / 
                                                             kWavelengthStepSizeNm;
      if(test2WavelengthIndex <0 ) {
        test2WavelengthIndex = 0;
      }
      else if (test2WavelengthIndex > grating.fTDef.fNumWavelengthSteps-2 ) {
        test2WavelengthIndex = grating.fTDef.fNumWavelengthSteps-2;
      }
    }
    if ( grating.fTDef.fAngleInMinRad >= 0.0) {
      test2AngleInIndex =0;
    }
    else { 
      test2AngleInIndex = -grating.fTDef.fAngleInMinRad /
                                            grating.fTDef.fAngleInStepSizeRad;
      // If Range doesn't cover 0. Get as close as I can
      if (test2AngleInIndex > grating.fTDef.fNumAngleInSteps-2) {
        test2AngleInIndex =  grating.fTDef.fNumAngleInSteps-2;
      }
    }
    test2AngleInIndex++;
    // Now generate this particular intensity distribution
    double angleInRad = (grating.fTDef.fAngleInMinRad + 
                      test2AngleInIndex * grating.fTDef.fAngleInStepSizeRad );
    double wavelengthNm = grating.fTDef.fWavelengthMinNm + 
                      test2WavelengthIndex * grating.fTDef.fWavelengthStepSizeNm;
    for (int iAngleOut = 1; iAngleOut < grating.fTDef.fNumAngleOutSteps+1; 
                                                            iAngleOut++) {
      double angleOutRad = grating.fTDef.fAngleOutMinRad + 
                         (iAngleOut * grating.fTDef.fAngleOutStepSizeRad);
      angleOutTest2Intensity.at(iAngleOut) = 
                 grating.calculateFunctionTransmissionGrating(angleInRad, 
                                                              angleOutRad, 
                                                              wavelengthNm);
      angleOutTest2Deg.at(iAngleOut) = angleOutRad * 180.0/PI;
    }

    // Setup the output text file for this test
    testFile << "#TestID:WavelengthNm:AngleInDeg:AangleOutDeg:Intensity" 
             << std::endl;
    angleInDeg = angleInRad * 180/PI;
    if ( abs(angleInDeg) < 1e-5 ) {
      angleInDeg = 0.0;
    }
    for (int i = 0; i <  angleOutSize; i++) {
      testFile << "2.1 "<< wavelengthNm << " " << angleInDeg << " "
               << angleOutTest2Deg.at(i) << " "<< angleOutTest2Intensity.at(i) 
               << std::endl;
    }
    
    // Now the cumulative probability 2.2
    testFile << "#TestID:WavelengthNm:AngleInDeg:"
      "AangleOutDeg:cumulative Probability" << std::endl;
    //Set pointers to within the table.
    std::vector< std::pair < float, float > >* pProbTableAngleOut = 
      &grating.fCumulativeProbTable.at(test2WavelengthIndex).
      at(test2AngleInIndex);
    int probSize = pProbTableAngleOut->size();
    
    for (int i = 0; i < probSize; i++) {
      testFile << "2.2 " << wavelengthNm << " " << angleInDeg << " "
               << pProbTableAngleOut->at(i).first * 180.0/PI << " " 
               << pProbTableAngleOut->at(i).second << std::endl;
    }
    // Test2 close up and exit if done.
    std::cout << "Test 2 complete!" << std::endl;
  }
  // **********************

  // Test 3
  // Tests grating.findAngleOutRadFromTables method.
  // **********************************************************************
  // Along table axis's. No interpolation
  // The table is made/loaded as normal and then random photons using fixed angleIn , 
  // and 5 different fixed wavelengths are read from the table. Results are 
  // written to the gratingtest.dat file for plotting
  // Approximate Center angleIn always used. Wavelengths:Min(test3.1), 
  // middle(test3.2), middle+1(test3.3), middle+2(test3.4)  and max(test3.5)
  // wavelengths are used. 
  // **********************************************************************
  if (testToEnable == 0 || testToEnable == 3) {
    testFile << "#Test 3:" << std::endl;
    testFile << "TestID:PhotonID:WavelengthNm:AngleInDeg:Probability:"
                "AngleOutDeg" << std::endl;    

    // Define middle input/
    //int testAngleInIndex    = grating.fTDef.fNumAngleInSteps/2; 
    int testAngleInIndex    = grating.fTDef.fNumAngleInSteps/2 + 1; 
    int testMinWavelengthIndex = 0;
    int testMidWavelengthIndex = (testWavelengthNM-minWavelengthNM) / 
                                                  grating.fTDef.fWavelengthStepSizeNm;
    int testMaxWavelengthIndex = grating.fTDef.fNumWavelengthSteps-1;

    int nPhotons = 1000000;
    testMidWavelengthIndex ++;
    for (int i = 0; i< nPhotons; i++) {
      // Tests  grating.findAngleOutRadFromTables method.

      lookupPhotToTestFile( "3.1", i, testFile, random, testMinWavelengthIndex, 
                            testAngleInIndex, grating);
      lookupPhotToTestFile( "3.2", i, testFile, random, testMidWavelengthIndex, 
                            testAngleInIndex, grating);
      lookupPhotToTestFile( "3.3", i, testFile, random, testMidWavelengthIndex+1, 
                            testAngleInIndex, grating);
      lookupPhotToTestFile( "3.4", i, testFile, random, testMidWavelengthIndex+2, 
                            testAngleInIndex, grating);
      lookupPhotToTestFile( "3.5", i, testFile, random, testMaxWavelengthIndex, 
                            testAngleInIndex, grating);
    }
    std::cout << "Test 3 Complete!" << std::endl;
  }

  //  Test 4
  // Tests grating.calculateAngleOutRadFromTable method. 
  // *************************
  // Fixed approxamate center angleIN used. 4 different wavelengths used that need 
  // interpolation within the table.
  // The table is made as normal and then random photons using fixed angleIn , 
  // and 5 different fixed wavelengths are read from the table. Results are 
  // written to a gratingtest.dat file for plotting
  if (testToEnable == 0 || testToEnable == 4) {
    // Test reconstruction
    grating.iTableFailCount1 = 0;
    grating.iTableFailCount2 = 0;
    grating.iTableFailCountAnd = 0;
    grating.iTableFailCountOr = 0;

    nPhotons = 100000;
    testFile << "#Test 4: " << std::endl;
    testFile << "#TestID:wavelengthNm:angleInDeg:probability:angleOutDeg"<< std::endl;
    int testAngleInIndex =  grating.fTDef.fNumAngleInSteps/2; 
    double angleInRad = (grating.fTDef.fAngleInMinRad + 
                         testAngleInIndex * grating.fTDef.fAngleInStepSizeRad );
    angleInRad = angleInRad + grating.fTDef.fAngleInStepSizeRad / 3; //Not on a line
    if ( abs(angleInRad/DEGREE) < 1e-5 ) {
      angleInRad = 0.0;
    }
    double testWavelengthTestStepSize = 1.+random.uniform()*.25;

    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                       angleInRad,
                                                       probability);
      testFile << "4.1 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }
    
    testWavelengthNM += testWavelengthTestStepSize;
    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "4.2 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }

    testWavelengthNM += testWavelengthTestStepSize;
    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "4.3 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }

    testWavelengthNM += testWavelengthTestStepSize;
    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "4.4 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }
   
    // Debug print
    std::cout << "Table Failure Counts: " << std::endl;
    std::cout << "iTableFailCount1 = " << grating.iTableFailCount1 << std::endl;
    std::cout << "iTableFailCount2 = " << grating.iTableFailCount2 << std::endl;
    std::cout << "iTableFailCountAnd = " << grating.iTableFailCountAnd 
              << std::endl;
    std::cout << "iTableFailCountOr = " << grating.iTableFailCountOr << std::endl;
    std::cout << " out of " << 4 * nPhotons << "photons" << std::endl;

    std::cout << "Test 4 Complete!" <<std::endl;
    
  } 

  // Test 5 
  // **************************
  // tests grating.calculateAngleOutRadFromTable
  // Fixed wavelength, not on axis(requires interpolation)
  // Various fixed angleIN , not on axis
  if (testToEnable == 0 || testToEnable == 5) {
    testFile << "#Test 5:" << std::endl;
    testFile << "#TestID:wavelengthNm:angleInDeg:probability:angleOutDeg"
             << std::endl;
    
    int nPhotons = 100000;
    int testAngleInIndex =  grating.fTDef.fNumAngleInSteps/2; 

    double angleInRad = (grating.fTDef.fAngleInMinRad + 
                         testAngleInIndex * grating.fTDef.fAngleInStepSizeRad );

    if (abs(angleInRad/DEGREE) < 1e-5 ) {
      angleInRad = 0.0;
    }
 

    // test 5.0a
    // tests calculateSingleSlitEnvelope() and 
    // calculateInterferenceIntensity()
    // Look around 14.0 deg < angleOut <18.0deg
    //Use angleIn as given
    testWavelengthNM=  503.0;
    for (int i=0; i < 10000; i ++ ) {
      double angleOutRad = ( 0.0 + i * 20.0 / 10000.0 ) *  DEGREE;

      double interference = 
        grating.calculateInterferenceIntensity(angleInRad, 
                                               angleOutRad,
                                               testWavelengthNM);
      double envelope = 
        grating.calculateSingleSlitEnvelope(angleInRad, 
                                            angleOutRad,
                                            testWavelengthNM);

      testFile << "5.0a " << testWavelengthNM << " " << angleInRad / DEGREE 
               << " " << angleOutRad / DEGREE<< " " 
               << interference * interference << " " << envelope * envelope 
               << std::endl;
    }


    //test 5.1
    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                       angleInRad,
                                                       probability);
      testFile << "5.1 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }
    
    // 1/8  step in angleIn
    angleInRad += .003142159/2;

    // test 5.0b
    // tests calculateSingleSlitEnvelope() and 
    // calculateInterferenceIntensity()
    // Look around 14.0 deg < angleOut <18.0deg
    //Use angleIn as given
    for (int i=0; i < 10000; i ++ ) {
      double angleOutRad = ( 0.0 + i * 20.0 / 10000.0 ) *  DEGREE;

      double interference = 
        grating.calculateInterferenceIntensity(angleInRad, 
                                               angleOutRad,
                                               testWavelengthNM);
      double envelope = 
        grating.calculateSingleSlitEnvelope(angleInRad, 
                                            angleOutRad,
                                            testWavelengthNM);

      testFile << "5.0b " << testWavelengthNM << " " << angleInRad / DEGREE 
               << " " << angleOutRad / DEGREE<< " " 
               << interference * interference << " " << envelope * envelope 
               << std::endl;
    }

    //Test 5.2
      for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      
      /*
      //This is just for gdb test. Delete when testing complete
      int testTemp;
      if (probability > .32 and probability < .34 ) {
        testTemp =1;
      }
      */

      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "5.2 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }

    // Half step in angleIn
    angleInRad += .003142159/2;

    // Test 5.3
    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "5.3 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }

    // Half step in angleIn
    angleInRad += .003142159/2;

    for (int i = 0; i< nPhotons; i++) {
      double probability = random.uniform();
      double angleOutRad = grating.calculateAngleOutRadFromTable(testWavelengthNM, 
                                                         angleInRad,
                                                         probability);
      testFile << "5.4 " << testWavelengthNM << " " << angleInRad / DEGREE << " " 
                << probability << " " << angleOutRad / DEGREE<< std::endl;
    }
    std::cout << "Test 5 Complete!" << std::endl;
  } 
 
 // **********************************
  // Test timing duration of grating.diffract  method
  if (testToEnable == 0 || testToEnable == 6) {
    // Test reconstruction timing
    std::cout << std::flush;
    int nPhotons = 1000000;

    std::clock_t start;
    double duration;
    start = std::clock();
    for (int i = 0; i< nPhotons; i++) {
      double wavelengthNm = minWavelengthNM + 
                            random.uniform() * (maxWavelengthNM - minWavelengthNM);
      double thetaRad = random.uniform()* 1.0 * DEGREE; 
      double phiRad = random.uniform()* kPi;
      vxIn = cos(phiRad) * sin(thetaRad);
      vyIn = sin(phiRad) * sin(thetaRad);
      vzIn = cos(thetaRad);

      grating.diffract(vxIn, vyIn, vzIn, vxOut, vyOut, vzOut,  wavelengthNm); 
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"Test 6 Complete! Duration "<< duration/60. << " min" << " for "
             << nPhotons << " calls to grating.diffract" << std::endl;
    std::cout << "By design Test 6 adds no output to gratingtest.dat file" << std::endl;
  }

// Test7.
  // Now test end to end call to grating.diffract
  if (testToEnable == 0 || testToEnable == 7) {
    testFile << "#Test 7:" << std::endl;
    testFile << "#TestID:wavelengthNm:angleXInDeg:angleYInDeg:"
                "angleXOutDeg:angleYOutDeg" << std::endl;
    
    // incoming photon is verticle
    double angleXInRad = 0.0;
    double angleYInRad = 0.0;
    grating.build3DVector(angleXInRad, angleYInRad, vxIn, vyIn, vzIn);

    //Test 7.1
    //Flat spectrum first
    int nPhotons = 100000;
    for (int i = 0; i< nPhotons; i++) {
        double wavelengthNm =  minWavelengthNM + 
                            random.uniform() * (maxWavelengthNM - minWavelengthNM);
        grating.diffract(vxIn, vyIn, vzIn, vxOut, vyOut, vzOut,  wavelengthNm); 
        
        double angleXOutRad = -asin( vxOut / sqrt(vxOut*vxOut + vzOut*vzOut) );
        double angleYOutRad = -asin( vyOut / sqrt(vyOut*vyOut + vzOut*vzOut) );

        testFile << "7.1 " << wavelengthNm << ' ' << angleXInRad/ DEGREE << ' ' 
                 << angleYInRad / DEGREE << ' ' << angleXOutRad/ DEGREE << ' ' 
                 << angleYOutRad / DEGREE << std::endl;
    }
  


    // step wavelengths by a constant step size
    // Random fraction from integer wavelength to start.
    double wavelengthStartNm = minWavelengthNM +  random.uniform();
    //Round down on purpose here 
    int numWavelengths = (maxWavelengthNM - wavelengthStartNm -1)/5;
    
    nPhotons = 1000;
   
    double wavelengthStep = 5.0;
    for (int j = 1; j < numWavelengths-1 ; j++ ) {
      double wavelengthNm = j*wavelengthStep + wavelengthStartNm;
      for (int i = 0; i< nPhotons; i++) {
        
        grating.diffract(vxIn, vyIn, vzIn, vxOut, vyOut, vzOut,  wavelengthNm); 
        
        double angleXOutRad = -asin( vxOut / sqrt(vxOut*vxOut + vzOut*vzOut) );
        double angleYOutRad = -asin( vyOut / sqrt(vyOut*vyOut + vzOut*vzOut) );
        testFile << "7.2 " << wavelengthNm << ' ' << angleXInRad/ DEGREE << ' ' 
                 << angleYInRad / DEGREE << " " <<angleXOutRad/ DEGREE << ' ' 
                 << angleYOutRad / DEGREE << std::endl;
      }
    }
    std::cout << "Test 7 Complete!" << std::endl;
    std::cout << std::flush;
  }
  testFile.close();
  std::cout << "gratingTest.dat file closed" << std::endl;
  return 0;
}
