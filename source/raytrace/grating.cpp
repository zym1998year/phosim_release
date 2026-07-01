
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///// @brief Modified by:
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include "grating.h"
  
Grating::Grating(std::string instrDir, Random* pRandm, double Minwavelength, 
                 double Maxwavelength, int TestToEnable)
{
  //TestToEnable can be 1-7 or 0 for all, -1 for no tests
  //for TestToEnable ==1 don't generate the grating lookup table
  testToEnable = TestToEnable;
  setupGrating(pRandm, Minwavelength, Maxwavelength);
  init(instrDir);
}
// **************************************************************************

Grating::Grating(std::string instrDir, Random* pRandm, double Minwavelength, 
                 double Maxwavelength)
{
 
  testToEnable = -1;
  setupGrating(pRandm, Minwavelength, Maxwavelength);
  init(instrDir);
}
// **************************************************************************

Grating::Grating(std::string instrDir, Random* pRandm, 
                 std::vector < double > &  gratingParams)
{
  testToEnable = -1;

  // The grating and table definition parameters normally come in from the 
  // optics_X.txt. If num iluminated lines not specified: (surface.three<1)
  // we use defualt values.) in which case this  constructor not called.
  pRandom = pRandm;
  
  fGDef.fAngleBlazeRad = gratingParams.at(0) * DEGREE;// Angle(radians) between 
                                        // face of groove and plane of grating
  fGDef.fNumSlits      = gratingParams.at(1);  // Number of illuminated grooves 
  fGDef.fDeltaNm       = gratingParams.at(2) * 1000.;  
                                               // Distance between grooves(uM). 
  fTDef.fWavelengthMinNm = (int) (gratingParams.at(3) * 1000.);
                                                            //Min wavelength (nm)
  fTDef.fWavelengthMaxNm = (int) ( gratingParams.at(4) * 1000.) + 1; 
                                                            //Max Wavelength (nm)
  fTDef.fAngleInMinRad =   gratingParams.at(5) * DEGREE;//Min Angle In (radians)
  fTDef.fAngleInMaxRad =   gratingParams.at(6) * DEGREE;//Max Angle In  (radians)
  fTDef.fAngleOutMinRad =  gratingParams.at(7) * DEGREE;//Min Angle Out (radians)
  fTDef.fAngleOutMaxRad =  gratingParams.at(8) * DEGREE;//Max Angle Out (radians)
  init(instrDir);
}
// ***********************************************************************

Grating::~Grating() 
{
  //Nothing here for now
}
// ***************************************************************************

void Grating::setupGrating(Random* pRandm, double Minwavelength, 
                           double Maxwavelength)
{
  pRandom = pRandm;

  // *************************************************************************
  // Wavelength range and number-of-steps are defined as input default 
  // parameters since we may want to specify them according to limits in SED 
  // file.
  // **************************************************************************
  fGDef.fAngleBlazeRad = kBlazeAngleDeg * DEGREE;//Angle(radians) between face 
                                           // of groove and plane of grating
  fGDef.fNumSlits      = kNumSlits;        // Number of illuminated grooves 
                                           // in the grating/object
  fGDef.fDeltaNm       = kDeltaNm;         // Distance between grooves. 
                                           // Inverse is groove density 
  // From site.txt
  fTDef.fWavelengthMinNm = (int)Minwavelength;    // Min wavelength (round down)
  fTDef.fWavelengthMaxNm = (int)Maxwavelength +1; // Max wavelength (round up) 
  fTDef.fAngleInMinRad   = kAngleInMinDeg * DEGREE;
  fTDef.fAngleInMaxRad   = kAngleInMaxDeg * DEGREE;
  fTDef.fAngleOutMinRad  = kAngleOutMinDeg * DEGREE;
  fTDef.fAngleOutMaxRad  = kAngleOutMaxDeg * DEGREE;
}
// **************************************************************************


void Grating::init(std::string InstrDir)
{
  fTDef.fTableVersion = kTableVersion;

  fTDef.fWavelengthStepSizeNm = kWavelengthStepSizeNm;
  fTDef.fAngleInStepSizeRad = kAngleInStepSizeDeg * DEGREE;
  fTDef.fAngleOutStepSizeRad = kAngleOutStepSizeDeg * DEGREE;
  fTDef.fAngleOutProbabilityStep = kAngleOutProbabilityStep;
  fCosBlazeAngle = cos(fGDef.fAngleBlazeRad);

  if (testToEnable != 1) {
    // Determine numberpof elements each axis, then check to make sure that 
    // the min/max table definition values result in a table 
    // with enough elements (3),  in each axis. If not increase max 
    // values until they do.
    generateCumulativeProbTableParameters();

    // At this point fGDef and fTDef have valid table parameters
    //If we already have a valid cumulative table, load it in.
    bool loadedCumProbTable = loadCumulativeProbTable(InstrDir);
    
    fGDef.print();
    fTDef.print();

   
    // If we don't have a valid table(file exists, same fGDef, same fTDef,
    // then make a new table and then save it.
     if ( !loadedCumProbTable ) {
       makeCumulativeProbTable();
      //Save the table to a file since its a new one.
      saveCumulativeProbTableToFile(InstrDir);
    }
  }
}
// ***************************************************************************

void Grating::generateCumulativeProbTableParameters()
// *************************************************************************
// Generate and check (and fix if necessarry) the Cumulative Probability 
// Table definition values.
// *************************************************************************
{
  //Check that the min/max values have enough range bewteen them to 
  // have at least 3 entries in the table of each axis.
  // Note that the generation of the number of step is a round down operation
  // (number of steps are all integer).
  fTDef.fNumWavelengthSteps = 
                       ( (fTDef.fWavelengthMaxNm - fTDef.fWavelengthMinNm) / 
                                           fTDef.fWavelengthStepSizeNm ) + 1;
  if (fTDef.fNumWavelengthSteps < 3) {
    // Min/max Wavelength check
    // The extra .5 of  fTDef.fWavelengthStepSizeNm is just for insurance.
    // Note Number of steps is integer so we round down.
    fTDef.fWavelengthMaxNm = fTDef.fWavelengthMinNm + 
                                         (2.5 * fTDef.fWavelengthStepSizeNm); 
    fTDef.fNumWavelengthSteps =  
                        ( (fTDef.fWavelengthMaxNm - fTDef.fWavelengthMinNm) / 
                                           fTDef.fWavelengthStepSizeNm ) + 1;
  }

  // Min max AngleIn check
  fTDef.fNumAngleInSteps = ( (fTDef.fAngleInMaxRad - fTDef.fAngleInMinRad) / 
                                             fTDef.fAngleInStepSizeRad ) + 1;
  if (fTDef.fNumAngleInSteps < 3) {
 
    // The extra .5 of  fTDef.fAngleInStepSizeRad is just for insurance.
    fTDef.fAngleInMaxRad = fTDef.fAngleInMinRad + 
                                           (2.5 * fTDef.fAngleInStepSizeRad); 
    fTDef.fNumAngleInSteps =  
                      ( (fTDef.fAngleInMaxRad - fTDef.fAngleInMinRad) / 
                                        fTDef.fAngleInStepSizeRad ) + 1;
  }

  // Min max AngleOut check
  fTDef.fNumAngleOutSteps = 
                        ( (fTDef.fAngleOutMaxRad - fTDef.fAngleOutMinRad) / 
                                          fTDef.fAngleOutStepSizeRad ) + 1;
  if (fTDef.fNumAngleOutSteps < 3) {
 
    // The extra .5 of  fTDef.fAngleOutStepSizeRad is just for insurance.
    fTDef.fAngleOutMaxRad = fTDef.fAngleOutMinRad + 
                                        (2.5 * fTDef.fAngleOutStepSizeRad); 
    fTDef.fNumAngleOutSteps =  
                      ( (fTDef.fAngleOutMaxRad - fTDef.fAngleOutMinRad) / 
                                        fTDef.fAngleOutStepSizeRad ) + 1;
  }
}
// *************************************************************************

double Grating::calculateFunctionTransmissionGrating(double AngleInRad, 
                                                     double AngleOutRad,
                                                     double WavelengthNm) 
// ****************************************************************************
// Calculate the value of the intensity function. AngleInRad is the angle to 
// the +z axis and  angleOutRad is the angle to the -z axis.
// ****************************************************************************
{
  // The intensity is given by the interference intensity expression
  // modulated by the single slit diffraction envelope for the slits 
  // which make up the grating.
  // *************************
  // The  interference intensity expression where fGDef.fDeltaNm is slit spacing
  //  and fGDef.fNumSlits is the number of slits of the grating illuminated is:
  // (see http://hyperphysics.phy-astr.gsu.edu/hbase/phyopt/gratint.html#c2) 
  // and  M. Born and E. Wolf, Principles of Optics (Macmillan, New York,964),
  // pg 401-405)
  // ****************************************************************
  // The plus sign on the sin(AngleInRad) terms seems to indidicate that
  // in the following (original You Dong equations),  angleOutRad is  relative 
  // to the -z axis  and AngleInRad is  relative to the +z axis.
  // ********************
  double interferenceIntensity =  calculateInterferenceIntensity(AngleInRad, 
                                                                 AngleOutRad,
                                                                 WavelengthNm);
  // Will be squared below.

  // ******************************************************************
  // I found a possible reference to the single slit (blaze grating) formula 
  // at:
  // https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-31-10-2179
  // pg 2181.
  // Though this is for a reflection grating a transmission grating will be 
  // almost the same. I would think the index of 
  // refraction would matter for a transmission grating? But,
  // I think all the index of refraction would do would be to shift the 
  // diffraction pattern and if the index of refraction  its nearly constant 
  // over wavelength it would just cause just a constant shift in the AngleOut 
  // which we can probably ignore.
  // Note here that the width of the face of the blaze (ie. width of the slot) 
  // is approximatly the grating spacing. ie. fGDef.fDeltaNm (see below)
  // ******************************************************************
  
  // ***********************************************************************
  // The single slit diffraction envelope  
  // for a transmission grating with blaze angle fGDef.fAngleBlazeRad and 
  // spacing fGDef.fDeltaNm:
  // *****************

  double singleSlitEnvelope = calculateSingleSlitEnvelope(AngleInRad, 
                                                          AngleOutRad,
                                                          WavelengthNm); 
  //WILL BE SQUARED BELOW
  //Return the relative intensity (with any constants pretty much ignored).

  double intensity = interferenceIntensity * singleSlitEnvelope;
  return intensity * intensity;
}
// ***************************************************************************

void Grating::defineCumulativeProbTable()
{
  std::vector < std::pair < float,float> >  angleOutTbl(1); 
  std::vector < std::vector < std::pair < float,float> > > 
                             angleInTbl(fTDef.fNumAngleInSteps,  angleOutTbl);
  fCumulativeProbTable.assign(fTDef.fNumWavelengthSteps, angleInTbl);
  return;
}

// ***************************************************************************

void Grating::makeCumulativeProbTable()
// ***************************************************************************
// We make a 3d Cumulative Probability table with steps in wavelength, angleIn 
// and angleOut (outer, middle,inner indices) as indices and with the cumulative 
// intensity along the angle out innermost vector as a value. Since this inner 
// vector is of cumulative values it is naturally sorted from 0 to 1.0 thus 
// allowing for use of the upper_bound vector function to find the angle out 
// index for a specific wavelength,angleIN, probability parameter set. 
// ****************************************************************************
// In order to make the table as small as possible we will be using variable 
// bin sizes for the angleOut axis. To do this we first for each wavelength 
// and angleIn generate the angleOut intensity array with very fine angleOut 
// bins. Then we will convert that array to a cumulative probability array. 
// Then we will go through that array finding the points where the cumulative 
// probability values change enough to warrant a bin step. This will give a 
// much smaller array that will be saved into the final table. For each bin in 
// the deepest angleOut vector we keep a pair. The first values is the value
// of angleOut where the cumulativeProbability step value (from the last bin) 
// equals or exceeds the allowed proability step size. The second value is the 
// cumulative probability at that angleOut.
{
  //vector of a vector of a vector:
  //The inner vector (for angleOut) will be of variable length
  //Note we don't know the size of this vector yet.
  // angleOut.at(n).first will be the angleOut value.
  // angleOut.at(n).second will be the coresponding cumulative probability at 
  // that angleOut
  std::cout << "  Starting build of new Grating Cumulative Probability table."
            << std::endl;
  
  defineCumulativeProbTable();

  std::pair < float, float > sumProb;
  // And now the working AngleOut vector
  std::vector <float> angleOutTemp(fTDef.fNumAngleOutSteps+1,0);

  // **********************************************************************
  int iAngleInPerCent;
  int iAngleInPerCentCount;
  //Test "#" table generation progress print 
  //using round downs here a lot!
  iAngleInPerCent = fTDef.fNumAngleInSteps/100.;
  if ( iAngleInPerCent == 0 ) {
    iAngleInPerCent = 1;
  } 
  iAngleInPerCentCount = 0;


  // ***********************************************************************
  // Now fill it up with the cumulative intensity. Orriginally this used 
  // Yu Dong's function but it was slow and npt clear as to how it worked. I 
  // replaced it with a better documented version. This 
  // will syill take a bit of time to do but it only has to happen once.
  double intensity;
  for (int iAngleIn = 0; iAngleIn < fTDef.fNumAngleInSteps; iAngleIn++ ) {
    //If we are in test mode print out the progress at each 10 %
    if ( ! (testToEnable == -1) ) {
      if ( iAngleIn % iAngleInPerCent == 0 ) {
        if( iAngleInPerCentCount % 10 ==0 ) {
          std::cout << iAngleInPerCentCount/10 << std::flush;
        } 
        else{
          std::cout << "#" << std::flush;
        }
        iAngleInPerCentCount++;
      }
    }

    double angleInRad = (fTDef.fAngleInMinRad + 
                                      iAngleIn * fTDef.fAngleInStepSizeRad );
    for (int iWavelength = 0; iWavelength < fTDef.fNumWavelengthSteps; 
                                                             iWavelength++ ) {
      double wavelengthNm = fTDef.fWavelengthMinNm + 
                                        iWavelength*fTDef.fWavelengthStepSizeNm;
 
      double previousIntensitySum=0;
      for (int iAngleOut = 1; iAngleOut < fTDef.fNumAngleOutSteps+1; 
                                                            iAngleOut++) {
        
        double angleOutRad = 
            fTDef.fAngleOutMinRad + (iAngleOut * fTDef.fAngleOutStepSizeRad);
 
        // *******************************************
        //Now we can fill the table using Yu Dong's function
        intensity = calculateFunctionTransmissionGrating(angleInRad, 
                                                           angleOutRad, 
                                                           wavelengthNm);
         // We are assuming that our step sizes in angleOutRad are small enough 
        // that the intensity is approximatly constant across the bin (or at 
        // least linear)
        double intensitySum = intensity + previousIntensitySum;
        previousIntensitySum = intensitySum;

  
        // Put this in the temporary vector. Its the cumulative intensity sum 
        // along the anlgeOut axis
        angleOutTemp.at(iAngleOut) = intensitySum;
      }
      // Note above sets intensity sum as 0 in bin 0 of angleOutTemp;

      // ************************
      // Now convert the temp angleOutTemp vector from a cumulative intensity 
      // sum to a cumulative probability
      // Divide all by last cumulative intensity sum value to get a cumulative 
      // probability. Then go though the temp Angle out vector which has fine 
      // bins and look for increases in the probability greater than 
      // kAngleOutProbabilityStep. Each time that happens, record the bin 
      // position and the probability in fCumulativeProbTable as a pair.. 
      // ************************

      // Use pointers to avoid lots of pointer arithmetic.
      std::vector< std::pair < float, float> >* pProbTableAngleOut = 
                      &fCumulativeProbTable.at(iWavelength).at(iAngleIn);

      // The element pProbTableAngleOut.at(0) already exists. See where 
      // fCumulativeProbTable was origially defined at the beginning of this method
      // Start at prob 0 in first bin which is at angleOut=fTDef.fAngleOutMinRad
      pProbTableAngleOut->at(0).first  =  fTDef.fAngleOutMinRad;
      pProbTableAngleOut->at(0).second = 0.0;

      double previousProbability =  pProbTableAngleOut->at(0).second;   //Should be 0
      double intensitySumMax    =  angleOutTemp.at(fTDef.fNumAngleOutSteps);

      for (int iAngleOut = 0; iAngleOut < fTDef.fNumAngleOutSteps+1; 
                                                              iAngleOut++ ) {
        
        double intensityProbability = angleOutTemp.at(iAngleOut) / intensitySumMax;
        

        //Search for probability changes > fTDef.fAngleOutProbabilityStep
        if ( (intensityProbability - previousProbability) >  
                                            fTDef.fAngleOutProbabilityStep ) {
          previousProbability =  intensityProbability;
          sumProb.second = intensityProbability;

          sumProb.first = fTDef.fAngleOutMinRad + 
                                       (iAngleOut * fTDef.fAngleOutStepSizeRad);
          //And save it to the table.
          pProbTableAngleOut->push_back(sumProb);
        }
 
      }
      // Finish up.
      if( pProbTableAngleOut->back().second < 1.0 ) {
        std::pair < float, float > lastPair { fTDef.fAngleOutMaxRad, 1.0 };
        pProbTableAngleOut->push_back(lastPair);
      }
      // **********************************************************************

    } //End angleIn loop
  } //End wavelength loop


 return;
}
// *****************************************************************************

void Grating::saveCumulativeProbTableToFile(std::string InstrDir)
{
  // ********************************************************************
  // Save the Cumulative Probability Table(CPT): "fCumulativeProbTable" 
  // to a file named "grating.tbl" 
  // which is in the directory specified in InstrDir(Set in the constructor).
  // Any exixting file by that name will be overwritten.
  // ********************************************************************
  // Structure of file will be first: 2 records that define things.  First is 
  // a record defining the grating itself (structure gratingDefinition). This
  // is follosed by the second record which defines the structure of the 
  // tables (structure tableDefinition). The first element of this rercod is 
  // an integer version number for use in identifying if the table structure 
  // gets changed.
  // ********************************************************************
  // After the defintion records will follow the data from the  table:
  // "fCumulativeProbTable". This table is structured as a vector of vectors 
  // of vectors of pairs.  The lowest  level (innermost) vector (indexed by 
  // the wavelength and angleIn indexes 
  // within the table is a vector of pairs where the first value in the pair 
  // is the angleOut (in Radians) and the second is the cumulative probability
  // at that angleOut summed along the angleOut axis. 
  // This vector was created as a variable length stepping by a minimum fixed
  // probability. 
  // *********************************************************************
  // Since this is a binary streaming file there are no indicators of record
  // breaks. So.. We add a byte count before each "record"

  // Open the file:
  std::string gratingFileName = InstrDir + fGratingTableFileName;
  std::ofstream tableFile(gratingFileName.c_str(), 
                                          std::ios::out | std::ios::binary);
   if(!tableFile) {
     std::cout << "Fatal: Cannot create file: " << gratingFileName <<std::endl;
     exit(1);
   }
 
   // Save the first 2 records. these structres were setup by the 
   // constructor, setup and init methods.
   int recordSize =   sizeof(gratingDefinition);
   tableFile.write((char*) &recordSize, sizeof(int) );
   tableFile.write((char*) &fGDef, recordSize );

   recordSize =   sizeof(tableDefinition);
   tableFile.write((char*) &recordSize, sizeof(int) );
   tableFile.write((char*) &fTDef, recordSize );
   if (!tableFile) {
     std::cout << "Fatal-gDef or tDef writes to " << gratingFileName 
               << " Failed!" << std::endl;
     exit(1);
   }

   // ************************************************************
   // Writing out a vector is problamatical. Its best to extract the data from 
   // the vector into an array and then write the array out.

   // Setup loops for the loops of the  tables 
   for (int iAngleIn = 0; iAngleIn < fTDef.fNumAngleInSteps; iAngleIn++ ) {
     for (int iWavelength = 0; iWavelength < fTDef.fNumWavelengthSteps; 
                                                        iWavelength++ ) {
       // set up vector pointers to each wavelength/angleIn vector. These 
       // vector is a vector of pairs, of which .first is the angleOut and 
       // .second is the cumulative probaility along angleOut for constant
       // wavelength and angleIn  of the vector elements.
       std::vector< std::pair < float, float > >* pProbTableAngleOut = 
                      &fCumulativeProbTable.at(iWavelength).at(iAngleIn);

       // Now create the array we will copy these vectors into
       // They will obviously have the same length. We could check that 
       // but why?

       std::vector<std::pair<float, float>> cumulativeArray(pProbTableAngleOut->size());

       std::copy(pProbTableAngleOut->begin(), pProbTableAngleOut->end(), 
		 cumulativeArray.begin());
       
       //Write out the the anglein and wavelegth indexes followed by the size of
       // the array of pairs, followed by the pair array contents
       recordSize =   static_cast<int>(cumulativeArray.size() * sizeof(std::pair<float, float>));

       tableFile.write((char*) &iAngleIn, sizeof(int) );
       tableFile.write((char*) &iWavelength, sizeof(int) );
       tableFile.write((char*) &recordSize, sizeof(int) );
       tableFile.write((char *) cumulativeArray.data(), recordSize);
       if (!tableFile) {
         std::cout << "Fatal-Data writes writes to " << gratingFileName 
                   << " Failed!" << std::endl;
         exit(1);
       }
     }
   }
   tableFile.close();
   std::cout << "  Grating Culmulative Probaility Table saved to : " 
             << gratingFileName << std::endl;
   return;
}
// **************************************************************************
 
bool Grating::loadCumulativeProbTable(std::string InstrDir)
{
  // *************************************************************************
  // Load in an existing InstDir/grating.tbl file.
  // *************************************************************************
  // Structure:
  // record1: # bytes in record; grating definition params: (gratingDefinition)
  // record2: # bytes in record; table definition params: (tableDefinition)
  // records 3,4....
  // #bytes per array :WvlngIndx,AngInIndx(0,0) Array of pair<angleOut,CumProb>
  // #bytes per array :WvlngIndx,AngInIndx(0,1) Array of pair<angleOut,CumProb>
  // #bytes per array :WvlngIndx,AngInIndx(0,2) Array of pair<angleOut,CumProb>
  //  etc...
  // *************************************************************************
  // read in the first grating definition record and check that its matches this 
  // gratings parameters.


  //Open the specified table file(if it exists)
  std::string tableFileName = InstrDir + fGratingTableFileName;
  std::ifstream tableFile( tableFileName.c_str(), 
                                            std::ios::in | std::ios::binary);
  if(!tableFile) {
    std::cout << "Unable to open exiating file" << tableFileName << std::endl;
    return false;
  }
  
  //There is a file. Read in the grating definition and check it against the 
  // required definition
  int recordSizeBytes;
  tableFile.read( (char*) &recordSizeBytes, sizeof(int) );
  tableFile.read( (char*) &fGTestDef, recordSizeBytes );
  if (!tableFile) {
    std::cout << "Fatal-Read of GDef from " << tableFileName 
              << " Failed!" << std::endl;
    exit(1);
  }
  if ( !fGDef.cmp( fGTestDef) ) {
    return false;
  }

  // Now the table definition record
  tableFile.read(  (char*) &recordSizeBytes, sizeof(int) );
  tableFile.read(  (char*) &fTTestDef, recordSizeBytes );
  if (!tableFile) {
    std::cout << "Fatal-Read of TDef from " << tableFileName 
              << " Failed!" << std::endl;
    exit(1);
  }
  if ( !fTDef.cmp( fTTestDef) ) {
    return false;
  }

  // **************************************************************
  // Now we can load up the table. See above for table structure.
  // We will read each wavelength,angleIn line of the table into an array
  // and then copy that array into the vector of vector of vector of pairs
  // ie CumulativeProbabilityTable.
  // **************************************************************
  defineCumulativeProbTable();

  // **********************************************************
  // Loop over reading in the angleIn,wavelegtn records until we reach the 
  // end of file.
  // use the iAngleIn adn iWavwelength indices to find where to put the
  // records
  int iAngleIn;
  int iWavelength;
  while(true) {
    tableFile.read( (char*) &iAngleIn,sizeof(int));
    tableFile.read( (char*) &iWavelength,sizeof(int));
    tableFile.read( (char*) &recordSizeBytes,sizeof(int));
    if(!tableFile) {
      break;
    }
  
    int numPairs = recordSizeBytes / sizeof( std::pair< float,float > );
    std::vector<std::pair<float, float>> tempCumulativeProbArray(numPairs);
    tableFile.read( (char*) tempCumulativeProbArray.data(), recordSizeBytes);
    if (!tableFile) {
      std::cout << "Fatal-Read of data from " << tableFileName 
                << " Failed!" << std::endl;
      exit(1);
    }

    //Now move this data into the vector table.
    std::vector< std::pair < float, float > >* pProbTableAngleOut = 
                      &fCumulativeProbTable.at(iWavelength).at(iAngleIn);
    pProbTableAngleOut->at(0) =  tempCumulativeProbArray[0];
    for (int iAngleOut = 1; iAngleOut < numPairs; iAngleOut++) {
      
      pProbTableAngleOut->push_back( tempCumulativeProbArray[ iAngleOut ] );
    }
  }
  tableFile.close();
  std::cout << "  Grating Culmulative Probability Table loaded from " 
            << tableFileName << std::endl;
  return true;
}
// ***************************************************************************


void Grating::diffract(double vxIn, double vyIn, double vzIn, double& vxOut,
                       double& vyOut, double& vzOut, double WavelengthNm) 
// ****************************************************************************
// Main method: Calculates the direction out of a photon after it interacts
// with the grating.
// vxIn,vyIn,vzIn is input direction of the photon
// vxOut,vyOut,vzOut is output direction of the photon.
// ****************************************************************************
{
  // ************************************************************************
  // To use the grating, we need to first transform the photons from the lab
  // frame to the frame of the grating (where the normal to the grating is 
  // (0,0,1). 
  // That's done before the call to this function using the transform() method.
  // In this frame the photons directions for a transmission grating should be
  // coming from abocve the grating. That is "vzIn" should be negative.  
  // Make assumption that in the grating frame the ruled lines are along the Y 
  // direction. Thus the photon retains its Y direction angle but after 
  // diffraction gets a new X direction angle.
  // This function will return a unit vector in the new direction of the 
  // photon in the gratings frame. After that the calling program will need to 
  // transform back to normal frame using the transformInverse() method.
  // If the normal is truly (0,0,1) and the lines are along the Y axis.
  // We want the angle between the Vin vector projected on the XZ plane and the
  // Z axis (angleXIn) and we also want the angle between the Vin vector 
  // projected on the YZ plane and the Z axis (angleYIn).
  // NOte: See notes in method  Grating::calculateFunctionTransmissionGrating:
  // " angleOutRad is  relative to the -z axis  and angleInRad is  relative to 
  // the +z axis. ( I think for both X and Y?)" Where the angle in vector is 
  // the reverse of the photon vector.

  //  
  // angleXOut will be calculated by the calculateAngleOutRadFromTable method
  // using a random number to  pick from the possible angleXOut values due to 
  // the diffraction of the grating.  //
  // We then use the angleYOut and angleXOut values to reconstruct the VOut 
  // unit vector.
  // *************************************************************************
  // (Sometime we should worry about vzIn == 0? situation. ie.the  photon is in 
  // the grating plane. No change of photon direction?)
  // I vzIn will always be negative, but we want angle of the photon to the +z 
  // axis. This is the convention.
  // may want to check for 0 divide here (photon misses the grating.

  //For the incoming the vector the calculateAngleOutRadFromTable method uses
  // the reverse of actual photon vector.
  double vxRIn = - vxIn;
  double vyRIn = - vyIn;
  double vzRIn = - vzIn;
  double angleXInRad  = asin( vxRIn / sqrt(vxRIn*vxRIn + vzRIn*vzRIn) );

  // Find a uniformly random probability for this photon
  double probability = pRandom->uniform();

  //Glenns method. Hopefully faster and very nerarly same results as Yu Dongs.
  //Note here that angleXInRad is relative to + z axis (reversed photon 
  //direction ) and angleXOutRad is final photon direction  relative to -z axis
  double angleXOutRad = calculateAngleOutRadFromTable(WavelengthNm, angleXInRad,
                                                       probability );
  // This angleXOutRad is the angle of the x axis component 
  // to the -z axis.

  // *****************************************************
  // Now generate the outgoing vector
  // angleYOutRad is relative to -z axis

  double angleYOutRad = -asin( vyRIn / sqrt(vyRIn*vyRIn + vzRIn*vzRIn) );

  build3DVector(angleXOutRad, angleYOutRad, vxOut, vyOut , vzOut);
  return;
}
// ***************************************************************************

inline double Grating::interpolate(double x1, double y1, double x2, double y2, 
                                   double x) 
{
  return y1 +  ( (x - x1) / (x2 - x1) )* (y2 - y1);
}
// *****************************************************************************

inline float Grating::interpolate(float x1, float y1, float x2, float y2, 
                                  float x) 
{
  return y1 +  ( (x - x1) / (x2 - x1) )* (y2 - y1);
}
// *****************************************************************************

// ****************************************************************************

double Grating::findAngleOutRadFromTables(int iWavelength,int iAngleIn, 
                                         double Probability) 
// *****************************************************************************
// Find the interpolated AngleOutRad for the supplied Probability in the 
// fCumulativeProbTable along the specifed wavelength and angleIn axis. 
// *****************************************************************************
// Using supplied Probability, Look it up in the table and return an 
// interpolated AngleOutRad
// *****************************************************************************
{
  // Use a pointer to get to the vector of pairs. 
  // Pairs are of the form: pair <float angleOut, float CumulativeProbability> 
  std::vector< std::pair < float, float > >* angleOutVector =
                            &fCumulativeProbTable.at(iWavelength).at(iAngleIn);

  // ***********************************************************************
  // We will do the search along this vector ourselves. It would have been nice 
  // to use the C++ upper_bound function but it doesn't seem to be working for 
  // use with pairs.
  // This search could be done as a binary search but we will use a linear search.
  // By definition both elements of the vector of pairs are ascendingly sorted.
  // We search for the index in the vector of the last pair that has a Cumulative 
  // Probability less than Probability. Note: Last pair has Probability == 1.
  int iAngleOutLow = 0;
  int iAngleOutSize = angleOutVector->size();
  for (int iAngleOut = 0; iAngleOut < iAngleOutSize; iAngleOut++ ) {
    if (angleOutVector->at(iAngleOut).second < Probability ) {
      continue;
    } 
    else {
      iAngleOutLow = iAngleOut -1;
      break;
    }
  }
  if (iAngleOutLow < 0 ) {
    iAngleOutLow = 0;
  }
  double angleOutLowRad =    (double) angleOutVector->at(iAngleOutLow).first;
  double cumulativeProbLow = (double) angleOutVector->at(iAngleOutLow).second;

  int iAngleOutHi = iAngleOutLow+1;
  double angleOutHiRad =    (double) angleOutVector->at(iAngleOutHi).first;
  double cumulativeProbHi = (double) angleOutVector->at(iAngleOutHi).second;

  // ***************************
  // Now interpolate
  // **************************
  double angleOut = interpolate(cumulativeProbLow, angleOutLowRad,
                                cumulativeProbHi,  angleOutHiRad, Probability);
  return angleOut;
}

// ****************************************************************************

double Grating::calculateAngleOutRadFromTable(double WavelengthNm, 
                                              double AngleInRad, 
                                              double Probability)
// ****************************************************************************
// Using WavelengthNm and AngleInRad determine for a random photon a grating 
// angleOutRad.  Interpolation over wavelength and AngleInRad  within 3D 
// fCumulativeProbTable.
// **********************************************************
// The cumuluativeProbTable which has Wavelength,AngleIn,AngleOut axis,  is not 
// used to get an interpolated content but rather the content is used to get an
// interpolated AngleOut given a wavelength,AngleIn and random 
// cumulativeProbability (content).
// ***********************************************************
// Procedure:
// 1: Find wavelength and AngleIn table indices values that bracket this request.
// 2: This gives 4 pairs of values that define 4 lines along the angleOut axis.
// 3: Get a random probability from A uniform distributed for this photon.
// 4: Generate the interpolated angleOutRad point on each of these lines where 
//    the interpolated cumulative probability matches the supplied random value. 
//    Note we are using the variable bins in angleOut here.
// 5: Take these 4 points as pairs with the same AngleIn index values and 
//    interpolate along the wavelength axis to get 2 new interpolated 
//    angleOut values. This gives us 2 points (with the same 
//    fractional wavelength value ) but differing AngleIn indices. 
//    Interpolate along the line between these 2 points to get a final 
//    angleOutRad point with the given wavelength and angleIn. 
// 5: Return the final angleOutRad.
// *************************************************************
{
  // ********************
  // 1: Find indicies. Limit wavelength to range
  double wavelengthFractionalIndex = 
        (WavelengthNm - fTDef.fWavelengthMinNm) / fTDef.fWavelengthStepSizeNm;
  int iWavelengthLo;
  if (WavelengthNm < fTDef.fWavelengthMinNm ) {
    iWavelengthLo = 0;
  }
  else if(WavelengthNm > fTDef.fWavelengthMaxNm ) {
    iWavelengthLo = fTDef.fNumWavelengthSteps-3;
  }
  else {

    iWavelengthLo = (int) wavelengthFractionalIndex;  //Rounds down
  }
  int iWavelengthHi = iWavelengthLo+1;

  // Now AngleInRad
  double angleInFractionalIndex = 
                 (AngleInRad-fTDef.fAngleInMinRad) / fTDef.fAngleInStepSizeRad;
  int iAngleInLo;
  if (AngleInRad < fTDef.fAngleInMinRad )  {
    iAngleInLo = 0;
  }
  else if(AngleInRad > fTDef.fAngleInMaxRad ) {
      iAngleInLo = fTDef.fNumAngleInSteps-3;
  }
  else {
    iAngleInLo = (int)angleInFractionalIndex; //Rounding down again.
  }
  int iAngleInHi = iAngleInLo + 1;

  // ******************************
  // 2: 4 Lines along AngleOut (A) are now defined defined: 
  // (iWavelengthLo, iAngleInLo,A) , (iWavelengthHi,iAngleInLo,A), 
  // (iWavelengthLo, iAngleInHi, A), (iWavelengthHi,iAngleInHi,A)
  // ******************************


  //4: Find intepolated AngleOutRad  values for these 4 lines , 
  //   giving us 4 points.

  // Pair up the poins to define 2 lines:
  // First two points have same lower angleIn indicies, different wavelength 
  // indices. (We interpolate over wavelength first!) 
  double angleOutRad1 =
    findAngleOutRadFromTables( iWavelengthLo, iAngleInLo, Probability);
  double angleOutRad2 =
    findAngleOutRadFromTables( iWavelengthHi, iAngleInLo, Probability);
  double angleOutRad3 =
    findAngleOutRadFromTables( iWavelengthLo, iAngleInHi, Probability);
  double angleOutRad4 =
    findAngleOutRadFromTables( iWavelengthHi, iAngleInHi, Probability);

  // Iterpolate over the 2 lines from the above 4 points in wavelength.
  // Group by same angleIn:

  double angleOutRad5 = interpolate((double) iWavelengthLo, angleOutRad1,
                                    (double) iWavelengthHi, angleOutRad2,
                                    wavelengthFractionalIndex);
  double angleOutRad6 = interpolate((double) iWavelengthLo, angleOutRad3,
                                    (double) iWavelengthHi, angleOutRad4,
                                    wavelengthFractionalIndex);
    

  //Now iterpolate over these the line these last 2 points make in AngleInRad
  double angleOutRad = interpolate( (double) iAngleInLo, angleOutRad5,
                              (double) iAngleInHi, angleOutRad6,
                              angleInFractionalIndex);

  // Test  for table failure
  // *****************************
  bool tableFail1 = false;
  bool tableFail2 = false;
  // ************************************************************************
  //Wide test of table use.
  // Test to see if the 4 lines above have different enough cumulative 
  // probabilty distributions to result in an erroneus angleOut after 
  // interpolation. Do this by checking that the 4 angleOuts are close enough
  // to each other.
  std::vector < double > angleOutRadVec;
  angleOutRadVec.push_back( angleOutRad1 );
  angleOutRadVec.push_back( angleOutRad2 );
  angleOutRadVec.push_back( angleOutRad3 );
  angleOutRadVec.push_back( angleOutRad4 );
  double angleOutRadMax = 
    *std::max_element(angleOutRadVec.begin(), angleOutRadVec.end());
  double angleOutRadMin = 
    *std::min_element(angleOutRadVec.begin(), angleOutRadVec.end());
  if ( ( abs( angleOutRadMax - angleOutRadMin ) / DEGREE ) > 
       kMaxAllowedAngleOutDiffDeg  ) {
    tableFail1 = true;   //Table use failed at "wide"selection.
    iTableFailCount1 ++;
  }
  // *************************************************************************
  
  //Narrow test:Not good enough
  if ( ( abs( angleOutRad5 - angleOutRad6 ) / DEGREE ) > 
         kMaxAllowedAngleOutDiffDeg ) {
      tableFail2 = true;  //Table failed at narrow selection
      iTableFailCount2 ++;
  }
   
  double angleOutRad7 = 0.0;
  //  if (tableFail1 || tableFail2 ) {
  if (tableFail1 ) {

    //*********
    // Table use failed,probably (it can happen!). 
    // Determine angleOutRad at explicitly at the priginal given parameters
    // This is a lot slower than using the table. Hopefully only about 2% of 
    // events will need this.
    // ****************************************
    angleOutRad7 = calculateSpecificAngleOutRad( WavelengthNm, AngleInRad, 
                                                 Probability );    
    
  }   //End of table use failure test 

  if(tableFail1 || tableFail2) {
    iTableFailCountOr++;
    /*std::cout << tableFail1 << " " <<tableFail2 << " " <<angleOutRad/DEGREE
              << " " << angleOutRad7/DEGREE << std::endl;
    */
    angleOutRad=angleOutRad7;
  }

  if(tableFail1 && tableFail2) {
    iTableFailCountAnd++;
  }

  return angleOutRad;
}
// **************************************************************************
 
void Grating::build3DVector(double angleXOutRad, double angleYOutRad, 
                            double& vxOut, double& vyOut, double& vzOut)
// **************************************************************************
// Using:
// angleXOutRad (projected vector in reversed XZ plane to  reversed z axis 
// angleYOutRad (projected vector in reversed YZ plane to  reversed z axis 
// : find the vector in the non-reversed  x,y,z axis system.
// **************************************************************************
// Following is a little tricky to derive but I think its correct. Some 
// special tests may have to be added to prevent infinities while taking the 
// tangents.
// ****************************************************************************
{
  // First find the direction of the vector in spherical coordinates. 
  // ie. find theta and phi.(of the reversed photon direction)
  
  //Worry here about :
  // angleXRad -> pi/2   :   tan(angleXoutRad) -> inf
  // angleYOutRad -> 0   :  tan(angleYoutRad) -> 0
  // angleYOutRad -> Pi/2: tan (angleYOutRad -> inf

  // phiRad -> 0         : thetaRad -> inf
  //(See above about pssibility we may have to change angleXOutRad sign here).

  double phiRad;
  double thetaRad;

  //Some special cases
  if (angleXOutRad == 0.0)  {
    vxOut = 0.0;
    vyOut = - sin(angleYOutRad);
    vzOut = - cos(angleYOutRad);
  }
  else if (angleYOutRad == 0.0) {
    vxOut = -sin(angleXOutRad);
    vyOut = 0.0;
    vzOut = - cos(angleXOutRad);
  }
  else {
    phiRad   = atan( tan(-angleYOutRad) / tan(-angleXOutRad) );
    thetaRad = atan( tan(-angleYOutRad) / sin(phiRad)  );

    //Generate VOut vector in the original photon coord system. Its the reverse
    // of the system that -angleyOutRad and -angleXOutRad are in.
    // (thus the - signs below)
    // Actualy I just added the sign to make things right! I can't figure out 
    // why I don't negate vxOut and vyOut also. But I give up for now.
    vxOut = -cos(phiRad) * sin(thetaRad);
    vyOut = -sin(phiRad) * sin(thetaRad);
    vzOut = -cos(thetaRad);
  }
  return;
}
// ****************************************************************************

double Grating::calculateInterferenceIntensity(double AngleInRad, 
                                               double AngleOutRad,
                                               double WavelengthNm)
{
  // *************************
  // The  interference intensity expression where fGDef.fDeltaNm is slit spacing
  //  and fGDef.fNumSlits is the number of slits of the grating illuminated is:
  // (see http://hyperphysics.phy-astr.gsu.edu/hbase/phyopt/gratint.html#c2) 
  // and  M. Born and E. Wolf, Principles of Optics (Macmillan, New York,964),
  // pg 401-405)
  // ****************************************************************
  // The plus sign on the sin(AngleInRad) terms seems to indidicate that
  // in the following (original You Dong equations),  angleOutRad is  relative 
  // to the -z axis  and AngleInRad is  relative to the +z axis.
  // ********************
  double vPrime = (PI * fGDef.fDeltaNm / WavelengthNm ) *
                                       (sin(AngleOutRad) - sin(AngleInRad)) ;
  //Avoid divide by 0: 
  double interferenceIntensity;
  if ( sin(vPrime) == 0 ) {
    interferenceIntensity = 1. / fGDef.fNumSlits;  //small angle approximation
  }
  else {
    interferenceIntensity = sin(fGDef.fNumSlits * vPrime) /  
                                               (fGDef.fNumSlits * sin(vPrime));
  }
  return interferenceIntensity;
  // Will be squared upon return
}
// *****************************************************************************

double Grating::calculateSingleSlitEnvelope(double AngleInRad, 
                                            double AngleOutRad,
                                            double WavelengthNm)
{
  // ******************************************************************
  // I found a possible reference to the single slit (blaze grating) formula 
  // at:
  // https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-31-10-2179
  // pg 2181.
  // Though this is for a reflection grating a transmission grating will be 
  // almost the same. I would think the index of 
  // refraction would matter for a transmission grating? But,
  // I think all the index of refraction would do would be to shift the 
  // diffraction pattern and if the index of refraction  its nearly constant 
  // over wavelength it would just cause just a constant shift in the AngleOut 
  // which we can probably ignore.
  // Note here that the width of the face of the blaze (ie. width of the slot) 
  // is approximatly the grating spacing. ie. fGDef.fDeltaNm (see below)
  // ******************************************************************
  
  // ***********************************************************************
  // The single slit diffraction envelope  
  // for a transmission grating with blaze angle fGDef.fAngleBlazeRad and 
  // spacing fGDef.fDeltaNm:
  // *****************
  // Shadowing: rho is the ratio of the "apparent" slit opening due to 
  // incoming photon angle. (There may be a sign issue here!)
  double rho;
  if(AngleInRad >= fGDef.fAngleBlazeRad ) {
    rho = cos(AngleInRad)/cos(AngleInRad - fGDef.fAngleBlazeRad); 
  }
  else{
    //This is not right yet. need to fix it. See above reference pg 2182
    rho = fCosBlazeAngle;   
  }
  
  double v = (PI * fGDef.fDeltaNm * rho / WavelengthNm ) *
                                 ( sin(AngleOutRad - fGDef.fAngleBlazeRad) + 
                                   sin(AngleInRad  - fGDef.fAngleBlazeRad) );
  //Avoid divide by 0
  double singleSlitEnvelope;
  if ( v == 0 ) {
    singleSlitEnvelope = 1.0;  //small angle approximation
  }
  else {
    singleSlitEnvelope = sin(v)/v;
  }
  
  return singleSlitEnvelope;
  //WILL BE SQUARED BELOW Upon return
}
// ****************************************************************************


double Grating::calculateSpecificAngleOutRad( double WavelengthNM, 
                                              double AngleInRad, 
                                              double Probability )
{
  // **************************************************************************
  // Generate the angleOut specifically for the given AngleInRad, WavelengthNM
  // and probability. This method is used when the use of the 
  // CumulativeProbabiltyTable fails. This failure should be down at the 3% 
  // level , depending on the various step sizes chosen. This method is comupute 
  // intensive, thus the perfered use of the table. See the 
  // makeCumulativeProbTable() method.
  // This method should be more exact since the interpolations are kept to a 
  // minimum.
  // **************************************************************************
  // Define the working AngleOut vector
  std::vector <float> angleOutTemp(fTDef.fNumAngleOutSteps+1,0);
  double previousIntensitySum = 0.0;

  //Fill this vector with the cumulative Intensity sum.
  for (int iAngleOut = 1; iAngleOut < fTDef.fNumAngleOutSteps+1; 
                                                            iAngleOut++) {
    double angleOutRad = 
            fTDef.fAngleOutMinRad + (iAngleOut * fTDef.fAngleOutStepSizeRad);
 
    // *******************************************
    //Now we generate the cumulative intensity vector
    double intensity = calculateFunctionTransmissionGrating(AngleInRad, 
                                                           angleOutRad, 
                                                           WavelengthNM);
    // We are assuming that our step sizes in angleOutRad are small enough 
    // that the intensity is approximatly constant across the bin (or at 
    // least linear)
    double intensitySum = intensity + previousIntensitySum;
    previousIntensitySum = intensitySum;

    // Put this in the temporary vector. Its the cumulative intensity sum 
    // along the anlgeOut axis
    angleOutTemp.at(iAngleOut) = intensitySum;
  }

  // ************************
  // Now convert the temp angleOutTemp vector from a cumulative intensity 
  // sum to a cumulative probability.Divide all by last cumulative intensity sum 
  // to get the cumulative probability. Then go though the temp Angle out vector
  // converting it to cumulative probability until the cumulative probability 
  // exceeds the requested Probability
  // ************************
  double previousProbability = 0.0;
  double intensitySumMax     = angleOutTemp.at(fTDef.fNumAngleOutSteps);

  // Search for value of iAngleOut which results in the cumulative Probabilty 
  // execceding the specified Probability. Thats as far as we need go.
  double angleOutRad=0.0;
  for (int iAngleOut = 0; iAngleOut < fTDef.fNumAngleOutSteps+1; 
       iAngleOut++ ) {

    double intensityProbability = angleOutTemp.at(iAngleOut) / intensitySumMax;

    if ( intensityProbability > Probability ) {
      // So angleOut is between indices iAngleOut-1 and iAngleOut
      // Interpolate between these 2 bins to get AngleOutRad
    double angleOutHiRad = 
            fTDef.fAngleOutMinRad + (iAngleOut * fTDef.fAngleOutStepSizeRad);
    double angleOutLowRad = 
      fTDef.fAngleOutMinRad + ( (iAngleOut-1) * fTDef.fAngleOutStepSizeRad);

    // ***************************
    // Now interpolate
    // **************************
    angleOutRad = interpolate(previousProbability,  angleOutLowRad,
                                     intensityProbability, angleOutHiRad, 
                                     Probability);
        break;
      }
      else {
        previousProbability =  intensityProbability;
      }
  }
  
  return angleOutRad;
}
// *************************************************************************
