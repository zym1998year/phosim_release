///
/// @package phosim
/// @file telescopesetup.cpp
/// @brief setup for telescope (part of image class)
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <vector>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include "ancillary/readtext.h"
#include "fea.h"

using namespace fea;

int Image::telSetup () {

    //        int ttsave=0; int nsurfsave=0;


    double runningz;
    double cumulativez;
    if (diagnostic) printf("$$%40s %22e %22.6lf\n","telSetup",(double)1,805.0);

    long seedchip = 0;
    for (size_t m = 0; m < chipid.size(); m++) {
        seedchip += static_cast<long>((static_cast<int>(chipid.c_str()[m]%10))*pow(10, m));
    }

    chip.nampx = (maxx - minx + 1);
    chip.nampy = (maxy - miny + 1);
    chip.buffer = 400;
      if (driftX == 0.0 && driftY == 0.0) {
        chip.extraX = chip.nampx;
        chip.extraY = chip.nampy;
        chip.drift = 0;
        state.focalPlane = static_cast<std::atomic<unsigned int>*>(calloc(chip.nampx*chip.nampy, sizeof(std::atomic<unsigned int>)));
        //        state.focalPlaneFloat = static_cast<float*>(calloc(chip.nampx*chip.nampy, sizeof(float)));
    } else {
        chip.extraX = floor(abs(driftX)*exptime)+1 + chip.nampx;
        chip.extraY = floor(abs(driftY)*exptime)+1 + chip.nampy;
        chip.drift = 1;
        state.focalPlane = static_cast<std::atomic<unsigned int>*>(calloc((chip.extraX)*(chip.extraY), sizeof(std::atomic<unsigned int>)));
        //state.focalPlaneFloat = static_cast<float*>(calloc((chip.extraX)*(chip.extraY), sizeof(float)));
    }
    if (exptime*driftX > 0) minxDrift = minx; else minxDrift=minx - (floor(exptime*abs(driftX)) +1);
    if (exptime*driftX > 0) maxxDrift = maxx + (floor(exptime*abs(driftX)) + 1); else maxxDrift=maxx;
    if (driftX == 0.0) {
        maxxDrift=maxx;
        minxDrift=minx;
    }
    if (exptime*driftY > 0) minyDrift = miny; else minyDrift=miny- (floor(exptime*abs(driftY))+1);
    if (exptime*driftY > 0) maxyDrift = maxy + (floor(exptime*abs(driftY)) + 1); else maxyDrift=maxy;
    if (driftY == 0.0) {
        maxyDrift=maxy;
        minyDrift=miny;
    }
      if (saturation == 1) {

          state.satmap = static_cast<std::atomic <int>*>(calloc(chip.extraX*chip.extraY, sizeof(int)));
          for (int i = minxDrift; i <= maxxDrift; i++) {
              for (int j = minyDrift; j <= maxyDrift; j++) {
                *(state.satmap + chip.extraX*(j - minyDrift) + (i - minxDrift)) = 0;
            }
        }

    }

      if (opdfile) {
        state.opd = static_cast<double*>(calloc(nsource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE, sizeof(double)));
        state.opdcount = static_cast<double*>(calloc(nsource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE, sizeof(double)));
        state.cx = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.cy = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.cz = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.r0 = static_cast<double*>(calloc(nsource, sizeof(double)));
        state.epR = static_cast<double*>(calloc(nsource, sizeof(double)));
    }


                //           FILE *outfile;
                // outfile = fopen("temp.txt","w");

    // OPTICS AND COATINGS
    fprintf(stdout, "Building Optics.\n");
    std::ostringstream opticsFile;
    opticsFile << instrdir << "/optics_" << filter << ".txt";
    readText opticsPars(opticsFile.str());
    int totSurf = opticsPars.getSize() + 1;  //detector, focalplane, exit pupil

    surface.setup(totSurf, SURFACE_POINTS);
    perturbation.scatteringsigma = new double[totSurf]();
    perturbation.scatteringb = new double[totSurf]();
    perturbation.scatteringc = new double[totSurf]();
    for (int i =0; i < totSurf; i++) {
        perturbation.scatteringsigma[i] = 0.0;
        perturbation.scatteringb[i] = 0.1;
        perturbation.scatteringc[i] = 2.5;
    }
    coating.setup(totSurf);
    contamination.surfaceDustTransmission = new float[totSurf]();
    contamination.surfaceCondensationTransmission = new float[totSurf]();
    for (int i =0; i < totSurf; i++) {
        contamination.surfaceDustTransmission[i] = 0.0;
        contamination.surfaceCondensationTransmission[i] = 0.0;
    }

    nsurf = 0;
    runningz = 0.0;
    cumulativez = 0.0;
    nmirror = 0;
    for (int t = 0; t < totSurf - 1; t++){
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName, surfacetype, coatingFile, mediumFile;
        iss >> surfaceName >> surfacetype;

        if (surfaceName != "optics") {

          if (surfacetype == "mirror") surface.surfacetype[nsurf] = MIRROR;
          else if (surfacetype == "lens") surface.surfacetype[nsurf] = LENS;
          else if (surfacetype == "filter") surface.surfacetype[nsurf] = FILTER;
          else if (surfacetype == "det") surface.surfacetype[nsurf] = DETECTOR;
          else if (surfacetype == "detector") surface.surfacetype[nsurf] = DETECTOR;
          else if (surfacetype == "grating") surface.surfacetype[nsurf] = GRATING;
          else if (surfacetype == "exitpupil") {
            nsurf++;
            surface.surfacetype[nsurf] = EXITPUPIL;
          }
        
          if (surfacetype == "mirror") {
            nmirror++;
          }
          
          iss >> surface.radiusCurvature[nsurf];
          double dz;
          iss >> dz;
          runningz += dz;
          cumulativez += fabs(dz);
          surface.height[nsurf] = runningz;
          surface.intz[nsurf] = cumulativez;
          iss >> surface.outerRadius[nsurf];
          iss >> surface.innerRadius[nsurf];
          iss >> surface.conic[nsurf];
          iss >> surface.two[nsurf];
          iss >> surface.three[nsurf];
          iss >> surface.four[nsurf];
          iss >> surface.five[nsurf];
          iss >> surface.six[nsurf];
          iss >> surface.seven[nsurf];
          iss >> surface.eight[nsurf];
          iss >> surface.nine[nsurf];
          iss >> surface.ten[nsurf];
          iss >> surface.eleven[nsurf];
          iss >> surface.twelve[nsurf];
          iss >> surface.thirteen[nsurf];
          iss >> surface.fourteen[nsurf];
          iss >> surface.fifteen[nsurf];
          iss >> surface.sixteen[nsurf];
          iss >> coatingFile;
          if (surfacetype == "det") {
            if (sensorCoatingFileFlag==1) {
              coatingFile = sensorCoatingFile;
            }
            if (telescopeMode == 0) {
              surface.radiusCurvature[nsurf] = 0.0;
              surface.conic[nsurf] = 0.0;
              surface.two[nsurf] = 0.0;
              surface.three[nsurf] = 0.0;
              surface.four[nsurf] = 0.0;
              surface.five[nsurf] = 0.0;
              surface.six[nsurf] = 0.0;
              surface.seven[nsurf] = 0.0;
              surface.eight[nsurf] = 0.0;
              surface.nine[nsurf] = 0.0;
              surface.ten[nsurf] = 0.0;
              surface.eleven[nsurf] = 0.0;
              surface.twelve[nsurf] = 0.0;
              surface.thirteen[nsurf] = 0.0;
              surface.fourteen[nsurf] = 0.0;
              surface.fifteen[nsurf] = 0.0;
              surface.sixteen[nsurf] = 0.0;
            }
          }
          iss >> mediumFile;
          surface.centerx[nsurf] = 0.0;
          surface.centery[nsurf] = 0.0;
          surface.rmax[nsurf] = surface.outerRadius[nsurf];
          surface.innerRadius0[nsurf] = surface.innerRadius[nsurf];
          

          // ************
          // Code to initalize the grating tables: This may take a bit of 
          // time (up to 30+  min for each grating) depending on wavelength, 
          // X angleIn and X angleOut allowences. Once made, the table is 
          // saved in instdir/grating.tbl. If the last grating table made 
          // (only latest made is saved) was  made using the identcle grating 
          // and table specifications then the existing table is loaded and 
          // a new table is not made. This saves loads of time.
          // **************
          
          if (surface.surfacetype[nsurf] == GRATING) {
            // This is where we instantiate  a Grating object.
            // surface.three is the number of ILUMINATED slits. If not defined 
            // ( < 1 ) then use all default values for the grating and table
            //  values which are specified in grating.h
            std::cout << "  Acquiring Grating Tables for optical surface " 
                      << nsurf << std::flush << std::endl;
            if ( int(surface.three[nsurf] ) < 1 ) {
              // Use default grating definition values and values from
              // grating.h: Note that the defualt values result is each axis of 
              // table having at least 3 elements.
              surface.ppGrating[nsurf] = new Grating(instrdir, random, 
                                                     minwavelength, 
                                                     maxwavelength);
            } 
            
            else {
              // Whenever the 2nd parameter argument in the optics file
              // (surface.three = number of illuminated slits)) is a positive 
              // non-zero value all 9 parameters must be specified. 
              // The min/max values for each axis of the table is  checked to 
              // result in at least 3 elements for each axis. In the case where 
              // they don't the max value is increased to result in at least 
              // 3 values for that axis. Note that the steps in the table for 
              // each axis is a fixed value set in grating.h. They have been 
              // carefully selected. Do not change.

              // Parameters are surface.:
              // two:   Blaze Angle Deg,        three: Number illuminated slits, 
              // four:  Slit spacing (nm),      five:  Min Wavelength nm,
              // six:   Maximum Wavelength,     seven: Minimum X AngleIn Deg,
              // eight: Maximum X AngleIn deg,  nine:  Minimum X AngleOut Deg
              // ten:   Maximum X AngleOut deg  
              // 
              //or else if params were given in the optics file we use:
              std::vector < double> gratingParameters;
              
              //Grating parameters
              gratingParameters.push_back(surface.two[nsurf]); // Blaze angle(deg)
              gratingParameters.push_back(surface.three[nsurf]);
              // Number of illuminated grooves
              gratingParameters.push_back(surface.four[nsurf]);  
              // Distance between grooves(uM).
            
              //Grating table parameters
              gratingParameters.push_back(surface.five[nsurf]);//Min wavelength(uM) 
              gratingParameters.push_back(surface.six[nsurf]); //Max wavelength(uM)
              gratingParameters.push_back(surface.seven[nsurf]);//Min X AngIn(deg)
              gratingParameters.push_back(surface.eight[nsurf]);//Max X AngIn(deg)
              gratingParameters.push_back(surface.nine[nsurf]);//Min X AngOut(deg)
              gratingParameters.push_back(surface.ten[nsurf]); //Max X AngOut(deg)
              
              // Using a vector here to aviod circular definitions of 
              // Gratin and Surface. The validity of the given values will be 
              // checked and if needed, modified,  to produce a usable table.
              surface.ppGrating[nsurf] = new Grating(instrdir, random, 
                                                     gratingParameters );
            }
          
            // Since grating must be absolutly flat zero out the surface parameters. 
            // so that the surface finding code sees the grating as flat.
            surface.two[nsurf]   = 0.0;
            surface.three[nsurf] = 0.0;
            surface.four[nsurf]  = 0.0;
            surface.five[nsurf]  = 0.0;
            surface.six[nsurf]   = 0.0;
            surface.seven[nsurf] = 0.0;
            surface.eight[nsurf] = 0.0;
            surface.nine[nsurf]  = 0.0;
            surface.ten[nsurf]   = 0.0;
          } 
          
          else {
            surface.ppGrating[nsurf] = NULL;  //This surface is not a grating
          }
          // ****************
          
          if (opdfile) surface.innerRadius[nsurf] = 0.0;

  
          if (surfacetype != "none") {
            surface.asphere(nsurf, SURFACE_POINTS);
            if (surface.surfacetype[nsurf] == MIRROR || surface.surfacetype[nsurf] == DETECTOR ||
                surface.surfacetype[nsurf] == LENS || surface.surfacetype[nsurf] == FILTER) {
              perturbation.zernikeflag.push_back(zernikemode);
            } else if (surface.surfacetype[nsurf] == EXITPUPIL) {
              perturbation.zernikeflag.push_back(0);
              perturbation.zernikeflag.push_back(0);
            } else {
              perturbation.zernikeflag.push_back(0);
            }
            
            if  (spaceMode < 0) {
              for (int i=0; i<nsurf; i++) {
                surface.height[i]=air.topAtmosphere*1e6-surface.height[i];
              }
            }


            // COATINGS
            int uncoatedflag=0;
            surface.surfacecoating[nsurf] = 0;
            if (coatingFile != "none" && coatingFile != "uncoated") {
                std::string tempFile;
                std::ifstream f(instrdir + "/" + coatingFile);
                if (f.good()) {
                    tempFile = instrdir + "/" + coatingFile;
                } else {
                    tempFile = instrdir + "/../standard/coating/" + coatingFile + ".txt";
                }
                readText coatingPars(tempFile);
                size_t nline = coatingPars.getSize();
                coating.allocate(nsurf, nline);
                int j = 0;
                int i = 0;
                int ttcc = 0;
                int flag = 0;
                double angle, angle0 = 0.0;
                for (size_t tt(0); tt<nline; tt++) {
                    std::istringstream issta(coatingPars[tt]);
                    std::string tempstring;
                    issta >> tempstring;
                    if (tempstring == "scattering") {
                        double a,b,c;
                        issta >> a >> b >> c;
                        perturbation.scatteringsigma[nsurf]=a;
                        perturbation.scatteringb[nsurf]=b;
                        perturbation.scatteringc[nsurf]=c;
                        flag++;
                    } else if (tempstring == "contamination") {
                        float a, b;
                        issta >> a >> b;
                        contamination.surfaceDustTransmission[nsurf]=a;
                        contamination.surfaceCondensationTransmission[nsurf]=b;
                        flag++;

                    } else if (tempstring == "uncoated") {
                        uncoatedflag=1;
                    } else if (tempstring == "monolayer") {
                        int nwave=0;
                        nwave = floor((maxwavelength - minwavelength)/(SED_TOLERANCE)) + 1;
                        int nangle = 100;
                        int nline = nwave*nangle;
                        coating.allocate(nsurf, nline);
                        coating.wavelengthNumber[nsurf]=nwave;
                        coating.angleNumber[nsurf]=nangle;
                        std::vector<double> refractionWavelength;
                        std::vector<double> refractionReal;
                        std::vector<double> refractionImag;
                        std::string monoFile;
                        std::string fileName;
                        float thickness;
                        issta >> fileName >> thickness;
                        //                        coating.nominalThickness[nsurf] = thickness;
                        monoFile = instrdir + "/../material/optical/"  + fileName + ".txt";
                        readText::readCol(monoFile,refractionWavelength, refractionReal, refractionImag);
                        int cc = 0;
                        for (int j =0; j < nangle; j++) {
                            double angle = PI/2.0*static_cast<double>(j)/static_cast<double>(nangle-1);
                            *(coating.angle[nsurf] + j) = angle/DEGREE;
                            for (int i = 0; i< nwave; i++) {
                                double wavelength;
                                wavelength = (minwavelength + SED_TOLERANCE*i)/1000.0;
                                *(coating.wavelength[nsurf] + i) = wavelength;
                                double index1 = 1.0;
                                double index1imag = 0.0;
                                if (nsurf != 0) {
                                    int index;
                                    if (surface.surfacemed[nsurf-1]==1) {
                                        find(medium.indexRefractionWavelength[nsurf-1],medium.indexRefractionNumber[nsurf-1],wavelength,&index);
                                        index1=interpolate(medium.indexRefraction[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                        index1imag=interpolate(medium.indexRefractionImag[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                    } else {
                                        index1=1.0;
                                        index1imag=0.0;
                                    }
                                }
                                //READ monolayer file
                                double index2 = 1.0;
                                double index2imag = 0.0;
                                int index;
                                find_v(refractionWavelength,wavelength,&index);
                                index2=interpolate_v(refractionReal,refractionWavelength,wavelength,index);
                                index2imag=interpolate_v(refractionImag,refractionWavelength,wavelength,index);
                                double index3 = 1.0;
                                double index3imag = 0.0;
                                if (surface.surfacemed[nsurf]==1) {
                                    find(medium.indexRefractionWavelength[nsurf],medium.indexRefractionNumber[nsurf],wavelength,&index);
                                    index3=interpolate(medium.indexRefraction[nsurf],medium.indexRefractionWavelength[nsurf],wavelength,index);
                                    index3imag=interpolate(medium.indexRefractionImag[nsurf],medium.indexRefractionWavelength[nsurf],wavelength,index);
                                } else {
                                    index3=1.0;
                                    index3imag=0.0;
                                }
                                double reflection = interfaceMonolayerAveragePolarization (angle, wavelength, index1, index1imag, index3, index3imag, thickness, index2, index2imag);
                                double transmission = 1.0 - reflection;
                                if (surfacetype == "mirror") {
                                    double swap;
                                    swap = transmission;
                                    transmission = reflection;
                                    reflection = swap;
                                }
                                *(coating.transmission[nsurf] + cc) = transmission;
                                *(coating.reflection[nsurf] + cc) = reflection;
                                cc++;
                            }
                        }
                        coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                        coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                        surface.surfacecoating[nsurf] = 1;
                        flag++;
                    } else if (tempstring == "protectedlayer") {
                        int nwave=0;
                        nwave = floor((maxwavelength - minwavelength)/(SED_TOLERANCE)) + 1;
                        int nangle = 100;
                        int nline = nwave*nangle;
                        coating.allocate(nsurf, nline);
                        coating.wavelengthNumber[nsurf]=nwave;
                        coating.angleNumber[nsurf]=nangle;
                        //READ monolayer file
                        std::vector<double> refractionWavelength;
                        std::vector<double> refractionReal;
                        std::vector<double> refractionImag;
                        std::vector<double> refraction2Wavelength;
                        std::vector<double> refraction2Real;
                        std::vector<double> refraction2Imag;
                        std::string monoFile;
                        std::string fileName;
                        std::string fileName2;
                        float thickness;
                        issta >> fileName >> fileName2 >>  thickness;
                        //                        coating.nominalThickness[nsurf] = thickness;
                        monoFile = instrdir + "/../material/optical/"  + fileName + ".txt";
                        readText::readCol(monoFile,refractionWavelength, refractionReal, refractionImag);
                        monoFile = instrdir + "/../material/optical/"  + fileName2 + ".txt";
                        readText::readCol(monoFile,refraction2Wavelength, refraction2Real, refraction2Imag);
                        int cc =0;
                        for (int j =0; j < nangle; j++) {
                            double angle = PI/2.0*static_cast<double>(j)/static_cast<double>(nangle-1);
                            *(coating.angle[nsurf] + j) = angle/DEGREE;
                            for (int i = 0; i< nwave; i++) {
                                double wavelength;
                                wavelength = (minwavelength + SED_TOLERANCE*i)/1000.0;
                                *(coating.wavelength[nsurf] + i) = wavelength;
                                double index1 = 1.0;
                                double index1imag = 0.0;
                                if (nsurf != 0) {
                                    int index;
                                    if (surface.surfacemed[nsurf-1]==1) {
                                        find(medium.indexRefractionWavelength[nsurf-1],medium.indexRefractionNumber[nsurf-1],wavelength,&index);
                                        index1=interpolate(medium.indexRefraction[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                        index1imag=interpolate(medium.indexRefractionImag[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                    } else {
                                        index1=1.0;
                                        index1imag=0.0;
                                    }
                                }
                                double index2 = 1.0;
                                double index2imag = 0.0;
                                int index;
                                find_v(refractionWavelength,wavelength,&index);
                                index2=interpolate_v(refractionReal,refractionWavelength,wavelength,index);
                                index2imag=interpolate_v(refractionImag,refractionWavelength,wavelength,index);

                                double index3 = 1.0;
                                double index3imag = 0.0;
                                find_v(refraction2Wavelength,wavelength,&index);
                                index3=interpolate_v(refraction2Real,refraction2Wavelength,wavelength,index);
                                index3imag=interpolate_v(refraction2Imag,refraction2Wavelength,wavelength,index);
                                double reflection = interfaceMonolayerAveragePolarization (angle, wavelength, index1, index1imag, index3, index3imag, thickness, index2, index2imag);
                                double transmission = 1.0 - reflection;
                                if (surfacetype == "mirror") {
                                    transmission = reflection;
                                    reflection = 0.0;
                                }
                                *(coating.transmission[nsurf] + cc) = transmission;
                                *(coating.reflection[nsurf] + cc) = reflection;
                                cc++;
                            }
                        }
                        coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                        coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                        surface.surfacecoating[nsurf] = 1;
                        flag++;
                    } else if (tempstring == "thicklayer") {
                        int nwave=0;
                        nwave = floor((maxwavelength - minwavelength)/(SED_TOLERANCE)) + 1;
                        int nangle = 100;
                        int nline = nwave*nangle;
                        coating.allocate(nsurf, nline);
                        coating.wavelengthNumber[nsurf]=nwave;
                        coating.angleNumber[nsurf]=nangle;
                        std::vector<double> refractionWavelength;
                        std::vector<double> refractionReal;
                        std::vector<double> refractionImag;
                        std::string monoFile;
                        std::string fileName;
                        issta >> fileName;
                        monoFile = instrdir + "/../material/optical/"  + fileName + ".txt";
                        readText::readCol(monoFile,refractionWavelength, refractionReal, refractionImag);
                        int cc =0;
                        for (int j =0; j < nangle; j++) {
                            double angle = PI/2.0*static_cast<double>(j)/static_cast<double>(nangle-1);
                            *(coating.angle[nsurf] + j) = angle/DEGREE;
                            for (int i = 0; i< nwave; i++) {
                                double wavelength;
                                wavelength = (minwavelength + SED_TOLERANCE*i)/1000.0;
                                *(coating.wavelength[nsurf] + i) = wavelength;
                                double index1 = 1.0;
                                double index1imag = 0.0;
                                if (nsurf != 0) {
                                    int index;
                                    if (surface.surfacemed[nsurf-1]==1) {
                                        find(medium.indexRefractionWavelength[nsurf-1],medium.indexRefractionNumber[nsurf-1],wavelength,&index);
                                        index1=interpolate(medium.indexRefraction[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                        index1imag=interpolate(medium.indexRefractionImag[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                    } else {
                                        index1=1.0;
                                        index1imag=0.0;
                                    }
                                }
                                double index2 = 1.0;
                                double index2imag = 0.0;
                                int index;
                                find_v(refractionWavelength,wavelength,&index);
                                index2=interpolate_v(refractionReal,refractionWavelength,wavelength,index);
                                index2imag=interpolate_v(refractionImag,refractionWavelength,wavelength,index);
                                double reflection = interfaceAveragePolarization(angle, index1, index1imag, index2, index2imag);
                                double transmission = 1.0 - reflection;
                                if (surfacetype == "mirror") {
                                    transmission = reflection;
                                    reflection = 0.0;
                                }
                                *(coating.transmission[nsurf] + cc) = transmission;
                                *(coating.reflection[nsurf] + cc) = reflection;
                                cc++;
                            }
                        }
                        coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                        coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                        surface.surfacecoating[nsurf] = 1;
                        flag++;
                    } else if (tempstring == "nonuniform") {
                        coating.nonUniformity[nsurf] = 1;
                        int nTerms = 0;
                        if ((surfacetype == "det") || (surfacetype == "detector")) {
                            coating.thickness[nsurf] = static_cast<float*>(calloc((floor((chip.nampx-1)/DET_RESAMPLE)+1)*(floor((chip.nampy-1)/DET_RESAMPLE)+1), sizeof(float)));
                            issta >> nTerms;
                            if (nTerms > NCHEB) nTerms=NCHEB;
                            double t[3];
                            float coeff;
                            for (int n=0; n<nTerms; n++) {
                                issta >> coeff;
                                if (n == 0) coating.nominalThickness[nsurf] = coeff;
                                for (int i=0; i< chip.nampy; i+= DET_RESAMPLE) {
                                    for (int j=0; j < chip.nampx; j+= DET_RESAMPLE) {
                                        double x = static_cast<double>(j-(chip.nampx-1)/2)/static_cast<double>((chip.nampx-1)/2);
                                        double y = static_cast<double>(i-(chip.nampy-1)/2)/static_cast<double>((chip.nampy-1)/2);
                                        chebyshevT2D(n,x,y,t);
                                        int location = (floor((chip.nampx - 1)/DET_RESAMPLE)+1)*
                                            floor((i)/DET_RESAMPLE) +
                                            floor((j)/DET_RESAMPLE);
                                        coating.thickness[nsurf][location] += t[0]*coeff;
                                    }
                                }
                            }
                            surface.surfacecoating[nsurf] = 1;
                        } else {
                            coating.thickness[nsurf] = static_cast<float*>(calloc(PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(float)));
                            issta >> nTerms;
                            if (nTerms > NZERN) nTerms=NZERN;
                            double t[4];
                            float coeff;
                            for (int n=0; n<nTerms; n++) {
                                issta >> coeff;
                                if (n == 0) coating.nominalThickness[nsurf] = coeff;
                                for (int i=0; i< PERTURBATION_POINTS; i++) {
                                    for (int j=0; j < PERTURBATION_POINTS; j++) {
                                        double x, y;
#ifdef PERTURBATION_COORDINATE
                                        double r, phi;
                                        r = (static_cast<double>(j))/(static_cast<double>(PERTURBATION_POINTS) - 1);
                                        phi= (static_cast<double>(i))/(static_cast<double>(PERTURBATION_POINTS))*2*PI;
                                        x = r*sin(phi);
                                        y= r*cos(phi);
#else
                                        x= -1.0 + (static_cast<double>(i))/(static_cast<double>(PERTURBATION_POINTS) - 1)*2.0;
                                        y= -1.0 + (static_cast<double>(j))/(static_cast<double>(PERTURBATION_POINTS) - 1)*2.0;
#endif
                                        zernikeValue(n,x,y,t);
                                        int location = i*PERTURBATION_POINTS + j;
                                        coating.thickness[nsurf][location] += coeff*t[0]*t[1];
                                    }
                                }
                            }
                            surface.surfacecoating[nsurf] = 1;
                        }
                        flag++;
                    } else {
                        std::istringstream isst(coatingPars[tt]);
                        isst >> angle;
                        if (i == 0 && j == 0) angle0 = angle;
                        if (angle > angle0) {
                            i++;
                            coating.wavelengthNumber[nsurf] = j;
                            if (j != coating.wavelengthNumber[nsurf] && i > 1) {
                                std::cout << "Error in format of " << coatingFile << std::endl;
                                exit(1);
                            }
                            j = 0;
                        }
                        *(coating.angle[nsurf] + i) = angle;
                        isst >> *(coating.wavelength[nsurf] + j);
                        isst >> *(coating.transmission[nsurf] + ttcc);
                        isst >> *(coating.reflection[nsurf] + ttcc);
                        ttcc++;
                        j++;
                        angle0 = angle;
                    }
                }
                i++;
                if ((nline-flag)>1) {
                    coating.wavelengthNumber[nsurf] = j;
                    coating.angleNumber[nsurf] = i;
                    //special approximate shift for incomplete numerical sets
                    if (i==1) {
                        int cc = 0;
                        double angleOrig = *(coating.angle[nsurf])/DEGREE;
                        // copy original
                        double *tWavelength;
                        double *tReflection;
                        double *tTransmission;
                        int nwave = coating.wavelengthNumber[nsurf];
                        tWavelength = static_cast<double*>(calloc(coating.wavelengthNumber[nsurf], sizeof(double)));
                        tReflection = static_cast<double*>(calloc(coating.wavelengthNumber[nsurf], sizeof(double)));
                        tTransmission = static_cast<double*>(calloc(coating.wavelengthNumber[nsurf], sizeof(double)));
                        for (int k=0;k<j;k++) {
                            tWavelength[k] = *(coating.wavelength[nsurf] + k);
                            tReflection[k] = *(coating.reflection[nsurf] + k);
                            tTransmission[k] = *(coating.transmission[nsurf] + k);
                        }
                        int nangle = 100;
                        nline = nwave*nangle;
                        coating.wavelengthNumber[nsurf]=nwave;
                        coating.angleNumber[nsurf]=nangle;
                        coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], nline*sizeof(double)));
                        coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], nline*sizeof(double)));
                        coating.reflection[nsurf] = static_cast<double*>(realloc(coating.reflection[nsurf], nline*sizeof(double)));
                        coating.transmission[nsurf] = static_cast<double*>(realloc(coating.transmission[nsurf], nline*sizeof(double)));
                        for (int l =0; l < nangle; l++) {
                            double angle = PI/2.0*static_cast<double>(l)/static_cast<double>(nangle-1);
                            *(coating.angle[nsurf] + l) = angle/DEGREE;
                            for (int k =0; k< nwave; k++) {
                                *(coating.wavelength[nsurf] + k) = tWavelength[k];
                                double neff = 1.7;
                                int index;
                                find(tWavelength,k,tWavelength[k]*sqrt(1+pow(sin(angleOrig)/neff, 2.0))/sqrt(1+pow(sin(angle)/neff, 2.0)),&index);
                                *(coating.reflection[nsurf] + cc) = tReflection[index];
                                *(coating.transmission[nsurf] + cc) = tTransmission[index];
                                cc++;
                            }
                        }
                    }
                    coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                    coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                    surface.surfacecoating[nsurf] = 1;
                }

                  
            }
            
            // MEDIUM
            surface.surfacemed[nsurf] = 0;
            if (mediumFile == "air") surface.surfacemed[nsurf] = 2;
            if (mediumFile != "vacuum" && mediumFile != "air") {
                std::ifstream f(instrdir + "/" + mediumFile);
                if (f.good()) {
                    medium.setup(nsurf, instrdir + "/" + mediumFile);
                } else {
                    medium.setup(nsurf, instrdir + "/../material/optical/" + mediumFile + ".txt");
                }
                surface.surfacemed[nsurf] = 1;
            }

            //Calculate interface for non-coatings
            // note:  for mirrors we don't need to do this since if its truly non-coated then we won't get the reflection and we
            // have no way of specifying what the medium of the actual mirror is (in many cases glass), so this is handled above
            if (((coatingFile == "uncoated") || (uncoatedflag==1))&& ((surfacetype == "lens") || (surfacetype == "filter"))) {
                int nwave=0;
                nwave = floor((maxwavelength - minwavelength)/(SED_TOLERANCE)) + 1;
                int nangle = 100;
                int nline = nwave*nangle;
                coating.allocate(nsurf, nline);
                coating.wavelengthNumber[nsurf]=nwave;
                coating.angleNumber[nsurf]=nangle;
                int tt =0;
                for (int j =0; j < nangle; j++) {
                    double angle = PI/2.0*static_cast<double>(j)/static_cast<double>(nangle-1);
                    *(coating.angle[nsurf] + j) = angle/DEGREE;
                    for (int i = 0; i< nwave; i++) {
                        double wavelength;
                        wavelength = (minwavelength + SED_TOLERANCE*i)/1000.0;
                        double index1 = 1.0;
                        double index1imag = 0.0;
                        if (nsurf != 0) {
                            int index;
                            if (surface.surfacemed[nsurf-1]==1) {
                                find(medium.indexRefractionWavelength[nsurf-1],medium.indexRefractionNumber[nsurf-1],wavelength,&index);
                                index1=interpolate(medium.indexRefraction[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                                index1imag=interpolate(medium.indexRefractionImag[nsurf-1],medium.indexRefractionWavelength[nsurf-1],wavelength,index);
                            } else {
                                index1=1.0;
                                index1imag=0.0;
                            }
                        }
                        int index;
                        double index2 = 1.0;
                        double index2imag = 0.0;
                        if (surface.surfacemed[nsurf]==1) {
                            find(medium.indexRefractionWavelength[nsurf],medium.indexRefractionNumber[nsurf],wavelength,&index);
                            index2=interpolate(medium.indexRefraction[nsurf],medium.indexRefractionWavelength[nsurf],wavelength,index);
                            index2imag=interpolate(medium.indexRefractionImag[nsurf],medium.indexRefractionWavelength[nsurf],wavelength,index);
                        } else {
                            index2=1.0;
                            index2imag=0.0;
                        }
                        double reflection = interfaceAveragePolarization(angle, index1, index1imag, index2, index2imag);
                        double transmission = 1.0 - reflection;
                        *(coating.wavelength[nsurf] + tt) = wavelength;
                        *(coating.transmission[nsurf] + tt) = transmission;
                        *(coating.reflection[nsurf] + tt) = reflection;
                        tt++;
                    }
                }
                    coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                    coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                    surface.surfacecoating[nsurf] = 1;
            }

                // for (int i=0; i < coating.angleNumber[nsurf]; i++) {
                //     for (int j=0; j < coating.wavelengthNumber[nsurf]; j++) {
                //         fprintf(outfile,"%d %f %f %f %f\n",nsurf,coating.angle[nsurf][i], coating.wavelength[nsurf][j], coating.reflection[nsurf][i*coating.wavelengthNumber[nsurf]+j],
                //                coating.transmission[nsurf][i*coating.wavelengthNumber[nsurf]+j]);
                //     }
                // }

            nsurf++;
          }
        }

    }

                    // fclose(outfile);

    // if (opdfile) nsurf -= 2; //exit pupil and focalplane

    if (telescopeMode == 1 || ((telescopeMode == 0) && (setApertureFlag==0))) {
        maxr = surface.outerRadius[0];
        minr = surface.innerRadius0[0];
    }
    maxr += 40.0;
    minr -= 40.0;
    if (minr < 0) minr=0.0;
    straylightcut += 2.5*log10((maxr*maxr-minr*minr)/(4180.0*4180.0-2558.0*2558.0));

    // OBSTRUCTIONS
    obstruction.nspid = 0;
    if (diffractionMode >= 1) {
        std::string sss;
        sss = instrdir + "/structure.txt";
        std::ifstream inStream(sss.c_str());
        if (inStream) {
            readText spiderPars(instrdir + "/structure.txt");
            std::cout << "Placing Obstructions." << std::endl;
            for (size_t t(0); t < spiderPars.getSize(); t++) {
                std::istringstream iss(spiderPars[t]);
                std::string spiderType;
                iss >> spiderType;
                if (spiderType == "pupilscreen"){
                    pupilscreenMode = 1;
                }
                if ((spiderType == "crossx") || (spiderType == "crossy")) {
                    iss >> obstruction.height[obstruction.nspid];
                    iss >> obstruction.width[obstruction.nspid];
                    iss >> obstruction.center[obstruction.nspid];
                    double tempf4;
                    iss >> tempf4;
                    obstruction.depth[obstruction.nspid] = tempf4;
                    iss >> tempf4;
                    obstruction.angle[obstruction.nspid] = tempf4;
                    iss >> tempf4;
                    obstruction.reference[obstruction.nspid] = tempf4;
                    if (spiderType == "crossx") obstruction.type[obstruction.nspid] = 1;
                    if (spiderType == "crossy") obstruction.type[obstruction.nspid] = 2;
                    obstruction.nspid++;
                }
            }
        }

   
    }


    // PERTURBATIONS
    fprintf(stdout, "Perturbing Design.\n");
    perturbation.zernike_coeff = static_cast<double*>(calloc(NTERM*nsurf, sizeof(double)));
    int perturbation_coordinate=1;
    #ifdef PERTURBATION_COORDINATE
    perturbation_coordinate=2;
    #endif
    // chuck's old model: 0.82312e-3*pow(zv, -1.2447)

    for (int i = 0; i < nsurf; i++) {
        for (int j = 0; j < NTERM; j++) {
            *(perturbation.zernike_coeff + i*NTERM + j) = izernike[i][j];
        }
    }

    perturbation.eulerPhi.reserve(totSurf);
    perturbation.eulerPsi.reserve(totSurf);
    perturbation.eulerTheta.reserve(totSurf);
    perturbation.decenterX.reserve(totSurf);
    perturbation.decenterY.reserve(totSurf);
    perturbation.defocus.reserve(totSurf);
    perturbation.rmax.reserve(totSurf);

    double *zernikeR, *zernikePhi, *zernikeNormalR, *zernikeNormalPhi;
    double *chebyshev, *chebyshevR, *chebyshevPhi;
    double *poly, *polyR, *polyPhi;
    perturbation.zernike_r_grid = static_cast<double*>(calloc(PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_phi_grid = static_cast<double*>(calloc(PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_x_grid = static_cast<double*>(calloc(PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_y_grid = static_cast<double*>(calloc(PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_a_grid = static_cast<double*>(calloc(PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_b_grid = static_cast<double*>(calloc(PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed_flag = static_cast<int*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(int)));
    perturbation.zernike_summed = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed_nr_p = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
    perturbation.zernike_summed_np_r = static_cast<double*>(calloc(nsurf*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));

    setupGrids(perturbation.zernike_r_grid,perturbation.zernike_phi_grid,perturbation.zernike_x_grid,perturbation.zernike_y_grid,
               perturbation.zernike_a_grid, perturbation.zernike_b_grid, PERTURBATION_POINTS, perturbation_coordinate);

    int initZern=0, initCheb=0, initPoly=0;

    for (int k = 0; k < nsurf; k++) {
        if (pertType[k] == "zern") {
            if (initZern==0) {
                zernikeR = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                zernikePhi = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                zernikeNormalR = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                zernikeNormalPhi = static_cast<double*>(calloc(NZERN*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                zernikes(zernikeR, zernikePhi, perturbation.zernike_r_grid,
                         perturbation.zernike_phi_grid, zernikeNormalR, zernikeNormalPhi,
                         PERTURBATION_POINTS, NZERN, perturbation_coordinate);
                initZern=1;
            }
            for (int j = 0; j < PERTURBATION_POINTS; j++) {
                for (int l = 0; l < PERTURBATION_POINTS; l++) {
                    *(perturbation.zernike_summed_flag + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    *(perturbation.zernike_summed + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    *(perturbation.zernike_summed_nr_p + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    *(perturbation.zernike_summed_np_r + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) = 0;
                    for (int m = 0; m < NZERN; m++) {
                        *(perturbation.zernike_summed + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeR + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j))*
                            (*(zernikePhi + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j));

                        *(perturbation.zernike_summed_nr_p + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeNormalR + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j))*
                            (*(zernikePhi + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j));

                        *(perturbation.zernike_summed_np_r + k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j) +=
                            *(perturbation.zernike_coeff + k*NTERM + m)*(*(zernikeR + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j))*
                            (*(zernikeNormalPhi + m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j));

                    }
                }
            }
        } else if (pertType[k] == "chebyshev") {
            if (initCheb==0) {
                chebyshev = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                chebyshevR = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                chebyshevPhi = static_cast<double*>(calloc(NCHEB*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                chebyshevs (perturbation.zernike_r_grid, perturbation.zernike_phi_grid,
                            chebyshev, chebyshevR, chebyshevPhi, PERTURBATION_POINTS, NCHEB, perturbation_coordinate);
                initCheb=1;
            }
            for (int j = 0; j < PERTURBATION_POINTS; j++) {
                for (int l = 0; l < PERTURBATION_POINTS; l++) {
                    int idx = k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                    perturbation.zernike_summed_flag[idx] = 0;
                    perturbation.zernike_summed[idx] = 0;
                    perturbation.zernike_summed_nr_p[idx] = 0;
                    perturbation.zernike_summed_np_r[idx] = 0;
                    for (int m = 0; m < NCHEB; m++) {
                        int idxi = m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                        perturbation.zernike_summed[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshev[idxi];
                        perturbation.zernike_summed_nr_p[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshevR[idxi];
                        perturbation.zernike_summed_np_r[idx] += perturbation.zernike_coeff[k*NTERM + m]*chebyshevPhi[idxi];
                    }
                }
            }
        } else if (pertType[k] == "poly") {
            if (initPoly==0) {
                poly = static_cast<double*>(calloc(NPOLY*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                polyR = static_cast<double*>(calloc(NPOLY*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                polyPhi = static_cast<double*>(calloc(NPOLY*PERTURBATION_POINTS*PERTURBATION_POINTS, sizeof(double)));
                polys (perturbation.zernike_r_grid, perturbation.zernike_phi_grid,
                       poly, polyR, polyPhi, PERTURBATION_POINTS, NPOLY, perturbation_coordinate);
                initPoly=1;
            }
            for (int j = 0; j < PERTURBATION_POINTS; j++) {
                for (int l = 0; l < PERTURBATION_POINTS; l++) {
                    int idx = k*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                    perturbation.zernike_summed_flag[idx] = 0;
                    perturbation.zernike_summed[idx] = 0;
                    perturbation.zernike_summed_nr_p[idx] = 0;
                    perturbation.zernike_summed_np_r[idx] = 0;
                    for (int m = 0; m < NPOLY; m++) {
                        int idxi = m*PERTURBATION_POINTS*PERTURBATION_POINTS + l*PERTURBATION_POINTS + j;
                        perturbation.zernike_summed[idx] += perturbation.zernike_coeff[k*NTERM + m]*poly[idxi];
                        perturbation.zernike_summed_nr_p[idx] += perturbation.zernike_coeff[k*NTERM + m]*polyR[idxi];
                        perturbation.zernike_summed_np_r[idx] += perturbation.zernike_coeff[k*NTERM + m]*polyPhi[idxi];
                    }
                }
            }
        }
    }
    if (initZern==1) {
        free(zernikeR);
        free(zernikePhi);
        free(zernikeNormalR);
        free(zernikeNormalPhi);
    }
    if (initCheb==1) {
        free(chebyshev);
        free(chebyshevR);
        free(chebyshevPhi);
    }
    if (initPoly==1) {
        free(poly);
        free(polyR);
        free(polyPhi);
    }
 
    perturbation.rotationmatrix = static_cast<double*>(calloc(totSurf*3*3, sizeof(double)));
    perturbation.inverserotationmatrix = static_cast<double*>(calloc(totSurf*3*3, sizeof(double)));
    for (int i = 0; i< nsurf + 1; i++) {
        perturbation.eulerPhi[i] = body[i][0];
        perturbation.eulerPsi[i] = body[i][1];
        perturbation.eulerTheta[i] = body[i][2];
        perturbation.decenterX[i] = body[i][3];
        perturbation.decenterY[i] = body[i][4];
        perturbation.defocus[i] = body[i][5];
    }
    satbuffer += static_cast<int>(2.0*fabs(perturbation.defocus[nsurf]/(pixsize)));
    perturbation.decenterX[nsurf] += centerx*1e-3;
    perturbation.decenterY[nsurf] += centery*1e-3;

    surface.centerx[nsurf] = 0.0;
    surface.centery[nsurf] = 0.0;
    surface.rmax[nsurf] = sqrt(pixelsx*pixelsx + pixelsy*pixelsy)*sqrt(2.0)/2.0*pixsize;
    surface.height[nsurf] =  surface.height[nsurf - 1];
    surface.intz[nsurf] =  surface.intz[nsurf - 1];
    for (int i = 0; i < totSurf; i++) {
        *(perturbation.rotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 1) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 0) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 2) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 1) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 1) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 0) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 2) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 1) = sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
    }
    for (int i = 0; i < totSurf; i++) {
        perturbation.rmax[i] = surface.rmax[surfaceLink[i]];
    }


    /* Update perturbation model from FEA file */
    size_t nNear = 4;
    // int degree=1;
    int leafSize = 10;
    for (int k = 0; k < nsurf; k++) {
        for (int m = 0; m < feaflag[k]; m++) {
            std::cout<<"Loading "<<feafile[k*MAX_SURF*2 + m]<<std::endl;
            std::vector<double> ix(PERTURBATION_POINTS*PERTURBATION_POINTS), iy(PERTURBATION_POINTS*PERTURBATION_POINTS);
            for (int l = 0; l < PERTURBATION_POINTS; l++) {
                for (int j = 0; j < PERTURBATION_POINTS; j++) {
// #ifdef PERTURBATION_COORDINATE
                    ix[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_r_grid[l*PERTURBATION_POINTS + j]*cos(perturbation.zernike_phi_grid[l*PERTURBATION_POINTS + j]);
                    iy[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_r_grid[l*PERTURBATION_POINTS + j]*sin(perturbation.zernike_phi_grid[l*PERTURBATION_POINTS + j]);
// #else
//                     ix[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_x_grid[l*PERTURBATION_POINTS + j];
//                     iy[l*PERTURBATION_POINTS + j] = perturbation.rmax[k]*perturbation.zernike_y_grid[l*PERTURBATION_POINTS + j];
// #endif
                    
                }
            }
            long sizeFea;
            Fea feaTree(feafile[k*MAX_SURF*2 + m], leafSize, surface, perturbation, k, feascaling[k*MAX_SURF*2 + m], &sizeFea);
            nNear = round(0.5*sqrt(static_cast<double>(sizeFea)));
            if (nNear < 4) nNear = 4;
                        feaTree.knnQueryFitDegree1(ix, iy, &perturbation.zernike_summed[k*PERTURBATION_POINTS*PERTURBATION_POINTS],
                                                   &perturbation.zernike_summed_nr_p[k*PERTURBATION_POINTS*PERTURBATION_POINTS],
                                                   &perturbation.zernike_summed_np_r[k*PERTURBATION_POINTS*PERTURBATION_POINTS],
                                                   &perturbation.zernike_summed_flag[k*PERTURBATION_POINTS*PERTURBATION_POINTS], nNear, perturbation_coordinate);

        }
    }

    /// Large Angle Scattering
    perturbation.scatteringprofile = static_cast<double*>(calloc(SCATTERING_POINTS*nsurf, sizeof(double)));
    for (int k =0; k < nsurf; k++) {
        for (int j = 0; j < SCATTERING_POINTS; j++) {
            double angle = static_cast<float>(j)/(static_cast<double>(SCATTERING_POINTS))*SCATTERING_LIMIT;
            perturbation.scatteringprofile[SCATTERING_POINTS*k+j] = 1.0/pow(1.0 + angle*angle,(perturbation.scatteringc[k]+1.0)/2.0)*angle;
        }
        double tempf1 = 0.0;
        for (int j = 0; j < SCATTERING_POINTS; j++) {
            tempf1+= perturbation.scatteringprofile[SCATTERING_POINTS*k + j];
        }
        for (int j = 0; j < SCATTERING_POINTS; j++) {
            perturbation.scatteringprofile[SCATTERING_POINTS*k+j] /= tempf1;
        }
        for (int j = 1; j < SCATTERING_POINTS; j++) {
            perturbation.scatteringprofile[SCATTERING_POINTS*k+j] +=perturbation.scatteringprofile[SCATTERING_POINTS*k + j - 1];
        }
    }

    // Tracking
    double angle = PI/2 - zenith;
    if (fabs(angle) < 0.5*PI/180.0) angle=0.5*PI/180.0;
    if (spaceMode <= 0) {
        rotationrate = 15.04*cos(latitude)*cos(azimuth)/cos(angle);
    } else {
        rotationrate = 0.0;
    }
    if (trackingfile == ".") trackingMode = 0;

    std::ifstream inStream(trackingfile.c_str());
    if (inStream) {
        readText trackingPars(trackingfile);
        size_t j = trackingPars.getSize();
        perturbation.jitterrot = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterele = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterazi = static_cast<double*>(calloc(j, sizeof(double)));
        screen.jitterwind = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jittertime = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.windshake = static_cast<double*>(calloc(j, sizeof(double)));
        trackinglines = j;
        for (size_t t(0); t < j; t++) {
            std::istringstream iss(trackingPars[t]);
            iss >> perturbation.jittertime[t];
            double f1, f2, f3, f4, f5;
            iss >> f1 >> f2 >> f3 >> f4 >> f5;
            perturbation.jitterele[t] = f1*elevationjitter*ARCSEC;
            perturbation.jitterazi[t] = f2*azimuthjitter*ARCSEC;
            perturbation.jitterrot[t] = f3*rotationjitter*ARCSEC;
            perturbation.windshake[t] = f5*windshake*ARCSEC;
            // printf("%e\n",perturbation.windshake[t]);
            if (trackingMode == 0) {
                perturbation.jitterele[t] = 0.0;
                perturbation.jitterazi[t] = 0.0;
                perturbation.jitterrot[t] = 0.0;
            }
            if (windShakeMode == 0) {
                perturbation.windshake[t] = 0.0;
            }
            screen.jitterwind[t] = f4*windjitter*DEGREE;
        }
    }

    /// Setup chip material

        std::string dir2;
        if (flatdir == 0) dir2 = datadir + "/material";
        else if (flatdir == 1) dir2 = ".";

    fprintf(stdout, "Electrifying Devices.\n");
    if (diagnostic) printf("$$%40s %22e %22.6lf\n","sensorSetup",(double)1,655.0);
    silicon.sensorSetup(sensorfilename, devmaterial, sensorTempNominal + sensorTempDelta, nbulk,
                  nf, nb, sf, sb, sensorthickness, overdepbias,
                  instrdir, dir2, chip.nampx, chip.nampy, pixsize, seedchip,
                  &impurityX, &impurityY, impurityvariation, minwavelength,
                  maxwavelength,activeBuffer, dopingflag, esatflag, vsatflag, detectorMode, siliconDebug, siliconSubSteps,obshistid);
    if (silicon.xTransferFlag == 1) {
        if (silicon.yTransferFlag ==1) {
            chip.midpoint = static_cast<int>(sqrt(pixelsx*pixelsx+pixelsy*pixelsy));
        } else {
            chip.midpoint = pixelsx/2;
        }
    } else {
        chip.midpoint = pixelsy/2;
    }

    if (contaminationmode == 1) {
        fprintf(stdout, "Contaminating Surfaces.\n");
        contamination.setup(surface.innerRadius, surface.outerRadius, nsurf - 1, PERTURBATION_POINTS,
                            pixsize*DET_RESAMPLE, floor((chip.nampx-1)/DET_RESAMPLE)+1, floor((chip.nampy-1)/DET_RESAMPLE)+1, qevariation, seedchip);
    }

    fprintf(stdout, "Diffracting.\n");
    pthread_mutex_init(&lock.lock1, NULL);
    pthread_mutex_init(&lock.lock3, NULL);
             {

        Photon photon;
        Vector position;
        Vector angle;
        //double shiftedAngle;
        double r,  phi;
        double rindex;
        int index;
        double radius;
        int atmtempdebug = 0;
        double bestvalue;
        double bestscale = 1.0;
        Vector largeAngle;
        double stdx = 0.0, stdy = 0.0;
        double radx = 0.0, rady = 0.0;
        double ncount = 0.0;
        screen.hffunc = static_cast<double*>(calloc(SCREEN_SIZE*2, sizeof(double)));
        screen.hffunc_n = static_cast<double*>(calloc(SCREEN_SIZE*2, sizeof(double)));
        photon.thread = 0;
        atmtempdebug = atmdebug;

        for (int l = 0; l <= atmtempdebug; l++) {
            atmdebug = l;

            for (int i = 0; i < SCREEN_SIZE; i++) {
                for (int j = 0; j < SCREEN_SIZE; j++) {
                    screen.tfocalscreen[i*SCREEN_SIZE + j] = 0;
                }
            }

            photon.prtime = -1.0;

            for (int k = 0; k < 5; k++) {


                photon.wavelength = 0.5;
                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                photon.time = random[0].uniform()*exptime;
                photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;
                photon.thread = 0;
                //shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                r = sqrt(random[0].uniform()*(maxr*maxr - minr*minr) + minr*minr);
                phi = random[0].uniform()*2*PI;
                position.x = r*cos(phi);
                position.y = r*sin(phi);
                index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                photon.xp = position.x;
                photon.yp = position.y;
                angle.x = 0.0;
                angle.y = 0.0;
                angle.z = -1.0;

                atmospherePropagate(&position, angle, -1, 2, &photon);
                for (int layer = 0; layer < natmospherefile; layer++) {
                    atmospherePropagate(&position, angle, layer, 2, &photon);
                    atmosphereIntercept(&position, layer, &photon);
                    atmosphereRefraction(&angle, layer, 3, &photon);
                }
                atmosphereDiffraction(&angle, &photon);
                for (int i = 0; i < SCREEN_SIZE; i++) {
                    for (int j = 0; j < SCREEN_SIZE; j++) {
                        screen.tfocalscreen[i*SCREEN_SIZE + j] += screen.focalscreen[i*SCREEN_SIZE + j];
                    }
                }


            }



            for (int i = 0; i < SCREEN_SIZE*2; i++) {
                *(screen.hffunc + i) = 0.0;
            }
            for (int i = 0; i < SCREEN_SIZE*2; i++) {
                *(screen.hffunc_n + i) = 0.0;
            }
            for (int i = 0; i < SCREEN_SIZE; i++) {
                for (int j = 0; j < SCREEN_SIZE; j++) {
                    radius = sqrt(pow((i - (SCREEN_SIZE/2 + 1)), 2) + pow((j - (SCREEN_SIZE/2 + 1)), 2));
                    *(screen.hffunc + ((int)(radius))) += screen.tfocalscreen[i*SCREEN_SIZE + j];
                    *(screen.hffunc_n + ((int)(radius))) += 1;
                }
            }

            for (int j = 0; j < SCREEN_SIZE*2; j++) {
                if (screen.hffunc_n[j] < 1) {
                    *(screen.hffunc + j) = 0;
                } else {
                    *(screen.hffunc + j) = *(screen.hffunc + j)/(screen.hffunc_n[j])*((double)(j));
                }
            }
            double tempf1 = 0.0;
            for (int j = 0; j < SCREEN_SIZE*2; j++) {
                tempf1 = tempf1 + *(screen.hffunc + j);
            }
            for (int j = 0; j < SCREEN_SIZE*2; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j)/tempf1;
            }
            for (int j = 1; j < SCREEN_SIZE*2; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j) + *(screen.hffunc + j - 1);
            }


            if (l == 0) {

                    bestscale = 1.0;
                    for (int run = 0; run < 2; run++) {
                        double lowscale, highscale, stepscale;
                        bestvalue = 1e10;
                        if (run == 0) {
                            lowscale=0.0;
                            highscale=2.0+0.1;
                            stepscale= 0.1;
                        } else {
                            lowscale=bestscale-0.1;
                            highscale=bestscale+0.1+0.01;
                            stepscale = 0.01;
                        }
                        for (double scales = lowscale; scales <= highscale; scales += stepscale) {
                            stdx = 0.0;
                            stdy = 0.0;
                            ncount = 0.0;

                            for (int k = 0; k < 10000; k++) {

                                screen.secondKickSize = scales;
                                photon.wavelength = 0.5;
                                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                                photon.time = random[0].uniform()*1000.0;
                                photon.absoluteTime = photon.time - 1000.0/2.0 + timeoffset;
                                //shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                                r = sqrt(random[0].uniform()*(maxr*maxr - minr*minr) + minr*minr);
                                phi = random[0].uniform()*2*PI;
                                position.x = r*cos(phi);
                                position.y = r*sin(phi);
                                index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                                position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                                photon.xp = position.x;
                                photon.yp = position.y;
                                angle.x = 0.0;
                                angle.y = 0.0;
                                angle.z = -1.0;
                                largeAngle.x = 0;
                                largeAngle.y = 0;

                                secondKick(&largeAngle, &photon);
                                atmospherePropagate(&position, angle, -1, 1, &photon);
                                for (int layer = 0; layer < natmospherefile; layer++) {
                                    atmospherePropagate(&position, angle, layer, 1, &photon);
                                    atmosphereIntercept(&position, layer, &photon);
                                    atmosphereRefraction(&angle, layer, 1, &photon);
                                }
                                angle.x = angle.x + largeAngle.x;
                                angle.y = angle.y + largeAngle.y;

                                radx = angle.x/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*photon.wavelengthFactor);
                                rady = angle.y/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*photon.wavelengthFactor);
                                stdx += radx*radx*exp(-(radx*radx + rady*rady)/2.0/1.0);
                                stdy += rady*rady*exp(-(radx*radx + rady*rady)/2.0/1.0);
                                ncount += exp(-(radx*radx + rady*rady)/2.0/1.0);

                            }

                            stdx /= ncount;
                            stdy /= ncount;
                            screen.secondKickSize = sqrt(stdx + stdy);
                            if (fabs(screen.secondKickSize - 1.0) < bestvalue) {
                                bestvalue = fabs(screen.secondKickSize - 1.0);
                                bestscale = scales;
                            }

                        }
                        screen.secondKickSize = bestscale;
                    }

            }
        }
        atmdebug = atmtempdebug;

             }
    if (areaExposureOverride == 0) {
        nphot = static_cast<long>(PI*(maxr/10.0*maxr/10.0 - minr/10.0*minr/10.0)*exptime*(totalnorm/H_CGS));
    } else {
        nphot = static_cast<long>(1e6*(totalnorm/H_CGS));
    }
    if (opdfile == 1) {
        nphot = totalnorm;
    }

    settings();

    fprintf(stdout, "Number of Sources: %d\n", nsource);
    fprintf(stdout, "Photons: %6.2e  Photons/cm^2/s: %6.2e\n", static_cast<double>(nphot), totalnorm/H_CGS);

    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Photon Raytrace\n");
    readText versionPars(bindir + "/version");
    fprintf(stdout, "%s\n", versionPars[0].c_str());
    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");

    return(0);

}
