///
/// @package phosim
/// @file silicon.cpp
/// @brief silicon class
///
/// @brief Created by:
/// @author Andy Rasmussen (SLAC)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <fitsio.h>
#include <fitsio2.h>

#include "silicon.h"
#include "constants.h"
#include "parameters.h"
#include "helpers.h"
#include "ancillary/readtext.h"
#include "raytrace/raytrace.h"

double Silicon::absorptionCoeffMCT (double lambda, double temperature, double x) {

    // Hg(1-x)Cd(x)Te follows Chu et al. J. Appl. Phys. 75 1234-5. (Kane)
    // and Finkman and Schacham J. Appl. Phys. 56 2896-900 (Urbach).
    double energy = 12398/lambda; // eV
    double alpha_o = exp(53.61*x - 18.88);
    double E_o = -0.3424 + 1.838*x + 0.148*pow(x, 2);
    double T_o = 81.9; // K
    double sigma = 3.267e4*(1 + x); // K eV-1
    double alpha_T = 100 + 5000*x;
    double E_T = (T_o + temperature)/sigma*log(alpha_T/alpha_o) + E_o; // eV
    double E_g = -0.302 + 1.93*x + 5.35e-4 * temperature * (1 - 2*x) - 0.81 * pow(x, 2) + 0.832*pow(x, 3); // eV  Hanson et al.
    double beta = alpha_T/pow(E_T - E_g, 0.5);
    double alpha = 0.0;

    if (energy < E_g) { // Urbach
        alpha = alpha_o*exp(sigma*(energy - E_o)/(temperature + T_o));
    } else { // Kane region (unsure about best model)
        alpha = beta*pow(energy - E_g, 0.5);
    }

    // alpha is in cm-1	
    return(alpha);

}


double Silicon::absorptionCoeffSi (double lambda, double temperature) {
    
    // follows Rajkanan et al. Solid-state electronics 22, pp793-795.
    double egT[2], eg0[] = {1.1557f, 2.5f}; // eV
    double edgT = 0.0, egd0 = 3.2f;         // eV
    double ep[] = {1.827e-2f, 5.773e-2f};   // eV
    double C[] = {5.5f, 4.0f};              // no dimension
    double A[] = {3.231e2f, 7.237e3f};      // cm-1 eV-2
    double ad = 1.052e6f;                   // cm-1 eV-2
    double kBoltzmann = 8.617e-5f;          // eV K-1
    double eV;                              // eV
    double beta = 7.021e-4f;                // eV K-1
    double gamma = 1108;                    // K
    double alpha = 0.0;                     // cm-1
    double deltaE0, deltaE1[2][2][2];       // eV

    edgT   = egd0   - (beta*temperature*temperature/(temperature + gamma));
    egT[0] = eg0[0] - (beta*temperature*temperature/(temperature + gamma));
    egT[1] = eg0[1] - (beta*temperature*temperature/(temperature + gamma));

    if (lambda > 3100) {
        eV = 12398/lambda;
        deltaE0 = eV - edgT;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            deltaE1[i][j][0] = eV - egT[j] + ep[i];
            deltaE1[i][j][1] = eV - egT[j] - ep[i];
        }
    }

    alpha = ad*sqrt((deltaE0>0)?deltaE0:0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            alpha += C[i]*A[j]*
                (((deltaE1[i][j][0]>0)?pow(deltaE1[i][j][0], 2)/(exp(ep[i]/(kBoltzmann*temperature)) - 1):0) +
                 ((deltaE1[i][j][1]>0)?pow(deltaE1[i][j][1], 2)/(1 - exp(-ep[i]/(kBoltzmann*temperature))):0));
        }
    }
    if (alpha==0.0) alpha=1e-7;
    // alpha is in cm-1.
    } else {
    // wavelength range is between 413 & 3100 Angstroms - not well characterized. use a constant.
        alpha = 1.35e6;
    }
    return(alpha);
}


double Silicon::indexRefractionMCT (double lambda, double temperature, double x) {

    // follows Liu et al. J. Appl. Phys. 75 4176.
    double A = 13.173 - 9.852*x + 2.909*pow(x, 2) + 0.001*(300 - temperature);
    double B = 0.83 - 0.246*x - 0.0961*pow(x, 2) + 8e-4 *(300 - temperature);
    double C = 6.706 - 14.437*x + 8.531*pow(x, 2) + 7e-4*(300 - temperature);
    double D = 1.953e-4 - 0.00128*x + 1.853e-4*pow(x, 2);
    double n = A + B/(1 - pow(C/lambda, 2)) + D * pow(lambda, 2);
    if (n < 1) n=1.0;
    n=pow(n, 0.5);

    return(n);
}


double Silicon::indexRefractionSi (double lambda) {

    double wave= lambda/10000.;
    int index;
    find(const_cast<double *>(&nsilicon_wavelength[0]), nsilicon_wavelength.size(), wave, &index);
    double n = interpolate(const_cast<double *>(&nsilicon_real[0]),const_cast<double *>(&nsilicon_wavelength[0]), wave,index);
    return(n);

}

double Silicon::coatingIndexRefractionReal (double lambda) {

    double wave= lambda/10000.;
    int index;
    find(const_cast<double *>(&ncoating_wavelength[0]), ncoating_wavelength.size(), wave, &index);
    double n = interpolate(const_cast<double *>(&ncoating_real[0]),const_cast<double *>(&ncoating_wavelength[0]), wave,index);

    return(n);

}

double Silicon::coatingIndexRefractionImag (double lambda) {

    double wave= lambda/10000.;
    int index;
    find(const_cast<double *>(&ncoating_wavelength[0]), ncoating_wavelength.size(), wave, &index);
    double n = interpolate(const_cast<double *>(&ncoating_imag[0]),const_cast<double *>(&ncoating_wavelength[0]), wave,index);

    return(n);

}


double Silicon::dopeProfile (double z, double nbulk, double nf, double nb, double sf, double sb, double tsi) {

    if ((sb > 0.0) && (sf > 0.0)) {
        return(nbulk + nb*exp(-(tsi - z)/sb) + nf*exp(-(z)/sf));
    }
    if ((sb <= 0.0) && (sf > 0.0)) {
        return(nbulk + nf*exp(-(z)/sf));
    }
    if ((sb > 0.0) && (sf <= 0.0)) {
        return(nbulk + nb*exp(-(tsi - z)/sb));
    }
    return(nbulk);

}

double Silicon::muMCT(double x, double temperature) {

    // follows Rosbeck et al. (1982)
    double s = pow(0.2/x, 7.5);
    double r = pow(0.2/x, 0.6);
    double z = temperature;
    if (temperature < 50) {
        z = 1.18e5/(2600.0-pow(fabs(temperature-35.0), 2.07));
    }
    double mu_e = 9e8*s/pow(z ,2*r);
    return(mu_e); //mu_h

}


double Silicon::muSi (double efield, double temperature, int polarity, int esat, int vsat) {

    // Jacobini et al. (1977) equation (9)
    double vm, ec, beta, mu;
    if (polarity == -1) {
        vm = 1.53e9 * pow(temperature, -0.87); // cm/s
        ec = 1.01 * pow(temperature, 1.55); // V/cm
        beta = 2.57e-2 * pow(temperature, 0.66); // index
    } else {
        vm = 1.62e8 * pow(temperature, -0.52); // cm/s
        ec = 1.24 * pow(temperature, 1.68); // V/cm
        beta = 0.46 * pow(temperature, 0.17); // index
    }
    if (esat == 1) {
        mu = (vm/ec)/pow(1 + pow(fabs(efield)/ec, beta), 1/beta);
    } else {
        mu = vm/ec;
    }
    if ((vsat == 1) && (polarity == -1) && (efield != 0.0)) {
        double vs = 2.4e7/(1 + 0.8*exp(temperature/600.0));
        if (mu > vs/fabs(efield)) {
            mu = vs/fabs(efield);
        }
    }
    return(mu);

}


double Silicon::epsilonMCT(double x) {

    // follows Dornhaus et al. (1983)
    // high frequency approximation
    double epsilon_inf = 15.2 - 15.6*x + 8.2*pow(x,2);
    return(epsilon_inf);

}


void Silicon::sensorSetup (std::string sensorfilename, std::string devmaterial, double sensorTemp, double nbulk, double nf, double nb,
                     double sf, double sb, double tsi, double overdepBias,
                     std::string instrdir, std::string mdir,
                     long nampx, long nampy, double pixsize, long seedchip,
                     float *impurityX, float *impurityY, int impurityvariation,
                     double minwavelength, double maxwavelength, long activeBuffer,
                     int dopingflag, int esatflag, int vsatflag, int detectorMode, int siliconDebug, long siliconSubSteps, std::string obshistid) {

    numWavelength = 1024;
    numTemperature = 16;
    numDopant = 16;
    numRho = numTemperature;
    numThickness = SILICON_STEPS;

        double compositionHgCd = 0.5;
    double period = 100.0;
    double amplitude = 0.0;
    double ratio = 10.0;
    double surfaceCharge=0.0;
    double edgeBorder = 0.0;
    long edgeBorderPixels = 0;
    double depth = 0.0;
    double width = 3000.0;
    double height = 500.0;
    double xoverlap = 0.0;
    double yoverlap = 0.0;
    double pixelVarX = 0.0;
    double pixelVarY = 0.0;
    double epsilon;
    double periodVariation = 0.0;
    double periodDecay = 100000.0;
    int polarity = -1;
    std::vector<double> xinit, yinit, angle;
    spaceChargeShield = 0.0;
    spaceChargeSpreadX = 0.0;
    spaceChargeSpreadY = 0.0;
    polarity = -1;
    chargeStopCharge = 0.0;
    double chebDead[NCHEB];
    xTransferFlag=1;
    yTransferFlag=0;
    midlineFlag=1;
    inactiveX = -activeBuffer;
    inactiveY = -activeBuffer;

    for (int i =0; i< NCHEB;i++) chebDead[i]=0.0;

    std::string sss;
    sss = instrdir + "/" + sensorfilename;
    std::ifstream f(instrdir + "/" + sensorfilename);
    if (f.good()) {
        sss = instrdir + "/" + sensorfilename;
    } else {
        sss = instrdir + "/../standard/sensor/" + sensorfilename + ".txt";
    }
    std::ifstream inStream(sss.c_str());
    if (inStream) {
        readText pars(sss);
        for (size_t t(0); t < pars.getSize(); t++) {
            std::string line(pars[t]);
            std::istringstream iss(line);
            std::string keyName;
            iss >> keyName;
            if (keyName == "deadLayerCheb") {
                int chebNumber;
                double chebValue;
                iss >> chebNumber;
                iss >> chebValue;
                chebDead[chebNumber]=chebValue;
            }
            readText::get(line, "compositionHgCd", compositionHgCd);
            readText::get(line, "treeRingPeriod", period);
            readText::get(line, "treeRingVariation", periodVariation);
            readText::get(line, "treeRingDecay", periodDecay);
            readText::get(line, "treeRingAmplitude", amplitude);
            readText::get(line, "treeRingRatio", ratio);
            readText::get(line, "edgeSurfaceCharge", surfaceCharge);
            readText::get(line, "edgeBorder", edgeBorder);
            readText::get(line, "deadLayerDepth", depth);
            readText::get(line, "deadLayerWidth", width);
            readText::get(line, "deadLayerHeight", height);
            readText::get(line, "deadLayerXoverlap", xoverlap);
            readText::get(line, "deadLayerYoverlap", yoverlap);
            readText::get(line, "pixelVarX", pixelVarX);
            readText::get(line, "pixelVarY", pixelVarY);
            readText::get(line, "deadLayerXinit", xinit);
            readText::get(line, "deadLayerYinit", yinit);
            readText::get(line, "deadLayerAngle", angle);
            readText::get(line, "spaceChargeShield", spaceChargeShield);
            readText::get(line, "spaceChargeSpreadX", spaceChargeSpreadX);
            readText::get(line, "spaceChargeSpreadY", spaceChargeSpreadY);
            readText::get(line, "chargeStopCharge", chargeStopCharge);
            readText::get(line, "siliconType", polarity);
            readText::get(line, "xTransfer", xTransferFlag);
            readText::get(line, "yTransfer", yTransferFlag);
            readText::get(line, "inactiveX", inactiveX);
            readText::get(line, "inactiveY", inactiveY);
            readText::get(line, "midline", midlineFlag);
        }
    }
    spaceChargeSpreadX = (spaceChargeSpreadX/(pixsize*1e3))*(spaceChargeSpreadX/(pixsize*1e3));
    spaceChargeSpreadY = (spaceChargeSpreadY/(pixsize*1e3))*(spaceChargeSpreadY/(pixsize*1e3));
    period = period/(pixsize*1e3);
    periodDecay = periodDecay/(pixsize*1e3);
    periodVariation = periodVariation/(pixsize*1e3);
    edgeBorderPixels = static_cast<long>(edgeBorder/(pixsize*1e3));

    // index of refraction
    readText::readCol(mdir + "/optical/si.txt",nsilicon_wavelength,nsilicon_real,nsilicon_imag);
    readText::readCol(mdir + "/optical/hfo2.txt",ncoating_wavelength,ncoating_real,ncoating_imag);

    sf = sf/1e4f;
    sb = sb/1e4f;
    tsi = tsi/10.0;

    rho = static_cast<float*>(calloc(numRho, sizeof(float)));
    if (rho == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    meanFreePath = static_cast<float*>(calloc(numWavelength*numTemperature, sizeof(float)));
    if (meanFreePath == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    wavelengthGrid = static_cast<float*>(calloc(numWavelength, sizeof(float)));
    if (wavelengthGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    temperatureGrid = static_cast<float*>(calloc(numTemperature, sizeof(float)));
    if (temperatureGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    indexRefraction = static_cast<float*>(calloc(numWavelength, sizeof(float)));
    if (indexRefraction == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    indexRefractionCoatingReal = static_cast<double*>(calloc(numWavelength, sizeof(double)));
    if (indexRefractionCoatingReal == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    indexRefractionCoatingImag = static_cast<double*>(calloc(numWavelength, sizeof(double)));
    if (indexRefractionCoatingImag == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }


        if (devmaterial == "MCT") { // HgCdTe
        for (long i = 0; i < numTemperature; i++) {
            for (long j = 0; j < numWavelength; j++) {
                temperatureGrid[i] = sensorTemp + (static_cast<double>(i))/(static_cast<double>(numTemperature))*2.0 - 1.0;
                rho[i] = tsi*(static_cast<double>(i + 1))/(static_cast<double>(numRho))/5.0;
                wavelengthGrid[j] = minwavelength/1000.0 +
                    (static_cast<double>(j))/(static_cast<double>(numWavelength - 1))*(maxwavelength - minwavelength)/1000.0;
                meanFreePath[i*numWavelength + j] = 10.0/absorptionCoeffMCT(wavelengthGrid[j]*10000.0, temperatureGrid[i], compositionHgCd);
            }
        }
        for (long i = 0; i < numWavelength; i++) {
            indexRefraction[i] = indexRefractionMCT(wavelengthGrid[i], sensorTemp, compositionHgCd);
            indexRefractionCoatingReal[i] = coatingIndexRefractionReal(wavelengthGrid[i]*10000.0);
            indexRefractionCoatingImag[i] = coatingIndexRefractionImag(wavelengthGrid[i]*10000.0);
        }
    } else { // Silicon
        for (long i = 0; i < numTemperature; i++) {
            for (long j = 0; j < numWavelength; j++) {
                temperatureGrid[i] = sensorTemp + (static_cast<double>(i))/(static_cast<double>(numTemperature))*2.0 - 1.0;
                rho[i] = tsi*(static_cast<double>(i + 1))/(static_cast<double>(numRho))/5.0;
                wavelengthGrid[j] = minwavelength/1000.0 +
                    (static_cast<double>(j))/(static_cast<double>(numWavelength - 1))*(maxwavelength - minwavelength)/1000.0;
                meanFreePath[i*numWavelength + j] = 10.0/absorptionCoeffSi(wavelengthGrid[j]*10000.0, temperatureGrid[i]);
            }
        }
        for (long i = 0; i < numWavelength; i++) {
            indexRefraction[i] = indexRefractionSi(wavelengthGrid[i]*10000.0);
            indexRefractionCoatingReal[i] = coatingIndexRefractionReal(wavelengthGrid[i]*10000.0);
            indexRefractionCoatingImag[i] = coatingIndexRefractionImag(wavelengthGrid[i]*10000.0);
        }
    }
        if (detectorMode > 0) {

    dopantGrid = static_cast<float*>(calloc(numDopant, sizeof(float)));
    if (dopantGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    thicknessGrid = static_cast<float*>(calloc(numThickness, sizeof(float)));
    if (thicknessGrid == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigma = static_cast<float*>(calloc(numTemperature*numThickness*numDopant, sizeof(float)));
    if (sigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    hsigma = static_cast<float*>(calloc(numRho*numThickness*numDopant, sizeof(float)));
    if (hsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    isigma = static_cast<float*>(calloc(numRho*numThickness*numDopant, sizeof(float)));
    if (isigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    fsigma = static_cast<float*>(calloc(numThickness*numDopant, sizeof(float)));
    if (fsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gsigma = static_cast<float*>(calloc(numThickness*numDopant, sizeof(float)));
    if (gsigma == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigmaX = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (sigmaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    sigmaY = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (sigmaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gammaX = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (gammaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    gammaY = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (gammaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deltaX = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (deltaX == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deltaY = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (deltaY == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    nbulkmap = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (nbulkmap == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    periodmap = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (periodmap == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }
    deadLayer = static_cast<float*>(calloc((floor((nampx+2*activeBuffer-1)/DET_RESAMPLE)+1)*(floor((nampy+2*activeBuffer-1)/DET_RESAMPLE)+1), sizeof(float)));
    if (deadLayer == NULL) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }

    // tree ring map
    random.setSeed32Fixed(1000 + seedchip);
    *impurityX = random.uniformFixed()*(ratio*nampx) - (ratio - 1)/2*nampx + activeBuffer;
    *impurityY = random.uniformFixed()*(ratio*nampy) - (ratio - 1)/2*nampy + activeBuffer;
    for (long k = 0; k < 10; k++) {
        double phase = random.uniformFixed()*2.0*PI;
        double localVariation = periodVariation*random.normalFixed();
        for (long i = 0; i < nampx + 2*activeBuffer; i+=DET_RESAMPLE) {
            for (long j = 0; j < nampy + 2*activeBuffer; j+=DET_RESAMPLE) {
                double radius = sqrt(pow(i - *impurityX, 2) + pow(j - *impurityY, 2));
                double localPeriod = (period + localVariation)*(exp(-radius/periodDecay));
                int index = floor((nampx + 2*activeBuffer-1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
                if (k == 0) periodmap[index] = localPeriod/10.0;
                if (k != 0) periodmap[index] += localPeriod/10.0;
                if (k == 0) nbulkmap[index] = 1.0 + amplitude/sqrt(10.0)*sin(2*PI*radius/localPeriod + phase);
                if (k != 0) nbulkmap[index] += amplitude/sqrt(10.0)*sin(2*PI*radius/localPeriod + phase);
            }
        }
    }

    if (devmaterial == "MCT") {
        epsilon = epsilonMCT(compositionHgCd);
        fano = FANO_MCT;
        pairenergy = PAIRENERGY_MCT;
        density = DENSITY_HGTE * (1-compositionHgCd) +DENSITY_CDTE * compositionHgCd;
        excitation = EXCITATION_HGTE * (1-compositionHgCd) + EXCITATION_CDTE * compositionHgCd;
        atomicA = (0.5*M_HGTE) *(1-compositionHgCd) + (0.5*M_CDTE) * compositionHgCd;
        atomicZ = (80*0.5+52*0.5)*(1-compositionHgCd) + (48*0.5+52*0.5)*compositionHgCd;
    } else {
        epsilon = EPSILON_SI;
        fano = FANO_SI;
        pairenergy = PAIRENERGY_SI;
        density = DENSITY_SI;
        excitation = EXCITATION_SI;
        atomicA = M_SI;
        atomicZ = 14;
    }
    maxWavelengthElectronPair = H_PLANCK*C_LIGHT/(2.0*pairenergy*E_CHARGE)*1e6;

    // lateral deflection map
    double scale = fabs(surfaceCharge)/nbulk;
    if (scale == 0.0) scale = 1e-1*pixsize;
    for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
        for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
            //            int index = floor((nampx + 2*activeBuffer)/DET_RESAMPLE)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            int index = floor((nampx + 2*activeBuffer- 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            int sign = 1;

            if  (((i - activeBuffer + edgeBorderPixels) < 0) && (surfaceCharge < 0)) sign = -1; else sign=1;
            double factor = exp(-sign*(std::abs(i - activeBuffer + edgeBorderPixels)*pixsize*1e-1)/scale);
            if (((i - activeBuffer + edgeBorderPixels) < 0) && (surfaceCharge > 0)) factor = -factor;

            if (((nampx - 1 - (i - activeBuffer) + edgeBorderPixels) < 0) && (surfaceCharge < 0)) sign = -1; else sign=1;
            double factor2 = exp(-sign*(std::abs(nampx - 1 - (i - activeBuffer) + edgeBorderPixels)*pixsize*1e-1)/scale);
            if (((nampx - 1 - (i - activeBuffer) + edgeBorderPixels) < 0) && (surfaceCharge > 0)) factor2 = -factor2;

            deltaX[index] = surfaceCharge/EPSILON_0_CM/epsilon*factor*E_CHARGE-
                surfaceCharge/EPSILON_0_CM/epsilon*factor2*E_CHARGE;

            if (((j - activeBuffer + edgeBorderPixels) < 0) && (surfaceCharge < 0)) sign = -1; else sign=1;
            factor = exp(-sign*(std::abs(j - activeBuffer + edgeBorderPixels)*pixsize*1e-1)/scale);
            if (((j - activeBuffer + edgeBorderPixels) < 0) && (surfaceCharge > 0)) factor = -factor;

            if (((nampy - 1 - (j - activeBuffer) + edgeBorderPixels) < 0) && (surfaceCharge < 0)) sign=-1; else sign = 1;
            factor2 = exp(-sign*(std::abs(nampy - 1 - (j - activeBuffer) + edgeBorderPixels)*pixsize*1e-1)/scale);
            if (((nampy - 1 - (j - activeBuffer) + edgeBorderPixels) < 0) && (surfaceCharge > 0)) factor2 = -factor2;

            deltaY[index] = surfaceCharge/EPSILON_0_CM/epsilon*factor*E_CHARGE-
                surfaceCharge/EPSILON_0_CM/epsilon*factor2*E_CHARGE;
        }
    }

    // add tree rings
    for (long i = 1; i < (nampx - 1 + 2*activeBuffer); i+= DET_RESAMPLE) {
        for (long j = 1; j < (nampy - 1 + 2*activeBuffer); j+= DET_RESAMPLE) {
            int index = floor((nampx + 2*activeBuffer - 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            int indexxp = floor((nampx + 2*activeBuffer -1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor((i+1)/DET_RESAMPLE);
            int indexxm = floor((nampx + 2*activeBuffer - 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor((i-1)/DET_RESAMPLE);
            int indexyp = floor((nampx + 2*activeBuffer -1)/DET_RESAMPLE+1)*floor((j+1)/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            int indexym = floor((nampx + 2*activeBuffer -1)/DET_RESAMPLE+1)*floor((j-1)/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            if (amplitude != 0.0) {
                sigmaX[index] = (nbulkmap[indexxp] - nbulkmap[indexxm])*
                    nbulk*(2.0*tsi)*(periodmap[index]*pixsize*1e-1)/(2.0*PI)/(2.0*PI)*E_CHARGE/
                    epsilon/EPSILON_0_CM/2.0/(pixsize*1e-1)/2.0;
                sigmaY[index] = (nbulkmap[indexyp] - nbulkmap[indexym])*
                    nbulk*(2.0*tsi)*(periodmap[index]*pixsize*1e-1)/(2.0*PI)/(2.0*PI)*E_CHARGE/
                    epsilon/EPSILON_0_CM/2.0/(pixsize*1e-1)/2.0;
            } else {
                sigmaX[index]=0.0;
                sigmaY[index]=0.0;
            }
        }
    }

    if (impurityvariation == 0) {
        for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
            for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
                int index = floor((nampx + 2*activeBuffer - 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
                nbulkmap[index] = 1.0;
            }
        }
    }

    // int status;
    // long naxesa[2];
    // fitsfile *faptr;
    // status = 0;
    // fits_create_file(&faptr,"!defl.fits",&status);
    // naxesa[0] = nampx; naxesa[1] = nampy;
    // fits_create_img(faptr,FLOAT_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TFLOAT,1,nampx*nampy,deltaX,&status);
    // fits_close_file(faptr,&status);

    // nonuniform pixels
    for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
        double lithographicShift = random.normalFixed()*pixelVarX*pixsize;
        for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
            int index = floor((nampx + 2*activeBuffer - 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            gammaX[index] = lithographicShift;
        }
    }
    for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
        double lithographicShift = random.normalFixed()*pixelVarY*pixsize;
        for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
            int index = floor((nampx + 2*activeBuffer- 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            gammaY[index] = lithographicShift;
        }
    }

    // dead layer
    for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
        for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
            int index = floor((nampx + 2*activeBuffer - 1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
            if ((i < activeBuffer + inactiveX) ||
                (i >= nampx + activeBuffer - inactiveX)) {
                deadLayer[index] += tsi*10.0*1e6;
            }
            if ((j < activeBuffer + inactiveY) ||
                (j >= nampy + activeBuffer - inactiveY)) {
                deadLayer[index] += tsi*10.0*1e6;
            }
        }
    }
    for (size_t k = 0; k < xinit.size(); k++) {
        for (long i = 0; i < nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
            for (long j = 0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
                int index = floor((nampx + 2*activeBuffer-1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
                double ip = i*cos(angle[k]*DEGREE) - j*sin(angle[k]*DEGREE);
                double jp = i*sin(angle[k]*DEGREE) + j*cos(angle[k]*DEGREE);
                double di = fmod(ip*pixsize*1e3 - xinit[k], width - xoverlap);
                double dj = fmod(jp*pixsize*1e3 - yinit[k], height - yoverlap);
                while (di<0) di += width - xoverlap;
                while (dj<0) dj += height - yoverlap;
                if (di < width && dj < height) {
                    deadLayer[index] += depth;
                    if (di < xoverlap) deadLayer[index] += depth;
                    if (dj < yoverlap) deadLayer[index] += depth;
                    if (di < xoverlap && dj < yoverlap) deadLayer[index] += depth;
                }
            }
        }
    }

    
    // dead layer chebyshev
    double t[3];
    for (long n=0; n<NCHEB; n++) {
        for (long i=0; i< nampx + 2*activeBuffer; i+= DET_RESAMPLE) {
            for (long j=0; j < nampy + 2*activeBuffer; j+= DET_RESAMPLE) {
                double x = static_cast<double>(i-(nampx-1)/2-activeBuffer)/static_cast<double>((nampx-1)/2 +activeBuffer);
                double y = static_cast<double>(j-(nampy-1)/2-activeBuffer)/static_cast<double>((nampy-1)/2 + activeBuffer);
                chebyshevT2D(n,x,y,t);
                int index = floor((nampx + 2*activeBuffer -1)/DET_RESAMPLE+1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
                deadLayer[index] += t[0]*chebDead[n];
            }
        }
    }



    
    for (long i = 0; i < numDopant; i++) {
        if (amplitude != 0.0) {
            dopantGrid[i] = (static_cast<double>(i))/(static_cast<double>(numDopant))*2*amplitude*nbulk +
                (1.0 - amplitude)*nbulk;
        } else {
            dopantGrid[i] = (static_cast<double>(i))/(static_cast<double>(numDopant))*2*0.01*nbulk +
                (1.0 - 0.01)*nbulk;
        }
    }

    double z[SILICON_STEPS], E[SILICON_STEPS], v[SILICON_STEPS],
        tcol[SILICON_STEPS], tp[SILICON_STEPS], gp[SILICON_STEPS], hp[SILICON_STEPS], ip[SILICON_STEPS],
        tcoli[SILICON_STEPS], tpi[SILICON_STEPS], gpi[SILICON_STEPS], hpi[SILICON_STEPS], ipi[SILICON_STEPS];

    double minE, dz;
    int i, minEindex;
    long n = SILICON_STEPS;

    if (devmaterial == "MCT") { // MCT

        double factor = E_CHARGE/EPSILON_0_CM/epsilonMCT(compositionHgCd)/(4*PI)*spaceChargeShield;

        for (long l = 0; l < numDopant; l++) {
            for (long k = 0; k < numTemperature; k++) {
                dz = (-tsi)/static_cast<double>(n - 1);
                for (i = 0; i < n; i++) {
                    z[i] = tsi + i*dz;
                    thicknessGrid[i] = -i*dz*10.0;
                    if (i == 0) {
                        E[0] = overdepBias/tsi;
                        v[0] = 0;
                    } else {
                        // sample the doping profile to evaluate E[i]
                        int nSample = 100;
                        double dp = 0;
                        int j;
                        dp = 0;
                        for (j = 0; j < nSample; j++) {
                            dp += dopeProfile(tsi + dz*(i + (j)/(static_cast<double>(nSample))), dopantGrid[l], nf, nb, sf, sb, tsi);
                        }
                        dp /= nSample;
                        if (dopingflag == 0) dp = 0.0;
                        E[i] = E[i - 1] + dp*E_CHARGE/(EPSILON_0_CM*epsilonMCT(compositionHgCd))*dz;
                        v[i] = v[i - 1] - 0.5f*(E[i] + E[i - 1])*dz;
                    }
                }
                minEindex = 0;
                minE = 0;
                for (i = 0; i < n; i++) {
                    // find minimum value for E and count from that index
                    if (E[i] < minE) {
                        minEindex = i;
                        minE = E[i];
                    }
                }
                i = n - 1;
                tcoli[i] = (dz/(minE*muMCT(compositionHgCd, temperatureGrid[k])));
                tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                    (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                gpi[i] = dz/minE;
                hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                while (i--) {
                    if (i > minEindex) {
                        tcoli[i] = (dz/(minE*muMCT(compositionHgCd, temperatureGrid[k])));
                        tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                            (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                        gpi[i] = dz/minE;
                        hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                    } else {
                        if (E[i] + E[i + 1] > 0) {
                            tcoli[i] = 0.1;
                            tpi[i] = 0.1*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 0.1;
                            hpi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        } else {
                            tcoli[i] = 2*dz/(E[i]*muMCT(compositionHgCd, temperatureGrid[k]) +
                                            E[i + 1]*muMCT(compositionHgCd, temperatureGrid[k]));
                            tpi[i] = 2*dz/(E[i] + E[i + 1])*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 2*dz/(E[i] + E[i + 1]);
                            hpi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        }
                    }
                }
                i = n;
                while (i--) {
                    tcol[i] = 0.0;
                    tp[i] = 0.0;
                    gp[i] = 0.0;
                    hp[i] = 0.0;
                    ip[i] = 0.0;
                    for (long j = i; j < i + siliconSubSteps; j++) {
                        if (j < n) {
                            tcol[i] += tcoli[j];
                            tp[i] += tpi[j];
                            gp[i] += gpi[j];
                            hp[i] += hpi[j];
                            ip[i] += ipi[j];
                        }
                    }
                }
                i = n;
                while (i--) {
                    sigma[l*numThickness*numTemperature + k*numThickness + i] = sqrt(2*K_BOLTZMANN*temperatureGrid[k]
                                                                                     /E_CHARGE*muMCT(compositionHgCd, temperatureGrid[k]) * tcol[i])*1e4f/1e3;
                    fsigma[l*numThickness + i] = tp[i]*1e4f/1e3;
                    gsigma[l*numThickness + i] = gp[i]*1e4f/1e3;
                    hsigma[l*numThickness*numRho + k*numThickness + i] = hp[i]*1e4f/1e3;
                    isigma[l*numThickness*numRho + k*numThickness + i] = ip[i]*1e4f/1e3;
                }

            }
        }
    } else { // Silicon
        
        double factor = E_CHARGE/EPSILON_0_CM/EPSILON_SI/(4*PI)*spaceChargeShield;

        for (long l = 0; l < numDopant; l++) {
            for (long k = 0; k < numTemperature; k++) {
                dz = (-tsi)/static_cast<double>(n - 1);
                for (i = 0; i < n; i++) {
                    z[i] = tsi + i*dz;
                    thicknessGrid[i] = -i*dz*10.0;
                    if (i == 0) {
                        E[0] = overdepBias/tsi;
                        v[0] = 0;
                    } else {
                        // sample the doping profile to evaluate E[i]
                        int nSample = 100;
                        double dp = 0;
                        int j;
                        dp = 0;
                        for (j = 0; j < nSample; j++) {
                            dp += dopeProfile(tsi + dz*(i + (j)/(static_cast<double>(nSample))), dopantGrid[l], nf, nb, sf, sb, tsi);
                        }
                        dp /= nSample;
                        if (dopingflag == 0) dp = 0.0;
                        E[i] = E[i - 1] + dp*E_CHARGE/(EPSILON_0_CM*EPSILON_SI)*dz;
                        v[i] = v[i - 1] - 0.5f*(E[i] + E[i - 1])*dz;
                    }
                }
                minEindex = 0;
                minE = 0;
                for (i = 0; i < n; i++) {
                    // find minimum value for E and count from that index
                    if (E[i] < minE) {
                        minEindex = i;
                        minE = E[i];
                    }
                }
                i = n - 1;
                tcoli[i] = (dz/(minE*muSi(minE, temperatureGrid[k], polarity, esatflag, vsatflag)));
                tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                    (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                gpi[i] = dz/minE;
                hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                while (i--) {
                    if (i > minEindex) {
                        tcoli[i] = (dz/(minE*muSi(minE, temperatureGrid[k], polarity, esatflag, vsatflag)));
                        tpi[i] = dz/minE*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                            (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                        gpi[i] = dz/minE;
                        hpi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        ipi[i] = dz/minE*factor/(z[i]*z[i] + rho[k]*rho[k]);
                    } else {
                        if (E[i] + E[i + 1] > 0) {
                            tcoli[i] = 0.1;
                            tpi[i] = 0.1*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 0.1;
                            hpi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 0.1*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        } else {
                            tcoli[i] = 2*dz/(E[i]*muSi(E[i], temperatureGrid[k], polarity, esatflag, vsatflag) +
                                             E[i + 1]*muSi(E[i + 1], temperatureGrid[k], polarity, esatflag, vsatflag));
                            tpi[i] = 2*dz/(E[i] + E[i + 1])*(static_cast<float>(i))/(static_cast<float>(n - 1))*
                                (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0;
                            gpi[i] = 2*dz/(E[i] + E[i + 1]);
                            hpi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                            ipi[i] = 2*dz/(E[i] + E[i + 1])*factor/(z[i]*z[i] + rho[k]*rho[k]);
                        }
                    }
                }
                i = n;
                while (i--) {
                    tcol[i] = 0.0;
                    tp[i] = 0.0;
                    gp[i] = 0.0;
                    hp[i] = 0.0;
                    ip[i] = 0.0;
                    for (long j = i; j < i + siliconSubSteps; j++) {
                        if (j < n) {
                            tcol[i] += tcoli[j];
                            tp[i] += tpi[j];
                            gp[i] += gpi[j];
                            hp[i] += hpi[j];
                            ip[i] += ipi[j];
                        }
                    }
                }
                i = n;
                while (i--) {
                    sigma[l*numThickness*numTemperature + k*numThickness + i] = sqrt(2*K_BOLTZMANN*temperatureGrid[k]
                                                                                     /E_CHARGE*muSi(0, temperatureGrid[k], polarity, esatflag, vsatflag) * tcol[i])*1e4f/1e3;
                    fsigma[l*numThickness + i] = tp[i]*1e4f/1e3;
                    gsigma[l*numThickness + i] = gp[i]*1e4f/1e3;
                    if (isnan(tp[i]) || isnan(gp[i])) printf("A %ld %d %e %e ",l,i,tp[i],gp[i]);
                    hsigma[l*numThickness*numRho + k*numThickness + i] = hp[i]*1e4f/1e3;
                    isigma[l*numThickness*numRho + k*numThickness + i] = ip[i]*1e4f/1e3;
                }

            }
        }
    }
        



        if (siliconDebug==1) {
            std::ostringstream silout;
            silout << "silicon_structure_" << obshistid << "_0.txt";
            std::ofstream silout1(silout.str().c_str());
            for (int i = 0; i < SILICON_STEPS; i++) {
                silout1 << std::scientific << std::setprecision(6) << i << " " << E[i] << " " << 1.0 << " " <<
                    (static_cast<float>(i))/(static_cast<float>(n - 1))*
                    (1.0 - (static_cast<float>(i))/(static_cast<float>(n - 1)))*4.0 << " " << std::endl;
            }
            silout1.close();
            std::ostringstream silout3;
            silout3 << "silicon_structure_" << obshistid << "_1.txt";
            std::ofstream silout2(silout3.str().c_str());
            for (long i = nampx/2+activeBuffer-100; i < nampx/2+activeBuffer + 200; i++) {
                for (long j = nampy/2+activeBuffer-100; j < nampy/2+activeBuffer+200; j++) {
                    int index = floor((nampx + 2*activeBuffer-1)/DET_RESAMPLE + 1)*floor(j/DET_RESAMPLE) + floor(i/DET_RESAMPLE);
                    silout2 << std::scientific << std::setprecision(6) << i << " " << j << " " <<
                        deltaX[index] << " " <<
                        deltaY[index] << " " <<
                        sigmaX[index] << " " <<
                        sigmaY[index]  << std::endl;
                }
            }
            silout2.close();
        }

        }

}
