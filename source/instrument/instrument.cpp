///
/// @package phosim
/// @file instrument.cpp
/// @brief instrument  application.
///
/// @brief Created by
/// @author Nathan Todd (Purdue)
///
/// @brief Modified by
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author Glenn Sembroski (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "instrument.h"
#include <vector>
#include <array>



int main(void) {

    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Instrument Configuration" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;

    readText pars(std::cin);
    int parsNumberObsExtra = pars.getSize();
    std::string instrumentDir = "../data/lsst";
    std::string workDir = ".";
    std::string obsID = "NO-obsid";
    double vistime = 15.0;
    long obsseed = -1;
    int filter = 0;
    double altitude = 90.0 * PI/180.0;
    double zenith = PI/2.0 - altitude;
    double jittertime = 0.1;
    double pressure = 520.0;
    double temperature = 20.0;
    double seeing = 0.67;
    double domeseeing = 0.1;
    double azimuth = 0.0;
    std::vector<int> pertSurf;
    std::vector<int> pertSecondSurf;
    std::vector<long> pertPoint;
    std::vector<double> pertTolerance;
    std::vector<int> pertControl;
    std::vector<double> pertWfe;
    std::vector<double> pertAct;
    std::vector<double> pertAngTol;
    std::vector<double> pertInitStep;
    std::vector<int> pertZernikeStart;
    std::vector<int> pertActuatorStart;
    std::vector<int> ipolySurfaces;
    std::vector<int> ipolyPolyIndex;
    std::vector<double> ipolyY;
    std::vector<double> ipolyZ;
    int pertDebug=0;
    double surfaceAccuracy=1.0;
    int camConfig;
    Random random;
    int perturbationMode = 3;
    int telescopeMode = 1;
    int control;
    double temperatureChange = 0.0;
    int overridePhysics = 0;
    int pertFlag = 0;
    int diagnostic = 0;
    int nthread=8;
    int burnin=200;
    int maxiter=100000;
    std::string site("none");

    // Read obsExtra file and get relevant commands
    for (size_t t(0); static_cast<int>(t) < parsNumberObsExtra; t++) {
        std::string line(pars[t]);
        readText::get(line, "thread", nthread);
        readText::get(line, "obshistid", obsID);
        readText::get(line, "instrdir", instrumentDir);
        readText::get(line, "workdir",  workDir);
        readText::get(line, "obsseed", obsseed);
        readText::get(line, "filter", filter);
        readText::get(line, "diagnostic", diagnostic);
        readText::get(line, "temperature", temperature);
        readText::get(line, "vistime",  vistime);
        readText::get(line, "jittertime",  jittertime);
        readText::get(line, "camconfig",   camConfig);
        readText::get(line, "constrainseeing", seeing);
        readText::get(line, "site", site);
        readText::get(line, "domeseeing", domeseeing);
        readText::get(line, "tempvar", temperatureChange);
        readText::get(line, "control", control);
        readText::get(line, "burnin", burnin);
        readText::get(line, "maxiter", maxiter);
        if (readText::getKey(line, "azimuth", azimuth)) {
            azimuth *= PI/180.0; // Degrees
        }
        if (readText::getKey(line, "altitude", altitude)) {
            altitude *= PI/180.0; // Degrees
            zenith = PI/2.0 - altitude; // Degrees
        }
        if (readText::getKey(line, "zenith", zenith)) {
            zenith *= PI/180.0; // Degrees
        }
        readText::get(line, "temperature", temperature);
        readText::get(line, "pressure", pressure);
        readText::get(line, "perturbationmode", perturbationMode);
        readText::get(line, "telescopemode", telescopeMode);
        readText::get(line, "overridesurfacephysics", overridePhysics);
        readText::get(line, "surfacephysicsdebug", pertDebug);
        readText::get(line, "surfacephysicsaccuracy", surfaceAccuracy);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if (keyName == "distortion" ) {
            int surfaceIndex;
            iss >> surfaceIndex;
            int secondSurfaceIndex;
            iss >> secondSurfaceIndex;
            long points;
            iss >> points;
            double tolerance;
            iss >> tolerance;
            int cont;
            iss >> cont;
            double wfe;
            iss >> wfe;
            double act;
            iss >> act;
            double angTol;
            iss >> angTol;
            double initStep;
            iss >> initStep;
            pertSurf.push_back(surfaceIndex);
            pertSecondSurf.push_back(secondSurfaceIndex);
            pertPoint.push_back(points);
            pertTolerance.push_back(tolerance);
            pertControl.push_back(cont);
            pertWfe.push_back(wfe);
            pertAct.push_back(act);
            pertAngTol.push_back(angTol);
            pertInitStep.push_back(initStep);
            int zernikeStart;
            iss >> zernikeStart;
            pertZernikeStart.push_back(zernikeStart);
            int actuatorStart;
            iss >> actuatorStart;
            pertActuatorStart.push_back(actuatorStart);
        }

    }
    if (telescopeMode == 0) perturbationMode = 0;

    int spaceMode = 0;
    std::string sitestring;
    if (site == "none") {
        sitestring = instrumentDir + "/site.txt";
    } else {
        sitestring = "../data/standard/location/" + site + ".txt";
    }
    std::ifstream inStreamLocation(sitestring.c_str());
    if (inStreamLocation) {
        readText locationPars(sitestring);
        for (size_t t(0); t < locationPars.getSize(); t++) {
            std::string line(locationPars[t]);
            readText::get(line,"spacemode",spaceMode);

        }
    }

    // Read optics_0.txt file to get map of surface names to surface ID
    // This makes optics_0.txt an ISC requirement
    std::string opticsFile = instrumentDir + "/optics_0.txt";
    instrumentFiles.makeSurfaceMap(opticsFile);

    int nsurf = 0;
    readText opticsPars(opticsFile);
    for (size_t t(0); t < opticsPars.getSize(); t++) {
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName;
        std::string surfaceType;
        iss >> surfaceName;
        iss >> surfaceType;
        if (surfaceType != "none" && surfaceName != "optics"){
            nsurf++;
        }
    }

    // Initialize random number generator
    if (obsseed == -1) {
        random.setSeedFromTime();
    } else {
        random.setSeed32(obsseed);
    }
    random.unwind(10000);


    // Create and write tracking file
    std::string trackingFileName = workDir + "/tracking_" + obsID + ".pars";
    instrumentFiles.makeTrackingFile(trackingFileName, vistime, jittertime);

    // Get current optics file
    std::string opticsFileName = workDir + "/optics_" + obsID + ".pars";

    std::ofstream ofs(opticsFileName.c_str());


    std::vector<double> minT(nsurf);
    std::vector<double> maxT(nsurf);
    std::vector<double> minA(nsurf);
    std::vector<double> maxA(nsurf);
    std::vector<double> minP(nsurf);
    std::vector<double> maxP(nsurf);
    std::vector<int> Tflag(nsurf, 0);
    std::vector<int> Aflag(nsurf, 0);
    std::vector<int> Pflag(nsurf, 0);
    std::vector<double> norm(nsurf);
    std::vector<std::array<double, NZERN>> zernike(nsurf);
    std::vector<std::array<double, 6>> body(nsurf);
    std::vector<std::array<double, 6>> bodyY(nsurf);

    for (long i=0; i < nsurf; i++) {
        for (long j = 0; j < 6; j++) {
            body[i][j] = 0.0;
            bodyY[i][j] = 0.0;
        }
    }
    for (long i=0; i < nsurf; i++) {
        for (long j = 0; j < NZERN; j++) {
            zernike[i][j] = 0.0;
        }
    }

    double actuatorValue[10000];
    double actuatorValueNormal[10000];
    double actuatorValueUniform[10000];
    double randomValueNormal[10000];
    double randomValueUniform[10000];
    double randomatmValueNormal[10000];
    double randomatmValueUniform[10000];
    double fabricationValueNormal[10000];
    double fabricationValueUniform[10000];

    for (long i = 0; i < nsurf*(NZERN + 6) - 1; i++) {
        actuatorValue[i]= 0.0;
        actuatorValueNormal[i] = random.normal();
        actuatorValueUniform[i] = (2.0*random.uniform() - 1.0);
        randomValueNormal[i] = random.normal();
        randomValueUniform[i] = (2.0*random.uniform() - 1.0);
        randomatmValueNormal[i] = random.normal();
        randomatmValueUniform[i] = (2.0*random.uniform() - 1.0);
        fabricationValueNormal[i] = random.normalFixed();
        fabricationValueUniform[i] = (2.0*random.uniformFixed() - 1.0);
    }
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        int counter = 0;
        int a = 0;
        double b;
        if (keyName == "move") {
            std::string value;
            while (iss >> value) {
                if (counter % 2 == 0) {
                    a = atoi(value.c_str());
                    std::string sss;
                    sss = instrumentDir + "/perturbation.txt";
                    std::ifstream inStream(sss.c_str());
                    if (inStream) {
                        readText perturbationPars4(sss.c_str());
                        for (size_t t(0); t < perturbationPars4.getSize(); t++) {
                            std::string line(perturbationPars4[t]);
                            std::istringstream iss(line);
                            std::string label;
                            std::string surfaceString;
                            std::string perturbationDof;
                            std::string perturbationValue;
                            std::string perturbationType;
                            std::string surface;
                            int actuatorNumber;
                            iss >> label >> surfaceString >> perturbationDof >> perturbationValue >> perturbationType >> actuatorNumber;
                            if (label == value) {
                                a = actuatorNumber;
                            }
                        }
                    }
                } else {
                    b = atof(value.c_str());
                    actuatorValue[a] = b;
                }
                counter++;
            }
        }
    }

    // Read perturbation file

    for (int surface = 0; surface < nsurf; surface++) {
        minT[surface] = 1e30;
        maxT[surface] = -1e30;
        minA[surface] = 1e30;
        maxA[surface] = -1e30;
        minP[surface] = 1e30;
        maxP[surface] = -1e30;
        norm[surface] = 0.0;
    }

    // Run through initial time
    std::string sss;
    sss = instrumentDir + "/perturbation.txt";
    std::ifstream inStream(sss.c_str());
    if (inStream) {
    readText perturbationPars(sss.c_str());
    for (size_t t(0); t < perturbationPars.getSize(); t++) {
        std::string line(perturbationPars[t]);
        std::istringstream iss(line);
        std::string label;
        std::string surfaceString;
        std::string perturbationDof;
        std::string perturbationValue;
        std::string perturbationType;
        std::string surface;
        iss >> label >> surfaceString >> perturbationDof >> perturbationValue >> perturbationType;
        std::istringstream ist(surfaceString);
        while (std::getline(ist, surface, '|')) {
            if (perturbationType == "environment") {
                double temperatureT, altitudeT, pressureT;
                int s = atoi(surface.c_str());
                iss >> temperatureT;
                iss >> altitudeT;
                iss >> pressureT;
                if (temperatureT < minT[s]) minT[s] = temperatureT;
                if (temperatureT > maxT[s]) maxT[s] = temperatureT;
                if (pressureT < minP[s]) minP[s] = pressureT;
                if (pressureT > maxP[s]) maxP[s] = pressureT;
                if (cos(altitudeT) < minA[s]) minA[s] = cos(altitudeT);
                if (cos(altitudeT) > maxA[s]) maxA[s] = cos(altitudeT);
            }
        }
    }

    for (int surface = 0; surface < nsurf; surface++) {
        if (minA[surface] < maxA[surface]) Aflag[surface] = 1;
        if (minT[surface] < maxT[surface]) Tflag[surface] = 1;
        if (minP[surface] < maxP[surface]) Pflag[surface] = 1;
    }

    // Run through second time
    readText perturbationPars2(sss.c_str());
    for (size_t t(0); t < perturbationPars2.getSize(); t++) {
        std::string line(perturbationPars2[t]);
        std::istringstream iss(line);
        std::string label;
        std::string surfaceString;
        std::string perturbationType;
        std::string surface;
        std::string perturbationDof;
        std::string perturbationValue;
        iss >> label >> surfaceString >> perturbationDof >> perturbationValue >> perturbationType;
        std::istringstream ist(surfaceString);
        while (getline(ist, surface, '|')) {
            if (perturbationType == "environment") {
                double temperatureT, altitudeT, pressureT;
                int s = atoi(surface.c_str());
                iss >> temperatureT;
                iss >> altitudeT;
                iss >> pressureT;
                double distance = 0.0;
                if (Tflag[s] == 1) distance += pow((temperature - temperatureT)/
                                                   (maxT[s] - minT[s]), 2.0);
                if (Aflag[s] == 1) distance += pow((cos(altitude) - cos(altitudeT))/
                                                   (maxA[s] - minA[s]), 2.0);
                if (Pflag[s] == 1) distance += pow((pressure - pressureT)/
                                                   (maxP[s] - minP[s]), 2.0);
                distance = sqrt(distance);
                norm[s] += sqrt(static_cast<double>(Tflag[s]*Tflag[s] +
                                Aflag[s]*Aflag[s] +
                                Pflag[s]*Pflag[s])) - distance;
            }
        }
    }


    // Now the third time
    readText perturbationPars3(sss.c_str());
    for (size_t t(0); t < perturbationPars3.getSize(); t++) {
        std::string line(perturbationPars3[t]);
        std::istringstream iss(line);
        std::string label;
        std::string surfaceString;
        std::string perturbationType;
        std::string surface;
        std::string perturbationDof;
        std::string perturbationValue;
        int actuatorNumber;
        double actuatorScaling, actuatorScalingConst, actuatorScalingCos, actuatorScalingSin;
        double actuatorScalingAngle, actuatorScalingGauss, actuatorScalingUnif, actuatorComponent, actuatorComponentScale;
        iss >> label >> surfaceString >> perturbationDof >> perturbationValue >> perturbationType;
        if (perturbationType == "actuator" || perturbationType == "randomatmscl" ||
            perturbationType == "random" || perturbationType == "fabrication") {
                    iss >> actuatorNumber;
                    iss >> actuatorScaling;
                    iss >> actuatorScalingConst;
                    iss >> actuatorScalingAngle;
                    iss >> actuatorScalingCos;
                    iss >> actuatorScalingSin;
                    iss >> actuatorScalingGauss;
                    iss >> actuatorScalingUnif;
                    iss >> actuatorComponent;
                    iss >> actuatorComponentScale;
        }
        if (perturbationType == "physics" && overridePhysics==0) {
            int secondSurfaceIndex=-1;
            int surfaceIndex=0;
            int surfaceCounter=0;
            std::istringstream ist(surfaceString);
            while (getline(ist, surface, '|')) {
                if (surfaceCounter==0) surfaceIndex = atoi(surface.c_str());
                if (surfaceCounter==1) secondSurfaceIndex = atoi(surface.c_str());
                surfaceCounter++;
            }
            long points = atoi(perturbationValue.c_str());;
            double tolerance;
            iss >> tolerance;
            int cont;
            iss >> cont;
            double wfe;
            iss >> wfe;
            double act;
            iss >> act;
            double angTol;
            iss >> angTol;
            double initStep;
            iss >> initStep;
            pertSurf.push_back(surfaceIndex);
            pertSecondSurf.push_back(secondSurfaceIndex);
            pertPoint.push_back(points);
            pertTolerance.push_back(tolerance);
            pertControl.push_back(cont);
            pertWfe.push_back(wfe);
            pertAct.push_back(act);
            pertAngTol.push_back(angTol);
            pertInitStep.push_back(initStep);
            int zernikeStart;
            iss >> zernikeStart;
            pertZernikeStart.push_back(zernikeStart);
            int actuatorStart;
            iss >> actuatorStart;
            pertActuatorStart.push_back(actuatorStart);
            pertFlag=1;
        }
        int v = atoi(perturbationValue.c_str());
        double value = 0.0;
        double valueX = 0.0;
        double valueY = 0.0;
        std::istringstream ist(surfaceString);
        int presurf = -1;
        if ((perturbationType == "actuator" && (perturbationMode==1 || perturbationMode==3)) ||
            ((perturbationType == "random" || perturbationType == "randomatmscl" ||
              perturbationType == "fabrication" || perturbationType == "environment") &&
             (perturbationMode == 2 || perturbationMode == 3))) {
        while (getline(ist, surface, '|')) {
            int s = atoi(surface.c_str());
            if (presurf != -1) {
                if (perturbationDof == "map" || perturbationDof ==  "zernike"  || perturbationDof == "zlist") {
                    ofs << "surfacelink " << presurf << " " << s << std::endl;
                }
            }
            presurf = s;
            if (perturbationDof == "map") {
                if (perturbationType == "actuator") {
                    double vvv;
                    vvv = actuatorValue[actuatorNumber] + actuatorScalingGauss*(actuatorValueNormal[actuatorNumber]) +
                        actuatorScalingUnif*(actuatorValueUniform[actuatorNumber]);
                    value = actuatorScaling*vvv + actuatorScalingConst +
                        actuatorScalingCos*cos(actuatorScalingAngle*vvv) +
                        actuatorScalingSin*sin(actuatorScalingAngle*vvv);
                    valueX = value*cos(actuatorComponent + actuatorComponentScale*vvv);
                    valueY = value*sin(actuatorComponent + actuatorComponentScale*vvv);
                    if (value != 0.0) {
                        ofs << "surfacemap " << surface << " " << instrumentDir << "/" << perturbationValue << " " << value << std::endl;
                    }
                }
                if (perturbationType == "environment") {
                    double temperatureT, altitudeT, pressureT;
                    iss >> temperatureT;
                    iss >> altitudeT;
                    iss >> pressureT;
                    double distance = 0.0;
                    if (Tflag[s] == 1) distance += pow((temperature - temperatureT)/
                                                       (maxT[s] - minT[s]), 2.0);
                    if (Aflag[s] == 1) distance += pow((cos(altitude) - cos(altitudeT))/
                                                       (maxA[s] - minA[s]), 2.0);
                    if (Pflag[s] == 1) distance += pow((pressure - pressureT)/
                                                       (maxP[s] - minP[s]), 2.0);
                    distance = sqrt(distance);
                    value = sqrt(static_cast<double>(Tflag[s]*Tflag[s] +
                                        Aflag[s]*Aflag[s] +
                                        Pflag[s]*Pflag[s])) - distance;
                    if (norm[s] != 0) {
                        value /= norm[s];
                    } else {
                        value = 1.0;
                    }
                    if (value != 0.0) {
                        ofs << "surfacemap " << surface << " " << instrumentDir << "/" << perturbationValue << " " << value << std::endl;
                    }
                }
            }

            if (perturbationDof == "zernike" || perturbationDof == "body" || perturbationDof == "zlist") {
                if (perturbationType == "actuator") {

                    double vvv;
                    vvv = actuatorValue[actuatorNumber] + actuatorScalingGauss*actuatorValueNormal[actuatorNumber] +
                        actuatorScalingUnif*actuatorValueUniform[actuatorNumber];
                    value = actuatorScaling*vvv + actuatorScalingConst +
                        actuatorScalingCos*cos(actuatorScalingAngle*vvv) +
                        actuatorScalingSin*sin(actuatorScalingAngle*vvv);
                    valueX = value*cos(actuatorComponent + actuatorComponentScale*vvv);
                    valueY = value*sin(actuatorComponent + actuatorComponentScale*vvv);
                }
                if (perturbationType == "random") {
                    double vvv;
                    vvv = actuatorScalingGauss*randomValueNormal[actuatorNumber] +
                        actuatorScalingUnif*randomValueUniform[actuatorNumber];
                    value = actuatorScaling*vvv + actuatorScalingConst +
                        actuatorScalingCos*cos(actuatorScalingAngle*vvv) +
                        actuatorScalingSin*sin(actuatorScalingAngle*vvv);
                    valueX = value*cos(actuatorComponent + actuatorComponentScale*vvv);
                    valueY = value*sin(actuatorComponent + actuatorComponentScale*vvv);

                }
                if (perturbationType == "randomatmscl") {
                    double vvv;
                    vvv = actuatorScalingGauss*randomatmValueNormal[actuatorNumber] +
                        actuatorScalingUnif*randomatmValueUniform[actuatorNumber];
                    value = actuatorScaling*vvv + actuatorScalingConst +
                        actuatorScalingCos*cos(actuatorScalingAngle*vvv) +
                        actuatorScalingSin*sin(actuatorScalingAngle*vvv);
                    value = value * sqrt(seeing*seeing + domeseeing*domeseeing)/0.67*pow(1/cos(zenith), 3.0/5.0);
                    valueX = value*cos(actuatorComponent + actuatorComponentScale*vvv);
                    valueY = value*sin(actuatorComponent + actuatorComponentScale*vvv);
                }
                if (perturbationType == "fabrication") {
                    double vvv;
                    vvv = actuatorScalingGauss*fabricationValueNormal[actuatorNumber] +
                        actuatorScalingUnif*fabricationValueUniform[actuatorNumber];
                    value = actuatorScaling*vvv + actuatorScalingConst +
                        actuatorScalingCos*cos(actuatorScalingAngle*vvv) +
                        actuatorScalingSin*sin(actuatorScalingAngle*vvv);
                    valueX = value*cos(actuatorComponent + actuatorComponentScale*vvv);
                    valueY = value*sin(actuatorComponent + actuatorComponentScale*vvv);
                }

                if (perturbationDof == "zernike") {
                    zernike[s][v] += value;
                }
                if (perturbationDof == "body") {
                    if (perturbationType == "ipoly") // For the special case of using the ipoly function
                    {
                        std::istringstream ist(surfaceString); // Keeping column 2 in the pertubration file as the surface (X value in ipoly command)
                        int surfaceIndex = -1;
                        if (getline(ist, surface, '|')) {
                            surfaceIndex = atoi(surface.c_str());
                        }

                        int polyIndex;
                        double yParam, zParam; //only using the column 6 (Y value) and column 7 (Z value)
                        iss >> polyIndex >> yParam >> zParam;

                        // Essentially skips next 8 columns, expected to be 0
                        double dummy;
                        for (int i = 0; i < 8; ++i) {
                            iss >> dummy;
                        }

                        // Store ipoly data and adds it to the end of the ipoly vectors
                        ipolySurfaces.push_back(surfaceIndex); 
                        ipolyPolyIndex.push_back(polyIndex);
                        ipolyY.push_back(yParam);
                        ipolyZ.push_back(zParam);

                        // Skip default body logic
                        continue;
                    }

                    if (v >= 3) body[s][v] += value;
                    if (v == 2) body[s][v] += value*value;
                    if (v < 2) {
                        body[s][v] += valueX;
                        bodyY[s][v] += valueY;
                    }
                }
                if (perturbationDof == "zlist") {
                    if (value != 0.0) {
                        readText zlist(instrumentDir + "/" + perturbationValue);
                        int vv = 0;
                        for (size_t t(0); t < zlist.getSize(); t++) {
                            std::string line(zlist[t]);
                            std::istringstream iss(line);
                            std::string zscale;
                            while (iss >> zscale) {
                                double zv = atof(zscale.c_str());
                                if (vv < NZERN) zernike[s][vv] += value*zv;
                                vv++;
                            }
                        }
                    }
                }
            }
        }
        }
    }

    for (long i=0; i < nsurf; i++) {
        for (long j = 0; j < NZERN; j++) {
            if (zernike[i][j] != 0.0) {
                ofs << "izernike " << i << " " << j << " " << std::scientific << std::setprecision(6) << zernike[i][j] << std::endl;
            }
        }
    }

    for (long i=0; i < nsurf; i++) {
        for (long j = 0; j < 6; j++) {
            if (body[i][j] != 0.0) {
                if (j < 2) body[i][j] = atan2(bodyY[i][j], body[i][j]);
                if (j == 2) body[i][j] = sqrt(body[i][j]);
                ofs << "body " << i << " " << j << " " << std::scientific << std::setprecision(6) << body[i][j] << std::endl;
            }
        }
    }
    }

    ofs << "trackingfile tracking_" << obsID << ".pars" << std::endl;
    ofs.close();
    // }

    
    std::string focalPlaneLayoutFileName = instrumentDir +
        "/focalplanelayout.txt";
    readText focalPlaneLayoutPars(focalPlaneLayoutFileName);

    std::string segmentationFileName = instrumentDir + "/segmentation.txt";
    readText segmentationPars(segmentationFileName);

    std::string readoutString = workDir + "/readout_" + obsID +"_";
    instrumentFiles.readoutPars(focalPlaneLayoutPars, segmentationPars,
                                readoutString, camConfig);



    // Set body, chipangle, izernike, qevariation
    std::string outChipString = workDir +"/chip_" + obsID + "_";
    instrumentFiles.focalPlanePars(focalPlaneLayoutPars, outChipString,
                                   camConfig, perturbationMode,pertSurf,
                                   pertSecondSurf,pertFlag,instrumentDir);


    // distortion
    if (overridePhysics==0) {
        long diagCall=0;
        double gravity;
        if (spaceMode <= 0) gravity=GRAV_ACCEL_EARTH; else gravity=0.0;
        
    for (size_t t(0); t < pertSurf.size(); t++){
        int flag = 0;
        int seedstep = 0;
        pertPoint[t]=round(pertPoint[t]*pow(surfaceAccuracy,0.33));
        pertTolerance[t]=pertTolerance[t]/surfaceAccuracy;
        pertAngTol[t]=pertAngTol[t]/surfaceAccuracy;
        do {
            //            flag = distortion(pertSurf[t], pertSecondSurf[t], pertPoint[t], pertTolerance[t], control*pertControl[t], pertWfe[t]*pow(seeing,5./3.0), pertAct[t], obsseed+seedstep, obsID, filter, instrumentDir,
            //                 zenith, temperature, azimuth, temperatureChange, pertAngTol[t], pertInitStep[t], pertZernikeStart[t], pertActuatorStart[t],actuatorValue,pertDebug,gravity,&diagCall,nthread,burnin,maxiter);
            flag = distortion(pertSurf[t], pertSecondSurf[t], pertPoint[t], pertTolerance[t], control*pertControl[t], pertWfe[t]*pow(seeing,5./6.),pertAct[t], obsseed+seedstep, obsID, filter, instrumentDir,
                              zenith, temperature, azimuth, temperatureChange, pertAngTol[t], pertInitStep[t], pertZernikeStart[t], pertActuatorStart[t],actuatorValue,pertDebug,gravity,&diagCall,nthread,burnin,maxiter);
            seedstep++;
        } while ((flag==1) && (seedstep<10));
        if (seedstep >= 10) {
            printf("Error in distortion\n");
            exit(1);
        }
    }
    if (diagnostic) printf("$$%40s %22e %22.6lf\n","elastic",(double)diagCall,323.0);
    }
    return 0;
}
