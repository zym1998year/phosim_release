///
/// @package phosim
/// @file main.cpp
/// @brief main for atmosphere program
///
/// @brief Created by
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <time.h>
#include "atmosphere/atmosphere.h"
#include "atmosphere/operator.cpp"
#include "raytrace/parameters.h"

int main(void) {

// Default values.
    std::string obshistid = "9999";
    double outerx = 500000.0;
    double pix = 256.0*ATM_RESAMPLE; // pixel size in cm
    int telconfig = 0;
    float constrainseeing = -1;

    double seeing;
    double constrainairglow = -1.0;
    char tempstring[4096];
    FILE *fptr, *faptr;
    char outputfilename[4096];
    double zenith, altitude;
    long monthnum;
    long seed = 0;
    double tai = 0.0;

    int atmospheremode = 1;
    int opticsonlymode = 0;
    int turbulencemode = 1;
    int opacitymode = 1;
    long atmosphericDispersion = 1;
    int skipAtmosphere = 0;
    int spaceMode = -1;

    int haveMonth = 0;
    int haveTai = 0;
    int haveSunAlt = 0;
    int haveSunAz = 0;
    double sunalt = -90.0;
    double sunaz = 0.0;
    int haveSunRa = 0;
    int haveSunDec = 0;
    double sunra = 0.0;
    double sundec = 0.0;

    int haveMoonRa = 0;
    double moonra = 0.0;
    int haveMoonDec = 0;
    double moondec = 0.0;
    int haveMoonAlt = 0;
    double moonalt = 0.0;
    int haveMoonAz = 0;
    double moonaz = 0.0;
    int haveMoonPhase = 0;
    double moonphase = 0.0;

    int havePointingRA = 0;
    double pointingra = 0.0;
    int havePointingDec = 0;
    double pointingdec = 0.0;
    int haveAltitude = 0;
    int haveAzimuth = 0;
    double azimuth = 0.0;
    int haveDistMoon = 0;
    double dist2moon = 0.0;

    int haveRotTelPos = 0;
    int haveRotSkyPos = 0;
    double rotTelPos = 0.0;
    double rotSkyPos = 0.0;

    int haveFilter = 0;
    int filter = 0;
    int control = 0;
    int threadFlag = 1;

    int haveTemperature = 0;
    double temperature = 20.0;
    double temperatureMaxCL[4];
    double temperatureMinCL[4];
    double temperatureCL[4];
    double predictedRainfallCL[4];
    int haveTemperatureVariation = 0;
    double temperatureVariation = 0.0;
    int haveHumidity = 0;
    double humidity = 0.0;
    int havePressure = 0;
    double pressure = 520.0;

    Atmosphere atmosphere;
    std::string mjdstring;
    int haveMjdString = 0;


    atmosphere.numlevel = 8;
    atmosphere.constrainclouddepth = -1.0;
    atmosphere.constraincloudcover = -1.0;
    static int maxlevel = 100;
    std::vector<int> cloudscreen(maxlevel, 0);
    std::vector<float> outerscale(maxlevel, -1);

    // Set some default values.
    atmosphere.datadir = "../data";
    atmosphere.instrdir = "../data/lsst";
    std::string workdir(".");
    std::string site("none");

    // Read parameters from stdin.
    readText pars(std::cin);
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        readText::get(line, "workdir", workdir);
        readText::get(line, "instrdir", atmosphere.instrdir);
        readText::get(line, "datadir", atmosphere.datadir);
        readText::get(line, "numlevel", atmosphere.numlevel);
        readText::get(line, "telconfig", telconfig);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "outerx", outerx);
        readText::get(line, "cloudscreen", cloudscreen);
        readText::get(line, "outerscale", outerscale);
        readText::get(line, "obsseed", seed);
        readText::get(line, "site",site);
        if (readText::getKey(line, "monthnum", monthnum)) haveMonth=1;
        if (readText::getKey(line, "solaralt", sunalt)) haveSunAlt=1;
        if (readText::getKey(line, "solaraz", sunaz)) haveSunAz=1;
        if (readText::getKey(line, "solarra", sunra)) haveSunRa=1;
        if (readText::getKey(line, "solardec", sundec)) haveSunDec=1;
        if (readText::getKey(line, "moonra", moonra)) haveMoonRa=1;
        if (readText::getKey(line, "moondec", moondec)) haveMoonDec=1;
        if (readText::getKey(line, "moonalt", moonalt)) haveMoonAlt=1;
        if (readText::getKey(line, "moonaz", moonaz)) haveMoonAz=1;
        if (readText::getKey(line, "phaseang", moonphase)) haveMoonPhase=1;
        if (readText::getKey(line, "tai", tai)) haveTai=1;
        if (readText::getKey(line, "pointingra", pointingra)) havePointingRA=1;
        if (readText::getKey(line, "pointingdec", pointingdec)) havePointingDec=1;
        if (readText::getKey(line, "azimuth", azimuth)) haveAzimuth=1;
        if (readText::getKey(line, "moondist", dist2moon)) haveDistMoon=1;
        if (readText::getKey(line, "tai", mjdstring)) haveMjdString=1;
        if (readText::getKey(line, "rotationangle", rotSkyPos)) haveRotSkyPos=1;
        if (readText::getKey(line, "spiderangle", rotTelPos)) haveRotTelPos=1;
        if (readText::getKey(line, "filter", filter)) haveFilter=1;
        if (readText::getKey(line, "temperature", temperature)) haveTemperature=1;
        if (readText::getKey(line, "humidity", humidity)) haveHumidity=1;
        if (readText::getKey(line, "tempvar", temperatureVariation)) haveTemperatureVariation=1;
        if (readText::getKey(line, "pressure", pressure)) havePressure=1;
        readText::get(line, "atmospheremode", atmospheremode);
        readText::get(line, "opticsonlymode", opticsonlymode);
        readText::get(line, "atmosphericdispersion", atmosphericDispersion);
        readText::get(line, "spacemode", spaceMode);
        readText::get(line, "control", control);
        readText::get(line, "thread", threadFlag);

        // check to see if atmosphere is off
        if (keyName == "cleareverything" || opticsonlymode == 1 || atmospheremode == 0) {
             skipAtmosphere = 1;
        }
        if (keyName == "clearturbulence") turbulencemode = 0;
        if (keyName == "clearopacity") opacitymode = 0;
        if (turbulencemode == 0 && atmosphericDispersion == 0.0 && opacitymode == 0) {
            skipAtmosphere = 1;
        }
        if (control == 1) skipAtmosphere=0;

        if (readText::getKey(line, "constrainseeing", seeing)) constrainseeing = seeing/2.35482;
        if (readText::getKey(line, "altitude", altitude)) {
            haveAltitude = 1;
            zenith = 90 - altitude;
        }
        readText::get(line, "zenith", zenith);
        readText::getKey(line, "constrainclouddepth", atmosphere.constrainclouddepth);
        readText::getKey(line, "constraincloudcover", atmosphere.constraincloudcover);
        readText::getKey(line, "constrainairglow", constrainairglow);

}

    atmosphere.random.setSeed32(seed);
    atmosphere.random.unwind(10000);
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl
    << "Operator" << std::endl
    << "----------------------------------------------------------------------------------------------------" << std::endl;

     // Get site location
    double latitude = 0.0, longitude = 0.0;
    atmosphere.groundlevel = 0.0;
    if (site == "none") {
        atmosphere.sitestring = atmosphere.instrdir + "/site.txt";
    } else {
        atmosphere.sitestring = "../data/standard/location/" + site + ".txt";
    }
    std::ifstream newStream(atmosphere.sitestring.c_str());
    if (newStream) {
        readText locationPars(atmosphere.sitestring);
        for (size_t t(0); t < locationPars.getSize(); t++) {
            std::string line(locationPars[t]);
            readText::getKey(line, "latitude", latitude);
            readText::getKey(line, "longitude", longitude);
            readText::get(line, "groundlevel", atmosphere.groundlevel);
            if (spaceMode == -1) readText::getKey(line, "spacemode", spaceMode);
        }
    }
    atmosphere.latitude = latitude;
    atmosphere.longitude = longitude;
    if (spaceMode == -1) spaceMode = 0;

    if (spaceMode > 0) {
        latitude = 0;
        longitude = 0;
        skipAtmosphere = 1;
        // say spaceMode 1 = low eath orbit
        // 2 = LagrangePt. L1
        // 3=L2, (JWST) ...etc
    }

    double rotationConstraint = 0.0;
    std:: string ss2;
    ss2 = atmosphere.instrdir + "/structure.txt";
    std::ifstream newStream2(ss2.c_str());
    if (newStream2) {
        readText trackingPars(atmosphere.instrdir + "/structure.txt");
        for (size_t t(0); t < trackingPars.getSize(); t++) {
            std::string line(trackingPars[t]);
            readText::getKey(line, "rotationconstraint", rotationConstraint);
        }
    }

    time_t now = time(0);
    struct tm tstruct;
    tstruct = *gmtime(&now);
    long curryr, currmo, currdy, currhr, currmn, currsc;
    curryr = tstruct.tm_year+1900;
    currmo = tstruct.tm_mon+1;
    currdy = tstruct.tm_mday;
    currhr = tstruct.tm_hour;
    currmn = tstruct.tm_min;
    currsc = tstruct.tm_sec;
    if ((currmo == 1) || (currmo ==2)) {
        curryr = curryr -1;
        currmo = currmo + 12;
    }
    long curra, currb;
    curra = floor ( curryr /100);
    currb= 2 - curra + floor(curra/4);
    double currjd = floor(365.25*(curryr + 4716)) + floor(30.6001*(currmo + 1)) + currdy + currb - 1524.5;
    double currmjd = currjd - 2400000.5;
    currmjd += currhr/24.0 + currmn/24.0/60.0 + currsc/24.0/60.0/60.0;

    // Get time
    int countTimes = 0;
    if (haveTai == 0) {
    redo2:;
        tai = 60676.0 + atmosphere.random.uniform()*(36525);
        //tai = currmjd + atmosphere.random.uniform()*(36525*2) - (36525);
        double lstt = localSiderealTime(tai, longitude);
        double tsunAlt, tsunAz, tsunRa, tsunDec;
        sunPos(tai,&tsunRa,&tsunDec);
        raDectoAltAz(lstt, latitude, tsunRa, tsunDec, &tsunAlt, &tsunAz);
        if ((tsunAlt > -18) && (spaceMode == 0)) goto redo2;
    }
    if (haveMjdString != 0) {
	if (mjdstring == "now" || mjdstring == "tonight") {
	    tai = currmjd;
	redo:;
	    if (mjdstring == "tonight") tai = currmjd + atmosphere.random.uniform();
	    double lstt = localSiderealTime(tai, longitude);
	    double tsunAlt, tsunAz, tsunRa, tsunDec;
	    sunPos(tai,&tsunRa,&tsunDec);
	    raDectoAltAz(lstt, latitude, tsunRa, tsunDec, &tsunAlt, &tsunAz);
	    if ((mjdstring == "tonight") && (tsunAlt > -18) && (spaceMode == 0)) goto redo;
	    if ((mjdstring == "now") && (tsunAlt > -18) && (spaceMode == 0)) {
		std::cout << "Warning:  This observation is in daytime or twilight." << std::endl;
	    }
	}
    }


    //Local Sidereal Time
    double lst = localSiderealTime(tai, longitude);

    double newra=0.0, newdec=0.0;

    if ((havePointingRA == 0) && (havePointingDec == 0) && (haveAltitude == 0) && (haveAzimuth == 0)) {
        do {
            pointingra = atmosphere.random.uniform()*360.0;
            pointingdec = (acos(2.0*atmosphere.random.uniform()-1.0)-PI/2.0)/DEGREE;
            raDectoAltAz(lst, latitude, pointingra, pointingdec, &altitude, &azimuth);
        } while (altitude < 30.0);
    } else {
        if ((havePointingRA ==1) && (havePointingDec==1)  && (haveAltitude == 0) && (haveAzimuth == 0)) {
            raDectoAltAz(lst, latitude, pointingra, pointingdec, &altitude, &azimuth);
            countTimes++;
            if ((altitude < 30.0) && (countTimes < 1000) && (spaceMode<=0)) goto redo2;
            if ((altitude < 30.0) && (spaceMode<=0)) std::cout << "Warning:  Could not find an altitude above 30 degrees." << std::endl;
        } else {
            if (haveAltitude == 0) {
                double temp =atmosphere.random.uniform();
                altitude = (acos(0.5*temp-1.0)-PI/2.0)/DEGREE;
            }
            if (haveAzimuth == 0) azimuth = atmosphere.random.uniform()*360.;
            altAztoRaDec(lst, latitude, altitude, azimuth, &newra, &newdec);
            if (havePointingRA == 0) {
                pointingra=newra;
            } else {
                if ((haveAltitude != 0) && (haveAzimuth !=0)) std::cout << "Warning Redundant Setting: Using the Right Ascension input as " << pointingra << " instead of " << newra << "." << std::endl;
            }
            if (havePointingDec == 0) {
                pointingdec=newdec;
            } else {
                if ((haveAltitude != 0) && (haveAzimuth !=0)) std::cout << "Warning Redundant Setting: Using the Declination input as " << pointingdec << " instead of " << newdec << "." << std::endl;
            }
        }
    }
    zenith = 90 - altitude;

    //Calendar Date
    long year, month, day;
    double frac;
    jdToCalendar(tai, &month, &day, &year, &frac);
    if (haveMonth == 1) {
        std::cout << "Warning Redundant Setting: Using the Month input as " << month << " instead of " << monthnum << "." << std::endl;
    } else {
        monthnum = month;
    }

    //Day of year
    double dayOfYear = calculateDayOfYear(month, day, year, frac);


    //Sun position
    double sunAlt, sunAz, sunRa, sunDec;
    sunPos(tai,&sunRa,&sunDec);
    raDectoAltAz(lst, latitude, sunRa, sunDec, &sunAlt, &sunAz);
    double distSun = distSphere(pointingra, pointingdec, sunRa, sunDec);
    if (haveSunAlt == 1) {
        std::cout << "Warning Redundant Setting: Using the Sun Altitude input as " << sunalt << " instead of " << sunAlt << "." << std::endl;
    } else {
        sunalt = sunAlt;
    }
    if (haveSunAz == 1) {
        std::cout << "Warning Redundant Setting: Using the Sun Azimuth input as " << sunaz << " instead of " << sunAz << "." << std::endl;
    } else {
        sunaz = sunAz;
    }
    if (haveSunRa == 1) {
        std::cout << "Warning Redundant Setting: Using the Sun Right Ascension input as " << sunra << " instead of " << sunRa << "." << std::endl;
    } else {
        sunra = sunRa;
    }
    if (haveSunDec == 1) {
        std::cout << "Warning Redundant Setting: Using the Sun Declination input as " << sunaz << " instead of " << sunDec << "." << std::endl;
    } else {
        sundec = sunDec;
    }


    //Moon position
    double moonAlt, moonAz, moonRa, moonDec, newMoonPhase;
    moonPos(tai,&moonRa,&moonDec);
    if ((haveMoonAlt==1) && (haveMoonAz==1)) {
        altAztoRaDec(lst, latitude, moonalt, moonaz, &moonRa, &moonDec);
    } else {
        raDectoAltAz(lst, latitude, moonRa, moonDec, &moonAlt, &moonAz);
    }
    newMoonPhase = moonPhase(sunRa,sunDec,moonRa,moonDec);
    double distMoon = distSphere(pointingra, pointingdec, moonRa, moonDec);
    if (haveMoonRa == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon RA input as " << moonra << " instead of " << moonRa << "." << std::endl;
    } else {
        moonra = moonRa;
    }
    if (haveMoonDec == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Dec input as " << moondec << " instead of " << moonDec << "." << std::endl;
    } else {
        moondec = moonDec;
    }
    if (haveMoonAlt == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Alt input as " << moonalt << " instead of " << moonAlt << "." << std::endl;
    } else {
        moonalt = moonAlt;
    }
    if (haveMoonAz == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Az input as " << moonaz << " instead of " << moonAz << "." << std::endl;
    } else {
        moonaz = moonAz;
    }
    if (haveMoonPhase == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Phase input as " << moonphase << " instead of " << newMoonPhase << "." << std::endl;
    } else {
        moonphase = newMoonPhase;
    }
    if (haveDistMoon == 1) {
        std::cout << "Warning Redundant Setting: Using the Moon Distance input as " << dist2moon << " instead of " << distMoon << "." << std::endl;
    } else {
        dist2moon = distMoon;
    }

    double ha = lst*15.0 - pointingra;

    // rottelpos = rotator angle
    // rotskypos + rotator = parallactic angle

    // for space instrument rotator should be 0 and also parallactic angle should be anything, so always make rottelpos =0 and 
    double newRotSkyPos = 0.0;
    if (haveRotTelPos == 1) {
        newRotSkyPos=calculateRotSkyPos(ha, latitude, pointingdec, rotTelPos);
        if (spaceMode > 0) newRotSkyPos = atmosphere.random.uniform()*360.0 - rotTelPos;
        if (haveRotSkyPos == 1) {
            if (spaceMode <= 0) std::cout << "Warning Redundant Setting: Using the Parallactic - Rotator input as " << rotSkyPos << " instead of " << newRotSkyPos << std::endl;
        } else {
            rotSkyPos = newRotSkyPos;
        }
    } else {
        if (haveRotSkyPos == 1) {
            rotTelPos=calculateRotTelPos(ha, latitude, pointingdec, rotSkyPos);
            if (spaceMode > 0) rotTelPos = 0.0;
        } else {
            rotTelPos = (atmosphere.random.uniform()*2.0 - 1.0)*rotationConstraint;
            rotSkyPos=calculateRotSkyPos(ha, latitude, pointingdec, rotTelPos);
            if (spaceMode > 0) {
                rotTelPos = 0.0;
                rotSkyPos = atmosphere.random.uniform()*360.0 - rotTelPos;
            }
        }
    }

    // Filter
    int numberFilter=0;
    char opticsFile[4096];
    int goodFilter=1;
    while (goodFilter) {
        snprintf(opticsFile, sizeof(opticsFile), "%s/optics_%d.txt", atmosphere.instrdir.c_str(), numberFilter);
        std::ifstream f(opticsFile);
        goodFilter = f.good();
        if (goodFilter) numberFilter++;
    }
    if (haveFilter == 1) {
        if ((filter < 0) || (filter > (numberFilter-1))) {
            std::cout << "Filter/Optics Configuration is out of range." << std::endl;
            exit(1);
        }
    } else {
        filter = floor(atmosphere.random.uniform()*numberFilter);
    }



    double artificial=0.0;
    char filename[4096];
    snprintf(filename, sizeof(filename), "%s/atmosphere/artificial.txt", atmosphere.datadir.c_str());
    FILE *f = fopen(filename, "r");
    long lines = 0;
    long c = 0;
    while (c != EOF) {
        c = fgetc(f);
        if (c == '\n') lines++;
    }
    fclose(f);
    f = fopen(filename, "r");
    double mindist=1e30;
    for (long i = 0; i < lines ; i++) {
        double phi = 0.0, theta = 0.0, p = 0.0;
        fscanf(f,"%lf %lf %lf\n",&phi,&theta,&p);
        double distance = acos(sin(latitude*DEGREE)*sin(theta) + cos(latitude*DEGREE)*cos(theta)*cos(longitude*DEGREE-phi));
        if (distance < mindist) {
            mindist=distance;
            artificial=p;
        }
    }
    fclose(f);
    if (artificial>0) artificial = -2.5*log10(artificial) + 27.6; else artificial=1000.0;

    // Weather & Climate
        double oceanDistance = 1000.0;
        double populationDensity = 0.0;
        double populationDensityCL[4];
        double temperatureMean=0.0;
        double temperatureRange=0.0;
        double temperatureDiurnal=0.0;
        double temperatureDiurnalRange=0.0;
        double pressureMean=0.0;
        double pressureRange=0.0;
        double pressureDiurnal=0.0;
        double pressureDiurnalRange=0.0;
        double humidityMean=0.0;
        double humidityRange=0.0;
        double humidityDiurnal=0.0;
        double humidityDiurnalRange=0.0;
        double temperatureMeanCL[4];
        double temperatureRangeCL[4];
        double temperatureDiurnalCL[4];
        double temperatureDiurnalRangeCL[4];
        double humidityMeanCL[4];
        double humidityRangeCL[4];
        double humidityDiurnalCL[4];
        double humidityDiurnalRangeCL[4];
        double altCL[4];
        double averageAlt=0.0;
        altCL[0] = 0.0;
        altCL[1] = 0.0;
        altCL[2] = 0.0;
        altCL[3] = 0.0;
        double hour = frac*24.0;
        double smoothLength = 1000.0/RADIUS_EARTH;
        double Tmax, Tmin, Tmean;
        populationDensityCL[0] = 0.0;
        populationDensityCL[1] = 0.0;
        populationDensityCL[2] = 0.0;
        populationDensityCL[3] = 0.0;
    if (spaceMode <= 0) {
        // double meanHigh;
        // if (std::abs(latitude) < 18) {
        //     meanHigh = 27.0;
        // } else {
        //     meanHigh = 27.0 - 0.75 * (std::abs(latitude) - 18.0);
        // }
        char filename[4096];
        snprintf(filename, sizeof(filename), "%s/atmosphere/earth.txt", atmosphere.datadir.c_str());
        FILE *f = fopen(filename, "r");
        long lines = 0, c = 0;
        while (c != EOF) {
            c = fgetc(f);
            if (c == '\n') lines++;
        }
        fclose(f);
        f = fopen(filename, "r");
        double mndist[5];
        mndist[0] = 1e30;
        mndist[1] = 1e30;
        mndist[2] = 1e30;
        mndist[3] = 1e30;
        mndist[4] = 1e30;
        for (long i = 0; i < lines; i++) {
            double phi = 0.0, theta = 0.0, o = 0.0, z = 0.0, l =0.0, of = 0.0, lf = 0.0;
            fscanf(f,"%lf %lf %lf %lf %lf %lf %lf\n",&phi,&theta,&o,&z,&l,&of,&lf);
            double distance = acos(sin(latitude*DEGREE)*sin(theta) + cos(latitude*DEGREE)*cos(theta)*cos(longitude*DEGREE-phi));
            double posang = atan2(sin(longitude*DEGREE-phi), cos(theta)*tan(latitude*DEGREE) - sin(theta)*cos(longitude*DEGREE-phi));
            double xx = distance*cos(posang);
            double yy = distance*sin(posang);
            if (distance < mndist[0]) {
                mndist[0]=distance;
                oceanDistance = o;
                averageAlt = z;
            }
            double x0, y0, rr;
            x0 = smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[1]) {
                mndist[1]=rr;
                altCL[0] = z;
            }
            x0 = -smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[2]) {
                mndist[2]=rr;
                altCL[1] = z;
            }
            x0 = 0.0;
            y0 = smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[3]) {
                mndist[3]=rr;
                altCL[2] = z;
            }
            x0 = 0.0;
            y0 = -smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[4]) {
                mndist[4]=rr;
                altCL[3] = z;
            }
        }
        fclose(f);

        snprintf(filename, sizeof(filename), "%s/atmosphere/population.txt", atmosphere.datadir.c_str());
        f = fopen(filename, "r");
        lines = 0;
        c = 0;
        while (c != EOF) {
            c = fgetc(f);
            if (c == '\n') lines++;
        }
        fclose(f);
        f = fopen(filename, "r");
        mndist[0] = 1e30;
        mndist[1] = 1e30;
        mndist[2] = 1e30;
        mndist[3] = 1e30;
        mndist[4] = 1e30;
        for (long i = 0; i < lines ; i++) {
            double phi = 0.0, theta = 0.0, p = 0.0;
            fscanf(f,"%lf %lf %lf\n",&phi,&theta,&p);
            double distance = acos(sin(latitude*DEGREE)*sin(theta) + cos(latitude*DEGREE)*cos(theta)*cos(longitude*DEGREE-phi));
            double posang = atan2(sin(longitude*DEGREE-phi), cos(theta)*tan(latitude*DEGREE) - sin(theta)*cos(longitude*DEGREE-phi));
            double xx = distance*cos(posang);
            double yy = distance*sin(posang);
                
            if (distance < mndist[0]) {
                mndist[0]=distance;
                populationDensity = p;
            }
            double x0, y0, rr;
            x0 = smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[1]) {
                mndist[1]=rr;
                populationDensityCL[0] = p;
            }
            x0 = -smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[2]) {
                mndist[2]=rr;
                populationDensityCL[1] = p;
            }
            x0 = 0.0;
            y0 = smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[3]) {
                mndist[3]=rr;
                populationDensityCL[2] = p;
            }
            x0 = 0.0;
            y0 = -smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[4]) {
                mndist[4]=rr;
                populationDensityCL[3] = p;
            }
        }
        fclose(f);

        snprintf(filename, sizeof(filename), "%s/atmosphere/climate.txt", atmosphere.datadir.c_str());
        f = fopen(filename, "r");
        lines = 0;
        c = 0;
        while (c != EOF) {
            c = fgetc(f);
            if (c == '\n') lines++;
        }
        fclose(f);
        f = fopen(filename, "r");
        mndist[0] = 1e30;
        mndist[1] = 1e30;
        mndist[2] = 1e30;
        mndist[3] = 1e30;
        mndist[4] = 1e30;
        for (long i = 0; i < lines ; i++) {
            double phi = 0.0, theta=0.0, t=0.0, ts=0.0, td=0.0, tds=0.0, p=0.0, ps=0.0, pd=0.0, pds=0.0, h=0.0, hs=0.0, hd=0.0, hds=0.0;
            int mo=0;
                fscanf(f,"%lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&phi,&theta,&mo,&t,&ts,&td,&tds,&p,&ps,&pd,&pds,&h,&hs,&hd,&hds);
            double distance = acos(sin(latitude*DEGREE)*sin(theta*DEGREE) + cos(latitude*DEGREE)*cos(theta*DEGREE)*cos(longitude*DEGREE-phi*DEGREE));
                        double posang = atan2(sin(longitude*DEGREE-phi*DEGREE), cos(theta*DEGREE)*tan(latitude*DEGREE) - sin(theta*DEGREE)*cos(longitude*DEGREE-phi*DEGREE));
            double xx = distance*cos(posang);
            double yy = distance*sin(posang);
                

            if ((distance < mndist[0]) && (mo == (month-1))) {
                mndist[0]=distance;
                temperatureMean = t;
                temperatureRange = ts;
                temperatureDiurnal = td;
                temperatureDiurnalRange = tds;
                pressureMean = p;
                pressureRange = ps;
                pressureDiurnal = pd;
                pressureDiurnalRange = pds;
                humidityMean = h;
                humidityRange = hs;
                humidityDiurnal = hd;
                humidityDiurnalRange = hds;
            }
            double x0, y0, rr;
            x0 = smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[1]) {
                mndist[1]=rr;
                temperatureMeanCL[0] = t;
                temperatureRangeCL[0] = ts;
                temperatureDiurnalCL[0] = td;
                temperatureDiurnalRangeCL[0] = tds;
                humidityMeanCL[0] = h;
                humidityRangeCL[0] = hs;
                humidityDiurnalCL[0] = hd;
                humidityDiurnalRangeCL[0] = hds;
            }
            x0 = -smoothLength;
            y0 = 0.0;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[2]) {
                mndist[2]=rr;
                temperatureMeanCL[1] = t;
                temperatureRangeCL[1] = ts;
                temperatureDiurnalCL[1] = td;
                temperatureDiurnalRangeCL[1] = tds;
                humidityMeanCL[1] = h;
                humidityRangeCL[1] = hs;
                humidityDiurnalCL[1] = hd;
                humidityDiurnalRangeCL[1] = hds;
            }
            x0 = 0.0;
            y0 = smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[3]) {
                mndist[3]=rr;
                temperatureMeanCL[2] = t;
                temperatureRangeCL[2] = ts;
                temperatureDiurnalCL[2] = td;
                temperatureDiurnalRangeCL[2] = tds;
                humidityMeanCL[2] = h;
                humidityRangeCL[2] = hs;
                humidityDiurnalCL[2] = hd;
                humidityDiurnalRangeCL[2] = hds;
            }
            x0 = 0.0;
            y0 = -smoothLength;
            rr = sqrt(pow(xx - x0, 2.0) + pow(yy- y0, 2.0));
            if (rr < mndist[4]) {
                mndist[4]=rr;
                temperatureMeanCL[3] = t;
                temperatureRangeCL[3] = ts;
                temperatureDiurnalCL[3] = td;
                temperatureDiurnalRangeCL[3] = tds;
                humidityMeanCL[3] = h;
                humidityRangeCL[3] = hs;
                humidityDiurnalCL[3] = hd;
                humidityDiurnalRangeCL[3] = hds;
            }
        }
        fclose(f);

        //mean should have no multiply
        //range should have no multiply
        //diurnal should multiply by 0.11*0.5
        // diurnal range should multiply by 0.11*0.5
        
        for (int kk = 0; kk < 4; kk++) {
            double TmeanCL = temperatureMeanCL[kk] - 273.15 + temperatureRangeCL[kk]*atmosphere.random.normalCorrel(tai + 70000 + (kk+1)*3000, 2.0);
            double dayNightDifferenceCL = temperatureDiurnalCL[kk] + temperatureDiurnalRangeCL[kk]*atmosphere.random.normalCorrel(tai + 71000 + (kk+1)*3000, 2.0);
        // sampling error
        dayNightDifferenceCL = dayNightDifferenceCL * 1.11;
        double TmaxCL = TmeanCL + 0.5*dayNightDifferenceCL;
        double TminCL = TmeanCL - 0.5*dayNightDifferenceCL;
        temperatureMaxCL[kk] = TmaxCL;
        temperatureMinCL[kk] = TminCL;
        temperatureCL[kk] = TmeanCL;
        }

        double bb = atmosphere.random.normalCorrel(tai + 70000, 2.0);
        Tmean = temperatureMean - 273.15 + temperatureRange*bb;
        double dayNightDifference = temperatureDiurnal + temperatureDiurnalRange*atmosphere.random.normalCorrel(tai + 71000, 2.0);
        // sampling error
        dayNightDifference = dayNightDifference * 1.11;
        Tmax = Tmean + 0.5*dayNightDifference;
        Tmin = Tmean - 0.5*dayNightDifference;

        double delta = 23.45*sin(360.0/365.0*(284.0+dayOfYear));
        double omega = 2.0/15.0*acos(-tan(latitude*DEGREE)*tan(delta*DEGREE))/DEGREE;
        double thermalMax = 12.0 + 2.2;
        double tSunset = 12.0 + omega/2.0;
        double localTemperature;
        if ((hour >= thermalMax - omega/2.0) && (hour <= tSunset)) {
            localTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(hour - thermalMax));
        } else {
            double k = omega/PI*atan(PI/omega*(tSunset-thermalMax));
            if (hour < thermalMax - omega/2.0) hour = hour + 24.0;
            localTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(tSunset - thermalMax))*exp(-(hour-tSunset)/k);
        }
        if (haveTemperature == 0) temperature = localTemperature;
        // previous hour temperature
        double prTemperature;
        double prHour = frac*24.0 - 1.0;
        if ((prHour >= thermalMax - omega/2.0) && (prHour <= tSunset)) {
            prTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(prHour - thermalMax));
        } else {
            double k = omega/PI*atan(PI/omega*(tSunset-thermalMax));
            if (prHour < thermalMax - omega/2.0) prHour = prHour + 24.0;
            prTemperature = Tmin + (Tmax - Tmin)*cos(PI/omega*(tSunset - thermalMax))*exp(-(prHour-tSunset)/k);
        }
        //        prTemperature += randomT*atmosphere.random.normalCorrel(tai + 12000 - 1.0/24.0, 2.0);
        //prTemperature -= lapseRate*atmosphere.groundlevel/1000.0;
        if (haveTemperatureVariation == 0) temperatureVariation = (localTemperature - prTemperature);

    } else {
        Tmean = 20.0;
        Tmax = 20.0;
        Tmin = 20.0;
        if (haveTemperature == 0) temperature = 20.0;

        if (spaceMode==3) {
            double prTemperature;
            prTemperature = -233.0 + 12.0*atmosphere.random.normalCorrel(tai + 13000 - 1.0/24.0, 2.0);
            if (haveTemperature ==0) {
                temperature = -233.0 + 12.0*atmosphere.random.normalCorrel(tai + 13000, 2.0);
            }
            if (haveTemperatureVariation == 0) temperatureVariation = (temperature - prTemperature);

        }
    }
        for (int kk = 0; kk < 4; kk++) {
            predictedRainfallCL[kk]=0;
            for (int ll =0; ll <= 1; ll++) {
                double waterVaporPressureCL = 0.0;
                double humidityMuCL = humidityMeanCL[kk]/100.0;
                if (ll == 0) {
                    humidityMuCL -= 0.5*1.11*humidityDiurnalCL[kk]/100.0;
                } else {
                    humidityMuCL += 0.5*1.11*humidityDiurnalCL[kk]/100.0;
                }
                double humiditySigmaCL = sqrt(pow(humidityRangeCL[kk], 2.0) + pow(0.5*1.11*humidityDiurnalRangeCL[kk], 2.0))/100.0;
                double localtemperature;
                 if (ll == 0) localtemperature= temperatureMaxCL[kk]; else localtemperature = temperatureMinCL[kk];
                double temperatureKelvinCL = 273.15 + localtemperature;
                if (localtemperature > 0) {
                    double Tc = 647.096;
                    double theta = 1.0 - temperatureKelvinCL/Tc;
                    double Pc=220640.0;
                    double c1=-7.85951783;
                    double c2=1.84408259;
                    double c3=-11.7866497;
                    double c4=22.6807411;
                    double c5=-15.9618719;
                    double c6=1.80122502;
                    waterVaporPressureCL = Pc*100.0/101325.0*760.0*
                        exp(Tc/temperatureKelvinCL*(c1*theta+c2*pow(theta,1.5)+c3*pow(theta,3.0)+
                                                    c4*pow(theta,3.5)+c5*pow(theta,4.0)+
                                                    c6*pow(theta,7.5)));
                } else {
                    double Tn = 273.16;
                    double theta = temperatureKelvinCL/Tn;
                    double Pn = 6.11657;
                    double a0 = -13.928169;
                    double a1 = 34.707823;
                    waterVaporPressureCL = Pn*100.0/101325.0*760.0*
                        exp(a0*(1 - pow(theta, -1.5)) + a1*(1 - pow(theta, -1.25)));
                }

                // predicted rainfall in centimeters if the humidity were to go above 100
                double humidFactorCL=0.0;
                for (long i = 0; i < 10000; i++) {
                    double humid = humidityMuCL + humiditySigmaCL*atmosphere.random.normal();
                    if (humid >= 1.0) {
                        humid = humid - 1.0;
                    } else {
                        humid = 0.0;
                    }
                    humidFactorCL += humid/10000.0;
                }
                predictedRainfallCL[kk] += 0.5*waterVaporPressureCL*PRESS_ATMOSPHERE/PRESS_ATMOSPHERE_MMHG/
                    GRAV_ACCEL_EARTH/(DENSITY_H2O*1000.0)*100.0*humidFactorCL;
            }
        }


        double waterVaporPressure = 0.0;
        double humidityBase = humidityMean/100.0 + humidityRange*atmosphere.random.normalCorrel(tai + 72000, 2.0)/100.0;
        double humidityAmplitude = 0.5*1.11*humidityDiurnal/100.0 + 0.5*1.11*humidityDiurnalRange*atmosphere.random.normalCorrel(tai + 73000, 2.0)/100.0;
        //double humidityMu = humidityMean/100.0 + 0.5*1.11*humidityDiurnal/100.0*cos(2*PI*(hour-3.0)/24.0);
        //double humiditySigma = sqrt(pow(humidityRange, 2.0) + pow(0.5*1.11*humidityDiurnalRange, 2.0))/100.0;
        double relativeHumidity = humidityBase + humidityAmplitude*cos(2*PI*((hour-3.0)/24.0));
        if (haveHumidity) relativeHumidity = humidity;
        if (relativeHumidity > 1.0) {
            if (temperature > 0) std::cout << "Warning:  It is raining." << std::endl;
            if (temperature <= 0) std::cout << "Warning:  It is snowing." << std::endl;
            relativeHumidity=1.0;
        }
        if (relativeHumidity < 0.0) relativeHumidity=0.0;

    double predictedRainfall = 0.0;
    for (int ll =0; ll <= 1; ll++) {
        double waterVaporPressure = 0.0;
        double humidityMu = humidityMean/100.0;
        if (ll == 0) {
            humidityMu -= 0.5*1.11*humidityDiurnal/100.0;
        } else {
            humidityMu += 0.5*1.11*humidityDiurnal/100.0;
        }
        double humiditySigma = sqrt(pow(humidityRange, 2.0) + pow(0.5*1.11*humidityDiurnalRange, 2.0))/100.0;
        double localtemperature;
        if (ll == 0) localtemperature= Tmax; else localtemperature = Tmin;
        double temperatureKelvin = localtemperature + 273.15;
        if (localtemperature > 0) {
            double Tc = 647.096;
            double theta = 1.0 - temperatureKelvin/Tc;
            double Pc=220640.0;
            double c1=-7.85951783;
            double c2=1.84408259;
            double c3=-11.7866497;
            double c4=22.6807411;
            double c5=-15.9618719;
            double c6=1.80122502;
            waterVaporPressure = Pc*100.0/101325.0*760.0*
                exp(Tc/temperatureKelvin*(c1*theta+c2*pow(theta,1.5)+c3*pow(theta,3.0)+
                                          c4*pow(theta,3.5)+c5*pow(theta,4.0)+
                                          c6*pow(theta,7.5)));
        } else {
            double Tn = 273.16;
            double theta = temperatureKelvin/Tn;
            double Pn = 6.11657;
            double a0 = -13.928169;
            double a1 = 34.707823;
            waterVaporPressure = Pn*100.0/101325.0*760.0*
                exp(a0*(1 - pow(theta, -1.5)) + a1*(1 - pow(theta, -1.25)));
        }

    // predicted rainfall in centimeters if the humidity were to go above 100
        double humidFactor=0.0;
        for (long i = 0; i < 10000; i++) {
            double humid = humidityMu + humiditySigma*atmosphere.random.normal();
            if (humid >= 1.0) {
                humid = humid - 1.0;
            } else {
                humid = 0.0;
            }
            humidFactor += humid/10000.0;
        }
        predictedRainfall += 0.5*waterVaporPressure*PRESS_ATMOSPHERE/PRESS_ATMOSPHERE_MMHG/
        GRAV_ACCEL_EARTH/(DENSITY_H2O*1000.0)*100.0*humidFactor;
    }



    double temperatureKelvin = temperature + 273.15;
    if (temperature > 0) {
        double Tc = 647.096;
        double theta = 1.0 - temperatureKelvin/Tc;
        double Pc=220640.0;
        double c1=-7.85951783;
        double c2=1.84408259;
        double c3=-11.7866497;
        double c4=22.6807411;
        double c5=-15.9618719;
        double c6=1.80122502;
        waterVaporPressure = Pc*100.0/101325.0*760.0*
            exp(Tc/temperatureKelvin*(c1*theta+c2*pow(theta,1.5)+c3*pow(theta,3.0)+
                                      c4*pow(theta,3.5)+c5*pow(theta,4.0)+
                                      c6*pow(theta,7.5)));
    } else {
        double Tn = 273.16;
        double theta = temperatureKelvin/Tn;
        double Pn = 6.11657;
        double a0 = -13.928169;
        double a1 = 34.707823;
        waterVaporPressure = Pn*100.0/101325.0*760.0*
            exp(a0*(1 - pow(theta, -1.5)) + a1*(1 - pow(theta, -1.25)));
    }
    waterVaporPressure=relativeHumidity*waterVaporPressure;


    if (havePressure == 0) {
        if (spaceMode <= 0) {
            double pressureMu = pressureMean  + pressureRange*atmosphere.random.normalCorrel(tai + 75000, 2.0);
            double pressureSigma = pressureDiurnal + pressureDiurnalRange*atmosphere.random.normalCorrel(tai + 76000, 2.0);
            //sampling correction
            pressureSigma = pressureSigma * 1.57 * 0.5;
            pressure = pressureMu + pressureSigma*cos(2*PI*((hour-10.5)/12.0));
            pressure = pressure * PRESS_ATMOSPHERE_MMHG / PRESS_ATMOSPHERE;
        } else {
            pressure = 760.0;
        }
    }
    
    
    // Aerosols
    for (int kk = 0; kk < 4; kk++) if (temperatureCL[kk]<0) temperatureCL[kk]=0.0;
    double seasaltOpticalDepth = 0.0, dustOpticalDepth = 0.0, smokeOpticalDepth = 0.0, pollutionOpticalDepth=0.0;
    double seasaltIndex = 0.0, dustIndex = 0.0, smokeIndex = 0.0, pollutionIndex = 0.0;
    double seasaltCurvature = 0.0, dustCurvature = 0.0, smokeCurvature = 0.0, pollutionCurvature = 0.0;
    if (spaceMode <= 0) {
        double scaledTemperature[5];
        double scaledRainfall[5];
        double scaledPopulation[5];
        double weight[5];

        for (int kk = 0; kk<5; kk++) {
            if (kk==0) scaledTemperature[kk] = Tmean/(49.0);
            if (kk!=0) scaledTemperature[kk] = temperatureCL[kk-1]/(49.0);
            if (scaledTemperature[kk]<0) scaledTemperature[kk]=0.0;
            if (kk==0) scaledRainfall[kk] = predictedRainfall/2.3;
            if (kk!=0) scaledRainfall[kk] = predictedRainfallCL[kk-1]/2.3;
            if (scaledRainfall[kk]<0) scaledRainfall[kk]=0.0;
            if (kk==0) scaledPopulation[kk] = populationDensity /16389.0;
            if (kk!=0) scaledPopulation[kk] = populationDensityCL[kk-1] /16389.0;
            if (scaledPopulation[kk]<0) scaledPopulation[kk]=0.0;
            if (kk==0) weight[kk]=1.0;
            if (kk!=0) weight[kk]=exp(-1.0);
        }
        for (int kk=0; kk<5; kk++) {
            weight[kk]=weight[kk]/(1.0+4.0*exp(-1.0));
        }

        double seasalt[5];
        double dust[5];
        double smoke[5];
        double pollution[5];
        double land[5];
        land[0]=0.0;
        land[1]=0.0;
        land[2]=0.0;
        land[3]=0.0;
        land[4]=0.0;
        if (averageAlt > 0.0) land[0] = 1.0;
        if (altCL[0] > 0.0) land[1] = 1.0;
        if (altCL[1] > 0.0) land[2] = 1.0;
        if (altCL[2] > 0.0) land[3] = 1.0;
        if (altCL[3] > 0.0) land[4] = 1.0;
        for (int kk =0 ; kk < 5 ; kk++) {
            seasalt[kk] = (1.0 - land[kk])*weight[kk];
            dust[kk] = -(scaledRainfall[kk] - 0.2*scaledTemperature[kk]);
            if (dust[kk] < 0.0) dust[kk]=0.0;
            dust[kk] = dust[kk]/0.12*land[kk]*weight[kk];
            smoke[kk] = (scaledRainfall[kk] - 0.2*scaledTemperature[kk]);
            if (smoke[kk] < 0.0) smoke[kk]=0.0;
            smoke[kk] = smoke[kk]/0.06*land[kk]*weight[kk];
            pollution[kk] = scaledPopulation[kk]/0.14*weight[kk];
        }

        double seasaltTotal=0.0;
        double dustTotal=0.0;
        double smokeTotal=0.0;
        double pollutionTotal=0.0;
        for (int kk=0; kk< 5; kk++) {
            seasaltTotal += seasalt[kk];
            dustTotal += dust[kk];
            smokeTotal += smoke[kk];
            pollutionTotal += pollution[kk];

        }

        seasaltOpticalDepth = seasaltTotal*exp(-2.56 + 0.29*atmosphere.random.normalCorrel(tai+55000, 2.0));
        dustOpticalDepth = dustTotal*exp(-0.80 + 0.46*atmosphere.random.normalCorrel(tai+56000, 2.0));
        smokeOpticalDepth = smokeTotal*exp(-1.83 + 0.80*atmosphere.random.normalCorrel(tai+57000, 2.0));
        pollutionOpticalDepth = pollutionTotal*exp(-1.28 + 0.40*atmosphere.random.normalCorrel(tai+58000, 2.0));

        seasaltIndex = -0.59 + 0.23*atmosphere.random.normalCorrel(tai+59000, 2.0);
        dustIndex = -0.20 + 0.20*atmosphere.random.normalCorrel(tai+60000, 2.0);
        smokeIndex = -1.33 + 0.26*atmosphere.random.normalCorrel(tai+61000, 2.0);
        pollutionIndex = -1.09 + 0.21*atmosphere.random.normalCorrel(tai+62000, 2.0);

        seasaltCurvature = 0.22 + 0.26*atmosphere.random.normalCorrel(tai+63000, 2.0);
        dustCurvature = -0.18 + 0.21*atmosphere.random.normalCorrel(tai+64000, 2.0);
        smokeCurvature = 0.17 + 0.40*atmosphere.random.normalCorrel(tai+65000, 2.0);
        pollutionCurvature = -0.03 + 0.23*atmosphere.random.normalCorrel(tai+66000, 2.0);

    }

    
    // Calculated inputs
    snprintf(outputfilename, sizeof(outputfilename), "%s/obs_%s.pars", workdir.c_str(), obshistid.c_str());
    fptr = fopen(outputfilename, "a+");
    fprintf(fptr, "pointingra %lf\n", pointingra);
    fprintf(fptr, "pointingdec %lf\n", pointingdec);
    fprintf(fptr, "altitude %lf\n", altitude);
    fprintf(fptr, "azimuth %lf\n", azimuth);
    fprintf(fptr, "tai %lf\n", tai);
    fprintf(fptr, "moonra %lf\n", moonra);
    fprintf(fptr, "moondec %lf\n", moondec);
    fprintf(fptr, "phaseang %lf\n", moonphase);
    fprintf(fptr, "moondist %lf\n", dist2moon);
    fprintf(fptr, "moonalt %lf\n", moonalt);
    fprintf(fptr, "moonaz %lf\n", moonaz);
    fprintf(fptr, "monthnum %ld\n", monthnum);
    fprintf(fptr, "solaralt %lf\n", sunalt);
    fprintf(fptr, "solaraz %lf\n", sunaz);
    fprintf(fptr, "solarra %lf\n", sunra);
    fprintf(fptr, "solardec %lf\n", sundec);
    fprintf(fptr, "solardist %lf\n", distSun);
    fprintf(fptr, "spiderangle %lf\n", rotTelPos);
    fprintf(fptr, "rotationangle %lf\n", rotSkyPos);
    fprintf(fptr, "filter %d\n", filter);
    fprintf(fptr, "dayofyear %lf\n", dayOfYear);
    fprintf(fptr, "temperature %lf\n", temperature);
    fprintf(fptr, "pressure %lf\n", pressure);
    fprintf(fptr, "tempvar %lf\n", temperatureVariation);
    fprintf(fptr, "waterpressure %lf\n", waterVaporPressure);
    fprintf(fptr, "hourangle %lf\n", ha);
    fprintf(fptr, "lst %lf\n", lst);

    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;

    printf("[rightascension] Pointing Right Ascension (degrees)        %lf\n", pointingra);
    printf("[declination] Pointing Declination (degrees)               %lf\n", pointingdec);
    if (spaceMode == 0) {
        printf("[altitude] Pointing Altitude (degrees)                     %lf\n", altitude);
        printf("[azimuth] Pointing Azimuth (degrees)                       %lf\n", azimuth);
    }
    printf("[mjd] Modified Julian Date (days)                          %lf\n", tai);
    printf("[date] Date (month/day/year/day-fraction)                  %02ld/%02ld/%2ld/%lf\n", monthnum,day,year,frac);
    printf("[rottelpos] Rotator Angle (degrees)                        %lf\n", rotTelPos);
    printf("[moonra] Moon Right Ascension (degrees)                    %lf\n", moonra);
    printf("[moondec] Moon Declination (degrees)                       %lf\n", moondec);
    printf("[moonphase] Moon Phase (percent)                           %lf\n", moonphase);
    printf("[moonalt] Altitude of Moon (degrees)                       %lf\n", moonalt);
    printf("[moonaz] Azimuth of Moon (degrees)                         %lf\n", moonaz);
    printf("[sunalt] Altitude of Sun (degrees)                         %lf\n", sunalt);
    printf("[sunaz] Azimuth of Sun (degrees)                           %lf\n", sunaz);
    printf("[sunra] Right Ascension of Sun (degrees)                   %lf\n", sunra);
    printf("[sundec] Declination of Sun (degrees)                      %lf\n", sundec);
    printf("[rotskypos] Parallactic - Rotator Angle (degrees)          %lf\n", rotSkyPos);
    printf("[pressure] Ambient Pressure (mmHg)                         %lf\n", pressure);
    printf("[waterpressure] Water Vapor Pressure (mmHg)                %lf\n", waterVaporPressure);
    printf("[humidity] Relative Humidity (%%)                           %lf\n", relativeHumidity*100.0);
    printf("Predicted Rainfall (cm)                                    %lf\n", predictedRainfall);
    printf("[temperature] Ambient Temperature (degrees C)              %lf\n", temperature);
    printf("[tempvar] Temperature Variation (degrees/hr)               %lf\n", temperatureVariation);

    double totalseeing=constrainseeing;
    double exosphereTemperature = 1050.0 + 550.0*sin(2*PI*(tai - 54466.0)/(365.24*11.0)- PI/2.0);

    if (skipAtmosphere == 0) {

        std::cout << "----------------------------------------------------------------------------------------------------" << std::endl
        << "Atmosphere Creator" << std::endl
        << "----------------------------------------------------------------------------------------------------" << std::endl;

        atmosphere.osests.reserve(atmosphere.numlevel);
        atmosphere.altitudes.reserve(atmosphere.numlevel);
        atmosphere.jests.reserve(atmosphere.numlevel);

        snprintf(outputfilename, sizeof(outputfilename), "%s/atmosphere_%s.pars", workdir.c_str(), obshistid.c_str());

        atmosphere.atmospheredir = atmosphere.datadir + "/atmosphere";
        atmosphere.createAtmosphere(monthnum, constrainseeing, outputfilename, cloudscreen, seed, tai, &totalseeing,latitude,oceanDistance,relativeHumidity,temperature, pressure, atmosphere.groundlevel, longitude, waterVaporPressure, exosphereTemperature);
        // overwrite outer scales if they have been provided.
        for (int m(0); m < atmosphere.numlevel; m++) {
            if (outerscale[m] >= 0) atmosphere.osests[m] = outerscale[m];
        }

        faptr = fopen(outputfilename, "a+");
        {
            pthread_mutex_init(&atmosphere.lock, NULL);
            pthread_t *thread;
            thread_args *args;
            long numthread = atmosphere.numlevel;
            thread = (pthread_t*)malloc(numthread*sizeof(pthread_t));
            args = (thread_args*)malloc(numthread*sizeof(thread_args));
            for (int i = 0; i < atmosphere.numlevel; i++) {
                char tempstring2[4096];
                printf("Creating layer %d.\n", i);
                double outers = atmosphere.osests[i]*100.0*100.0;
                snprintf(tempstring2, sizeof(tempstring2), "%s/atmospherescreen_%s_%d", workdir.c_str(), obshistid.c_str(), i);
                //                atmosphere.turb2d(seed*10 + i, seeing, outerx, outers, zenith, 0.5, tempstring2);
                args[i].instance = &atmosphere;
                args[i].seed = seed*10 + i;
                args[i].see5 = seeing;
                args[i].outerx = outerx;
                args[i].outers = outers;
                args[i].zenith = zenith;
                args[i].wavelength = 0.5;
                args[i].resample = ATM_RESAMPLE;
                strcpy(args[i].name,tempstring2);
                args[i].N_size = SCREEN_SIZE;
                pthread_create(&thread[i], NULL, &Atmosphere::threadFunction, &args[i]);
                fprintf(faptr, "atmospherefile %d atmospherescreen_%s_%d\n", i, obshistid.c_str(), i);
                if (threadFlag==0) pthread_join(thread[i], NULL);
            }
            if (threadFlag==1) {
                for (int i = 0; i < atmosphere.numlevel; i++) {
                    pthread_join(thread[i], NULL);
                }
            }
        }
        for (int i = 0; i < atmosphere.numlevel; i++) {
            if (cloudscreen[i]) {
                double height = (atmosphere.altitudes[i] - atmosphere.groundlevel)/1000.0;
                snprintf(tempstring, sizeof(tempstring), "%s/cloudscreen_%s_%d", workdir.c_str(), obshistid.c_str(), i);
                atmosphere.cloud(seed*10 + i, height, pix, tempstring, SCREEN_SIZE);
                fprintf(faptr, "cloudfile %d cloudscreen_%s_%d\n", i, obshistid.c_str(), i);
            }
        }

        snprintf(tempstring, sizeof(tempstring), "%s/airglowscreen_%s", workdir.c_str(), obshistid.c_str());
        atmosphere.airglow(seed*10, tempstring, SCREEN_SIZE);

        atmosphere.random.setSeed32(seed);
        atmosphere.random.unwind(10000);
        if (constrainairglow > 0.0) {
            fprintf(faptr, "airglowpintensity %.3f\n", constrainairglow);
            fprintf(faptr, "airglowcintensity %.3f\n", constrainairglow);
        } else {
            double temp;
            temp = 21.6 + 0.5*atmosphere.random.normalCorrel(tai, 2700.0/(24.0*3600.0))
                + 2.5*log10(cos(zenith*PI/180.0))
                     + (0.50 - 0.50*sin(2*PI*(tai - 54466.0)/(365.24*11.0)
                                        - PI/2.0));
            fprintf(faptr, "airglowpintensity %.3f\n", temp);
            temp = (25.3 + 0.6*atmosphere.random.normalCorrel(tai+10000, 1200.0/(24.0*3600.0))
                    + 2.5*log10(cos(zenith*PI/180.0)));
            fprintf(faptr, "airglowcintensity %.3f\n", temp);
        }
        fclose(faptr);

    }


    fprintf(fptr,"exospheretemperature %lf\n", exosphereTemperature);
    fprintf(fptr, "constrainseeing %lf\n", totalseeing);
    fprintf(fptr, "artificial %lf\n", artificial);

    fprintf(fptr, "aerosolpar 0 %lf\n",seasaltOpticalDepth);
    fprintf(fptr, "aerosolpar 1 %lf\n", dustOpticalDepth);
    fprintf(fptr, "aerosolpar 2 %lf\n", smokeOpticalDepth);
    fprintf(fptr, "aerosolpar 3 %lf\n",pollutionOpticalDepth);
    fprintf(fptr, "aerosolpar 4 %lf\n",seasaltIndex);
    fprintf(fptr, "aerosolpar 5 %lf\n", dustIndex);
    fprintf(fptr, "aerosolpar 6 %lf\n", smokeIndex);
    fprintf(fptr, "aerosolpar 7 %lf\n",pollutionIndex);
    fprintf(fptr, "aerosolpar 8 %lf\n",seasaltCurvature);
    fprintf(fptr, "aerosolpar 9 %lf\n", dustCurvature);
    fprintf(fptr, "aerosolpar 10 %lf\n", smokeCurvature);
    fprintf(fptr, "aerosolpar 11 %lf\n",pollutionCurvature);
    fclose(fptr);

    return(0);
}
