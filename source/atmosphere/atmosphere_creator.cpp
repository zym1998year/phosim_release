///
/// @package phosim
/// @file turb2d.cpp
/// @brief turbulence generator
///
/// @brief Created by:
/// @author Mallory Young (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.

void Atmosphere::createAtmosphere(float monthnum, float constrainseeing,
                                  const std::string & outputfilename, std::vector<int> & cloudscreen, long seed, double tai, double *tseeing,
                                  double latitude, double oceanDistance, double humidity, double temperature, double pressure, double groundlevel,
                                  double longitude, double waterPressure, double exosphereTemperature) {

    random.setSeed32(seed);
    random.unwind(10000);

    float totalj;
    float totalseeing;
    //float bestseeing = 0.0;
    std::vector<float> relseeing(numlevel, 0);

    altitudes[numlevel - 1] = 10.0 + groundlevel;
    altitudes[numlevel - 2] = 20.0 + groundlevel;
    for (int i = 0; i < numlevel - 2; i++) {
        altitudes[i] = pow(2.0, (numlevel - 2) - i/((static_cast<float>(numlevel) - 2)/6.0))
            *(16000.0/pow(2.0, numlevel - 2)) + groundlevel;
    }

    std::vector<float> magests(numlevel, 0);
    std::vector<float> dirests(numlevel, 0);
    for (int i = 0; i < numlevel; i++) {
        if (i < (numlevel-1)) osests[i] = outerscale(altitudes[i], tai, i);
        if (i == (numlevel-1)) osests[i] = 1.0;
        if (i == 0 || i == 1) {
            osests[i] *= 1.5;
        }
        if (osests[i] > (altitudes[i] - groundlevel)/2) {
            osests[i] = (altitudes[i] - groundlevel)/2;
        }
    }

    int haveTurbulence = 0;
    double groundLayerMean = 2.7;
    double groundLayerVariation = 0.3;
    std::ifstream inStream(sitestring);
    if (inStream) {
        readText pars2(sitestring);
        for (size_t t(0); t < pars2.getSize(); t++) {
            std::string line(pars2[t]);
            std::istringstream iss(line);
            std::string keyName;
            iss >> keyName;
            if (keyName == "turbulence") haveTurbulence = 1;
            readText::get(pars2[t], "groundlayermean", groundLayerMean);
            readText::get(pars2[t], "groundlayervar", groundLayerVariation);
        }
    }
    if (haveTurbulence == 0) setupStructureFunction(temperature, groundlevel, latitude, pressure, exosphereTemperature, longitude, waterPressure);
    if (constrainseeing < 0.0) {
        if (haveTurbulence == 1) ccalc(tai, 0);
        if (haveTurbulence == 0) alternativeStructureFunction(tai, 0, groundLayerMean, groundLayerVariation);
        totalj = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            totalj += jests[b];
        }
        totalseeing = 5.25/(pow(5000e-10, 0.2))/(2.0*sqrt(2.0*log(2.0)))
            /(PI/180/3600)*pow(totalj, 0.6);
        for (int b = 0; b < numlevel-1; b++) {
            relseeing[b] = pow(jests[b]/totalj, 0.6);
        }

        float renorm = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            renorm += relseeing[b]*relseeing[b];
        }
        renorm = sqrt(renorm);
        for (int b=0; b < numlevel-1; b++) {
            relseeing[b] = relseeing[b]*totalseeing/renorm;
        }
    } else if (constrainseeing > 0.0) {
        int bestindex = 0;
        float seeingdiff = 0.0, minseeingdiff = 0.0;
        std::vector<std::vector<float> > jestarray(100);
        float besttotalj = 0.0, totalj = 0.0;
        for (int c = 0; c < 100; c++) {
            // 100 different versions of cests are created
            jestarray[c].resize(numlevel-1, 0);
            if (haveTurbulence == 1) ccalc(tai, c);
            if (haveTurbulence == 0) alternativeStructureFunction(tai, c, groundLayerMean, groundLayerVariation);
            for (int b = 0; b < numlevel-1; b++){
                jestarray[c][b] = jests[b];
            }
            totalj = 0.0;
            for (int b = 0; b < numlevel-1; b++) {
                totalj += jestarray[c][b];
            }
            totalseeing = 5.25/(pow(5000e-10, 0.2))/(2.0*sqrt(2.0*log(2.0)))
                /(PI/180/3600)*pow(totalj, 0.6);
            // here the totalseeing value closest to constrainseeing is
            // found by comparing differences between the two
            seeingdiff = fabs(constrainseeing - totalseeing);
            if (c == 0 || seeingdiff < minseeingdiff) {
                minseeingdiff = seeingdiff;
                //bestseeing = totalseeing;
                bestindex = c;
                besttotalj = totalj;
            }
        }
        totalseeing = constrainseeing;
        // the totalseeing value that is closest to constrainseeing is
        // used to generate the relseeing at each level
        for (int b = 0; b < numlevel-1; b++){
            relseeing[b] = pow(jestarray[bestindex][b]/besttotalj, 0.6);
        }
        float renorm = 0.0;
        for (int b = 0; b < numlevel-1; b++) {
            renorm += relseeing[b]*relseeing[b];
        }
        renorm = sqrt(renorm);
        for (int b = 0; b < numlevel-1; b++) {
            relseeing[b] = relseeing[b]*constrainseeing/renorm;
        }
    } else {
        for (int c = 0; c < numlevel-1; c++) {
            relseeing[c] = 0.0;
            jests[c] = 0.0;
        }
        totalseeing = 0.0;
    }

    double cloudmean1=0.0, cloudsigma1=0.0, cloudmean2=0.0, cloudsigma2=0.0;
    double cloudvariationscale=10.0, domeseeingmean=0.05, domeseeingmedian=0.05;
    std::ifstream inStream2(sitestring);
    if (inStream2) {
        readText pars(sitestring);
        for (size_t t(0); t < pars.getSize(); t++) {
            readText::get(pars[t], "cloudmean1", cloudmean1);
            readText::get(pars[t], "cloudsigma1", cloudsigma1);
            readText::get(pars[t], "cloudmean2", cloudmean2);
            readText::get(pars[t], "cloudsigma2", cloudsigma2);
            readText::get(pars[t], "cloudvariationscale", cloudvariationscale);
            readText::get(pars[t], "domeseeingmean", domeseeingmean);
            readText::get(pars[t], "domeseeingmedian", domeseeingmedian);
        }
    }

    //domeseeing
    double domeseeing=0.0;
    if ((domeseeingmean!=0.0) && (domeseeingmedian != 0.0)) {
        float mu = log(domeseeingmedian);
        float sigma = pow(2.0*(log(domeseeingmean) - mu), 0.5);
        mu = log(domeseeingmean) - (1/2)*(pow(sigma, 2.0));
        domeseeing = exp(mu + sigma*random.normalCorrel(tai + 1234, 600.0/(24.0*3600.0)));
    }
    relseeing[numlevel -1]=domeseeing/2.35;

    FILE *afile;
    afile = fopen(outputfilename.c_str(), "wt");
    float rvalue[11];
    rvalue[0] = 1.0;
    rvalue[1] = 1.0;
    rvalue[2] = 1.0;
    rvalue[3] = 1.0;
    rvalue[4] = 1.0;
    rvalue[5] = 1.0;
    rvalue[6] = 1.0;
    rvalue[7] = 1.0;
    rvalue[8] = 1.0;
    rvalue[9] = 1.0;
    rvalue[10] = 1.0;
    fprintf(afile, "natmospherefile %i\n", numlevel);
    fprintf(afile, "totalseeing %f\n", totalseeing*(2.0*sqrt(2.0*log(2.0))));
    *tseeing=totalseeing*(2.0*sqrt(2.0*log(2.0)));
    fprintf(afile, "reldensity %f\n", rvalue[0]);
    fprintf(afile, "relo2 %f\n", rvalue[1]);
    fprintf(afile, "relh2o %f\n", rvalue[2]);
    fprintf(afile, "relo3 %f\n\n", rvalue[3]);
    fprintf(afile, "reln2o %f\n\n", rvalue[4]);
    fprintf(afile, "relco %f\n\n", rvalue[5]);
    fprintf(afile, "relch4 %f\n\n", rvalue[6]);
    rvalue[0] = 0.0004*random.normalCorrel(tai+98000, 1.0/24.0);
    rvalue[1] = 0.0004*random.normalCorrel(tai+98100, 1.0/24.0);
    rvalue[2] = 0.005*random.normalCorrel(tai+98200, 1.0/24.0);
    rvalue[3] = 0.002*random.normalCorrel(tai+98300, 1.0/24.0);
    rvalue[4] = 0.002*random.normalCorrel(tai+98400, 1.0/24.0);
    rvalue[5] = 0.002*random.normalCorrel(tai+98500, 1.0/24.0);
    rvalue[6] = 0.002*random.normalCorrel(tai+98600, 1.0/24.0);
    rvalue[7] = 0.002*random.normalCorrel(tai+98700, 1.0/24.0);
    rvalue[8] = 0.002*random.normalCorrel(tai+98800, 1.0/24.0);
    rvalue[9] = 0.002*random.normalCorrel(tai+98900, 1.0/24.0);
    rvalue[10] = 0.002*random.normalCorrel(tai+99000, 1.0/24.0);
    fprintf(afile, "raygradient %f\n\n", rvalue[0]);
    fprintf(afile, "o2gradient %f\n\n", rvalue[1]);
    fprintf(afile, "o3gradient %f\n\n", rvalue[2]);
    fprintf(afile, "h2ogradient %f\n\n", rvalue[3]);
    fprintf(afile, "n2ogradient %f\n\n", rvalue[4]);
    fprintf(afile, "cogradient %f\n\n", rvalue[5]);
    fprintf(afile, "ch4gradient %f\n\n", rvalue[6]);
    fprintf(afile, "aerosolseasaltgradient %f\n\n", rvalue[7]);
    fprintf(afile, "aerosoldustgradient %f\n\n", rvalue[8]);
    fprintf(afile, "aerosolpollutiongradient %f\n\n", rvalue[9]);
    fprintf(afile, "aerosolsmokegradient %f\n\n", rvalue[10]);
    fprintf(afile, "rayangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "o2angle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "o3angle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "h2oangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "n2oangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "coangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "ch4angle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "aerosolseasaltangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "aerosoldustangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "aerosolpollutionangle %f\n\n", random.uniform()*2.0*PI);
    fprintf(afile, "aerosolsmokeangle %f\n\n", random.uniform()*2.0*PI);
    double rr1 = random.normalCorrel(tai, 600.0/(24.0*3600.0));
    double rr2 = random.normalCorrel(tai + 1000, 600.0/(24.0*3600.0));
    fprintf(afile, "windrandom1 %lf\n", rr1);
    fprintf(afile, "windrandom2 %lf\n", rr2);

    //    double cloudCover = humidity*0.65 + 0.1*random.normalCorrel(tai + 12000, 2.0);
    //    double cloudCover = pow(humidity, 3.0) +  0.1*random.normalCorrel(tai + 12000, 2.0);
    double cloudCover = 0.65*exp(-(1.0-humidity)/0.3) +  0.1*random.normalCorrel(tai + 12000, 2.0);
    if (constraincloudcover >= 0.0) {
        cloudCover = constraincloudcover;
    }
    if (cloudCover > 1.0) cloudCover = 1.0;
    if (cloudCover < 0.0) cloudCover = 0.0;
    double cloudCover1 = 1.0 - sqrt(1.0 - cloudCover);
    double cloudCover2 = 1.0 - sqrt(1.0 - cloudCover);
    double sigma = 0.5;
    double cloudDepth = sigma*random.uniformCorrel(tai+15000, 8.0/24.0, 0);
    if (constrainclouddepth >= 0.0) {
        cloudDepth = constrainclouddepth;
    }
    double cloudDepth1 = cloudDepth /2.0;
    double cloudDepth2 = cloudDepth /2.0;
    double q1 = ((cloudCover1 - 0.5) + PI/24.0*pow(2*cloudCover1 - 1.0, 3.0)) + PI*PI*7.0/960.0*pow(2*cloudCover1 - 1.0, 5.0);
    double q2 = ((cloudCover2 - 0.5) + PI/24.0*pow(2*cloudCover2 - 1.0, 3.0)) + PI*PI*7.0/960.0*pow(2*cloudCover2 - 1.0, 5.0);
    cloudsigma1 = cloudDepth1/(sqrt(2*PI)*q1*cloudCover1 + 1.0/sqrt(2*PI)*exp(-PI*q1));
    cloudsigma2 = cloudDepth2/(sqrt(2*PI)*q2*cloudCover2 + 1.0/sqrt(2*PI)*exp(-PI*q2));
    cloudmean1 = sqrt(2*PI)*cloudsigma1*q1;
    cloudmean2 = sqrt(2*PI)*cloudsigma2*q2;

    // determine cloud heights
    // int maxcloud = 0;
    // double closest = 1e30;
    // for (int i = 0; i < numlevel - 2; i++) {
    //     double difference = abs(6000.0 - altitudes[i]);
    //     if (difference < closest) {
    //         closest = difference;
    //         maxcloud = i;
    //     }
    // }
    int mincloud = 0;
    int maxcloud = 0;
    double closest = 1e30;
    int cloudclose1 = 0;
    int cloudclose2 = 0;
    for (int i =0; i < numlevel - 2; i++) {
        double difference = abs(-8500.0*log(humidity) - altitudes[i]);
        if (difference < closest) {
            closest = difference;
            cloudclose1 = i;
        }
    }
    closest=1e30;
    for (int i =0; i < numlevel - 2; i++) {
        double difference = abs(-8500.0*log(humidity) - altitudes[i]);
        if ((difference < closest) && (i != cloudclose1)) {
            closest = difference;
            cloudclose2 = i;
        }
    }
    if (cloudclose1 < cloudclose2) {
        maxcloud = cloudclose1;
        mincloud = cloudclose2;
    } else {
        maxcloud = cloudclose2;
        mincloud = cloudclose1;
    }
    //    maxcloud = mincloud - 1;
    //    if (mincloud <= maxcloud) mincloud=maxcloud + 1;
    cloudscreen[mincloud]=1;
    cloudscreen[maxcloud]=1;

    for (int i = 0; i < numlevel; i++) {
        fprintf(afile, "height %d %f\n", i, (altitudes[i] - groundlevel)/1000.0);
        fprintf(afile, "outerscale %d %f\n", i, osests[i]);
        fprintf(afile, "seeing %d %f\n", i, relseeing[i]);

        if (cloudscreen[i]) {
            float cloudmean = 0.0, cloudvary = 0.0;
            if (i == maxcloud) cloudmean = cloudmean1;
            if (i == mincloud) cloudmean = cloudmean2;
            if (i == maxcloud) cloudvary = cloudsigma1;
            if (i == mincloud) cloudvary = cloudsigma2;
            fprintf(afile, "cloudmean %d %f\n\n", i, cloudmean);
            fprintf(afile, "cloudvary %d %f\n\n", i, cloudvary);
        }
    }
    fclose(afile);

}

void Atmosphere::ccalc(double tai, int call) {

    float data[100][11];
    int numTurbulence = 0;
    std::ifstream inStream(sitestring);
    if (inStream) {
        readText pars(sitestring);
        for (size_t t(0); t < pars.getSize(); t++) {
            std::string line(pars[t]);
            std::istringstream iss(line);
            std::string keyName;
            iss >> keyName;
            if (keyName == "turbulence") {
                for (int i = 0; i < 11; i++) {
                    iss >> data[numTurbulence][i];
                }
                numTurbulence++;
                if (numTurbulence > 100) {
                    printf("Error:  Too many turbulence lines\n");
                    exit(1);
                }
            }
        }
    }
    std::vector<float> j1grid(numTurbulence, 0);
    std::vector<float> j2grid(numTurbulence, 0);
    std::vector<float> jtgrid(numTurbulence, 0);
    std::vector<float> hlow(numlevel - 1, 0);
    std::vector<float> hhigh(numlevel - 1, 0);
    float r1, r2, xl, xh, overlap;

    r1 = random.normalCorrel(tai + 20000 + 1000*call, 600.0/(3600.0*24.0));
    r2 = random.normalCorrel(tai + 20000 + 1000*call + 500, 600.0/(3600.0*24.0));

    for (int i = 0; i < numTurbulence; i++) {


        if (data[i][4] != 0.0) {
            j1grid[i] = exp(log(data[i][4]) + data[i][6]*r1);
        } else {
            j1grid[i] = 0.0;
        }
        if (j1grid[i] < 0) {
            j1grid[i] = 0.0;
        }
        if (data[i][8] != 0.0) {
            j2grid[i] = exp(log(data[i][8]) + data[i][10]*r2);
        } else {
            j2grid[i] = 0.0;
        }
        if (j2grid[i] < 0) {
            j2grid[i] = 0.0;
        }
        jtgrid[i] = (j1grid[i] + j2grid[i])*1e-13;
    }
    for (int i = 0; i < numlevel-1; i++) {
        if (i > 0) {
            hhigh[i] = 0.5*(altitudes[i - 1] + altitudes[i]);
        }
        if (i < numlevel - 2) {
            hlow[i] = 0.5*(altitudes[i] + altitudes[i + 1]);
        }
    }
    hlow[numlevel - 2] = 0.0;
    hhigh[0] = 20000 + groundlevel;

    for (int i=0; i < numlevel-1; i++) {
        jests[i] = 0.0;
        for (int j = 0; j < numTurbulence; j++) {
            if (data[j][0] > hlow[i]) {
                xl = data[j][0];
            } else {
                xl = hlow[i];
            }
            if (data[j][2] < hhigh[i]) {
                xh = data[j][2];
            } else {
                xh = hhigh[i];
            }
            overlap = (xh - xl)/(data[j][2] - data[j][0]);
            if (data[j][2] < hlow[i]) {
                overlap = 0.0;
            }
            if (data[j][0] > hhigh[i]) {
                overlap = 0.0;
            }
            jests[i] = jests[i] + overlap*jtgrid[j];
        }
    }
}

float Atmosphere::outerscale(float altitude, double tai, int n) {
    /// @fn float AtmosphereCreator::outerscale(float altitude)
    //  @brief creates an outerscale estimate based on a lognormal distribution.

    float mean, median;
    mean = 38.0;
    median = 32.0;

    std::ifstream inStream(sitestring);
    if (inStream) {
        readText pars(sitestring);
        for (size_t t(0); t < pars.getSize(); t++) {
            readText::get(pars[t], "outerscalemean", mean);
            readText::get(pars[t], "outerscalemedian", median);
        }
    }

    float mu = log(median);
    float sigma = pow(2.0*(log(mean) - mu), 0.5);
    mu = log(mean) - (1/2)*(pow(sigma, 2.0));
    float osest = exp(mu + sigma*random.normalCorrel(tai + n*1000 + 2000, 600.0/(24.0*3600.0)));
    return osest;
}

void Atmosphere::setupStructureFunction(double temperature, double groundlevel, double latitude, double pressure, double exosphereTemperature, double longitude, double waterPressure) {

    air.readReferenceAtmosphere (atmospheredir, temperature + 4.65*groundlevel/1000.0, latitude, pressure/exp(-groundlevel/7640.0), waterPressure/exp(-groundlevel/7640), exosphereTemperature, longitude);
    air.potentialTemperatureProfile.resize(air.temperatureProfile.size());
    for (unsigned long i = 0 ; i < air.potentialTemperatureProfile.size(); i++) {
        //0.286 = R/cp/(mu*mp) generalize this later
        air.potentialTemperatureProfile[i] = air.temperatureProfile[i]*pow(air.pressureProfile[0]/air.pressureProfile[i], 0.286);
    }

    air.potentialTemperatureDerivative.resize(air.temperatureProfile.size());
    for (unsigned long i = 1;  i < air.potentialTemperatureDerivative.size() - 1; i++) {
        air.potentialTemperatureDerivative[i] = (air.potentialTemperatureProfile[i+1] - air.potentialTemperatureProfile[i-1])/
            (air.altitudeProfile[i + 1] - air.altitudeProfile[i - 1]);
    }
    air.potentialTemperatureDerivative[0] =  (air.potentialTemperatureProfile[1] - air.potentialTemperatureProfile[0])/
        (air.altitudeProfile[1] - air.altitudeProfile[0]);
    air.potentialTemperatureDerivative[air.potentialTemperatureDerivative.size() - 1] = (air.potentialTemperatureProfile[air.potentialTemperatureDerivative.size() - 1] - air.potentialTemperatureProfile[air.potentialTemperatureDerivative.size() - 2])/
        (air.altitudeProfile[air.potentialTemperatureDerivative.size() - 1] - air.altitudeProfile[air.potentialTemperatureDerivative.size() - 2]);

    air.structureFunction.resize(air.temperatureProfile.size());
    air.structureFunctionVar.resize(air.temperatureProfile.size());
    air.structureFunctionCurrent.resize(air.temperatureProfile.size());
    for (unsigned long i = 0 ; i < air.potentialTemperatureProfile.size(); i++) {

        air.structureFunction[i] = air.potentialTemperatureDerivative[i] * air.numberDensityProfile[i] * K_BOLTZMANN * 50.0/1e6;
        air.structureFunctionVar[i] = air.structureFunction[i] * air.numberDensityProfile[i] * 1e6/1e6 * 2e-24;
        // printf("XXX %e %e %e %e %e %e %e\n",air.altitudeProfile[i],air.structureFunction[i],air.structureFunctionVar[i],air.numberDensityProfile[i],air.temperatureProfile[i],air.potentialTemperatureProfile[i],air.potentialTemperatureDerivative[i]);
    }

}


void Atmosphere::alternativeStructureFunction(double tai, int call, double groundLayerMean, double groundLayerVariation) {

    // choose from log-normal for each layer with correlation
    //  add in something for ground layer

    double r1 = random.normalCorrel(tai + 20000 + 1000*call, 600.0/(3600.0*24.0));
    double r2 = random.normalCorrel(tai + 20000 + 1000*call + 500, 600.0/(3600.0*24.0));
    for (unsigned long i = 0; i < air.potentialTemperatureProfile.size(); i++) {
        air.structureFunctionCurrent[i] = 1e-13 * (exp(log(air.structureFunction[i]) + r1*air.structureFunctionVar[i]));
    }


    for (int j = 0; j < numlevel - 1; j++) {
        jests[j] = 0;
    }
    for (unsigned long i = 0; i < air.potentialTemperatureProfile.size(); i++) {
        if (air.altitudeProfile[i]*1000.0 >= groundlevel) {
            double distance = 1e30;
            int bestvalue = 0;
            int up = 0;
            for (int j = 0; j < numlevel - 1; j++) {
                double localdistance = abs(altitudes[j] - air.altitudeProfile[i]*1000.0);
                if (localdistance < distance) {
                    bestvalue = j;
                    distance = localdistance;
                    if (altitudes[j] < air.altitudeProfile[i]*1000.0) up =1; else up = 0;
                }
            }
            jests[bestvalue] += air.structureFunctionCurrent[i]*0.5;
            if (up==1) jests[bestvalue + 1] += air.structureFunctionCurrent[i]*0.5;
            if (up==0) jests[bestvalue] += air.structureFunctionCurrent[i]*0.5;
        }
    }

    jests[numlevel -2] = 1e-13 * (exp(log(groundLayerMean) + r2*groundLayerVariation));

}
