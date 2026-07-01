///
/// @package phosim
/// @file atmosphere.h
/// @brief header for atmosphere class
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
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <mutex>

#include "raytrace/helpers.h"
#include "ancillary/random.h"
#include "ancillary/fft.h"
#include "ancillary/fits.h"
#include "ancillary/readtext.h"
#include "raytrace/constants.h"
#include "raytrace/air.h"

class Atmosphere {

 public:

    int numlevel;
    float groundlevel;
    double latitude;
    double longitude;
    double constrainclouddepth;
    double constraincloudcover;
    std::string datadir;
    std::string instrdir;
    std::string atmospheredir;
    std::string sitestring;
    std::vector<float> osests;
    std::vector<float> altitudes;
    std::vector<float> jests;
    Random random;
    std::mutex lock;
    Air air;

    void createAtmosphere(float monthnum, float constrainseeing, const std::string & outputfilename, std::vector<int> & cloudscreen, long seed, double tai, double *tseeing, double latitude, double oceanDistance, double humidity, double temperature, double pressure, double groundlevel, double longitude, double waterPressure, double exosphereTemperature);

    void turb2d(long seed, double see5, double outerx, double outers,
                double zenith, double wavelength, char *name,
                long N_size=1024, double resample=1.0);

    void cloud(long seed, double cloheight, double pixsz,
               const std::string & name, long N_size = 1024);

    void airglow(long seed, const std::string & name, long screenSize=1024);

    void ccalc(double tai, int call);
    void alternativeStructureFunction(double tai, int call, double groundLayerMean, double groundLayerVariation);
    void setupStructureFunction(double temperature, double groundlevel, double latitude, double pressure, double exosphereTemperature, double longitude, double waterPressure);

    float outerscale(float altitude, double tai, int n);
    static void* threadFunction(void *voidArgs);

};


struct thread_args {
    Atmosphere *instance;
    long seed;
    double see5;
    double outerx;
    double outers;
    double zenith;
    double wavelength;
    char name[4096];
    long N_size;
    double resample;
};
