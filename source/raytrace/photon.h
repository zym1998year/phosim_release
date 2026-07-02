///
/// @package phosim
/// @file photon.h
/// @brief header for photon class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

struct Photon {
    int direction;
    int polarization;
    double shiftedAngle;
    double airRefraction;
    double airRefractionPrevious;
    double airRefractionADC;
    double airRefractionPreviousADC;
    double adcx, adcy, adcz;
    double nprev, ncurr;
    double wavelength;
    double wavelengthFactor;
    int indexlx0, indexly0, indexlx1, indexly1;
    int indexcx0, indexcy0, indexcx1, indexcy1;
    int indexmx0, indexmy0, indexmx1, indexmy1;
    int indexfx0, indexfy0, indexfx1, indexfy1;
    float dlx, dly, dcx, dcy, dmx, dmy, dfx, dfy;
    int uuint, vvint, wwint;
    double windx, windy;
    int lindex;
    double xp, yp;
    double xporig, yporig, zporig;
    double opdx, opdy;
    TwoVectorInt pos;
    double xPosR, yPosR;
    double xpos, ypos;
    float time, prtime;
    float absoluteTime;
    double dvr;
    double op;
    int oindex;
    long counter;
    long maxcounter;
    int ghostFlag;
    int preGhostFlag;
    long sourceOver_m;
    int sourceSaturationRadius;
    int sourceSaturationPreviousRadius;
    float sourceSaturationRadiusCorrected;
    float sourceSaturationFraction1;
    float sourceSaturationFraction2;
    float sourceSaturationFraction3;
    float sourceSaturationFraction4;
    float sourceSaturationRadius1;
    float sourceSaturationRadius2;
    float sourceSaturationRadius3;
    float sourceSaturationRadius4;
    float sourceSaturationIndex;
    int sourceSaturationFirst;
    int saturationFlag;
    double z0;
    double collect_z;
    int xindex;
    double rxindex;
    int dopantTemperatureThicknessIndex;
    int dopantThicknessIndex;
    int dopantRhoThicknessIndex;
    int zindex;
    double sensorDirection;
    int location;
    int location2;
    float saveRand[MAX_BOUNCE];
    double refractionInvariant;
    double refractionInvariantCenter;
    int thread;
    float transmissionSave[5];
    float transmissionCheckRand;
    double sat_xpos;
    double sat_ypos;
    double sat_n;
    double sat_xx;
    double sat_yy;
    double sat_sigma;
    double sat_maxr;
    double sat_r;
};
