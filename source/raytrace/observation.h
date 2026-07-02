///
/// @package phosim
/// @file observation.h
/// @brief observation header file
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

#include <string>
#include <vector>
#include "parameters.h"
#include "vector_type.h"
#include "source.h"

class Observation {

public:
    double rotationjitter;
    double zenith;
    double altitude;
    double azimuth;
    double rotationrate;
    double elevationjitter;
    double azimuthjitter;
    double windjitter;
    double windshake;
    double groundlevel;
    double xtelloc;
    double ytelloc;
    double latitude;
    double longitude;
    std::vector<double> seefactor;
    std::vector<double> wind;
    std::vector<double> winddir;
    std::vector<int> windset;
    double windrandom1;
    double windrandom2;
    std::vector<double> outerscale;
    std::vector<double> height;
    std::vector<double> aerosol;
    std::vector<double> cloudmean;
    std::vector<double> cloudvary;
    std::vector<std::string> extraCommandString;
    double *dtau;
    double pressure;
    double waterPressure;
    double temperature;
    double tempvar;
    double raynorm;
    double o2norm;
    double o3norm;
    double h2onorm;
    double n2onorm;
    double conorm;
    double ch4norm;
    double raygradient;
    double o2gradient;
    double o3gradient;
    double h2ogradient;
    double cogradient;
    double n2ogradient;
    double ch4gradient;
    double aerosolSeasaltgradient;
    double aerosolDustgradient;
    double aerosolPollutiongradient;
    double aerosolSmokegradient;
    double rayAngle;
    double o2Angle;
    double o3Angle;
    double h2oAngle;
    double n2oAngle;
    double coAngle;
    double ch4Angle;
    double aerosolSeasaltAngle;
    double aerosolDustAngle;
    double aerosolPollutionAngle;
    double aerosolSmokeAngle;
    std::vector<std::string> pertType;
    int NTERM;
    std::vector<std::vector<double> > izernike;
    std::vector<int> surfaceLink;
    std::vector<std::vector<double> > body;
    double minr;
    int setApertureFlag;
    double shuttererror;
    double zodiacalLightMagnitude;
    double overrideZodiacalLightMagnitude;
    double artificialLightMagnitude;
    double bortle;
    float transtol;
    double backAlpha, backBeta, backRadius, backGamma, backDelta,
        backBuffer, backEpsilon, np, finiteDistance;
    int satbuffer;
    int satbuffermin;
    int satindex;
    int activeBuffer;
    double screentol;
    double maxr;
    double domeseeing;
    double toypsf;
    double toyg1;
    double toyg2;
    double exosphereTemperature;
    double pixsize;
    double pra;
    double pdec;
    double spiderangle;
    double platescale;
    double centerx;
    double centery;
    double tai;
    double day;
    double exptime;
    double vistime;
    double timeoffset;
    double rotatex;
    double rotatey;
    double rotatez;
    double standardwavelength;
    int flipX;
    int flipY;
    int swap;
    std::vector<double> sedW;
    std::vector<double> sedC;
    double largeScale;
    double coarseScale;
    double mediumScale;
    double fineScale;
    double totalnorm;
    double totalseeing;
    double moonalt;
    double moonaz;
    double moondist;
    double moonra;
    double moondec;
    double phaseang;
    double solarzen;
    double solaralt;
    double solaraz;
    double solardist;
    double solarra;
    double solardec;
    double solarfactor;
    double lst;
    double hourangle;
    double airglowcintensity;
    double airglowpintensity;
    double watervar;
    double minwavelength;
    double maxwavelength;
    double centralwavelength;
    double domelight;
    double domewave;
    double chipangle;
    double decenterx;
    double decentery;
    int telconfig;
    int checkpointcount;
    int checkpointtotal;
    int impurityvariation;
    int fieldanisotropy;
    int fringeflag;
    int detectorcollimate;
    int deadlayer;
    int chargesharing;
    int pixelerror;
    int chargediffusion;
    int photoelectric;
    int sensorrefraction;
    int dopingflag;
    int esatflag;
    int vsatflag;
    int airrefraction;
    int aberration;
    int nutation;
    int precession;
    int onlyadjacent;
    int splitsources;
    int focusSteps;
    int focusFlag;
    int backgroundTypes;
    int backgroundProbability;
    int siliconDebug;
    int rotationTracking;
    int rotationPoints;
    double bpVxMin;
    double bpVxMax;
    double bpVyMin;
    double bpVyMax;
    int bpNx;
    int bpNy;
    double bpStep;
    double bpBeta;
    double bpIndex;
    double raydensity;
    double devvalue;
    double airglowvariation;
    double focusStepZ;
    double focusStepX;
    float nbulk;
    float nf;
    float nb;
    float sf;
    float sb;
    float sensorthickness;
    float impurityX;
    float impurityY;
    float overdepbias;
    float siliconreflectivity;
    float sensorTempNominal;
    float sensorTempDelta;
    float qevariation;
    float *airglow;
    long nphot;
    int microroughness;
    int airglowScreenSize;
    int telescopeMode;
    int backgroundMode;
    int coatingmode;
    int poissonMode;
    int contaminationmode;
    int trackingMode;
    int windShakeMode;
    int ranseed;
    int obsseed;
    int zernikemode;
    int atmospheric_dispersion;
    int atmosphericdispcenter;
    int natmospherefile;
    int straylight;
    double straylightcut;
    int detectorMode;
    int opticsonlymode;
    int additivemode;
    int diffractionMode;
    int spiderMode;
    int pupilscreenMode;
    int spaceMode;
    int sequenceMode;
    int aperturemode;
    int areaExposureOverride;
    int filter;
    std::string filterName;
    int saturation;
    int eventfile;
    int diagnostic;
    int opdfile;
    int opdsize;
    int opdsampling;
    int centroidfile;
    int throughputfile;
    int pixelsx, pixelsy, minx, maxx, miny, maxy;
    int minxDrift, maxxDrift, minyDrift, maxyDrift;
    float driftX, driftY;
    double halfpixelsx, halfpixelsy, inversePixsize;
    std::string obshistid;
    int pairid;
    int blooming;
    int siliconSubSteps;
    int siliconSegments;
    unsigned long well_depth;
    int nsedptr;
    int sedptr;
    int nsource;
    int nsourcerotation;
    int nimage;
    int ghostonly;
    int atmdebug;
    int largeGrid;
    int coarseGrid;
    int mediumGrid;
    int fineGrid;
    int nsnap;
    int nframes;
    int nskip;
    int ngroups; // not exposed
    std::vector<int> ghost;
    int flatdir;
    int tarfile;
    int date;
    int numthread;
    int sunsafety;
    std::string devmaterial;
    std::string devtype;
    std::string devmode;
    std::vector<int> feaflag;
    std::vector<std::string> feafile;
    std::vector<double> feascaling;
    double normwave;
    int atmospheremode;
    int opacitymode;
    int sourceperthread;

    std::vector<long> sourceXpos;
    std::vector<long> sourceYpos;
    std::vector<long> sourcePhoton;
    std::vector<int> sedN;
    std::vector<int> sedPtr;
    std::vector<int> sedMode;
    std::vector<double> sedCorr;
    std::vector<double> sedDwdp;

    // should be part of image but have conflict with settings.c
    int nsurf;
    int nmirror;
    double airmass;

    // remainder should be ok
    std::vector<std::string> atmospherefile;
    std::vector<std::string> cloudfile;
    std::string trackingfile;
    int trackinglines;
    std::string outputfilename;
    std::string sensorfilename;
    std::string chipid;
    std::string focalplanefile;
    std::string sensorCoatingFile;
    int sensorCoatingFileFlag;
    std::string eventFitsFileName;
    std::string outputdir, workdir, seddir, imagedir, datadir, instrdir, bindir;
    Vector tpx, tpy, tpz;

    int naxesb[MAX_IMAGE][2];
    float *tempptr[MAX_IMAGE];
    double *cumulativex[MAX_IMAGE];
    double cumulative[MAX_IMAGE];
    double cumulativef[MAX_IMAGE];

    Source sources;
    int nspatial = 0;
    int nsed = 0;
    std::vector<std::string> spatiallist;
    std::vector<std::string> sedlist;

    // functions
    int parser();
    int background();
    int addSource(const std::string & object, int sourcetype);
    int addOpd(const std::string & opd);
    int header(fitsfile *faptr);
    int settings();
    int filterTruncateSources();
    int splitSources();
    void readSed(const std::string & filename, int mode);

};
