///
/// @package phosim
/// @file silicon.h
/// @brief silicon header file
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

#include "ancillary/random.h"

class Silicon {

public:

    float *meanFreePath;
    float *sigma, *sigmaX, *sigmaY, *fsigma, *gsigma, *hsigma, *isigma;
    float *gammaX, *gammaY;
    float *deltaX, *deltaY;
    float *nbulkmap, *deadLayer, *periodmap;
    float *indexRefraction;
    double *indexRefractionCoatingReal;
    double *indexRefractionCoatingImag;
    float *temperatureGrid;
    float *rho;
    float *wavelengthGrid;
    float *thicknessGrid;
    float *dopantGrid;
    std::vector<double> nsilicon_real;
    std::vector<double> nsilicon_imag;
    std::vector<double> nsilicon_wavelength;
    std::vector<double> ncoating_real;
    std::vector<double> ncoating_imag;
    std::vector<double> ncoating_wavelength;
    double spaceChargeShield;
    double spaceChargeSpreadX;
    double spaceChargeSpreadY;
    double chargeStopCharge;
    double x0;
    double y0;
    int numWavelength;
    int numTemperature;
    int numRho;
    int numThickness;
    int numDopant;
    Random random;
    double fano;
    double pairenergy;
    double excitation;
    double atomicA;
    double atomicZ;
    double density;
    double maxWavelengthElectronPair;
    int xTransferFlag;
    int yTransferFlag;
    int midlineFlag;
    int inactiveX;
    int inactiveY;

    double absorptionCoeffMCT(double lambda, double temperature, double x);
    double absorptionCoeffSi(double lambda, double temperature);
    double indexRefractionMCT (double lambda, double temperature, double x);
    double indexRefractionSi(double lambda);
    double coatingIndexRefractionReal(double lambda);
    double coatingIndexRefractionImag(double lambda);
    void sensorSetup(std::string sensorfilename, std::string devmaterial, double ccdtemp, double nBulk, double nF, double nB, double sF,
               double sB, double tSi, double overdepBias, std::string instrdir, std::string mdir,
               long nampx, long nampy, double pixsize, long seedchip,
               float *impurityX, float *impurityY, int impurityvariation,
               double minwavelength, double maxwavelength, long activeBuffer,
               int dopingflag, int esatflag, int vsatflag, int detectorMode, int siliconDebug, long siliconSubSteps, std::string obshistid);
    double dopeProfile(double z, double nbulk, double nf, double nb, double sf, double sb, double tsi);
    double muMCT(double x, double temperature);
    double muSi (double efield, double temperature, int polarity, int esatflag, int vsatflag);
    double epsilonMCT(double x);

};
