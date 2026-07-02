///
/// @package phosim
/// @file mie.h
/// @brief header file for air class
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

#include <complex>

class Mie {

 public:

    std::vector<double> refractionWavelengthWater;
    std::vector<double> refractionRealWater;
    std::vector<double> refractionImagWater;
    double *muValueWater;
    double *cumulativeMuWater;
    double *muWavelengthWater;
    double *muXWater;
    double *absorptionWater;
    double *scatteringWater;
    double *totalExtinctionWater;
    double *totalProfileWater;
    double *tauArrayWater;
    double dropMaxWater;
    double dropMinWater;
    double dropMeanWater;
    double dropSigmaWater;
    double qmaxWater;
    double *effsWater;
    std::vector<double> refractionWavelengthSeasalt;
    std::vector<double> refractionRealSeasalt;
    std::vector<double> refractionImagSeasalt;
    double *muValueSeasalt;
    double *cumulativeMuSeasalt;
    double *muWavelengthSeasalt;
    double *muXSeasalt;
    double *absorptionSeasalt;
    double *scatteringSeasalt;
    double *totalExtinctionSeasalt;
    double *totalProfileSeasalt;
    double *tauArraySeasalt;
    double dropMaxSeasalt;
    double dropMinSeasalt;
    double dropMeanSeasalt;
    double dropSigmaSeasalt;
    double qmaxSeasalt;
    double *effsSeasalt;
    std::vector<double> refractionWavelengthDust;
    std::vector<double> refractionRealDust;
    std::vector<double> refractionImagDust;
    double *muValueDust;
    double *cumulativeMuDust;
    double *muWavelengthDust;
    double *muXDust;
    double *absorptionDust;
    double *scatteringDust;
    double *totalExtinctionDust;
    double *totalProfileDust;
    double *tauArrayDust;
    double dropMaxDust;
    double dropMinDust;
    double dropMeanDust;
    double dropSigmaDust;
    double qmaxDust;
    double *effsDust;
    std::vector<double> refractionWavelengthPollution;
    std::vector<double> refractionRealPollution;
    std::vector<double> refractionImagPollution;
    double *muValuePollution;
    double *cumulativeMuPollution;
    double *muWavelengthPollution;
    double *muXPollution;
    double *absorptionPollution;
    double *scatteringPollution;
    double *totalExtinctionPollution;
    double *totalProfilePollution;
    double *tauArrayPollution;
    double dropMaxPollution;
    double dropMinPollution;
    double dropMeanPollution;
    double dropSigmaPollution;
    double qmaxPollution;
    double *effsPollution;
    std::vector<double> refractionWavelengthSmoke;
    std::vector<double> refractionRealSmoke;
    std::vector<double> refractionImagSmoke;
    double *muValueSmoke;
    double *cumulativeMuSmoke;
    double *muWavelengthSmoke;
    double *muXSmoke;
    double *absorptionSmoke;
    double *scatteringSmoke;
    double *totalExtinctionSmoke;
    double *totalProfileSmoke;
    double *tauArraySmoke;
    double dropMaxSmoke;
    double dropMinSmoke;
    double dropMeanSmoke;
    double dropSigmaSmoke;
    double qmaxSmoke;
    double *effsSmoke;


    void mieTyndall(double radius, double wavelength, double re, double im, double *extinction, double *scattering, int angles, std::vector<double> & phi, std::vector<double> & psf, std::vector<std::complex<double> > & s1, std::vector<std::complex<double> > & s2,
                     std::vector<double> & mu, std::vector<double> & piN,
                     std::vector<double> & piNm2,std::vector<double> & piNm1,
                    std::vector<double> & tauN);
    void setupMieGrid(int elementsWavelength, int elementsX, double minWavelength, double maxWavelength, double minX, double maxX, double *muValue, double *cumulativeMu, int elementsMu, double *absorption, double *scattering, std::vector<double> & refractionWavelength, std::vector<double> & refractionReal, std::vector<double> & refractionImag, double *muWavelength, double *muX, double *qmax, std::vector<double> & phi, std::vector<double> & psf,
                     std::vector<std::complex<double> > & s1, std::vector<std::complex<double> > & s2,
                     std::vector<double> & mu, std::vector<double> & piN,
                     std::vector<double> & piNm2,std::vector<double> & piNm1,
                     std::vector<double> & tauN);

    void setupMieIndexRefraction (std::string dir, double minwavelength, double maxwavelength,double meanWater, double sigmaWater,double meanSeasalt, double sigmaSeasalt,double meanDust, double sigmaDust,double meanPollution, double sigmaPollution,double meanSmoke, double sigmaSmoke);
    void setupMieScale (double dropmin, double dropmax, double dropmean, double dropsigma, int elementsTau, int elementsWavelength, double *effs, double *tauArray, double *muWavelength, double *muX, double *absorption, double *scattering, double *totalExtinction, double *psfCumulative, double *totalProfile);
};
