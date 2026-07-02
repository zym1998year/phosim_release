///
/// @package phosim
/// @file air.h
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

class Air {

 public:
    double *tauWavelength;
    double **tau;
    double *airmassLayer;
    double **tauMoon;
    double **tauMoonScatter;
    double **tauSunScatter;
    double *airmassLayerMoon;
    double air_refraction_adc;
    double topAtmosphere;



    void opacitySetup(double zenith, double moonalt, double solaralt, std::vector<double> height, double groundlevel, double raynorm, double o2norm, double h2onorm, double o3norm, double n2onorm, double conorm, double ch4norm, std::vector<double> aerosol, int layers, const std::string & dir, double *airmass, double pressure, double waterPresssure, double temperature, double latitude, double exosphereTemperature, double longitude, double minwavelength, double maxwavelength, std::vector<double> & cloudmean, std::vector<double> & cloudvary);

    void opticalDepth(double wavelength, double height, double cheight, double prheight, double *transmission, double *transmission2, double *transmission3, double *transmission4, double *transmission5, double *transmission6, double *transmission7, double *transmission8, double *transmission9, double *transmission10, double *transmission11, double *rayprofile,  std::vector<double> & o3cs, std::vector<double> & o3cs_wavelength, std::vector<std::vector<double> > & o2cs, std::vector<double> & o2cs_wavelength, std::vector<std::vector<double> > & h2ocs, std::vector<double> & h2ocs_wavelength, std::vector<double> & aerosol, int flag);

    std::vector<double> temperatureProfile;
    std::vector<double> altitudeProfile;
    std::vector<double> densityProfile;
    std::vector<double> numberDensityProfile;
    std::vector<double> pressureProfile;
    std::vector<double> molecularMassProfile;
    std::vector<double> windUMuProfile;
    std::vector<double> windVMuProfile;
    std::vector<double> windUSigmaProfile;
    std::vector<double> windVSigmaProfile;
    std::vector<double> potentialTemperatureProfile;
    std::vector<double> potentialTemperatureDerivative;
    std::vector<double> structureFunction;
    std::vector<double> structureFunctionVar;
    std::vector<double> structureFunctionCurrent;

    std::vector<double> h2oProfile;
    std::vector<double> co2Profile;
    std::vector<double> n2Profile;
    std::vector<double> o2Profile;
    std::vector<double> o3Profile;
   std::vector<double> n2oProfile;
   std::vector<double> coProfile;
   std::vector<double> ch4Profile;
    std::vector<double> oProfile;
    std::vector<double> arProfile;
    std::vector<double> heProfile;
    std::vector<double> hProfile;
    std::vector<double> noProfile;
    std::vector<double> so2Profile;
    std::vector<double> no2Profile;
    std::vector<double> nh3Profile;
    std::vector<double> hno3Profile;
    std::vector<double> ohProfile;
    std::vector<double> hfProfile;
    std::vector<double> hclProfile;
    std::vector<double> hbrProfile;
    std::vector<double> hiProfile;
    std::vector<double> cloProfile;
    std::vector<double> ocsProfile;
    std::vector<double> h2coProfile;
    std::vector<double> hoclProfile;
    std::vector<double> hcnProfile;
    std::vector<double> ch3clProfile;
    std::vector<double> h2o2Profile;
    std::vector<double> c2h2Profile;
    std::vector<double> c2h6Profile;
    std::vector<double> ph3Profile;
    std::vector<double> neProfile;
    std::vector<double> krProfile;


    void readReferenceAtmosphere(std::string dir, double groundTemperature, double latitude, double pressure,  double waterPressure, double exosphereTemperature, double longitude);

   void readWindPattern(std::string dir, double latitude, double longitude, int currentMonth);

    void setWindLayers(std::vector<double> height, double groundlevel, std::vector<double> & wind, std::vector<double> & winddir, int layers, double r1, double r2);

};
