///
/// @package phosim
/// @file air.cpp
/// @brief air class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Emily Grace (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author Michael Wood-Vasey (Pitt)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <math.h>
#include <vector>

#include "ancillary/readtext.h"
#include "air.h"
#include "constants.h"
#include "helpers.h"
#include "parameters.h"


void Air::opticalDepth(double wavelength, double height, double cheight, double prheight, double *transmission, double *transmission2, double *transmission3, double *transmission4, double *transmission5, double *transmission6, double *transmission7, double *transmission8, double *transmission9, double *transmission10, double *transmission11,  double *rayprofile,  std::vector<double> & o3cs, std::vector<double> & o3cs_wavelength, std::vector<std::vector<double> > & o2cs, std::vector<double>  & o2cs_wavelength, std::vector<std::vector<double> > & h2ocs, std::vector<double> & h2ocs_wavelength, std::vector<double> & aerosol, int flag) {

    int index, index2, index3;
    double rindex;

    *transmission = 0.0;
    *transmission2 = 0.0;
    *transmission3 = 0.0;
    *transmission4 = 0.0;
    *transmission5 = 0.0;
    *transmission6 = 0.0;
    *transmission7 = 0.0;
    *transmission8 = 0.0;
    *transmission9 = 0.0;
    *transmission10 = 0.0;
    *transmission11 = 0.0;

            // rayleigh scattering
            find_v(altitudeProfile, height, &index3);
            index = find_linear(tauWavelength, AIR_OPACITY_POINT, wavelength, &rindex);

            double temperature = temperatureProfile[index3];
            double pressure = pressureProfile[index3]/PRESS_ATMOSPHERE*PRESS_ATMOSPHERE_MMHG;
            double waterPressure = h2oProfile[index3]*pressure;
            double sigma=1.0/wavelength;
            double ps=pressure/760.00*1013.25;
            double temp=(273.15+temperature);
            double pw=waterPressure/760.00*1013.25;
            double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
            double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
            double n=(2371.34+683939.7/(130.0-pow(sigma,2.0))+4547.3/(38.9-pow(sigma,2.0)))*ds;
            n=n+(6478.31+58.058*pow(sigma,2.0)-0.71150*pow(sigma,4.0)+0.08851*pow(sigma,6.0))*dw;
            double dn = 1e-8*n;
            n = dn + 1.0;
            double dens=pressure*133.32/(1.38e-23*temp);
            double cross=24.0*PI*PI*PI/pow(wavelength*1e-6,4.0)*1.016/dens/dens*pow((n*n-1.0),2.0)/pow((n*n+2.0),2.0);
            *transmission = cross*1e4*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
            //            printf("%lf %lf %lf %lf %e\n",wavelength,*transmission,height,altitudeProfile[index3],numberDensityProfile[index3]);

            //*transmission = rayprofile[index]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;

            // O3 absorption
            index = find_linear_v(o3cs_wavelength, wavelength, &rindex);
            *transmission4 = o3cs[index]*o3Profile[index3]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;

            // O2 absorption
            index = find_linear_v(o2cs_wavelength, wavelength, &rindex);
            index2 = ((int)((height)/7.0));
            if (index2 <= 8) {
                if (index2 < 8) {
                    double tempf1 = ((height) - index2*7.0)/7.0;
                    *transmission2 = (1 - tempf1)*o2cs[index][index2]*o2Profile[index3]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                    *transmission2 += tempf1*o2cs[index][index2 + 1]*o2Profile[index3]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    *transmission3 = h2oProfile[index3]*(1 - tempf1)*
                        h2ocs[index][index2]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                    *transmission3 += h2oProfile[index3]*tempf1*
                        h2ocs[index][index2 + 1]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                } else {
                    *transmission2 = o2cs[index][index2]*o2Profile[index3]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                    // H2O absorption
                    index = find_linear_v(h2ocs_wavelength, wavelength, &rindex);
                    *transmission3 = h2oProfile[index3]*h2ocs[index][index2]*numberDensityProfile[index3]/1e6*fabs(prheight - cheight)*1e5;
                }
            }

            // 5-7 nothing currently
            *transmission5 = 0.0;
            *transmission6 = 0.0;
            *transmission7 = 0.0;

            //*transmission8 += (aerosol[0])*pow(wavelength/0.5, aerosol[4] + aerosol[8]*log(wavelength/0.5))*
            //   exp(-height/2.20)/2.20*fabs(prheight - cheight);
            //*transmission9 += (aerosol[1])*pow(wavelength/0.5, aerosol[5] + aerosol[9]*log(wavelength/0.5))*
            //    exp(-height/3.60)/3.60*fabs(prheight - cheight);
            //*transmission10 += (aerosol[2])*pow(wavelength/0.5, aerosol[6] + aerosol[10]*log(wavelength/0.5))*
            //    exp(-height/6.60)/6.60*fabs(prheight - cheight);
            //*transmission11 += (aerosol[3])*pow(wavelength/0.5, aerosol[7] + aerosol[11]*log(wavelength/0.5))*
            //   exp(-height/1.70)/1.70*fabs(prheight - cheight);

            *transmission8 += (aerosol[0])*exp(-height/2.20)/2.20*fabs(prheight - cheight);
            *transmission9 += (aerosol[1])*exp(-height/3.60)/3.60*fabs(prheight - cheight);
            *transmission10 += (aerosol[2])*exp(-height/6.60)/6.60*fabs(prheight - cheight);
            *transmission11 += (aerosol[3])*exp(-height/1.70)/1.70*fabs(prheight - cheight);

            // 1)   rayleigh
            // 2)   non-variable molecules:  of the above only o2
            // 3)   water
            // 4)   o3
            // 5)   n2o:   probably not much until 2 microns
            // 6)   co:   unknown
            // 7)   methane:  probably not much until 2 microns
            // 8)   aerosol:  seasalt
            // 9) aerosol:  dust
            // 10) aerosol:  smoke
            // 11) aerosol:  pollution
            // 12)   clouds: mie scattering
}


void Air::opacitySetup(double zenith, double moonalt, double solaralt, std::vector<double> height, double groundlevel, double raynorm, double o2norm, double h2onorm,
                       double o3norm, double n2onorm, double conorm, double ch4norm, std::vector<double> aerosol, int layers, const std::string & dir, double *airmass, double pressure, double waterPressure, double temperature, double latitude, double exosphereTemperature, double longitude, double minwavelength, double maxwavelength, std::vector<double> & cloudmean, std::vector<double> & cloudvary) {

    double *rayprofile;
    std::vector<double> o3cs;
    std::vector<double> o3cs_wavelength;
    std::vector<std::vector<double> > o2cs;
    std::vector<double> o2cs_wavelength;
    std::vector<std::vector<double> > h2ocs;
    std::vector<double> h2ocs_wavelength;
    std::vector<double> airmass_height;
    std::vector<double> airmass_heightMoon;
    std::vector<std::vector<double> > airmass_heightScatterMoon;
    std::vector<std::vector<double> > airmass_heightScatterSun;
    std::vector<double> cloud;
    std::vector<double> cloudSmall;

    double chi, minchi, cheight, cheightnext, prheight, cheightprime, prheightprime;
    int closestlayer;
    double transmission, transmission2, transmission3, transmission4, transmission5,transmission6;
    double ttransmission, ttransmission2, ttransmission3, ttransmission4, ttransmission5, ttransmission6;
    double transmissionMoon, transmissionMoon2, transmissionMoon3, transmissionMoon4, transmissionMoon5, transmissionMoon6;
    double transmission7, transmission8, transmission9, transmission10, transmission11,transmission12;
    double ttransmission7, ttransmission8, ttransmission9, ttransmission10, ttransmission11, ttransmission12;
    double transmissionMoon7, transmissionMoon8, transmissionMoon9, transmissionMoon10, transmissionMoon11, transmissionMoon12;
    double wavelength;
    double moonZenith;
    double sunZenith;

    int points = AIR_OPACITY_POINT;

    topAtmosphere = 2000.0;
    moonZenith = PI/2.0 - moonalt;
    sunZenith = PI/2.0 - solaralt;

    readText::readCol(dir + "/o3cs.txt", o3cs_wavelength, o3cs);
    readText::readMultiCol(dir + "/o2cs.txt", o2cs_wavelength, o2cs);
    readText::readMultiCol(dir + "/h2ocs.txt", h2ocs_wavelength, h2ocs);

    // correct for ground level

    readReferenceAtmosphere (dir, temperature + 4.65*groundlevel/1000.0, latitude, pressure/exp(-groundlevel/7640.0), waterPressure/exp(-groundlevel/7640)/exp(-groundlevel/2200), exosphereTemperature,longitude);
    airmassLayer = new double[layers + 1]();
    airmassLayerMoon = new double[layers + 1]();
    tauWavelength = new double[points]();
    rayprofile = new double[points]();
    tau = new double*[12*(layers + 1)]();
    tauMoon = new double*[12*(layers + 1)]();
    tauMoonScatter = new double*[12*(layers + 1)]();
    tauSunScatter = new double*[12*(layers + 1)]();
    for (int i = 0; i < (12*(layers + 1)); i++) {
        tau[i] = new double[points]();
        tauMoon[i] = new double[points]();
        tauMoonScatter[i] = new double[points]();
        tauSunScatter[i] = new double[points]();
    }
    for (int i = 0; i < points; i++) {
        tauWavelength[i] = minwavelength/1000.0 + (maxwavelength/1000.0-minwavelength/1000.0)*((static_cast<double>(i))/(static_cast<double>(points)-1.0));
    }

    for (int i =0; i < points; i++) {
        double sigma=1.0/tauWavelength[i];
        double ps=pressure/760.00*1013.25;
        double temp=(273.15+temperature);
        double pw=waterPressure/760.00*1013.25;
        double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
        double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
        double n=(2371.34+683939.7/(130.0-pow(sigma,2.0))+4547.3/(38.9-pow(sigma,2.0)))*ds;
        n=n+(6478.31+58.058*pow(sigma,2.0)-0.71150*pow(sigma,4.0)+0.08851*pow(sigma,6.0))*dw;
        double dn = 1e-8*n;
        n = dn + 1.0;
        double dens=pressure*133.32/(1.38e-23*temp);
        double cross=24.0*PI*PI*PI/pow(tauWavelength[i]*1e-6,4.0)*1.016/dens/dens*pow((n*n-1.0),2.0)/pow((n*n+2.0),2.0);
        rayprofile[i]=cross*1e4;
    }
    prheight = 2000.0; // km

    // Airmass calculations

    // approximate airmass (for fits header)
    *airmass = 1.0/ (cos(zenith) + 0.025*exp(-11.0*cos(zenith)));

    airmass_height.resize(AIR_HEIGHT, 0);
    airmass_heightMoon.resize(AIR_HEIGHT, 0);

    // cloud opacity
    cloud.resize(AIR_HEIGHT, 0);
    cloudSmall.resize(AIR_HEIGHT_SCATTER, 0);
    for (int l = 0; l < layers; l++) {
        minchi=1e30;
        int bestk = 0;
        for (int k = 0; k < AIR_HEIGHT; k++) {
            cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT-1));
            chi = fabs(log10(height[l] + groundlevel/1000.0) - log10(cheight));
            if (chi < minchi) {
                minchi = chi;
                bestk = k;
            }
        }
        double v = cloudmean[l]/2.0*(1-erf(-cloudmean[l]/sqrt(2.0)/cloudvary[l])) +
            1/sqrt(2*PI)*cloudvary[l]*exp(-cloudmean[l]*cloudmean[l]/2.0/cloudvary[l]/cloudvary[l]);
        if (cloudvary[l] == 0) v=cloudmean[l];
        if (v >= 0) cloud[bestk] =  -log(pow(10.0, -0.4*v));
    }

    for (int l = 0; l < layers; l++) {
        minchi=1e30;
        int bestk = 0;
        for (int k = 0; k < AIR_HEIGHT_SCATTER; k++) {
            cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT_SCATTER-1));
            chi = fabs(log10(height[l] + groundlevel/1000.0) - log10(cheight));
            if (chi < minchi) {
                minchi = chi;
                bestk = k;
            }
        }
        double v = cloudmean[l]/2.0*(1-erf(-cloudmean[l]/sqrt(2.0)/cloudvary[l])) +
            1/sqrt(2*PI)*cloudvary[l]*exp(-cloudmean[l]*cloudmean[l]/2.0/cloudvary[l]/cloudvary[l]);
        if (cloudvary[l] == 0) v=cloudmean[l];
        if (v >= 0) cloudSmall[bestk] =  -log(pow(10.0, -0.4*v));
    }

    airmass_heightScatterMoon.resize(AIR_HEIGHT_SCATTER);
    for (int i = 0; i < AIR_HEIGHT_SCATTER; i++) {
        airmass_heightScatterMoon[i].resize(AIR_HEIGHT_SCATTER);
    }
    for (int k = 0; k < AIR_HEIGHT_SCATTER; k++) {
        cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT_SCATTER-1));
        for (int l = 0; l < AIR_HEIGHT_SCATTER; l++) {
            double cheightprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l))/(AIR_HEIGHT_SCATTER-1));
            double cheightnextprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l) + 1)/(AIR_HEIGHT_SCATTER-1));
            airmass_heightScatterMoon[k][l] = (sqrt(pow(cos(moonZenith), 2.0)*pow(RADIUS_EARTH + cheight, 2.0) + 2*(cheightprime - cheight)*(RADIUS_EARTH + cheight) + pow(cheightprime - cheight, 2.0)) -
                                               sqrt(pow(cos(moonZenith), 2.0)*pow(RADIUS_EARTH + cheight, 2.0) + 2*(cheightnextprime - cheight)*(RADIUS_EARTH + cheight) + pow(cheightnextprime - cheight, 2.0)))
                /(cheightprime - cheightnextprime);
        }
    }
    airmass_heightScatterSun.resize(AIR_HEIGHT_SCATTER);
    for (int i = 0; i < AIR_HEIGHT_SCATTER; i++) {
        airmass_heightScatterSun[i].resize(AIR_HEIGHT_SCATTER);
    }
    for (int k = 0; k < AIR_HEIGHT_SCATTER; k++) {
        cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT_SCATTER-1));
        for (int l = 0; l < AIR_HEIGHT_SCATTER; l++) {
            double cheightprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l))/(AIR_HEIGHT_SCATTER-1));
            double cheightnextprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l) + 1)/(AIR_HEIGHT_SCATTER-1));
            airmass_heightScatterSun[k][l] = (sqrt(pow(cos(sunZenith), 2.0)*pow(RADIUS_EARTH + cheight, 2.0) + 2*(cheightprime - cheight)*(RADIUS_EARTH + cheight) + pow(cheightprime - cheight, 2.0)) -
                                              sqrt(pow(cos(sunZenith), 2.0)*pow(RADIUS_EARTH + cheight, 2.0) + 2*(cheightnextprime - cheight)*(RADIUS_EARTH + cheight) + pow(cheightnextprime- cheight, 2.0)))
                /(cheightprime - cheightnextprime);
        }
    }

    for (int k = 0; k < AIR_HEIGHT; k++) {
        cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT-1));
        cheightnext = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k) + 1)/(AIR_HEIGHT-1));
        if (zenith == 0.0) {
            airmass_height[k] = 1.0;
        } else {
            if (zenith <= PI/2.0) {
                airmass_height[k] = RADIUS_EARTH*(sqrt(pow(1 + cheight/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                                  sqrt(pow(1 + cheightnext/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                    sin(zenith)/(cheight - cheightnext);
            } else {
                airmass_height[k]=1e6;
            }
        }
        if (moonZenith == 0.0) {
            airmass_heightMoon[k] = 1.0;
        } else {
            if (moonZenith <= PI/2.0) {
                airmass_heightMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + cheight/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                                      sqrt(pow(1 + cheightnext/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                    sin(moonZenith)/(cheight - cheightnext);
            } else {
                airmass_heightMoon[k] = 1e6;
            }
        }
    }

    for (int k = 0; k < layers + 1; k++) {
        if (zenith == 0.0) {
            airmassLayer[k] = 1.0;
        } else {
            if (zenith <= PI/2.0) {
                if (k == 0) {
                    airmassLayer[k] = RADIUS_EARTH*(sqrt(pow(1 + topAtmosphere/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                                    sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                        sin(zenith)/(topAtmosphere - height[k]);
                } else {
                    airmassLayer[k] = RADIUS_EARTH*(sqrt(pow(1 + height[k - 1]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1)-
                                                    sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(zenith)/sin(zenith) - 1))*
                        sin(zenith)/(height[k - 1] - height[k]);
                }
            } else {
                airmassLayer[k]=1e6;
            }
        }
        if (moonZenith == 0.0) {
            airmassLayerMoon[k] = 1.0;
        } else {
            if (moonZenith <= PI/2.0) {
                if (k == 0) {
                    airmassLayerMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + topAtmosphere/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                                        sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                        sin(moonZenith)/(topAtmosphere - height[k]);
                } else {
                    airmassLayerMoon[k] = RADIUS_EARTH*(sqrt(pow(1 + height[k - 1]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1)-
                                                        sqrt(pow(1 + height[k]/RADIUS_EARTH, 2.0)/sin(moonZenith)/sin(moonZenith) - 1))*
                        sin(moonZenith)/(height[k - 1] - height[k]);
                }
            } else {
                airmassLayerMoon[k]=1e6;
            }
        }
    }

    for (int i = 0; i < points; i++) {
        for (int k = 0; k < AIR_HEIGHT; k++) {
            wavelength = tauWavelength[i];
            cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT-1));

            if (cheight > groundlevel/1000.0) {

                opticalDepth(wavelength, cheight, cheight, prheight, &ttransmission, &ttransmission2, &ttransmission3, &ttransmission4, &ttransmission5, &ttransmission6, &ttransmission7, &ttransmission8, &ttransmission9, &ttransmission10, &ttransmission11, rayprofile,  o3cs, o3cs_wavelength, o2cs, o2cs_wavelength, h2ocs, h2ocs_wavelength, aerosol, 0);

                ttransmission12 = cloud[k];

                transmission = raynorm*airmass_height[k]*ttransmission;
                transmission2 = o2norm*airmass_height[k]*ttransmission2;
                transmission3 = h2onorm*airmass_height[k]*ttransmission3;
                transmission4 = o3norm*airmass_height[k]*ttransmission4;
                transmission5 = n2onorm*airmass_height[k]*ttransmission5;
                transmission6 = conorm*airmass_height[k]*ttransmission6;
                transmission7 = ch4norm*airmass_height[k]*ttransmission7;
                transmission8 = airmass_height[k]*ttransmission8;
                transmission9 = airmass_height[k]*ttransmission9;
                transmission10 = airmass_height[k]*ttransmission10;
                transmission11 = airmass_height[k]*ttransmission11;
                transmission12 = airmass_height[k]*ttransmission12;

                transmissionMoon = raynorm*airmass_heightMoon[k]*ttransmission;
                transmissionMoon2 = o3norm*airmass_heightMoon[k]*ttransmission2;
                transmissionMoon3 = h2onorm*airmass_heightMoon[k]*ttransmission3;
                transmissionMoon4 = o3norm*airmass_heightMoon[k]*ttransmission4;
                transmissionMoon5 = n2onorm*airmass_heightMoon[k]*ttransmission5;
                transmissionMoon6 = conorm*airmass_heightMoon[k]*ttransmission6;
                transmissionMoon7 = ch4norm*airmass_heightMoon[k]*ttransmission7;
                transmissionMoon8 = airmass_heightMoon[k]*ttransmission8;
                transmissionMoon9 = airmass_heightMoon[k]*ttransmission9;
                transmissionMoon10 = airmass_heightMoon[k]*ttransmission10;
                transmissionMoon11 = airmass_heightMoon[k]*ttransmission11;
                transmissionMoon12 = airmass_heightMoon[k]*ttransmission12;

                minchi = fabs(log10(2000.0) - log10(cheight));
                closestlayer = -1;
                for (int l = 0; l < layers; l++) {
                    chi = fabs(log10(height[l] + groundlevel/1000.0) - log10(cheight));
                    if (chi < minchi) {
                        closestlayer = l;
                        minchi = chi;
                    }
                }

                *(tau[12*(closestlayer + 1) + 0] + i) += transmission;
                *(tau[12*(closestlayer + 1) + 1] + i) += transmission2;
                *(tau[12*(closestlayer + 1) + 2] + i) += transmission3;
                *(tau[12*(closestlayer + 1) + 3] + i) += transmission4;
                *(tau[12*(closestlayer + 1) + 4] + i) += transmission5;
                *(tau[12*(closestlayer + 1) + 5] + i) += transmission6;
                *(tau[12*(closestlayer + 1) + 6] + i) += transmission7;
                *(tau[12*(closestlayer + 1) + 7] + i) += transmission8;
                *(tau[12*(closestlayer + 1) + 8] + i) += transmission9;
                *(tau[12*(closestlayer + 1) + 9] + i) += transmission10;
                *(tau[12*(closestlayer + 1) + 10] + i) += transmission11;
                *(tau[12*(closestlayer + 1) + 11] + i) += transmission12;
                *(tauMoon[12*(closestlayer + 1) + 0] + i) += transmissionMoon;
                *(tauMoon[12*(closestlayer + 1) + 1] + i) += transmissionMoon2;
                *(tauMoon[12*(closestlayer + 1) + 2] + i) += transmissionMoon3;
                *(tauMoon[12*(closestlayer + 1) + 3] + i) += transmissionMoon4;
                *(tauMoon[12*(closestlayer + 1) + 4] + i) += transmissionMoon5;
                *(tauMoon[12*(closestlayer + 1) + 5] + i) += transmissionMoon6;
                *(tauMoon[12*(closestlayer + 1) + 6] + i) += transmissionMoon7;
                *(tauMoon[12*(closestlayer + 1) + 7] + i) += transmissionMoon8;
                *(tauMoon[12*(closestlayer + 1) + 8] + i) += transmissionMoon9;
                *(tauMoon[12*(closestlayer + 1) + 9] + i) += transmissionMoon10;
                *(tauMoon[12*(closestlayer + 1) + 10] + i) += transmissionMoon11;
                *(tauMoon[12*(closestlayer + 1) + 11] + i) += transmissionMoon12;
                prheight = cheight;
            }
        }
    }

    // new moon loop
    std::vector<double> taumin(layers+1);
    for (int l = -1 ; l < layers; l++) {
        taumin[l+1]=1e30;
    }
    for (int i = 0; i < points; i++) {
        std::vector<int> counter(layers+1);
        for (int k = 0; k < AIR_HEIGHT_SCATTER; k++) {
            wavelength = tauWavelength[i];
            cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT_SCATTER-1));

            int l;
            l = k;
            prheightprime = cheight;
            double hmin = (RADIUS_EARTH + cheight)*sin(moonZenith) - RADIUS_EARTH;
            int direction;
            if (moonZenith < PI/2.0) {
                direction=-1;
            } else {
                direction=1;
            }
            double transmissionMoonScatter=0.0;
            double transmissionMoonScatter2=0.0;
            double transmissionMoonScatter3=0.0;
            double transmissionMoonScatter4=0.0;
            double transmissionMoonScatter5=0.0;
            double transmissionMoonScatter6=0.0;
            double transmissionMoonScatter7=0.0;
            double transmissionMoonScatter8=0.0;
            double transmissionMoonScatter9=0.0;
            double transmissionMoonScatter10=0.0;
            double transmissionMoonScatter11=0.0;
            double transmissionMoonScatter12=0.0;

            while (l > 0) {
                cheightprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l))/(AIR_HEIGHT_SCATTER-1));
                if (cheightprime <= hmin) direction=-1;
                if (direction == 1) {
                    if (l < AIR_HEIGHT_SCATTER - 1) {
                        if (isnan(airmass_heightScatterMoon[k][l+1])) direction=-1;
                    }
                }
                if (direction == 1) {
                    l++;
                } else {
                    l--;
                }
                if (l > AIR_HEIGHT_SCATTER-1) {
                    transmissionMoonScatter = 1e12;
                    transmissionMoonScatter2 = 1e12;
                    transmissionMoonScatter3 = 1e12;
                    transmissionMoonScatter4 = 1e12;
                    transmissionMoonScatter5 = 1e12;
                    transmissionMoonScatter6 = 1e12;
                    transmissionMoonScatter7 = 1e12;
                    transmissionMoonScatter8 = 1e12;
                    transmissionMoonScatter9 = 1e12;
                    transmissionMoonScatter10 = 1e12;
                    transmissionMoonScatter11 = 1e12;
                    transmissionMoonScatter12 = 1e12;
                    goto abort;
                }

                opticalDepth(wavelength, cheightprime, cheightprime, prheightprime, &ttransmission, &ttransmission2, &ttransmission3, &ttransmission4, &ttransmission5,&ttransmission6,&ttransmission7,&ttransmission8,&ttransmission9,&ttransmission10,&ttransmission11, rayprofile,  o3cs, o3cs_wavelength, o2cs, o2cs_wavelength, h2ocs, h2ocs_wavelength, aerosol, 0);


                ttransmission12 = cloudSmall[l];

                transmissionMoonScatter += raynorm*airmass_heightScatterMoon[k][l]*ttransmission;
                transmissionMoonScatter2 += o2norm*airmass_heightScatterMoon[k][l]*ttransmission2;
                transmissionMoonScatter3 += h2onorm*airmass_heightScatterMoon[k][l]*ttransmission3;
                transmissionMoonScatter4 += o3norm*airmass_heightScatterMoon[k][l]*ttransmission4;
                transmissionMoonScatter5 += n2onorm*airmass_heightScatterMoon[k][l]*ttransmission5;
                transmissionMoonScatter6 += conorm*airmass_heightScatterMoon[k][l]*ttransmission6;
                transmissionMoonScatter7 += ch4norm*airmass_heightScatterMoon[k][l]*ttransmission7;
                transmissionMoonScatter8 += airmass_heightScatterMoon[k][l]*ttransmission8;
                transmissionMoonScatter9 += airmass_heightScatterMoon[k][l]*ttransmission9;
                transmissionMoonScatter10 += airmass_heightScatterMoon[k][l]*ttransmission10;
                transmissionMoonScatter11 += airmass_heightScatterMoon[k][l]*ttransmission11;
                transmissionMoonScatter12 += airmass_heightScatterMoon[k][l]*ttransmission12;

                prheightprime=cheightprime;
            }

            abort:;
            minchi = fabs(log10(2000.0) - log10(cheight));
            closestlayer = -1;
            for (int l = 0; l < layers; l++) {
                chi = fabs(log10(height[l] + groundlevel/1000.0) - log10(cheight));
                if (chi < minchi) {
                    closestlayer = l;
                    minchi = chi;
                }
            }

            if (isnan(transmissionMoonScatter)) printf("error 1 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter2)) printf("error 2 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter3)) printf("error 3 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter4)) printf("error 4 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter5)) printf("error 5 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter6)) printf("error 6 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter7)) printf("error 7 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter8)) printf("error 8 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter9)) printf("error 9 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter10)) printf("error 10 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter11)) printf("error 11 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionMoonScatter12)) printf("error 12 %d %d %e %d\n",i,k,cheight,l);

            *(tauMoonScatter[12*(closestlayer + 1) + 0] + i) += transmissionMoonScatter;
            *(tauMoonScatter[12*(closestlayer + 1) + 1] + i) += transmissionMoonScatter2;
            *(tauMoonScatter[12*(closestlayer + 1) + 2] + i) += transmissionMoonScatter3;
            *(tauMoonScatter[12*(closestlayer + 1) + 3] + i) += transmissionMoonScatter4;
            *(tauMoonScatter[12*(closestlayer + 1) + 4] + i) += transmissionMoonScatter5;
            *(tauMoonScatter[12*(closestlayer + 1) + 5] + i) += transmissionMoonScatter6;
            *(tauMoonScatter[12*(closestlayer + 1) + 6] + i) += transmissionMoonScatter7;
            *(tauMoonScatter[12*(closestlayer + 1) + 7] + i) += transmissionMoonScatter8;
            *(tauMoonScatter[12*(closestlayer + 1) + 8] + i) += transmissionMoonScatter9;
            *(tauMoonScatter[12*(closestlayer + 1) + 9] + i) += transmissionMoonScatter10;
            *(tauMoonScatter[12*(closestlayer + 1) + 10] + i) += transmissionMoonScatter11;
            *(tauMoonScatter[12*(closestlayer + 1) + 11] + i) += transmissionMoonScatter12;
            counter[(closestlayer+1)]=counter[(closestlayer+1)]+1;
            prheight = cheight;

            // then average it over the closest layers
        }
        for (int l = -1; l< layers; l++) {
            if (counter[l+1] > 0) {
                *(tauMoonScatter[12*(l + 1) + 0] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 1] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 2] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 3] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 4] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 5] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 6] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 7] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 8] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 9] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 10] + i) /= counter[l+1];
                *(tauMoonScatter[12*(l + 1) + 11] + i) /= counter[l+1];
            }
        }
        for (int l = -1; l< layers; l++) {
            if (counter[l+1] == 0) {
                minchi = 1e30;
                int bestk = 0;
                for (int k = -1; k < layers; k++) {
                    if (counter[k+1] > 0) {
                        chi = fabs(log10(height[l]) - log10(height[k]));
                        if (chi < minchi) {
                            minchi=chi;
                            bestk =k;
                        }
                    }
                }
                *(tauMoonScatter[12*(l + 1) + 0] + i) = *(tauMoonScatter[12*(bestk + 1) + 0] + i);
                *(tauMoonScatter[12*(l + 1) + 1] + i) = *(tauMoonScatter[12*(bestk + 1) + 1] + i);
                *(tauMoonScatter[12*(l + 1) + 2] + i) = *(tauMoonScatter[12*(bestk + 1) + 2] + i);
                *(tauMoonScatter[12*(l + 1) + 3] + i) = *(tauMoonScatter[12*(bestk + 1) + 3] + i);
                *(tauMoonScatter[12*(l + 1) + 4] + i) = *(tauMoonScatter[12*(bestk + 1) + 4] + i);
                *(tauMoonScatter[12*(l + 1) + 5] + i) = *(tauMoonScatter[12*(bestk + 1) + 5] + i);
                *(tauMoonScatter[12*(l + 1) + 6] + i) = *(tauMoonScatter[12*(bestk + 1) + 6] + i);
                *(tauMoonScatter[12*(l + 1) + 7] + i) = *(tauMoonScatter[12*(bestk + 1) + 7] + i);
                *(tauMoonScatter[12*(l + 1) + 8] + i) = *(tauMoonScatter[12*(bestk + 1) + 8] + i);
                *(tauMoonScatter[12*(l + 1) + 9] + i) = *(tauMoonScatter[12*(bestk + 1) + 9] + i);
                *(tauMoonScatter[12*(l + 1) + 10] + i) = *(tauMoonScatter[12*(bestk + 1) + 10] + i);
                *(tauMoonScatter[12*(l + 1) + 11] + i) = *(tauMoonScatter[12*(bestk + 1) + 11] + i);
            }
        }
        for (int l = -1; l< layers; l++) {
            if (isnan(*(tauMoonScatter[12*(l + 1) + 0] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 0] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 1] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 1] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 2] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 2] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 3] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 3] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 4] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 4] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 5] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 5] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 6] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 6] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 7] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 7] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 8] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 8] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 9] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 9] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 10] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 10] + i));
            if (isnan(*(tauMoonScatter[12*(l + 1) + 11] + i))) printf("Error %e\n",*(tauMoonScatter[12*(l + 1) + 11] + i));
        }

    }

    // sun loop
   for (int i = 0; i < points; i++) {
        std::vector<int> counter(layers+1);
        for (int k = 0; k < AIR_HEIGHT_SCATTER; k++) {
            wavelength = tauWavelength[i];
            cheight = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(k))/(AIR_HEIGHT_SCATTER-1));

            int l;
            l = k;
            prheightprime = cheight;
            double hmin = (RADIUS_EARTH + cheight)*sin(sunZenith) - RADIUS_EARTH;
            int direction;
            if (sunZenith < PI/2.0) {
                direction=-1;
            } else {
                direction=1;
            }
            double transmissionSunScatter=0.0;
            double transmissionSunScatter2=0.0;
            double transmissionSunScatter3=0.0;
            double transmissionSunScatter4=0.0;
            double transmissionSunScatter5=0.0;
            double transmissionSunScatter6=0.0;
            double transmissionSunScatter7=0.0;
            double transmissionSunScatter8=0.0;
            double transmissionSunScatter9=0.0;
            double transmissionSunScatter10=0.0;
            double transmissionSunScatter11=0.0;
            double transmissionSunScatter12=0.0;

            while (l > 0) {
                cheightprime = pow(10.0, log10(2000.0) - 6.0*(static_cast<double>(l))/(AIR_HEIGHT_SCATTER-1));
                if (cheightprime <= hmin) direction=-1;
                if (direction == 1) {
                    if (l < AIR_HEIGHT_SCATTER - 1) {
                        if (isnan(airmass_heightScatterSun[k][l+1])) direction=-1;
                    }
                }
                if (direction == 1) {
                    l++;
                } else {
                    l--;
                }
                if (l > AIR_HEIGHT_SCATTER-1) {
                    transmissionSunScatter = 1e12;
                    transmissionSunScatter2 = 1e12;
                    transmissionSunScatter3 = 1e12;
                    transmissionSunScatter4 = 1e12;
                    transmissionSunScatter5 = 1e12;
                    transmissionSunScatter6 = 1e12;
                    transmissionSunScatter7 = 1e12;
                    transmissionSunScatter8 = 1e12;
                    transmissionSunScatter9 = 1e12;
                    transmissionSunScatter10 = 1e12;
                    transmissionSunScatter11 = 1e12;
                    transmissionSunScatter12 = 1e12;
                    goto abort2;
                }

                opticalDepth(wavelength, cheightprime, cheightprime, prheightprime, &ttransmission, &ttransmission2, &ttransmission3, &ttransmission4, &ttransmission5, &ttransmission6, &ttransmission7, &ttransmission8,&ttransmission9,&ttransmission10,&ttransmission11, rayprofile,  o3cs, o3cs_wavelength, o2cs, o2cs_wavelength, h2ocs, h2ocs_wavelength, aerosol, 0);

                ttransmission12 = cloudSmall[l];

                transmissionSunScatter += raynorm*airmass_heightScatterSun[k][l]*ttransmission;
                transmissionSunScatter2 += o2norm*airmass_heightScatterSun[k][l]*ttransmission2;
                transmissionSunScatter3 += h2onorm*airmass_heightScatterSun[k][l]*ttransmission3;
                transmissionSunScatter4 += o3norm*airmass_heightScatterSun[k][l]*ttransmission4;
                transmissionSunScatter5 += n2onorm*airmass_heightScatterSun[k][l]*ttransmission5;
                transmissionSunScatter6 += conorm*airmass_heightScatterSun[k][l]*ttransmission6;
                transmissionSunScatter7 += ch4norm*airmass_heightScatterSun[k][l]*ttransmission7;
                transmissionSunScatter8 += airmass_heightScatterSun[k][l]*ttransmission8;
                transmissionSunScatter9 += airmass_heightScatterSun[k][l]*ttransmission9;
                transmissionSunScatter10 += airmass_heightScatterSun[k][l]*ttransmission10;
                transmissionSunScatter11 += airmass_heightScatterSun[k][l]*ttransmission11;
                transmissionSunScatter12 += airmass_heightScatterSun[k][l]*ttransmission12;

                prheightprime=cheightprime;
            }

            abort2:;
            minchi = fabs(log10(2000.0) - log10(cheight));
            closestlayer = -1;
            for (int l = 0; l < layers; l++) {
                chi = fabs(log10(height[l] + groundlevel/1000.0) - log10(cheight));
                if (chi < minchi) {
                    closestlayer = l;
                    minchi = chi;
                }
            }

            if (isnan(transmissionSunScatter)) printf("error 1 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter2)) printf("error 2 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter3)) printf("error 3 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter4)) printf("error 4 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter5)) printf("error 5 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter6)) printf("error 6 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter7)) printf("error 7 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter8)) printf("error 8 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter9)) printf("error 9 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter10)) printf("error 10 %d %d %e %d\n",i,k,cheight,l);
            if (isnan(transmissionSunScatter11)) printf("error 11 %d %d %e %d\n",i,k,cheight,l);


            *(tauSunScatter[12*(closestlayer + 1) + 0] + i) += transmissionSunScatter;
            *(tauSunScatter[12*(closestlayer + 1) + 1] + i) += transmissionSunScatter2;
            *(tauSunScatter[12*(closestlayer + 1) + 2] + i) += transmissionSunScatter3;
            *(tauSunScatter[12*(closestlayer + 1) + 3] + i) += transmissionSunScatter4;
            *(tauSunScatter[12*(closestlayer + 1) + 4] + i) += transmissionSunScatter5;
            *(tauSunScatter[12*(closestlayer + 1) + 5] + i) += transmissionSunScatter6;
            *(tauSunScatter[12*(closestlayer + 1) + 6] + i) += transmissionSunScatter7;
            *(tauSunScatter[12*(closestlayer + 1) + 7] + i) += transmissionSunScatter8;
            *(tauSunScatter[12*(closestlayer + 1) + 8] + i) += transmissionSunScatter9;
            *(tauSunScatter[12*(closestlayer + 1) + 9] + i) += transmissionSunScatter10;
            *(tauSunScatter[12*(closestlayer + 1) + 10] + i) += transmissionSunScatter11;
            *(tauSunScatter[12*(closestlayer + 1) + 11] + i) += transmissionSunScatter12;
            counter[(closestlayer+1)]=counter[(closestlayer+1)]+1;
            prheight = cheight;

            // then average it over the closest layers
        }

        for (int l = -1; l< layers; l++) {
            if (counter[l+1] > 0) {
                *(tauSunScatter[12*(l + 1) + 0] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 1] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 2] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 3] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 4] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 5] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 6] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 7] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 8] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 9] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 10] + i) /= counter[l+1];
                *(tauSunScatter[12*(l + 1) + 11] + i) /= counter[l+1];
            }
        }
       for (int l = -1; l< layers; l++) {
            if (counter[l+1] == 0) {
                minchi = 1e30;
                int bestk = 0;
                for (int k = -1; k < layers; k++) {
                    if (counter[k+1] > 0) {
                        chi = fabs(log10(height[l]) - log10(height[k]));
                        if (chi < minchi) {
                            minchi=chi;
                            bestk =k;
                        }
                    }
                }
                *(tauSunScatter[12*(l + 1) + 0] + i) = *(tauSunScatter[12*(bestk + 1) + 0] + i);
                *(tauSunScatter[12*(l + 1) + 1] + i) = *(tauSunScatter[12*(bestk + 1) + 1] + i);
                *(tauSunScatter[12*(l + 1) + 2] + i) = *(tauSunScatter[12*(bestk + 1) + 2] + i);
                *(tauSunScatter[12*(l + 1) + 3] + i) = *(tauSunScatter[12*(bestk + 1) + 3] + i);
                *(tauSunScatter[12*(l + 1) + 4] + i) = *(tauSunScatter[12*(bestk + 1) + 4] + i);
                *(tauSunScatter[12*(l + 1) + 5] + i) = *(tauSunScatter[12*(bestk + 1) + 5] + i);
                *(tauSunScatter[12*(l + 1) + 6] + i) = *(tauSunScatter[12*(bestk + 1) + 6] + i);
                *(tauSunScatter[12*(l + 1) + 7] + i) = *(tauSunScatter[12*(bestk + 1) + 7] + i);
                *(tauSunScatter[12*(l + 1) + 8] + i) = *(tauSunScatter[12*(bestk + 1) + 8] + i);
                *(tauSunScatter[12*(l + 1) + 9] + i) = *(tauSunScatter[12*(bestk + 1) + 9] + i);
                *(tauSunScatter[12*(l + 1) + 10] + i) = *(tauSunScatter[12*(bestk + 1) + 10] + i);
                *(tauSunScatter[12*(l + 1) + 11] + i) = *(tauSunScatter[12*(bestk + 1) + 11] + i);
            }
       }
      for (int l = -1; l< layers; l++) {
           if (isnan(*(tauSunScatter[12*(l + 1) + 0] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 0] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 1] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 1] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 2] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 2] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 3] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 3] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 4] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 4] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 5] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 5] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 6] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 6] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 7] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 7] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 8] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 8] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 9] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 9] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 10] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 10] + i));
            if (isnan(*(tauSunScatter[12*(l + 1) + 11] + i))) printf("Error %e\n",*(tauSunScatter[12*(l + 1) + 11] + i));
       }

        for (int l = -1; l< layers; l++) {
            if (*(tauSunScatter[12*(l + 1) + 0] + i) < taumin[l+1]) {
                taumin[l+1]=*(tauSunScatter[12*(l + 1) + 0] + i);
            }
        }


   }


}



void Air::readReferenceAtmosphere (std::string dir, double groundTemperature, double latitude, double pressure,  double waterPressure, double exosphereTemperature, double longitude) {

    double temperatureKelvin = groundTemperature + CELCIUS_KELVIN;
    double latitudeDegree = latitude/DEGREE;
    double waterFraction = waterPressure/pressure;
    double pressurePascal =  pressure/PRESS_ATMOSPHERE_MMHG*PRESS_ATMOSPHERE;

    std::vector<double> co2TRA;
    std::vector<double> n2TRA;
    std::vector<double> o2TRA;
    std::vector<double> oTRA;
    std::vector<double> arTRA;
    std::vector<double> noTRA;
    std::vector<double> so2TRA;
    std::vector<double> no2TRA;
    std::vector<double> nh3TRA;
    std::vector<double> hno3TRA;
    std::vector<double> ohTRA;
    std::vector<double> hfTRA;

    std::vector<double> hclTRA;
    std::vector<double> hbrTRA;
    std::vector<double> hiTRA;
    std::vector<double> cloTRA;
    std::vector<double> ocsTRA;
    std::vector<double> h2coTRA;
    std::vector<double> hoclTRA;

    std::vector<double> hcnTRA;
    std::vector<double> ch3clTRA;
    std::vector<double> h2o2TRA;
    std::vector<double> c2h2TRA;
    std::vector<double> c2h6TRA;
    std::vector<double> ph3TRA;

    std::vector<double> neTRA;
    std::vector<double> krTRA;

    std::vector<double> h2oSTA;
    std::vector<double> o3STA;
    std::vector<double> n2oSTA;
    std::vector<double> coSTA;
    std::vector<double> ch4STA;

    //    {
        std::vector<std::vector<double> > tempArr;
        std::vector<double> tempCol;
        readText::readMultiCol(dir + "/us_standard.txt", tempCol, tempArr);
        for (unsigned int i = 0; i < tempCol.size() ; i++) h2oSTA.push_back(tempArr[i][3]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) o3STA.push_back(tempArr[i][4]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) n2oSTA.push_back(tempArr[i][5]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) coSTA.push_back(tempArr[i][6]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) ch4STA.push_back(tempArr[i][7]/1e6);

        std::vector<std::vector<double> > tempArr1;
        std::vector<double> tempCol1;
        readText::readMultiCol(dir + "/us_standard_trace_1.txt", tempCol1, tempArr1);
        for (unsigned int i = 0; i < tempCol1.size() ; i++) co2TRA.push_back(tempArr1[i][1]/1e6);
        for (unsigned int i = 0; i < tempCol1.size() ; i++) o2TRA.push_back(tempArr1[i][6]/1e6);

        std::vector<std::vector<double> > tempArr2;
        std::vector<double> tempCol2;
        readText::readMultiCol(dir + "/us_standard_trace_2.txt", tempCol2, tempArr2);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) noTRA.push_back(tempArr2[i][0]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) so2TRA.push_back(tempArr2[i][1]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) no2TRA.push_back(tempArr2[i][2]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) nh3TRA.push_back(tempArr2[i][3]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) hno3TRA.push_back(tempArr2[i][4]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) ohTRA.push_back(tempArr2[i][5]/1e6);
        for (unsigned int i = 0; i < tempCol2.size() ; i++) hfTRA.push_back(tempArr2[i][6]/1e6);

        std::vector<std::vector<double> > tempArr3;
        std::vector<double> tempCol3;
        readText::readMultiCol(dir + "/us_standard_trace_3.txt", tempCol3, tempArr3);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) hclTRA.push_back(tempArr2[i][0]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) hbrTRA.push_back(tempArr2[i][1]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) hiTRA.push_back(tempArr2[i][2]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) cloTRA.push_back(tempArr2[i][3]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) ocsTRA.push_back(tempArr2[i][4]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) h2coTRA.push_back(tempArr2[i][5]/1e6);
        for (unsigned int i = 0; i < tempCol3.size() ; i++) hoclTRA.push_back(tempArr2[i][6]/1e6);

        std::vector<std::vector<double> > tempArr4;
        std::vector<double> tempCol4;
        readText::readMultiCol(dir + "/us_standard_trace_4.txt", tempCol4, tempArr4);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) n2TRA.push_back(tempArr4[i][0]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) hcnTRA.push_back(tempArr4[i][1]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) ch3clTRA.push_back(tempArr4[i][2]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) h2o2TRA.push_back(tempArr4[i][3]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) c2h2TRA.push_back(tempArr4[i][4]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) c2h6TRA.push_back(tempArr4[i][5]/1e6);
        for (unsigned int i = 0; i < tempCol4.size() ; i++) ph3TRA.push_back(tempArr4[i][6]/1e6);

        std::vector<std::vector<double> > tempArr5;
        std::vector<double> tempCol5;
        readText::readMultiCol(dir + "/us_standard_trace_5.txt", tempCol5, tempArr5);
        for (unsigned int i = 0; i < tempCol5.size() ; i++) arTRA.push_back(tempArr5[i][0]/1e6);
        for (unsigned int i = 0; i < tempCol5.size() ; i++) oTRA.push_back(tempArr5[i][1]/1e6);
        for (unsigned int i = 0; i < tempCol5.size() ; i++) neTRA.push_back(tempArr5[i][2]/1e6);
        for (unsigned int i = 0; i < tempCol5.size() ; i++) krTRA.push_back(tempArr5[i][3]/1e6);


     for (unsigned int i = 0; i < tempCol.size(); i++) {
         double total = h2oSTA[i] + o3STA[i] -n2oSTA[i] + coSTA[i] + ch4STA[i] + co2TRA[i] +
              o2TRA[i] + noTRA[i] +so2TRA[i] + no2TRA[i] + nh3TRA[i] + hno3TRA[i] + ohTRA[i] + hfTRA[i] +
              hclTRA[i] + hbrTRA[i] + hiTRA[i] + cloTRA[i] + ocsTRA[i] + h2coTRA[i] + hoclTRA[i] +
              n2TRA[i] + hcnTRA[i] + ch3clTRA[i] + h2o2TRA[i] + c2h2TRA[i] + c2h6TRA[i] + ph3TRA[i] +
              arTRA[i] + oTRA[i] + neTRA[i] + krTRA[i];
     total = 1 - total;

     }
     //  }

    std::vector<std::vector<double> > altitudeRA;
    std::vector<std::vector<double> > densityRA;
    std::vector<std::vector<double> > temperatureRA;
    std::vector<std::vector<double> > pressureRA;
    std::vector<std::vector<double> > waterRA;
    std::vector<std::vector<double> > o3RA;
    std::vector<std::vector<double> > n2oRA;
    std::vector<std::vector<double> > coRA;
    std::vector<std::vector<double> > ch4RA;
    std::vector<double> referenceTemperatureRA;
    std::vector<double> referenceLatitudeRA;
    std::vector<double> referencePressureRA;
    std::vector<double> referenceWaterRA;

    altitudeRA.resize(5);
    temperatureRA.resize(5);
    pressureRA.resize(5);
    referenceTemperatureRA.resize(5);
    referenceLatitudeRA.resize(5);
    referencePressureRA.resize(5);
    waterRA.resize(5);
    o3RA.resize(5);
    n2oRA.resize(5);
    coRA.resize(5);
    ch4RA.resize(5);
    referenceWaterRA.resize(5);

    for (int j = 0; j < 5; j++) {
        std::vector<std::vector<double> > tempArr;
        std::vector<double> tempCol;
        if (j==0) readText::readMultiCol(dir + "/tropical.txt", tempCol, tempArr);
        if (j==1) readText::readMultiCol(dir + "/midlatitude_summer.txt", tempCol, tempArr);
        if (j==2) readText::readMultiCol(dir + "/midlatitude_winter.txt", tempCol, tempArr);
        if (j==3) readText::readMultiCol(dir + "/subarctic_summer.txt", tempCol, tempArr);
        if (j==4) readText::readMultiCol(dir + "/subarctic_winter.txt", tempCol, tempArr);
        for (unsigned int i = 0; i < tempCol.size(); i++) altitudeRA[j].push_back(tempCol[i]);
        for (unsigned int i = 0; i < tempCol.size() ; i++) temperatureRA[j].push_back(tempArr[i][1]);
        for (unsigned int i = 0; i < tempCol.size() ; i++) pressureRA[j].push_back(tempArr[i][0]*100.0);
        for (unsigned int i = 0; i < tempCol.size() ; i++) waterRA[j].push_back(tempArr[i][3]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) o3RA[j].push_back(tempArr[i][4]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) n2oRA[j].push_back(tempArr[i][5]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) coRA[j].push_back(tempArr[i][6]/1e6);
        for (unsigned int i = 0; i < tempCol.size() ; i++) ch4RA[j].push_back(tempArr[i][7]/1e6);
        referenceTemperatureRA[j] = temperatureRA[j][0];
        referencePressureRA[j] = pressureRA[j][0];
        referenceWaterRA[j] = waterRA[j][0];
        if (j==0) referenceLatitudeRA[j] = 0.0;
        if (j==1) referenceLatitudeRA[j] = 45.0;
        if (j==2) referenceLatitudeRA[j] = 45.0;
        if (j==3) referenceLatitudeRA[j] = 60.0;
        if (j==4) referenceLatitudeRA[j] = 60.0;
    }

    std::vector<double> weight(5,0);
    std::vector<double> distance(5,0);

    double scale = 10.0;
    for (int i =0 ; i < 5 ; i++) {
        distance[i] = sqrt(pow((temperatureKelvin - referenceTemperatureRA[i] + 1.0)/scale, 2.0) + pow(fabs(latitudeDegree) - referenceLatitudeRA[i] + 0.1, 2.0));
        weight[i] = 1/pow(distance[i], 2.0);
    }
    double sumWeight = 0.0;
    for (int i = 0; i < 5; i++) {
        sumWeight += weight[i];
    }

    for (unsigned int j = 0; j < altitudeRA[0].size()-1 ; j++) {
        double temp = 0.0;
        double press = 0.0;
        double water = 0.0;
        double o3 = 0.0;
        double n2o = 0.0;
        double co = 0.0;
        double ch4 = 0.0;
        for (int i =0; i < 5; i++ ) {
            temp+= temperatureRA[i][j] * weight[i]/ sumWeight *
                (temperatureKelvin/referenceTemperatureRA[i] + altitudeRA[i][j]/altitudeRA[i][altitudeRA[0].size()-2] * (1-temperatureKelvin/referenceTemperatureRA[i]));
            press+= pressureRA[i][j] * weight[i]/ sumWeight *
                (pressurePascal/referencePressureRA[i] + altitudeRA[i][j]/altitudeRA[i][altitudeRA[0].size()-2]*(1-pressurePascal/referencePressureRA[i]));
            water += waterRA[i][j] * weight[i]/ sumWeight *
                (waterFraction/referenceWaterRA[i] + altitudeRA[i][j]/altitudeRA[i][altitudeRA[0].size()-2]*(1-waterFraction/referenceWaterRA[i]));
            o3 += o3RA[i][j] * weight[i]/ sumWeight;
            n2o += n2oRA[i][j] * weight[i]/ sumWeight;
            co += coRA[i][j] * weight[i]/ sumWeight;
            ch4 += ch4RA[i][j] * weight[i]/ sumWeight;
        }
        altitudeProfile.push_back(altitudeRA[0][j]);
        temperatureProfile.push_back(temp);
        pressureProfile.push_back(press);
        //        if (j==0) printf("%e %e\n",press,pressurePascal);
        //if (j==0) printf("%e %e\n",temp,temperatureKelvin);
        //if (j==0) printf("%e %e\n",water,waterFraction);
        //        printf("%lf %lf\n",altitudeRA[0][j],press*PRESS_ATMOSPHERE_MMHG/PRESS_ATMOSPHERE);
        // variable components
        double waterfactor = 1 - water;
        h2oProfile.push_back(water);
        o3Profile.push_back(o3*waterfactor);
        n2oProfile.push_back(n2o*waterfactor);
        coProfile.push_back(co*waterfactor);
        ch4Profile.push_back(ch4*waterfactor);
        //from trace
        //renormalize
        double total = o3 + n2o + co + ch4;
        double totalSTA = o3STA[j] + n2oSTA[j] +coSTA[j] + ch4STA[j];
        double renorm = (1-total)/(1-totalSTA)*waterfactor;

        co2Profile.push_back(co2TRA[j]*renorm);
        n2Profile.push_back(n2TRA[j]*renorm);
        o2Profile.push_back(o2TRA[j]*renorm);
        hProfile.push_back(0.0);
        heProfile.push_back(0.0);

        noProfile.push_back(noTRA[j]*renorm);
        so2Profile.push_back(so2TRA[j]*renorm);
        no2Profile.push_back(no2TRA[j]*renorm);
        nh3Profile.push_back(nh3TRA[j]*renorm);
        hno3Profile.push_back(hno3TRA[j]*renorm);
        ohProfile.push_back(ohTRA[j]*renorm);
        hfProfile.push_back(hfTRA[j]*renorm);

        hclProfile.push_back(hclTRA[j]*renorm);
        hbrProfile.push_back(hbrTRA[j]*renorm);
        hiProfile.push_back(hiTRA[j]*renorm);
        cloProfile.push_back(cloTRA[j]*renorm);
        ocsProfile.push_back(ocsTRA[j]*renorm);
        h2coProfile.push_back(h2coTRA[j]*renorm);
        hoclProfile.push_back(hoclTRA[j]*renorm);

        hcnProfile.push_back(hcnTRA[j]*renorm);
        ch3clProfile.push_back(ch3clTRA[j]*renorm);
        h2o2Profile.push_back(h2o2TRA[j]*renorm);
        c2h2Profile.push_back(c2h2TRA[j]*renorm);
        c2h6Profile.push_back(c2h6TRA[j]*renorm);
        ph3Profile.push_back(ph3TRA[j]*renorm);

        arProfile.push_back(arTRA[j]*renorm);
        oProfile.push_back(oTRA[j]*renorm);
        neProfile.push_back(neTRA[j]*renorm);
        krProfile.push_back(krTRA[j]*renorm);
    }


    std::vector<std::vector<double> > altitudeREA;
    std::vector<std::vector<double> > densityREA;
    std::vector<std::vector<double> > temperatureREA;
    std::vector<std::vector<double> > pressureREA;
    std::vector<std::vector<double> > n2REA;
    std::vector<std::vector<double> > o2REA;
    std::vector<std::vector<double> > oREA;
    std::vector<std::vector<double> > arREA;
    std::vector<std::vector<double> > heREA;
    std::vector<std::vector<double> > hREA;
    std::vector<std::vector<double> > muREA;
    std::vector<double> referenceTemperatureREA;

    altitudeREA.resize(9);
    temperatureREA.resize(9);
    pressureREA.resize(9);
    referenceTemperatureREA.resize(9);
    n2REA.resize(9);
    o2REA.resize(9);
    oREA.resize(9);
    arREA.resize(9);
    heREA.resize(9);
    hREA.resize(9);
    muREA.resize(9);

for (int j = 0; j < 9; j++) {
        std::vector<std::vector<double> > tempArr;
        std::vector<double> tempCol;
        if (j == 0) readText::readMultiCol(dir + "/exosphere_500.txt", tempCol, tempArr);
        if (j == 1) readText::readMultiCol(dir + "/exosphere_600.txt", tempCol, tempArr);
        if (j == 2) readText::readMultiCol(dir + "/exosphere_700.txt", tempCol, tempArr);
        if (j == 3) readText::readMultiCol(dir + "/exosphere_800.txt", tempCol, tempArr);
        if (j == 4) readText::readMultiCol(dir + "/exosphere_900.txt", tempCol, tempArr);
        if (j == 5) readText::readMultiCol(dir + "/exosphere_1000.txt", tempCol, tempArr);
        if (j == 6) readText::readMultiCol(dir + "/exosphere_1200.txt", tempCol, tempArr);
        if (j == 7) readText::readMultiCol(dir + "/exosphere_1400.txt", tempCol, tempArr);
        if (j == 8) readText::readMultiCol(dir + "/exosphere_1600.txt", tempCol, tempArr);
        for (unsigned int i = 0; i < tempCol.size(); i++) altitudeREA[j].push_back(tempCol[i]);
        for (unsigned int i = 0; i < tempCol.size() ; i++) temperatureREA[j].push_back(tempArr[i][0]);
        for (unsigned int i = 0; i < tempCol.size() ; i++) pressureREA[j].push_back(K_BOLTZMANN*tempArr[i][0]/(tempArr[i][7]/N_AVOGADRO/1000.0)*tempArr[i][8]);
        for (unsigned int i = 0; i < tempCol.size(); i++) {
            if (tempArr[i][1] != 0.0) {
                n2REA[j].push_back(pow(10.0, tempArr[i][1])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
            } else {
                n2REA[j].push_back(0.0);
            }
        }
        for (unsigned int i = 0; i < tempCol.size(); i++) o2REA[j].push_back(pow(10.0, tempArr[i][2])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
        for (unsigned int i = 0; i < tempCol.size(); i++) oREA[j].push_back(pow(10.0, tempArr[i][3])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
        for (unsigned int i = 0; i < tempCol.size(); i++) arREA[j].push_back(pow(10.0, tempArr[i][4])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
        for (unsigned int i = 0; i < tempCol.size(); i++) heREA[j].push_back(pow(10.0, tempArr[i][5])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
        for (unsigned int i = 0; i < tempCol.size(); i++) {
            if (tempArr[i][6] != 0.0) {
                hREA[j].push_back(pow(10.0, tempArr[i][6])/(tempArr[i][8]*N_AVOGADRO*1000.0/tempArr[i][7]));
            } else {
                hREA[j].push_back(0.0);
            }
        }
        for (unsigned int i = 0; i < tempCol.size(); i++) muREA[j].push_back(tempArr[i][7]);
referenceTemperatureREA[j] = temperatureREA[j][tempCol.size() - 1];
 for (unsigned int i = 0; i < tempCol.size(); i++) {
     double total = n2REA[j][i] + o2REA[j][i] + oREA[j][i] + arREA[j][i] + heREA[j][i] + hREA[j][i];
     total = 1 - total;
     if (total > 5e-3) {
         if (hREA[j][i]==0.0) {
             hREA[j][i]=total;
         }
     }
 }
 }


    double low = 0;
    double high = 1e30;
    int highv = 0;
    int lowv = 0;

    for (int j = 0; j < 9; j++) {
        if ((exosphereTemperature < referenceTemperatureREA[j]) && (referenceTemperatureREA[j] < high)) {
            high = referenceTemperatureREA[j];
            highv = j;
        }
        if ((exosphereTemperature >= referenceTemperatureREA[j]) && (referenceTemperatureREA[j] > low)) {
            low = referenceTemperatureREA[j];
            lowv = j;
        }
    }

    for (unsigned int j = 0; j < altitudeREA[0].size() ; j++) {
        double f1 = (high - exosphereTemperature)/(high - low);
        double f2= (exosphereTemperature - low)/(high - low);
        temperatureProfile.push_back(f1 * temperatureREA[highv][j] +f2 * temperatureREA[lowv][j]);
        altitudeProfile.push_back(altitudeREA[0][j]);
        pressureProfile.push_back(f1 * pressureREA[highv][j] + f2 * pressureREA[lowv][j]);

        h2oProfile.push_back(0.0);
        co2Profile.push_back(0.0);
        o3Profile.push_back(0.0);
        n2oProfile.push_back(0.0);
        coProfile.push_back(0.0);
        ch4Profile.push_back(0.0);
        o2Profile.push_back(f1*o2REA[lowv][j] + f2*o2REA[highv][j]);


        noProfile.push_back(0.0);
        so2Profile.push_back(0.0);
        no2Profile.push_back(0.0);
        nh3Profile.push_back(0.0);
        hno3Profile.push_back(0.0);
        ohProfile.push_back(0.0);
        hfProfile.push_back(0.0);

        hclProfile.push_back(0.0);
        hbrProfile.push_back(0.0);
        hiProfile.push_back(0.0);
        cloProfile.push_back(0.0);
        ocsProfile.push_back(0.0);
        h2coProfile.push_back(0.0);
        hoclProfile.push_back(0.0);

        n2Profile.push_back(f1*n2REA[lowv][j] + f2*n2REA[highv][j]);
        hcnProfile.push_back(0.0);
        ch3clProfile.push_back(0.0);
        h2o2Profile.push_back(0.0);
        c2h2Profile.push_back(0.0);
        c2h6Profile.push_back(0.0);
        ph3Profile.push_back(0.0);

        oProfile.push_back(f1*oREA[lowv][j] + f2*oREA[highv][j]);
        arProfile.push_back(f1*arREA[lowv][j] + f2*arREA[highv][j]);
        neProfile.push_back(0.0);
        krProfile.push_back(0.0);

        heProfile.push_back(f1*heREA[lowv][j] + f2*heREA[highv][j]);
        hProfile.push_back(f1*hREA[lowv][j] + f2*hREA[highv][j]);

    }

    //adjust values between 100 km and 120 km to avoid discontinuities
    unsigned int highpoint=0;
    unsigned int lowpoint=0;
    for (unsigned int i = 0; i < altitudeProfile.size(); i++) {
        if (altitudeProfile[i] <= 120.0)  highpoint=i;
        if (altitudeProfile[i] <= 100.0)  lowpoint=i;
    }
    for (unsigned int i = lowpoint; i < highpoint; i++) {
        temperatureProfile[i] = temperatureProfile[i] * (120.0 - altitudeProfile[i])/20.0 + temperatureProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        pressureProfile[i] = pressureProfile[i] * (120.0 - altitudeProfile[i])/20.0 + pressureProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;

        h2oProfile[i] = h2oProfile[i] * (120.0 - altitudeProfile[i])/20.0 + h2oProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        co2Profile[i] = co2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + co2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        o3Profile[i] = o3Profile[i] * (120.0 - altitudeProfile[i])/20.0 + o3Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        n2oProfile[i] = n2oProfile[i] * (120.0 - altitudeProfile[i])/20.0 + n2oProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        coProfile[i] = coProfile[i] * (120.0 - altitudeProfile[i])/20.0 + coProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        ch4Profile[i] = ch4Profile[i] * (120.0 - altitudeProfile[i])/20.0 + ch4Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        o2Profile[i] = o2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + o2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;

        noProfile[i] = noProfile[i] * (120.0 - altitudeProfile[i])/20.0 + noProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        so2Profile[i] = so2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + so2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        no2Profile[i] = no2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + no2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        nh3Profile[i] = nh3Profile[i] * (120.0 - altitudeProfile[i])/20.0 + nh3Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hno3Profile[i] = hno3Profile[i] * (120.0 - altitudeProfile[i])/20.0 + hno3Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        ohProfile[i] = ohProfile[i] * (120.0 - altitudeProfile[i])/20.0 + ohProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hfProfile[i] = hfProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hfProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;

        hclProfile[i] = hclProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hclProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hbrProfile[i] = hbrProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hbrProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hiProfile[i] = hiProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hiProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        cloProfile[i] = cloProfile[i] * (120.0 - altitudeProfile[i])/20.0 + cloProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        ocsProfile[i] = ocsProfile[i] * (120.0 - altitudeProfile[i])/20.0 + ocsProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        h2coProfile[i] = h2coProfile[i] * (120.0 - altitudeProfile[i])/20.0 + h2coProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hoclProfile[i] = hoclProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hoclProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;

        n2Profile[i] = n2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + n2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hcnProfile[i] = hcnProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hcnProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        ch3clProfile[i] = ch3clProfile[i] * (120.0 - altitudeProfile[i])/20.0 + ch3clProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        h2o2Profile[i] = h2o2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + h2o2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        c2h2Profile[i] = c2h2Profile[i] * (120.0 - altitudeProfile[i])/20.0 + c2h2Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        c2h6Profile[i] = c2h6Profile[i] * (120.0 - altitudeProfile[i])/20.0 + c2h6Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        ph3Profile[i] = ph3Profile[i] * (120.0 - altitudeProfile[i])/20.0 + ph3Profile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;

        oProfile[i] = oProfile[i] * (120.0 - altitudeProfile[i])/20.0 + oProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        arProfile[i] = arProfile[i] * (120.0 - altitudeProfile[i])/20.0 + arProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        neProfile[i] = neProfile[i] * (120.0 - altitudeProfile[i])/20.0 + neProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        krProfile[i] = krProfile[i] * (120.0 - altitudeProfile[i])/20.0 + krProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        heProfile[i] = heProfile[i] * (120.0 - altitudeProfile[i])/20.0 + heProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
        hProfile[i] = hProfile[i] * (120.0 - altitudeProfile[i])/20.0 + hProfile[highpoint] * (altitudeProfile[i] - 100.0)/20.0;
    }

    for (unsigned int i = 0; i < altitudeProfile.size() ; i++) {
        double mu=0.0;
        mu += h2oProfile[i]*M_H2O + co2Profile[i]*M_CO2 + n2Profile[i]*M_N2+o2Profile[i]*M_O2+o3Profile[i]*M_O3;
        mu += n2oProfile[i]*M_N2O + coProfile[i]*M_CO + ch4Profile[i]*M_CH4 + oProfile[i]*M_O + arProfile[i]*M_AR;
        mu += heProfile[i]*M_HE + hProfile[i]*M_H + noProfile[i]*M_NO + so2Profile[i]*M_SO2 + no2Profile[i]*M_NO2 + nh3Profile[i]*M_NH3;
        mu += hno3Profile[i]*M_HNO3 + ohProfile[i]*M_OH + hfProfile[i]*M_HF + hclProfile[i]*M_HCL + hbrProfile[i]*M_HBR + hiProfile[i]*M_HI;
        mu += cloProfile[i]*M_CLO + ocsProfile[i]*M_OCS + h2coProfile[i] + hoclProfile[i]*M_HOCL + hcnProfile[i]*M_HCN + ch3clProfile[i]*M_CH3CL;
        mu += h2o2Profile[i]*M_H2O2 + c2h2Profile[i]*M_C2H2 + c2h6Profile[i]*M_C2H6 + ph3Profile[i]*M_PH3 + neProfile[i]*M_NE + krProfile[i]*M_KR;
        molecularMassProfile.push_back(mu);
        densityProfile.push_back(pressureProfile[i]*mu/N_AVOGADRO/1000.0/K_BOLTZMANN/temperatureProfile[i]);
        numberDensityProfile.push_back(pressureProfile[i]/K_BOLTZMANN/temperatureProfile[i]);
        //        if (i<=10) printf("%lf %e\n",altitudeProfile[i],numberDensityProfile[i]);
    }

    readWindPattern(dir, latitude, longitude, 0);

}

void Air::readWindPattern (std::string dir, double latitude, double longitude, int currentMonth) {


        std::vector<std::vector<double> > tempArr;
        std::vector<double> tempCol;

        readText::readMultiCol(dir + "/wind.txt", tempCol, tempArr);

        std::vector<double> longitudeA;
        std::vector<double> latitudeA;
        std::vector<double> level;
        std::vector<double> month;
        std::vector<double> umean;
        std::vector<double> usigma;
        std::vector<double> vmean;
        std::vector<double> vsigma;

        for (unsigned int i = 0; i < tempCol.size(); i++) longitudeA.push_back(tempCol[i]);
        for (unsigned int i = 0; i < tempCol.size(); i++) latitudeA.push_back(tempArr[i][0]);
        for (unsigned int i = 0; i < tempCol.size(); i++) level.push_back(tempArr[i][1]);
        for (unsigned int i = 0; i < tempCol.size(); i++) month.push_back(tempArr[i][2]);
        for (unsigned int i = 0; i < tempCol.size(); i++) umean.push_back(tempArr[i][3]);
        for (unsigned int i = 0; i < tempCol.size(); i++) usigma.push_back(tempArr[i][4]);
        for (unsigned int i = 0; i < tempCol.size(); i++) vmean.push_back(tempArr[i][5]);
        for (unsigned int i = 0; i < tempCol.size(); i++) vsigma.push_back(tempArr[i][6]);


        double minDis = 1e30;
        unsigned int best=0;
        for (unsigned int i = 0; i < longitudeA.size(); i++) {
            if ((level[i] == 1000.0) && (month[i]==currentMonth)) {
                double dis = acos(sin(latitude)*sin(latitudeA[i]) + cos(latitudeA[i])*cos(latitude)*cos(longitudeA[i] - longitude));
                if (dis < minDis) {
                    best = i;
                    minDis = dis;
                }
            }
        }

        std::vector<double> bestumean;
        std::vector<double> bestusigma;
        std::vector<double> bestvmean;
        std::vector<double> bestvsigma;
        std::vector<double> bestlevel;
        for (unsigned int i = 0; i < longitudeA.size(); i++) {
            if ((latitudeA[i] == latitudeA[best])  && (longitudeA[i] == longitudeA[best]) && (month[i]==currentMonth)) {
                bestlevel.push_back(level[i]);
                bestumean.push_back(umean[i]);
                bestusigma.push_back(usigma[i]);
                bestvmean.push_back(vmean[i]);
                bestvsigma.push_back(vsigma[i]);
            }
        }

        for (unsigned int i =0; i < altitudeProfile.size(); i++) {
            int index = 0;
            double press = pressureProfile[i]/PRESS_ATMOSPHERE*PRESS_ATMOSPHERE_MB;
            find_v_reverse(bestlevel, press, &index);
            windUMuProfile.push_back(interpolate_v(bestumean, bestlevel, press, index));
            windVMuProfile.push_back(interpolate_v(bestvmean, bestlevel, press, index));
            windUSigmaProfile.push_back(interpolate_v(bestusigma, bestlevel, press, index));
            windVSigmaProfile.push_back(interpolate_v(bestvsigma, bestlevel, press, index));


        }

}


void Air::setWindLayers (  std::vector<double> height, double groundlevel, std::vector<double> & wind, std::vector<double> & winddir, int layers, double r1, double r2) {

    for (int i = 0; i < layers; i++) {
        int index;
        double altitude = height[i] + groundlevel/1000.0;
        find_v(altitudeProfile, altitude, &index);
        double uwind = windUMuProfile[index] + windUSigmaProfile[index]*r1;
        double vwind = windVMuProfile[index] + windVSigmaProfile[index]*r2;
        wind[i] = sqrt(uwind*uwind + vwind*vwind);
        winddir[i] = atan2(uwind, vwind)/DEGREE;
    }



}
