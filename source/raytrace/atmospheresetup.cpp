///
/// @package phosim
/// @file atmospheresetup.cpp
/// @brief setup for atmosphere (part of image class)
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

#include <stdlib.h>
#include "ancillary/readtext.h"
#include "constants.h"

int Image::atmSetup () {

    char tempstring[4096];
    std::string dir;
    std::string dir2;

    if (diagnostic) printf("$$%40s %22e %22.6lf\n","atmSetup",(double)1,236.0);

    if (flatdir == 0) dir = datadir + "/atmosphere";
    else if (flatdir == 1) dir = ".";

    if (flatdir == 0) dir2 = datadir + "/material";
    else if (flatdir == 1) dir2 = ".";

    ranseed = obsseed + pairid;

    posix_memalign((void **)&random, PHOSIM_RANDOM_ALIGN, numthread*sizeof(Random));
    posix_memalign((void **)&galaxy.random, PHOSIM_RANDOM_ALIGN, numthread*sizeof(Random));

    for (int i = 0; i < numthread; i++) {
        new (&random[i]) Random();
        new (&galaxy.random[i]) Random();
    }
    if (obsseed == -1) {
        for (int i = 0; i < numthread; i++) {
            random[i].setSeedFromTime();
            galaxy.random[i].setSeedFromTime();
        }
    } else {
        for (int i = 0; i < numthread; i++) {
            random[i].setSeed32(ranseed+i*100);
            galaxy.random[i].setSeed32(ranseed+9999+i*100);
        }
    }

    for (int i = 0; i < numthread; i++) {
        random[i].unwind(100);
        galaxy.random[i].unwind(100);
    }

    if (devtype == "CMOS") {
        exptime = devvalue;
        timeoffset = 0.5*exptime + pairid*(vistime/nsnap) - 0.5*vistime;
    } else {
        exptime = (vistime - (nsnap - 1)*devvalue)/nsnap;
        timeoffset = 0.5*exptime + pairid*(devvalue + exptime) - 0.5*vistime;
    }

    exptime = exptime + shuttererror*random[0].normal();

    // SKY SETUP
    char sersicdata[4096];
    if (flatdir == 0) snprintf(sersicdata, sizeof(sersicdata), "%s/sky/sersic_const.txt", datadir.c_str());
    else if (flatdir == 1) snprintf(sersicdata, sizeof(sersicdata), "sersic_const.txt");
    galaxy.sampleSersic(sersicdata);
    galaxy.sampleSersic2d(sersicdata, numthread);
    dust.setup(maxwavelength);
    filterTruncateSources();
    if (splitsources == 1) splitSources();
    totalnorm = 0.0;
    for (long i = 0; i < nsource; i++) {
        totalnorm += sources.norm[i];
    }

    // BP SETUP
    bpVxMin -= 3.0*backRadius*ARCSEC;
    bpVxMax += 3.0*backRadius*ARCSEC;
    bpVyMin -= 3.0*backRadius*ARCSEC;
    bpVyMax += 3.0*backRadius*ARCSEC;
    bpStep = (bpBeta*pixsize)/platescale*DEGREE;
    bpNx = static_cast<int>((bpVxMax - bpVxMin)/bpStep);
    bpNy = static_cast<int>((bpVyMax - bpVyMin)/bpStep);

    // ATMOSPHERIC SETUP
    for (int i = 0; i < natmospherefile; i++) {
        seefactor[i] = pow(1/cos(zenith), 0.6)*seefactor[i];
    }

    screen.large_sizeperpixel = 5120.0*ATM_RESAMPLE;
    screen.coarse_sizeperpixel = 640.0*ATM_RESAMPLE;
    screen.medium_sizeperpixel = 80.0*ATM_RESAMPLE;
    screen.fine_sizeperpixel = 10.0*ATM_RESAMPLE;

    screen.wavelengthfactor_nom = pow(0.5, -0.2);

    // Air Setup
    if (atmospheremode) {
        fprintf(stdout, "Creating Air.\n");
        air.opacitySetup(zenith, moonalt, solaralt, height, groundlevel, raynorm,
                         o2norm, h2onorm, o3norm, n2onorm, conorm, ch4norm, aerosol,
                         natmospherefile, dir, &airmass, pressure, waterPressure, temperature, latitude, exosphereTemperature,
                         longitude, minwavelength, maxwavelength, cloudmean, cloudvary);
    } else {
        air.topAtmosphere = 0.010;
    }

    //Wind setup
    air.setWindLayers(height, groundlevel, wind, winddir, natmospherefile, windrandom1, windrandom2);

    // Diffraction Screens
    int over=1;
    if (atmospheremode>=2) over=numthread;
    //over= numthread;
    screen.phasescreen = new double[over*SCREEN_SIZE*SCREEN_SIZE]();
    screen.focalscreen = new double[over*SCREEN_SIZE*SCREEN_SIZE]();
    screen.tfocalscreen = new double[over*SCREEN_SIZE*SCREEN_SIZE]();
    screen.outscreen = (fftw_complex*)fftw_malloc(over*sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.inscreen = (fftw_complex*)fftw_malloc(over*sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.pupil_values = new float[over*SCREEN_SIZE*SCREEN_SIZE]();
    screen.pupil_temp = new float[over*SCREEN_SIZE*ATM_RESAMPLE*SCREEN_SIZE*ATM_RESAMPLE]();
    screen.focalscreencum = new double[over*SCREEN_SIZE*SCREEN_SIZE]();

    // Atmospheric Turbulence & Clouds
        screen.turbulenceLargeX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceLargeY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceCoarseX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceCoarseY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceMediumX = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.turbulenceMediumY = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.phaseLarge = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseCoarse = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseMedium = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseFine = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.phaseMediumH = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
        screen.phaseFineH = new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

        screen.see_norm = new float[natmospherefile]();
        screen.phase_norm = new float[natmospherefile]();

        if (atmospheremode) {
            fprintf(stdout, "Generating Turbulence.\n");

        for (int layer = 0; layer < natmospherefile; layer++) {

            if (cloudvary[layer] != 0 || cloudmean[layer] != 0) {
                screen.cloud[layer] = static_cast<float*>(calloc(SCREEN_SIZE*SCREEN_SIZE, sizeof(float)));
                snprintf(tempstring, sizeof(tempstring), "%s.fits.gz", cloudfile[layer].c_str());
                screen.readScreen(-1, screen.cloud[layer], tempstring);
             }

            snprintf(tempstring, sizeof(tempstring), "%s_largep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseLarge + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring,  sizeof(tempstring), "%s_coarsep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseCoarse + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_mediump.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseMedium + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_finep.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseFine + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            snprintf(tempstring, sizeof(tempstring), "%s_mediumh.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseMediumH + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_fineh.fits.gz", atmospherefile[layer].c_str());
            screen.phase_norm[layer] = screen.readScreen(9, screen.phaseFineH + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            snprintf(tempstring, sizeof(tempstring), "%s_largex.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_largey.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_coarsex.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_coarsey.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_mediumx.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);
            snprintf(tempstring, sizeof(tempstring), "%s_mediumy.fits.gz", atmospherefile[layer].c_str());
            screen.see_norm[layer] = screen.readScreen(9, screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE, tempstring);

            for (int i = 0; i < SCREEN_SIZE; i++) {
                for (int j = 0; j < SCREEN_SIZE; j++) {
                    *(screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceLargeX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceLargeY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceCoarseX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceCoarseY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceMediumX + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                    *(screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j) = (float)
                        (*(screen.turbulenceMediumY + layer*SCREEN_SIZE*SCREEN_SIZE + SCREEN_SIZE*i + j)*seefactor[layer]);
                }
            }

        }
        }

    double sigma=1.0/centralwavelength;
    double temp=temperature+273.15;
    double ps=pressure/760.00*1013.25;
    double pw=waterPressure/760.00*1013.25;
    double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
    double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
    double n=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
    n=n+(6478.31-58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
    air.air_refraction_adc = 1e-8*n;

    // air.air_refraction_adc = 64.328 + 29498.1/(146 - 1/centralwavelength/centralwavelength) + 255.4/(41 - 1/centralwavelength/centralwavelength);
    // air.air_refraction_adc = air.air_refraction_adc*pressure*(1 + (1.049 - 0.0157*temperature)*1e-6*pressure)/720.883/(1 + 0.003661*temperature);
    // air.air_refraction_adc = air.air_refraction_adc - ((0.0624 - 0.000680/centralwavelength/centralwavelength)/
    //                                                    (1 + 0.003661*temperature)*waterPressure);
    // air.air_refraction_adc = air.air_refraction_adc/1e6;


    // MIE setup
    if (atmospheremode) {
        fprintf(stdout, "Sprinkling Aerosols.\n");
        double x1 = 0.7 + 2.6*random[0].normalCorrel(tai + 16000,1.0);
        double y1 = (x1 - 0.7)*(-0.08) + 0.93 + 0.40*random[0].normalCorrel(tai + 11000, 1.0);
        double x2 = -4.2 + 1.3*random[0].normalCorrel(tai + 17000,1.0);
        double y2 = (x2 + 4.2)*(-0.20) + 1.13 + 0.11*random[0].normalCorrel(tai + 12000, 1.0);
        double x3 = -3.6  +  1.2*random[0].normalCorrel(tai + 18000,1.0);
        double y3 = (x3 + 3.6)*(-0.22) + 1.16 + 0.10*random[0].normalCorrel(tai + 13000, 1.0);
        double x4 = -4.3  + 1.3*random[0].normalCorrel(tai + 19000,1.0);
        double y4 = (x4 + 4.3 )*(-0.21) + 0.99 + 0.06*random[0].normalCorrel(tai + 14000, 1.0);
        double x5 = -4.9  + 1.2*random[0].normalCorrel(tai + 20000,1.0);
        double y5 = (x5 + 4.9)*(-0.22) + 0.95 + 0.06*random[0].normalCorrel(tai + 15000, 1.0);
        x1 = exp(x1);
        x2 = exp(x2);
        x3 = exp(x3);
        x4 = exp(x4);
        x5 = exp(x5);
        mie.setupMieIndexRefraction(dir2, minwavelength, maxwavelength,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
    }

    return(0);

}
