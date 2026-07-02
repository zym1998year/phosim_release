///
/// @package phosim
/// @file observation.cpp
/// @brief observation class for raytrace code
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Alan Meert (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "raytrace.h"
#include "observation.h"
#include "constants.h"
#include "helpers.h"
#include "constants.h"
#include "ancillary/readtext.h"
#include "coordinate.cpp"
#include "../atmosphere/operator.cpp"
#include "settings.cpp"
#include "sed.cpp"

int Observation::addSource (const std::string & object, int sourcetype) {

    char tempstring[4096];
    int nspatialpar, ndustpar, ndustparz;
    int status;
    char *ffptr;
    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    double x, y;
    double nn;
    double mag;
    double ra, dec, redshift, gamma1, gamma2, kappa, magnification,
        deltara, deltadec;
    double cd11 = 0.0, cd21 = 0.0, cd22 = 0.0, cd12 = 0.0;
    int tstat = 0;
    std::string id;


    nspatialpar = 0;

    std::istringstream iss(object);
    iss >> id >> ra >> dec;
    sources.id.push_back(id);
    sources.split.push_back(0);
    sources.ra.push_back(ra);
    sources.dec.push_back(dec);
    sources.backgroundAngle.push_back(0.0);
    sources.ra[nsource] *= DEGREE;
    sources.dec[nsource] *= DEGREE;
    std::string sedfilename;
    iss >> mag >> sedfilename >> redshift >> gamma1 >> gamma2 >> kappa >> deltara >> deltadec;
    sources.redshift.push_back(redshift);
    sources.gamma1.push_back(gamma1);
    sources.gamma2.push_back(gamma2);
    sources.kappa.push_back(kappa);
    deltara *= DEGREE;
    deltadec *= DEGREE;
    sources.type.push_back(sourcetype);
    magnification = 2.5*log10((1 - sources.kappa[nsource])*(1 - sources.kappa[nsource])-
                              sources.gamma1[nsource]*sources.gamma1[nsource]-
                              sources.gamma2[nsource]*sources.gamma2[nsource]);
    sources.norm.push_back(pow(10.0, ((mag + magnification + 48.6)/(-2.5))));
    sources.mag.push_back(mag + magnification);
    sedPtr.push_back(0);
    sedN.push_back(0);
    sedMode.push_back(0);
    sedCorr.push_back(0.0);
    sedDwdp.push_back(0.0);
    sourceXpos.push_back(0);
    sourceYpos.push_back(0);
    sourcePhoton.push_back(0);
    sources.spatialpar.push_back(std::vector<double>(15, 0));
    sources.dustpar.push_back(std::vector<double>(2, 0));
    sources.dustparz.push_back(std::vector<double>(2, 0));

    // read SED file
    if (nsource > 0) {
        for (int i = 0; i < nsed; i++) {
            if (sedlist[i] == sedfilename) {
                sources.sedptr.push_back(sources.sedptr[i]);
                goto skipsedread;
            }
        }
    }
    sources.sedptr.push_back(nsedptr);
    sedlist.push_back(sedfilename);
    nsed++;
    readSed(sedfilename, 0);

skipsedread:;
    sources.norm[nsource] = sources.norm[nsource]/(normwave)*(1 + sources.redshift[nsource])*sedDwdp[sources.sedptr[nsource]];

    std::string spatialname;
    iss >> spatialname;

    if (spatialname == "point") {
        sources.spatialtype.push_back(POINT);
        nspatialpar = 0;
    } else if (spatialname.find("fit") != std::string::npos) {
        sources.spatialtype.push_back(IMAGE);
        nspatialpar = 2;
    } else if (spatialname == "gauss") {
        sources.spatialtype.push_back(GAUSSIAN);
        nspatialpar = 1;
    } else if (spatialname == "sersic") {
        sources.spatialtype.push_back(SERSIC);
        nspatialpar = 7;
    } else if (spatialname == "sersic2d") {
        sources.spatialtype.push_back(SERSIC2D);
        nspatialpar = 4;
    } else if (spatialname == "sersic2D") {
        sources.spatialtype.push_back(SERSIC2D);
        nspatialpar = 4;
    } else if (spatialname == "sersicComplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersiccomplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersicDiskComplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersicdiskcomplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersicDiscComplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersicdisccomplex") {
        sources.spatialtype.push_back(SERSICCOMPLEX);
        nspatialpar = 15;
    } else if (spatialname == "sersicDisk") {
        sources.spatialtype.push_back(SERSICDISK);
        nspatialpar = 7;
    } else if (spatialname == "sersicdisc") {
        sources.spatialtype.push_back(SERSICDISK);
        nspatialpar = 7;
    } else if (spatialname == "sersicDisc") {
        sources.spatialtype.push_back(SERSICDISK);
        nspatialpar = 7;
    } else if (spatialname == "sersicdisk") {
        sources.spatialtype.push_back(SERSICDISK);
        nspatialpar = 7;
    } else if (spatialname == "distortedSphere") {
        sources.spatialtype.push_back(DISTORTEDSPHERE);
        nspatialpar = 10;
    } else if (spatialname == "distortedsphere") {
        sources.spatialtype.push_back(DISTORTEDSPHERE);
        nspatialpar = 10;
    } else if (spatialname == "movingpoint") {
        sources.spatialtype.push_back(MOVINGPOINT);
        nspatialpar = 2;
    } else if (spatialname == "pinhole") {
        sources.spatialtype.push_back(PINHOLE);
        nspatialpar = 4;
    } else if (spatialname == "sphere") {
        sources.spatialtype.push_back(SPHERE);
        nspatialpar = 1;
    } else if (spatialname == "reflectedSphere") {
        sources.spatialtype.push_back(REFLECTEDSPHERE);
        nspatialpar = 3;
    } else if (spatialname == "reflectedsphere") {
        sources.spatialtype.push_back(REFLECTEDSPHERE);
        nspatialpar = 3;
    } else if (spatialname == "disk") {
        sources.spatialtype.push_back(DISK);
        nspatialpar = 1;
    } else if (spatialname == "movingsphere") {
        sources.spatialtype.push_back(MOVINGSPHERE);
        nspatialpar = 3;
    } else if (spatialname == "movingSphere") {
        sources.spatialtype.push_back(MOVINGSPHERE);
        nspatialpar = 3;
    } else if (spatialname == "movingdisk") {
        sources.spatialtype.push_back(MOVINGDISK);
        nspatialpar = 3;
    } else if (spatialname == "movingDisk") {
        sources.spatialtype.push_back(MOVINGDISK);
        nspatialpar = 3;
    } else if (spatialname == "movingreflectedsphere") {
        sources.spatialtype.push_back(MOVINGREFLECTEDSPHERE);
        nspatialpar = 5;
    } else if (spatialname == "movingReflectedSphere") {
        sources.spatialtype.push_back(MOVINGREFLECTEDSPHERE);
        nspatialpar = 5;
    } else {
        sources.spatialtype.push_back(POINT);
        nspatialpar = 0;
    }
        
    if (nspatialpar > 0) {
        for (int i = 0; i < nspatialpar; i++) {
            iss >> sources.spatialpar[nsource][i];
        }
    }

    sources.skysameas.push_back(-1);

    if (sources.spatialtype[nsource] == 1) {
        if (nsource > 0) {
            sources.skysameas[nsource] = -1;
            for (int i = 0; i < nspatial; i++) {
                if (spatialname == spatiallist[i]) sources.skysameas[nsource] = i;
            }
            if (sources.skysameas[nsource] == -1) {
                spatiallist.push_back(spatialname);
                nspatial++;
            }
        } else {
            sources.skysameas[nsource] = -1;
            spatiallist.push_back(spatialname);
            nspatial++;
        }


        // read image file

        if (sources.skysameas[nsource] == -1) {
            snprintf(tempstring, sizeof(tempstring), "%s/%s+0", imagedir.c_str(), spatialname.c_str());

            ffptr = tempstring;
            status = 0;

            if (fits_open_file(&faptr, ffptr, READONLY, &status)) printf("Error opening %s\n", ffptr);
            fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
	    ffgkyd(faptr, "CD1_1", &cd11, NULL, &tstat);
	    ffgkyd(faptr, "CD2_1", &cd21, NULL, &tstat);
	    ffgkyd(faptr, "CD1_2", &cd12, NULL, &tstat);
	    ffgkyd(faptr, "CD2_2", &cd22, NULL, &tstat);
            tempptr[nimage] = static_cast<float*>(malloc(naxes[0]*naxes[1]*sizeof(float)));
            naxesb[nimage][0] = naxes[1];
            naxesb[nimage][1] = naxes[0];
            fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval,
                          tempptr[nimage], &anynull, &status);
	    
            fits_close_file(faptr, &status);


            cumulative[nimage] = 0;
            cumulativex[nimage] = static_cast<double*>(malloc(naxesb[nimage][0]*sizeof(double)));
            for (int i = 0; i < naxesb[nimage][0]; i++) {
                cumulativex[nimage][i] = cumulative[nimage];
                for (int j = 0; j < naxesb[nimage][1]; j++) {
                    if (*(tempptr[nimage] + i*naxesb[nimage][1] + j) < 0) {
                        *(tempptr[nimage] + i*naxesb[nimage][1] + j) = 0;
                    }
                    cumulative[nimage] += *(tempptr[nimage] + i*naxesb[nimage][1] + j);
                }
            }

            sources.spatialpar[nsource][2] = nimage;
            nimage++;
	    //If CD keywords are present and scale is -1 use CD info to set scale and rotation
	    if(tstat == 0 && sources.spatialpar[nsource][0] == -1)
	    {
		
		sources.spatialpar[nsource][0] = sqrt( cd11*cd11 + cd21*cd21) * 3600;
		sources.spatialpar[nsource][1] = 180 - atan2(cd21, cd11)*180/PI;
		

	    }
	    //If CD keywords are absent and scale is -1 use default i.e scale = 1 and rot = 0
	    if(tstat != 0 && sources.spatialpar[nsource][0] == -1)
	    {
		
		sources.spatialpar[nsource][0] = 1;
		sources.spatialpar[nsource][1] = 0;
		printf("Error in reading CD Matrix. Using default values \n");
		
	    }
        }

    }



    std::string dustnamez;
    iss >> dustnamez;

    ndustparz = 0;
    if (dustnamez == "ccm") {
        sources.dusttypez.push_back(1);
        ndustparz = 2;
    } else if (dustnamez == "calzetti") {
        sources.dusttypez.push_back(2);
        ndustparz = 2;
    } else if (dustnamez == "CCM") {
        sources.dusttypez.push_back(1);
        ndustparz = 2;
    } else if (dustnamez == "CALZETTI") {
        sources.dusttypez.push_back(2);
        ndustparz = 2;
    } else {
        sources.dusttypez.push_back(0);
    }

    if (ndustparz > 0) {
        for (int i = 0; i < ndustparz; i++) {
            iss >> sources.dustparz[nsource][i];
        }
    }

    std::string dustname;
    iss >> dustname;

    ndustpar = 0;
    if (dustname == "ccm") {
        sources.dusttype.push_back(1);
        ndustpar = 2;
    } else if (dustname == "calzetti") {
        sources.dusttype.push_back(2);
        ndustpar = 2;
    } else if (dustname== "CCM") {
        sources.dusttype.push_back(1);
        ndustpar = 2;
    } else if (dustname == "CALZETTI") {
        sources.dusttype.push_back(2);
        ndustpar = 2;
    } else {
        sources.dusttype.push_back(0);
    }

    if (ndustpar > 0) {
        for (int i = 0; i < ndustpar; i++) {
            iss >> sources.dustpar[nsource][i];
        }
    }

    int stepRotation = rotationPoints;
    if (rotationTracking >= 1) stepRotation=1;
    for (int i = 0; i < stepRotation; i++) {
        double localtai = tai - vistime/2.0/3600.0/24.0 +
            vistime/3600./24.0*(static_cast<double>(i))/(static_cast<double>(stepRotation - 1));
        if (stepRotation == 1) localtai = tai;

        double spra, spdec;
        spra=pra;
        spdec=pdec;
        double newrotatez=rotatez;
        if (rotationTracking == 0) {
            double lstLocal = localSiderealTime(localtai, longitude/DEGREE);
            double newra, newdec;
            altAztoRaDec(lstLocal, latitude/DEGREE, altitude/DEGREE, azimuth/DEGREE, &newra, &newdec);
            double ha = lstLocal*15.0 - newra;
            newrotatez=calculateRotSkyPos(ha, latitude/DEGREE, newdec/DEGREE, spiderangle/DEGREE);
            spra=newra*DEGREE;
            spdec=newdec*DEGREE;
            newrotatez=-newrotatez*DEGREE;
        }
        if (aberration == 1) {
            double xprime, yprime, zprime;
            aberrationAxes(localtai, &xprime, &yprime, &zprime);
            aberrationShift(&spra, &spdec, xprime, yprime, zprime);
        }
        if (precession == 1) {
            precessionShift(&spra, &spdec, localtai, 51544.5);
        }
        if (nutation ==1) {
            double dlong, dobl, obl;
            nutationValues(localtai,&dlong,&dobl);
            obl=obliquity(localtai);
            nutationShift(&spra,&spdec,obl,dobl,dlong);
        }
        setup_tangent(spra, spdec, &tpx, &tpy, &tpz);

        ra = sources.ra[nsource] + deltara;
        dec = sources.dec[nsource] + deltadec;
        if (aberration == 1) {
            double xprime, yprime, zprime;
            aberrationAxes(localtai, &xprime, &yprime, &zprime);
            aberrationShift(&ra, &dec, xprime, yprime, zprime);
        }
        if (precession == 1) {
            precessionShift(&ra, &dec, localtai, 51544.5);
        }
        if (nutation ==1) {
            double dlong, dobl, obl;
            nutationValues(localtai,&dlong,&dobl);
            obl=obliquity(localtai);
            nutationShift(&ra,&dec,obl,dobl,dlong);
        }
        tangent(ra, dec, &x, &y, &tpx, &tpy, &tpz);

        sources.vx.push_back(x*cos(newrotatez) - y*sin(newrotatez) + rotatex);
        sources.vy.push_back(x*sin(newrotatez) + y*cos(newrotatez) + rotatey);
        sources.vz.push_back(-1.0);
        nn = sqrt((sources.vx[nsourcerotation])*(sources.vx[nsourcerotation]) +
                  (sources.vy[nsourcerotation])*(sources.vy[nsourcerotation]) + 1);
        sources.vx[nsourcerotation] = sources.vx[nsourcerotation]/nn;
        sources.vy[nsourcerotation] = sources.vy[nsourcerotation]/nn;
        sources.vz[nsourcerotation] = sources.vz[nsourcerotation]/nn;
        nsourcerotation++;
    }
        
    nsource++;
    return(0);
}

int Observation::addOpd (const std::string & opd) {

    double x, y;
    double nn;
    double ra, dec;
    double wave;
    std::string id;


    std::istringstream iss(opd);
    iss >> id >> ra >> dec;
    sources.id.push_back(id);
    sources.split.push_back(0);
    sources.ra.push_back(ra);
    sources.dec.push_back(dec);
    sources.backgroundAngle.push_back(0.0);
    sources.ra[nsource] *= DEGREE;
    sources.dec[nsource] *= DEGREE;
    iss >> wave;
    sources.redshift.push_back(0.0);
    sources.gamma1.push_back(wave);
    sources.gamma2.push_back(0.0);
    sources.kappa.push_back(0.0);
    sources.type.push_back(9);
    sources.norm.push_back(opdsize*opdsize*opdsampling*opdsampling);
    sources.mag.push_back(opdsize*opdsize*opdsampling*opdsampling);
    sources.spatialtype.push_back(OPD);
    sources.dusttypez.push_back(0);
    sources.dusttype.push_back(0);
    sources.sedptr.push_back(0);
    sources.skysameas.push_back(-1);
    sedPtr.push_back(0);
    sedN.push_back(0);
    sedMode.push_back(0);
    sedCorr.push_back(0.0);
    sedDwdp.push_back(0.0);
    sourceXpos.push_back(0);
    sourceYpos.push_back(0);
    sourcePhoton.push_back(0);
    sources.spatialpar.push_back(std::vector<double>(15, 0));
    sources.dustpar.push_back(std::vector<double>(2, 0));
    sources.dustparz.push_back(std::vector<double>(2, 0));

    setup_tangent(0.0, 0.0, &tpx, &tpy, &tpz);

    tangent(sources.ra[nsource], sources.dec[nsource], &x, &y, &tpx, &tpy, &tpz);

    int stepRotation = rotationPoints;
    if (rotationTracking >= 1) stepRotation=1;
    for (int i = 0; i < stepRotation; i++) {
        sources.vx.push_back(x);
        sources.vy.push_back(y);
        sources.vz.push_back(-1.0);
        nn = sqrt((sources.vx[nsourcerotation])*(sources.vx[nsourcerotation]) +
                  (sources.vy[nsourcerotation])*(sources.vy[nsourcerotation]) + 1);
        sources.vx[nsourcerotation] = sources.vx[nsourcerotation]/nn;
        sources.vy[nsourcerotation] = sources.vy[nsourcerotation]/nn;
        sources.vz[nsourcerotation] = sources.vz[nsourcerotation]/nn;
        nsourcerotation++;
    }
    nsource++;
    return(0);
}


int Observation::background () {

    // Bookkeeping setup
    char tempstring[4096];
    double focalLength = platescale*1e3/DEGREE;
    double aa, f, xp, yp, ra, dec, xv, yv, maxDistanceY, currbuffer, sourceDistanceX, sourceDistanceY, maxDistanceX;
    double dx , dy;
    long ndeci,  nrai, deci, rai;
    long over;
    double dra, dis, cosdis;
    std::string dir;
    double x, y;
    double nn;
    double mag = 100;
    double angularSepRadians, moonMagnitude, moonApparentMagnitude,
        scatter_function;
    double sunApparentMagnitude;
    int btype[9];

    for (int i = 0; i < 9; i++) btype[i]=0;

    if (backgroundMode == 1 || domelight < 100) {

    if (flatdir == 0) {
        dir = datadir + "/sky";
    } else if (flatdir == 1) {
        dir = ".";
    }

    if (spaceMode <= 0) {
        airglow = static_cast<float*>(calloc((airglowScreenSize)*(airglowScreenSize), sizeof(float)));
        snprintf(tempstring, sizeof(tempstring), "airglowscreen_%s.fits.gz", obshistid.c_str());
        fitsReadImage(tempstring, airglow);
    }
    backgroundTypes=0;

    // Calculations setup
    sunApparentMagnitude = -26.74;
    moonApparentMagnitude = -12.73 + 0.026 * (fabs(phaseang/DEGREE)) + (4E-9) * pow(phaseang/DEGREE,  4);

    // Zodiacal light
    double zl; // Zodiacal light luminosity
    double xx;
    double yy;
    double zz;
    double xEcliptic;
    double yEcliptic;
    double zEcliptic;
    double lambdaEcliptic; // Ecliptic longitude
    double betaEcliptic; // Ecliptic latitude
    double eclipticObliquity = 23.4*DEGREE;
    double lambdaSun; // Sun longitude
    double deltaLambda; // Ecliptic longitude - Sun longitude 
    std::string sedfilename;


    // Populate sources
    double cpdec=cos(pdec);
    double spdec=sin(pdec);
    double cca=cos(chipangle);
    double sca=sin(chipangle);
    double crz=cos(-rotatez);
    double srz=sin(-rotatez);
    currbuffer = backBuffer + (2.5*backBeta*backRadius*sqrt(1+1/backDelta/backDelta)
                               + 120)*ARCSEC*focalLength/(pixsize*1e3);
    if (backAlpha<=0.001)  currbuffer = backBuffer + (2.5*1.0*backRadius*sqrt(1+1/backDelta/backDelta) +
                                                      + 120)*ARCSEC*focalLength/(pixsize*1e3);
    maxDistanceX = pixelsx * (pixsize*1e3)/2.0 + currbuffer * (pixsize*1e3);
    maxDistanceY = pixelsy * (pixsize*1e3)/2.0 + currbuffer * (pixsize*1e3);
    double maxDistanceXs = pixelsx * pixsize * 1e3/2.0 + currbuffer * pixsize * 1e3 - 120.0*ARCSEC*focalLength;
    double maxDistanceYs = pixelsy * pixsize * 1e3/2.0 + currbuffer * pixsize * 1e3 - 120.0*ARCSEC*focalLength;

    ndeci = (long)(180*3600/backRadius);
    over = (long)(60/backRadius);
    if (over == 0) over=1;
    for (deci = 0; deci < ndeci; deci += over) {
        dec = ((deci + 0.5 - ndeci/2)/(ndeci))*PI;
        nrai = ((long)(ndeci*cos(dec)*2));

        // Coarse placement
        for (rai = 0; rai < nrai; rai += over) {
            ra = 2*PI*((rai + 0.5 - nrai/2)/(nrai));

            aa = cos(dec)*cos(ra - pra);
            f = focalLength/(spdec*sin(dec) + aa*cpdec);
            yp = flipY*f*(cpdec*sin(dec) - aa*spdec);
            xp = flipX*f*cos(dec)*sin(ra - pra);
            if (swap==1) {
                double tt;
                tt = xp;
                xp = yp;
                yp = tt;
            }
            xv = xp*crz + yp*srz;
            yv = -xp*srz + yp*crz;
            dx = xv - centerx - decenterx;
            dy = yv - centery - decentery;
            sourceDistanceX = fabs(cca*dx + sca*dy);
            sourceDistanceY = fabs(-sca*dx + cca*dy);

            dis = acos(focalLength/f);

            if ((sourceDistanceX <= maxDistanceX) && (sourceDistanceY <= maxDistanceY) &&
                (dis < PI/2) ){

                // Fine placement
                for (int i = 0; i < over; i++) {
                    for (int j = 0; j < over; j++) {
                        dec = ((deci + i + 0.5 - ndeci/2)/ndeci)*PI;
                        ra = 2*PI*((rai + j + 0.5 - nrai/2)/(nrai));
                        aa = cos(dec)*cos(ra - pra);
                        f = focalLength/(spdec*sin(dec) + aa*cpdec);
                        yp = flipY*f*(cpdec*sin(dec) - aa*spdec);
                        xp = flipX*f*cos(dec)*sin(ra - pra);
                        if (swap==1) {
                            double tt;
                            tt = xp;
                            xp = yp;
                            yp = tt;
                        }
                        xv = xp * crz + yp * srz;
                        yv = -xp * srz + yp * crz;
                        dx = xv - centerx - decenterx;
                        dy = yv - centery - decentery;
                        sourceDistanceX = fabs(cca*dx + sca*dy);
                        sourceDistanceY = fabs(-sca*dx + cca*dy);

                        dis = acos(focalLength/f);

                        if ((sourceDistanceX <= maxDistanceXs) && (sourceDistanceY <= maxDistanceYs) &&
                            (dis < PI/2)){

                            // Source type loop
                            for (int diffusetype = 0; diffusetype < 9; diffusetype++) {
                                //                                if (diffusetype!=5) goto skipdiff;
                                // Conditional checks for observational parameters and modes
                                if ((diffusetype == 0 && domelight < 100 && spaceMode == 0) ||
                                    (diffusetype == 1 && airglowcintensity < 100 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 2 && airglowpintensity < 100 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 3 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 4 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 5 && backgroundMode == 1 && telconfig == 0) ||
                                    (diffusetype == 6 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 7 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0) ||
                                    (diffusetype == 8 && backgroundMode == 1 && telconfig == 0 && spaceMode == 0)) {

                                    double backAngle = 0.0;
                                    if (diffusetype == 0) mag = domelight - 2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 1) {
                                        double airglowv;
                                        // Airglow variation
                                        int ax0, ax1, ay0, ay1;
                                        double dx, dy;
                                        find_linear_wrap(xv, platescale*1e3*15.0/3600, airglowScreenSize, &ax0, &ax1, &dx);
                                        find_linear_wrap(yv, platescale*1e3*15.0/3600, airglowScreenSize, &ay0, &ay1, &dy);
                                        airglowv = airglowvariation*
                                            (static_cast<double>(interpolate_bilinear_float_wrap(airglow,airglowScreenSize, ax0, ax1, dx, ay0, ay1, dy))); 
                                        mag = airglowcintensity + airglowv - 2.5*log10(backRadius*backRadius);

                                    }
                                    if (diffusetype == 2) {

                                        mag = airglowpintensity - 2.5*log10(backRadius*backRadius);
                                    }
                                    if ((diffusetype == 3) || (diffusetype == 4)) {
                                        dra = fabs(ra - moonra);
                                        if (dra > PI) dra = 2*PI - dra;
                                        cosdis = sin(dec)*sin(moondec) + cos(dec)*cos(moondec)*cos(dra);
                                        if (cosdis > 1) cosdis = 1.0;
                                        if (cosdis < -1) cosdis = -1.0;
                                        angularSepRadians = acos(cosdis);
                                        backAngle = angularSepRadians;
                                        if (diffusetype == 3) {
                                            scatter_function = - 1000.0;
                                        } else {
                                            scatter_function = - 2.5*log10(4.0*PI) +
                                                2.5*log10(PI/180.0/3600.0*PI/180.0/3600.0) +
                                                2.5*log10(4.0);
                                            //Due factor of 4 oversample
                                        }
                                        moonMagnitude = moonApparentMagnitude - scatter_function;
                                        mag = moonMagnitude - 2.5*log10(backRadius*backRadius);

                                        // effectively no light below 20 degrees
                                        if (moonalt < -20.0*DEGREE) mag=10000.0;

                                    }
                                    if (diffusetype == 5) {

                                        if (overrideZodiacalLightMagnitude == 0) {
                                            // Calculate ecliptic longitude and latitude
                                            // Source position in Cartesian coordinates 
                                            xx = cos(ra) * cos(dec);
                                            yy = sin(ra) * cos(dec);
                                            zz = sin(dec);
                                            // Transform to ecliptic coordinates
                                            xEcliptic = xx;
                                            yEcliptic = yy*cos(eclipticObliquity) + zz*sin(eclipticObliquity);
                                            zEcliptic = -yy*sin(eclipticObliquity) + zz*cos(eclipticObliquity);
                                            // Ecliptic longitue and latitude
                                            lambdaEcliptic = atan2(yEcliptic, xEcliptic)/DEGREE; // Degrees
                                            betaEcliptic = atan(zEcliptic/sqrt(pow(xEcliptic, 2) + pow(yEcliptic, 2)))/DEGREE; // Degrees

                                            // Zodiacal light variation
                                            // Follows Kwon et al. New Astronomy 10-2 (2004) pp91-107
                                            if (lambdaEcliptic < 0) {
                                                lambdaEcliptic += 360.0;
                                            }
                                            lambdaSun = fmod((640.466 - 0.985607*(day - 51544.0)), 360.0); // Degrees
                                            if (lambdaSun < 0) {
                                                lambdaSun += 360.0;
                                            }
                                            deltaLambda = fmod(lambdaEcliptic - lambdaSun, 360.0); // Degrees
                                            if (deltaLambda < 0) {
                                                deltaLambda += 360.0;
                                            }
                                            zl = (140.0 + 20.0/(1 + pow(fabs(180.0 - deltaLambda)/20.0, 2.0)))*
                                                (5.0/14.0 + 9.0/14.0/(1 + pow((betaEcliptic)/20.0, 2.0)))*0.7;
                                            zodiacalLightMagnitude = -2.5*(log10(zl)) + 27.78;

                                        }

                                        mag = zodiacalLightMagnitude - 2.5*log10(backRadius*backRadius);

                                    }
                                    if ((diffusetype == 6) || (diffusetype == 7)) {
                                        dra = fabs(ra - solarra);
                                        if (dra > PI) dra = 2*PI - dra;
                                        cosdis = sin(dec)*sin(solardec) + cos(dec)*cos(solardec)*cos(dra);
                                        if (cosdis > 1) cosdis = 1.0;
                                        if (cosdis < -1) cosdis = -1.0;
                                        angularSepRadians = acos(cosdis);
                                        backAngle = angularSepRadians;

                                        if (diffusetype == 6) {
                                            scatter_function = - 1000.0;
                                        } else {
                                            scatter_function = - 2.5*log10(4.0*PI) +
                                                2.5*log10(PI/180.0/3600.0*PI/180.0/3600.0) +
                                                2.5*log10(4.0);
                                            //factor of 4 oversamples
                                        }
                                        mag = sunApparentMagnitude - scatter_function - 2.5*log10(backRadius*backRadius);

                                        // safety switch
                                        if (sunsafety==1) {
                                            if ((solaralt >= 0) || (solaralt < -20.0*DEGREE)) {
                                                mag=10000.0;
                                                solarfactor=10.0;
                                            }
                                        }

                                        if ((solaralt < 0) && (solaralt >= -5.0*DEGREE)) solarfactor=10.0;
                                        if ((solaralt < -5.0*DEGREE) && (solaralt >= -20.0*DEGREE)) solarfactor=1e6;
                                        mag = mag + 2.5*log10(solarfactor);

                                    }

                                    if (diffusetype==8) {
                                        mag = artificialLightMagnitude - 2.5*log10(backRadius*backRadius) - 1.05;
                                        // 1.05 is approx correction for reflected light
                                    }


                                    if (mag < 100) {

                                        btype[diffusetype]=1;

                                        //nspatialpar = 0;


                                        sources.id.push_back("0.0");
                                        sources.split.push_back(0);
                                        sources.ra.push_back(ra);
                                        sources.dec.push_back(dec);
                                        sources.backgroundAngle.push_back(backAngle);
                                        if (diffusetype == 0 && domewave == 0.0) sedfilename = dir + "/sed_dome.txt";
                                        if (diffusetype == 0 && domewave != 0.0) sedfilename = "laser";
                                        if (diffusetype == 1) sedfilename = dir + "/airglowc_sed.txt";
                                        if (diffusetype == 2) sedfilename = dir + "/airglowp_sed.txt";
                                        if (diffusetype == 3) sedfilename = dir + "/lunar_sed.txt";
                                        if (diffusetype == 4) sedfilename = dir + "/lunar_sed.txt";
                                        if (diffusetype == 5) sedfilename = dir + "/zodiacal_sed.txt";
                                        if (diffusetype == 6) sedfilename = dir + "/solar_sed.txt";
                                        if (diffusetype == 7) sedfilename = dir + "/solar_sed.txt";
                                        if (diffusetype == 8) sedfilename = dir + "/artificial.txt";
                                        sources.redshift.push_back(0.0);
                                        sources.gamma1.push_back(0.0);
                                        sources.gamma2.push_back(0.0);
                                        sources.kappa.push_back(0.0);
                                        sources.type.push_back(diffusetype);
                                        sources.norm.push_back(pow(10.0, ((mag + 48.6)/(-2.5))));
                                        sources.mag.push_back(mag);
                                        sources.skysameas.push_back(-1);
                                        sedPtr.push_back(0);
                                        sedN.push_back(0);
                                        sedMode.push_back(0);
                                        sedCorr.push_back(0.0);
                                        sedDwdp.push_back(0.0);
                                        sourceXpos.push_back(0);
                                        sourceYpos.push_back(0);
                                        sourcePhoton.push_back(0);
                                        sources.spatialpar.push_back(std::vector<double>(15, 0));
                                        sources.dustpar.push_back(std::vector<double>(2, 0));
                                        sources.dustparz.push_back(std::vector<double>(2, 0));

                                        // Read SED file
                                        if (nsource > 0) {
                                            for (int ii = 0; ii < nsed; ii++) {
                                                if (sedlist[ii] == sedfilename) {
                                                    sources.sedptr.push_back(sources.sedptr[ii]);
                                                    goto skipsedread;
                                                }
                                            }
                                        }

                                        sedlist.push_back(sedfilename);
                                        nsed++;
                                        sources.sedptr.push_back(nsedptr);
                                        readSed(sedfilename, 1);

                                    skipsedread:;
                                        sources.norm[nsource] = sources.norm[nsource]/(normwave)*(1 + sources.redshift[nsource])*
                                            sedDwdp[sources.sedptr[nsource]];

                                        sources.spatialtype.push_back(2);
                                        sources.spatialpar[nsource][0] = backRadius*backEpsilon;
                                        sources.dusttypez.push_back(0);
                                        sources.dusttype.push_back(0);

                                        setup_tangent(pra, pdec, &tpx, &tpy, &tpz);

                                        tangent(sources.ra[nsource], sources.dec[nsource],
                                                &x, &y, &tpx, &tpy, &tpz);

                                        int stepRotation = rotationPoints;
                                        if (rotationTracking >= 1) stepRotation=1;
                                        for (int k = 0; k < stepRotation; k++) {
                                            sources.vx.push_back(x*cos(rotatez) - y*sin(rotatez) + rotatex);
                                            sources.vy.push_back(x*sin(rotatez) + y*cos(rotatez) + rotatey);
                                            sources.vz.push_back(-1.0);
                                            nn = sqrt((sources.vx[nsourcerotation])*(sources.vx[nsourcerotation]) +
                                                      (sources.vy[nsourcerotation])*(sources.vy[nsourcerotation]) + 1);
                                            sources.vx[nsourcerotation] = sources.vx[nsourcerotation]/nn;
                                            sources.vy[nsourcerotation] = sources.vy[nsourcerotation]/nn;
                                            sources.vz[nsourcerotation] = sources.vz[nsourcerotation]/nn;

                                        //                                        if (backgroundProbability > 0) {
                                            if (sources.vx[nsourcerotation] < bpVxMin) bpVxMin = sources.vx[nsourcerotation];
                                            if (sources.vx[nsourcerotation] > bpVxMax) bpVxMax = sources.vx[nsourcerotation];
                                            if (sources.vy[nsourcerotation] < bpVyMin) bpVyMin = sources.vy[nsourcerotation];
                                            if (sources.vy[nsourcerotation] > bpVyMax) bpVyMax = sources.vy[nsourcerotation];
                                            //}
                                            nsourcerotation++;
                                        }
                                        nsource++;

                                    }
                                    //                                                                   skipdiff:;
                                }


                            }
                        }
                    }

                }
            }
        }
    }

   }
    for (int i = 0; i < 9; i++) if (btype[i]==1) backgroundTypes++;
    return(0);
}


int Observation::filterTruncateSources () {

    double filterlow,  filterhigh;
    long lowptr, highptr;
    double lowvalue, highvalue;

    filterlow = 0;
    filterhigh = maxwavelength;

    for (int i = 0; i < nsedptr; i++) {
        if (sedMode[i] == 1) {
            filterlow = minwavelength;
        } else {
            filterlow = 0;
        }
        highvalue = 0;
        lowvalue = 0;
        lowptr = 0;
        highptr = 0;
        for (int j = 0; j < sedN[i]; j++) {
            if (sedW[sedPtr[i] + j] < filterlow) {
                lowptr = j;
                lowvalue = sedC[sedPtr[i] + j];
            }
            if (sedW[sedPtr[i] + j] < filterhigh) {
                highptr = j;
                highvalue = sedC[sedPtr[i] + j];
            }

        }
        sedCorr[i]=(sedC[sedPtr[i] + highptr] - sedC[sedPtr[i] + lowptr]);

        for (int j = lowptr; j < highptr + 1; j++) {
            sedC[sedPtr[i] + j]=(sedC[sedPtr[i] + j] - lowvalue)/
                (highvalue - lowvalue);
        }
        sedPtr[i] = sedPtr[i] + lowptr;
        sedN[i] = highptr - lowptr;
    }

    for (int i = 0; i < nsource; i++) {
        if (sources.spatialtype[i] != OPD) {
            sources.norm[i] = sources.norm[i]*sedCorr[sources.sedptr[i]];
        }
    }

    return(0);

}


int Observation::splitSources () {

    int count = 0;

    if (numthread > 1) {

        while (1) {
            int brightest = -1;
            double bright = 0.0;
            for (int i = 0; i < nsource; i++) {
                if ((sources.norm[i] > bright) && (sources.type[i] == 9) && (sources.spatialtype[i] != OPD)  && (sources.split[i]==0)) {
                    bright = sources.norm[i];
                    brightest = i;
                }
            }

            if (brightest == -1) break;
            int totalbrightsplit = 0;
            for (int i = 0; i < nsource; i++) {
                if ((sources.norm[i] > bright) && (sources.split[i]==1)) totalbrightsplit++;
            }
            if (totalbrightsplit >= numthread*2) break;

            sources.split[brightest]=1;
            sources.norm[brightest]=sources.norm[brightest]/(static_cast<double>(numthread));
            sources.mag[brightest]=sources.mag[brightest]+2.5*log10(numthread);

            for (int j = 0; j < (numthread-1); j++) sources.ra.push_back(sources.ra[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.dec.push_back(sources.dec[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.backgroundAngle.push_back(sources.backgroundAngle[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.redshift.push_back(sources.redshift[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.gamma1.push_back(sources.gamma1[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.gamma2.push_back(sources.gamma2[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.kappa.push_back(sources.kappa[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.split.push_back(sources.split[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.id.push_back(sources.id[brightest]);
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k < rotationPoints; k++) sources.vx.push_back(sources.vx[brightest*rotationPoints + k]);
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k < rotationPoints; k++) sources.vy.push_back(sources.vy[brightest*rotationPoints + k]);
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k < rotationPoints; k++) sources.vz.push_back(sources.vz[brightest*rotationPoints + k]);
            for (int j = 0; j < (numthread-1); j++) sources.mag.push_back(sources.mag[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.norm.push_back(sources.norm[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.spatialtype.push_back(sources.spatialtype[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.dusttype.push_back(sources.dusttype[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.dusttypez.push_back(sources.dusttypez[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.type.push_back(sources.type[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.skysameas.push_back(sources.skysameas[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.sedptr.push_back(sources.sedptr[brightest]);
            for (int j = 0; j < (numthread-1); j++) sedPtr.push_back(sedPtr[brightest]);
            for (int j = 0; j < (numthread-1); j++) sedN.push_back(sedN[brightest]);
            for (int j = 0; j < (numthread-1); j++) sedMode.push_back(sedMode[brightest]);
            for (int j = 0; j < (numthread-1); j++) sedCorr.push_back(sedCorr[brightest]);
            for (int j = 0; j < (numthread-1); j++) sedDwdp.push_back(sedDwdp[brightest]);
            for (int j = 0; j < (numthread-1); j++) sourceXpos.push_back(sourceXpos[brightest]);
            for (int j = 0; j < (numthread-1); j++) sourceYpos.push_back(sourceYpos[brightest]);
            for (int j = 0; j < (numthread-1); j++) sourcePhoton.push_back(sourcePhoton[brightest]);
            for (int j = 0; j < (numthread-1); j++) sources.spatialpar.push_back(std::vector<double>(15, 0));
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k< 15; k++) {
                    sources.spatialpar[nsource+j][k]=sources.spatialpar[brightest][k];
                }
            for (int j = 0; j < (numthread-1); j++) sources.dustpar.push_back(std::vector<double>(2, 0));
            for (int  j = 0; j < (numthread-1); j++) for (int k = 0; k< 2; k++) sources.dustpar[nsource+j][k]=sources.dustpar[brightest][k];
            for (int j = 0; j < (numthread-1); j++) sources.dustparz.push_back(std::vector<double>(2, 0));
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k< 2; k++) sources.dustparz[nsource+j][k]=sources.dustparz[brightest][k];

            for (int j = 0; j < (numthread-1); j++) nsource++;
            for (int j = 0; j < (numthread-1); j++) for (int k = 0; k < rotationPoints; k++) nsourcerotation++;

            count++;
            if (count > 100000) {
                std::cout << "Warning:  Large number of split sources." << std::endl;
            }
        }

    }
    return(0);

}
