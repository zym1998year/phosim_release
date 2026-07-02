///
/// @package phosim
/// @file trim.cpp
/// @brief trim program: removes sources that have no chance of producing photons on a chip
///
/// @brief Created by:
/// @author Alan Meert (Purdue)
///
/// @brief Modified by:
/// @author Justin Bankert (Purdue)
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <fstream>

#include "trim/trim.h"
#include "raytrace/constants.h"
#include "ancillary/readtext.h"

/// @addtogroup trim_group trim
/// @{

void Trim::xyPosition(double alpha, double delta, double *x, double *y, double *dis) {

    double a = cos(delta)*cos(alpha - pointingRa);
    double f = focalLength/(sin(pointingDec)*sin(delta) + a*cos(pointingDec));
    double yp = flipY*f*(cos(pointingDec)*sin(delta) - a*sin(pointingDec));
    double xp = flipX*f*cos(delta)*sin(alpha - pointingRa);
    if (swap==1) {
        double tt;
        tt = xp;
        xp = yp;
        yp = tt;
    }
    *x = (xp*cos(rotatez) + yp*sin(rotatez));
    *y = (-xp*sin(rotatez) + yp*cos(rotatez));
    double cosdis = focalLength/f;
    if (cosdis > 1) cosdis = 1.0;
    if (cosdis < -1) cosdis = -1.0;
    *dis = acos(cosdis);
}

void Trim::readOpdCatalog() {

    FILE *fp;
    FILE *fp2;
    char line[4096];
 
    snprintf(outputOpdFilename, sizeof(outputOpdFilename), "trimcatalog_%s_opd.pars", obshistid.c_str());
    fp2 = fopen(outputOpdFilename, "wt");
    if (fp2 == NULL) {
        std::cout << "Cannot open for writing output file " << outputOpdFilename << "\n.";
        exit(1);
    }
    fprintf(fp2, "\n");
    fprintf(fp2, "\n");

    fp = fopen(catalogOpd.c_str(), "rt");
    if (fp == NULL){
        std::cout << "Cannot open catalog file " << catalogOpd << "." << std::endl;
        exit(1);
    }
    long totalSourceCount = 0;
    while (fgets(line, 4096, fp)) {
        fprintf(fp2, "%s\n", line);
        totalSourceCount++;
    }
    fclose(fp);
    fclose(fp2);

    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Trim Catalog" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Found " << totalSourceCount << " opd source(s)." << std::endl;


}

void Trim::readCatalog() {

    gzFile gzfp[MAX_CATALOG];
    FILE *fp[MAX_CATALOG];
    FILE *fp2[MAX_CHIP];
    char line[4096];
    long sourceCount[MAX_CHIP];
    char ignore[3][4096];
    double ra, dec, mag;
    double x, y, dis;
    double dx, dy, xp, yp;
    double currentBuffer;

    for (int d = 0; d < nChip; d++) {
        snprintf(outputFilename[d], sizeof(outputFilename[d]), "trimcatalog_%s_%s.pars", obshistid.c_str(), chipid[d].c_str());
    }

    for (int d = 0; d < nChip; d++) {
        fp2[d] = fopen(outputFilename[d], "wt");
        if (fp2[d]==NULL) {
            std::cout << "Cannot open for writing output file " << outputFilename[d] << "\n.";
            exit(1);
        }
        fprintf(fp2[d], "\n");
        fprintf(fp2[d], "\n");
        sourceCount[d] = 0;
    }

    for (int c = 0; c < nCatalog; c++) {
        if (strstr(catalog[c].c_str(), ".gz")==NULL) {

            fp[c] = fopen(catalog[c].c_str(), "rt");
            if (fp[c] == NULL){
                std::cout << "Cannot open catalog file " << catalog[c] << "." << std::endl;
                exit(1);
            }
            while (fgets( line, 4096, fp[c] )) {
                sscanf(line, "%s %s %lf %lf %lf %s", ignore[0], ignore[2], &ra, &dec, &mag, ignore[1]);
                ra = ra*DEGREE;
                dec = dec*DEGREE;
                for (int d = 0; d < nChip; d++) {
                    currentBuffer = (buffer*pixelSize[d] + scale*pow(2.5, 17 - mag));
                    if (strayLight == 0 && currentBuffer > extendedBuffer + buffer*pixelSize[d]) currentBuffer = extendedBuffer + buffer*pixelSize[d];
                    xyPosition(ra, dec, &x, &y, &dis);
                    dx = x - xPosition[d] - deltaX[d];
                    dy = y - yPosition[d] - deltaY[d];
                    xp = cos(angle[d])*dx + sin(angle[d])*dy;
                    yp = -sin(angle[d])*dx + cos(angle[d])*dy;
                    if ((fabs(xp) <= xDimension[d] + currentBuffer) &&
                        (fabs(yp) <= yDimension[d] + currentBuffer) &&
                        (dis < PI/2)){
                        fprintf( fp2[d], "%s\n", line);
                        sourceCount[d]++;
                    }
                }
            }
            fclose(fp[c]);

        } else {

            gzfp[c] = gzopen( catalog[c].c_str(), "r" );
            if (gzfp[c] == NULL) {
                std::cout << "Cannot open catalog file " << catalog[c] << "." << std::endl;
                exit(1);
            }
            while (gzgets( gzfp[c] , line, 4096)) {
                sscanf(line, "%s %s %lf %lf %lf %s", ignore[0], ignore[2], &ra, &dec, &mag, ignore[1]);
                ra = ra*DEGREE;
                dec = dec*DEGREE;
                for (int d = 0; d < nChip; d++) {
                    currentBuffer = (buffer*pixelSize[d] + scale*pow(2.5, 17 - mag));
                    if (strayLight == 0 && currentBuffer > extendedBuffer + buffer*pixelSize[d]) currentBuffer = extendedBuffer + buffer*pixelSize[d];
                    xyPosition(ra, dec, &x, &y, &dis);
                    dx = x - xPosition[d] - deltaX[d];
                    dy = y - yPosition[d] - deltaY[d];
                    xp = cos(angle[d])*dx + sin(angle[d])*dy;
                    yp = -sin(angle[d])*dx + cos(angle[d])*dy;
                    if ((fabs(xp) <= xDimension[d] + currentBuffer) &&
                        (fabs(yp) <= yDimension[d] + currentBuffer) &&
                        (dis < PI/2)){
                        if (sqrt(xp*xp+yp*yp) < sqrt(xDimension[d]*xDimension[d] +
                                                     yDimension[d]*yDimension[d]) + currentBuffer) {
                            fprintf(fp2[d], "%s\n", line);
                            sourceCount[d]++;
                        }
                    }
                }
            }
            gzclose(gzfp[c]);

        }
    }

    for (int d = 0; d < nChip; d++) {
        fclose(fp2[d]);
    }

    long totalSourceCount = 0;
    for (int d = 0; d < nChip; d++) {
        totalSourceCount += sourceCount[d];
    }
    if (totalSourceCount >= minSource) {
        std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "Trim Catalog" << std::endl;
        std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
        for (int d = 0; d < nChip; d++) {
            if (sourceCount[d] >= minSource) {
                std::cout << "Found " << sourceCount[d] << " source(s) for sensor " << chipid[d] << "." << std::endl;
            }
        }
    }

}

void Trim::getDetectorProperties( int d ){

    std::istringstream focalplanePars(readText::get(instrdir + "/focalplanelayout.txt", chipid[d]));
    focalplanePars >> xPosition[d] >> yPosition[d] >> pixelSize[d] >> xDimension[d] >> yDimension[d];
    xDimension[d] *= pixelSize[d]/2.0;
    yDimension[d] *= pixelSize[d]/2.0;
    std::string tempstring;
    double temp;
    focalplanePars >> tempstring;
    focalplanePars >> tempstring;
    focalplanePars >> tempstring;
    focalplanePars >> temp;
    focalplanePars >> temp;
    focalplanePars >> tempstring;
    focalplanePars >> angle[d];
    angle[d] *= DEGREE;
    focalplanePars >> temp;
    focalplanePars >> temp;
    focalplanePars >> deltaX[d];
    focalplanePars >> deltaY[d];
    deltaX[d] *= 1000.0;
    deltaY[d] *= 1000.0;
}

void Trim::setup() {

    instrdir = "../data/lsst";
    flatDirectory = 0;

    readText pars(std::cin);

    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "pointingra", pointingRa);
        readText::get(line, "pointingdec", pointingDec);
        readText::get(line, "rotatez", rotatez);
        readText::get(line, "rotationangle", rotatez); //deprecated
        readText::get(line, "filter", filter);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "flatdir", flatDirectory);
        readText::get(line, "minsource", minSource);
        readText::get(line, "chipid", chipid);
        readText::get(line, "catalog", catalog);
        readText::get(line, "catalogopd", catalogOpd);
    }
    pointingRa *= DEGREE;
    pointingDec *= DEGREE;
    rotatez *= DEGREE;

    nCatalog = catalog.size();
    nChip = chipid.size();
    if (flatDirectory == 1) instrdir = ".";
    if (chipid[0] == "opd") {
        opdMode = 1;
        if (catalogOpd.size()==0) {
            FILE *fp2;
            snprintf(outputOpdFilename, sizeof(outputOpdFilename), "trimcatalog_%s_opd.pars", obshistid.c_str());
            fp2 = fopen(outputOpdFilename, "wt");
            if (fp2 == NULL) {
                std::cout << "Cannot open for writing output file " << outputOpdFilename << "\n.";
                exit(1);
            }
            fprintf(fp2, "\n");
            fprintf(fp2, "\n");
            exit(0);
        }
    }

    flipX = 1;
    flipY = 1;
    swap = 0;

    readText wavelengthPars(instrdir + "/optics_" + std::to_string(filter) + ".txt");
    int haveOptics = 0;
    double maxr, minr;
    int first = 0;
    for (size_t t(0); t < wavelengthPars.getSize(); t++) {
        std::string line(wavelengthPars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;
        if ((haveOptics==1) && (first==0)) {
            std::string nameT;
            double rcurv, dz;
            iss >> nameT  >> rcurv >> dz >> maxr >> minr;
            first++;
        }
        if (keyName == "optics") {
            haveOptics = 1;
            std::string nameT;
            double wavelength, minW, maxW;
            iss >> nameT >> minW >> maxW >> wavelength >> plateScale >> flipX >> flipY >> swap;
        }
    }
    if (haveOptics == 0) {
        std::cout << "Error:  Need optics approximation information on first line of optics file." << std::endl;
        exit(1);
    }

    buffer = 20;
    strayLight = 1;
    focalLength = plateScale/DEGREE;
    scale = (0.2*ARCSEC/DEGREE)*(focalLength*DEGREE)*0.0;
    magconst = pow(100.0, 0.2);
    maglimit = 17.0 + 2.5*log10((maxr*maxr-minr*minr)/(4180.0*4180.0-2558.0*2558.0));
    extendedBuffer = (60.0*ARCSEC/DEGREE)*(focalLength*DEGREE)/6.0;

    if (nChip > MAX_CHIP || nCatalog > MAX_CATALOG ){
        std::cout << "Cannot split into that many files." << std::endl;
        exit(1);
    }

}


/// @}
