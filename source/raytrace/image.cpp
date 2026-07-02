///
/// @package phosim
/// @file image.cpp
/// @brief image class
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

#include "constants.h"
#include "parameters.h"
#include "ancillary/random.h"
#include "helpers.h"
#include "image.h"
#include "ext.h"

#include "atmospheresetup.cpp"
#include "telescopesetup.cpp"
#include "photonmanipulate.cpp"
#include "photonoptimization.cpp"
#include "photonloop.cpp"
#include "particleloop.cpp"
#include "sourceloop.cpp"
#include "../ancillary/fits.h"

void Image::writeImageFile () {

    int status = 0;
	double cd11, cd12, cd21, cd22, temp_cd;
	cd11 = pixsize/platescale*cos(rotatez-chipangle)* flipX;
	cd12 = pixsize/platescale*sin(rotatez-chipangle)*flipY;
	cd21 = -pixsize/platescale*sin(rotatez-chipangle)*flipX;
	cd22 = pixsize/platescale*cos(rotatez-chipangle)*flipY;
    char tempstring[4096], tempstring2[4096], line[4096];
    fitsfile *fptr = NULL;
    FILE *indafile;
	if(swap == 1)
	{
		temp_cd = cd12;
		cd12 = cd11;
		cd11 = temp_cd;
		temp_cd=cd21;
		cd21 = cd22;
		cd22 = temp_cd;
	}
    

    std::string filename = "!" + workdir + "/" + outputfilename + ".fits.gz";
    fitsCreateImage(&fptr, filename.c_str());
	
    fitsWriteKey(fptr, "CTYPE1", "RA---TAN", "");
    fitsWriteKey(fptr, "CRPIX1", -(centerx+decenterx)/(pixsize*1e3)*cos(chipangle)
                 -(centery+decentery)/(pixsize*1e3)*sin(chipangle) + pixelsx/2, "");
    fitsWriteKey(fptr, "CRVAL1", pra/DEGREE, "");
    fitsWriteKey(fptr, "CTYPE2", "DEC--TAN", "");
    fitsWriteKey(fptr, "CRPIX2", -(centerx+decenterx)/(pixsize*1e3)*(-sin(chipangle))-
                 (centery+decentery)/(pixsize*1e3)*(cos(chipangle))+ pixelsy/2, "");
    fitsWriteKey(fptr, "CRVAL2", pdec/DEGREE, "");
    fitsWriteKey(fptr, "CD1_1", cd11, "");
    fitsWriteKey(fptr, "CD1_2", cd12, "");
    fitsWriteKey(fptr, "CD2_1", cd21, "");
    fitsWriteKey(fptr, "CD2_2", cd22, "");
    fitsWriteKey(fptr, "RADESYS", "ICRS", "");
    fitsWriteKey(fptr, "EQUINOX", 2000.0, "");
    header(fptr);
    fitsWriteKey(fptr, "CREATOR", "PHOSIM", "");

    snprintf(tempstring, sizeof(tempstring), "%s/version", bindir.c_str());
    indafile = fopen(tempstring, "r");
    fgets(line, 4096, indafile);
    sscanf(line, "%s %s", tempstring, tempstring2);

    fitsWriteKey(fptr, "VERSION", tempstring2, "");

    fgets(line, 4096, indafile);
    sscanf(line, "%s", tempstring2);
    fclose(indafile);

    fitsWriteKey(fptr, "BRANCH", tempstring2, "");

    if (date) fits_write_date(fptr, &status);

    state.focalPlaneFloat = static_cast<float*>(calloc((chip.extraX)*(chip.extraY), sizeof(float)));
    for (long i = 0; i < chip.extraX; i++) {
        for (long j = 0; j < chip.extraY; j++) {
            *(state.focalPlaneFloat + chip.extraX*j + i) =
                static_cast<float>(*(state.focalPlane + chip.extraX*j + i));
        }
    }
    fitsWriteImage(fptr, chip.extraX, chip.extraY, state.focalPlaneFloat);

}

void Image::writeOPD () {

#ifdef EXT
#include "ext.cpp"
#endif

}

void Image::writeCheckpoint(int checkpointcount) {

    fitsfile *faptr;
    long naxes[2];
    int status;
    double *tempDynamicTransmission;
    Uint32 z = 0, w = 0, zg = 0, wg = 0;
    char tempstring[4096];

    tempDynamicTransmission = static_cast<double*>(calloc(((natmospherefile+1)*6 + nsurf*2 + 3)*(maxwavelength - minwavelength + 1), sizeof(double)));
    for (int i = 0; i < maxwavelength - minwavelength + 1; i++) {
        for (int j = 0; j < ((natmospherefile+1)*6 + nsurf*2 + 3); j++) {
            tempDynamicTransmission[i*((natmospherefile+1)*6 + nsurf*2 + 3) + j] =
                state.dynamicTransmission[i*((natmospherefile+1)*6 + nsurf*2 + 3) + j];
        }
    }

    status = 0;
    std::ostringstream filename;
    filename << "!"<< workdir << "/" << outputfilename << "_ckptdt_" << checkpointcount << ".fits.gz";
    fits_create_file(&faptr, filename.str().c_str(), &status);
    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_img(faptr, DOUBLE_IMG, 2, naxes, &status);
    naxes[0] = maxwavelength - minwavelength + 1;
    naxes[1] = (natmospherefile+1)*6 + nsurf*2 + 3;
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    for (int i = 0; i < numthread; i++) {
        random[i].getSeed(&z, &w);
        galaxy.random[i].getSeed(&zg, &wg);
        snprintf(tempstring, sizeof(tempstring), "M_Z%4d", i);
        fits_update_key(faptr, TUINT, tempstring, &z, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "M_W%4d", i);
        fits_update_key(faptr, TUINT, tempstring, &w, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "G_Z%4d", i);
        fits_update_key(faptr, TUINT, tempstring, &zg, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "G_W%4d", i);
        fits_update_key(faptr, TUINT, tempstring, &wg, NULL, &status);
    }
    fits_write_img(faptr, TDOUBLE, 1, ((natmospherefile+1)*6 + nsurf*2 + 3)*(maxwavelength - minwavelength + 1), tempDynamicTransmission, &status);
    fits_close_file(faptr, &status);

    free(tempDynamicTransmission);

    filename.str("");
    filename << "!"<< workdir << "/" << outputfilename << "_ckptfp_" << checkpointcount << ".fits.gz";
    fits_create_file(&faptr, filename.str().c_str(), &status);
    naxes[0] = 1;
    naxes[1] = 1;
    fits_create_img(faptr, FLOAT_IMG, 2, naxes, &status);
    naxes[0] = chip.extraX;
    naxes[1] = chip.extraY;
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    fits_write_img(faptr, TFLOAT, 1, chip.extraX*chip.extraY, state.focalPlane, &status);
    fits_close_file(faptr, &status);

}

void Image::readCheckpoint(int checkpointcount) {

    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    int status;
    double *tempDynamicTransmission;
    Uint32 z = 0, w = 0, zg = 0, wg = 0;
    char tempstring[4096];

    tempDynamicTransmission = static_cast<double*>(calloc(((natmospherefile+1)*6 + nsurf*2 + 3)*(maxwavelength - minwavelength + 1), sizeof(double)));

    std::ostringstream filename;
    filename << workdir << "/" << outputfilename << "_ckptdt_" << checkpointcount - 1 << ".fits.gz";
    status = 0;
    if (fits_open_file(&faptr, filename.str().c_str(), READONLY, &status)) {
        printf("Error opening %s\n", filename.str().c_str());
        exit(1);
    }
    fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
    for (int i = 0; i < numthread; i++) {
        snprintf(tempstring, sizeof(tempstring), "M_Z%4d", i);
        fits_read_key(faptr, TUINT, tempstring, &z, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "M_W%4d", i);
        fits_read_key(faptr, TUINT, tempstring, &w, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "G_Z%4d", i);
        fits_read_key(faptr, TUINT, tempstring, &zg, NULL, &status);
        snprintf(tempstring, sizeof(tempstring), "G_W%4d", i);
        fits_read_key(faptr, TUINT, tempstring, &wg, NULL, &status);
        random[i].setSeed(z, w);
        galaxy.random[i].setSeed(zg, wg);
    }
    fits_read_img(faptr, TDOUBLE, 1, naxes[0]*naxes[1], &nullval, tempDynamicTransmission, &anynull, &status);
    fits_close_file(faptr, &status);

    for (int i = 0; i < maxwavelength - minwavelength + 1; i++) {
        for (int j = 0; j < ((natmospherefile+1)*6 + nsurf*2 + 3); j++) {
            state.dynamicTransmission[i*((natmospherefile+1)*6 + nsurf*2 + 3) + j] =
                tempDynamicTransmission[i*((natmospherefile+1)*6 + nsurf*2 + 3) + j];
        }
    }

    free(tempDynamicTransmission);


    filename.str("");
    filename << workdir << "/" << outputfilename << "_ckptfp_" << checkpointcount - 1 << ".fits.gz";
    status = 0;
    if (fits_open_file(&faptr, filename.str().c_str(), READONLY, &status)) {
        printf("Error opening %s\n", filename.str().c_str());
        exit(1);
    }
    fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
    fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval, state.focalPlane, &anynull, &status);
    fits_close_file(faptr, &status);

}


void Image::cleanup () {


}
