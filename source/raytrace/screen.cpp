///
/// @package phosim
/// @file screen.cpp
/// @brief screen class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "screen.h"

double Screen::readScreen (int keynum, float *array, char *tempstring) {


        char *ffptr;
        fitsfile *faptr;
        long naxes[2];
        int nfound;
        int anynull;
        float nullval=0;  //must be defined for call to fits_read_img
        char keyname[4096];
        char comment[4096];
        char value[4096];
        int status;

        ffptr=tempstring;
        status=0;
        if (fits_open_file(&faptr,ffptr,READONLY,&status)) {
            printf("Error opening %s\n",ffptr);
            exit(1);
        }
        fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
        if (keynum >= 0) fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
        fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,array,&anynull,&status);
        fits_close_file(faptr,&status);

        return(atof(value));

}

void Screen::downsample (float *newScreen, float *oldScreen) {

    float frac = 1.0/(float(ATM_RESAMPLE))/(float(ATM_RESAMPLE));

    for (int i = 0; i < SCREEN_SIZE;i++) {
        for (int j = 0; j < SCREEN_SIZE; j++) {
            for (int k=0; k < ATM_RESAMPLE; k++) {
                for (int l=0; l < ATM_RESAMPLE; l++) {
                    newScreen[i*SCREEN_SIZE + j] += frac*oldScreen[(ATM_RESAMPLE*i + k)*SCREEN_SIZE + (ATM_RESAMPLE*j + l)];
                }
            }
        }
    }

}
