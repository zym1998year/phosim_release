///
/// @package phosim
/// @file sed.cpp
/// @brief sed reader
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

void Observation::readSed(const std::string & filename, int mode) {

    long lsedptr;
    char *sptr;
    char *sptr2;
    char *sptr3;
    char tempstring[4096];
    double tempf1;
    double tempf2;
    double dw;
    double cdw;
    long closestw = 0;
    long j;
    FILE *indafile;
    gzFile ingzfile;
    char line[4096];
    char tempstring1[512];
    char tempstring2[512];
    double monochromaticFactor;

    normwave = standardwavelength;
    monochromaticFactor = 1.0;

    lsedptr = 0;
    if (mode == 0) snprintf(tempstring, sizeof(tempstring), "%s/%s", seddir.c_str(), filename.c_str());
    if (mode == 1) snprintf(tempstring, sizeof(tempstring), "%s", filename.c_str());
    if (flatdir == 1) {
        sptr = strtok_r(tempstring, "/", &sptr2);
        do {
            sptr3 = sptr;
            sptr = strtok_r(NULL, "/", &sptr2);
        } while (sptr != NULL);
        snprintf(tempstring, sizeof(tempstring), "%s", sptr3);
    }

    if (strstr(tempstring, "laser") != NULL) {
        normwave = domewave;
        closestw = 0;
        cdw = 1e30;

        tempf1 = 0.0;
        tempf2 = 0.0;
        for (int k = 0; k < 2; k++) {

            if (k == 0) {
                tempf1 = domewave - 0.51e-6;
                tempf2 = 0.0;
            }
            if (k == 1) {
                tempf1 = domewave + 0.49e-6;
                tempf2 = 1.0/tempf1;
            }
            sedW.push_back(tempf1);
            sedC.push_back(tempf2*tempf1);

            dw = fabs(tempf1 - domewave);

            if (dw < cdw) {
                cdw = dw;
                closestw = lsedptr;
            }
            lsedptr++;
        }
        monochromaticFactor = normwave/(0.5e-6)*(0.01*H_CGS/pow(10.0, (20.0 + 48.6)/(-2.5)));

    } else {


        if (strstr(tempstring, ".gz") == NULL) {
            //printf("%s\n",tempstring);
            indafile = fopen(tempstring, "r");
            if (indafile == NULL) {
                fprintf(stderr, "Can't find SED file: %s\n", tempstring);
                exit(1);
            }

            closestw = 0;
            cdw = 1e30;
            double lastwavelength = 0;
            double runningflux = 0.0;
            double runningwavelength = 0.0;
            int count = 0;
            int init = 0;
            while (fgets(line, 4096, indafile)) {
                sscanf(line, "%s %s", tempstring1, tempstring2);
                if (tempstring1[0] != '#') {
                tempf1 = strtod(tempstring1, NULL);
                tempf2 = strtod(tempstring2, NULL);
                if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                    fprintf(stderr, "Warning:   SED file: %s contains a NaN!\n", tempstring);
                } else {
                        runningwavelength += tempf1;
                        runningflux += tempf2;
                        count++;
                        if (init == 0)  lastwavelength = tempf1;
                        //init = 1;
                        if (fabs(lastwavelength - tempf1) >= SED_TOLERANCE || fabs(tempf1 - standardwavelength) <= SED_TOLERANCE || init==0) {
                            sedW.push_back(runningwavelength/count);
                            sedC.push_back(runningflux/count*sedW[sedptr + lsedptr]);
                          lastwavelength = tempf1;
                          init=1;
                            runningwavelength = 0.0;
                            runningflux = 0.0;
                            count = 0;
                            dw = fabs(sedW[sedptr+lsedptr] - standardwavelength);
                            if (dw < cdw && sedC[sedptr+lsedptr] > 0.0) {
                                cdw = dw;
                                closestw = lsedptr;
                            }
                            lsedptr++;
                        }
                }
                }
            }
            fclose(indafile);
        } else {

            ingzfile = gzopen(tempstring, "r");
            if (ingzfile == NULL) {
                fprintf(stderr, "Can't find SED file: %s\n", tempstring);
                exit(1);
            }

            closestw = 0;
            cdw = 1e30;
            double lastwavelength = 0;
            double runningflux = 0.0;
            double runningwavelength = 0.0;
            int count = 0;
            int init = 0;
            while (gzgets(ingzfile, line, 4096)) {
                sscanf(line, "%s %s", tempstring1, tempstring2);
                if (tempstring1[0] != '#') {
                tempf1 = strtod(tempstring1, NULL);
                tempf2 = strtod(tempstring2, NULL);
                if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                    fprintf(stderr, "Warning:   SED file: %s contains a NaN!\n", tempstring);
                } else {
                        runningwavelength += tempf1;
                        runningflux += tempf2;
                        count++;
                        if (init == 0)  lastwavelength = tempf1;
                        //init = 1;
                        if (fabs(lastwavelength - tempf1) >= SED_TOLERANCE || fabs(tempf1 - standardwavelength) <= SED_TOLERANCE || (init==0)) {
                            sedW.push_back(runningwavelength/count);
                            sedC.push_back(runningflux/count*sedW[sedptr + lsedptr]);
                            lastwavelength = tempf1;
                            runningwavelength = 0.0;
                            runningflux = 0.0;
                            count = 0;
                            init = 1;
                            dw = fabs(sedW[sedptr+lsedptr] - standardwavelength);
                            if (dw < cdw && sedC[sedptr+lsedptr] > 0.0) {
                                cdw = dw;
                                closestw = lsedptr;
                            }
                            lsedptr++;
                        }
                }
                }
            }
            gzclose(ingzfile);

        }


        // monochromatic exception
        if (lsedptr == 1) {
            lsedptr--;
            sedptr++;
            normwave = tempf1;
            closestw = 0;
            cdw = 1e30;

            tempf1 = 0.0;
            tempf2 = 0.0;
            for (int k = 0; k < 2; k++) {

                if (k == 0) {
                    tempf1 = normwave - 0.51e-6;
                    tempf2 = 0.0;
                }
                if (k == 1) {
                    tempf1 = normwave + 0.49e-6;
                    tempf2 = 1.0/normwave;
                }

                sedW.push_back(tempf1);
                sedC.push_back(tempf2*tempf1);

                dw = fabs(tempf1 - normwave);

                if (dw < cdw) {
                    cdw = dw;
                    closestw = lsedptr;
                }
                lsedptr++;
            }
            monochromaticFactor = normwave/(0.5e-6)*(0.01*H_CGS/pow(10.0, (20.0 + 48.6)/(-2.5)));
        }


    }

    for (j = 0; j < lsedptr; j++) {
        if (j != 0 && j != (lsedptr - 1)) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j + 1] - sedW[sedptr + j - 1])/2.0;
        if (j == 0) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j + 1] - sedW[sedptr + j]);
        if (j == (lsedptr - 1)) sedC[sedptr + j] = sedC[sedptr + j]*(sedW[sedptr + j] - sedW[sedptr + j - 1]);
    }

    tempf1 = 0;
    for (j = 0; j < lsedptr; j++) {
        tempf1 += sedC[sedptr + j];
    }
    for (j = 0; j < lsedptr; j++) {
        sedC[sedptr + j] = sedC[sedptr + j]/tempf1;
    }

    if (closestw == 0) {
        sedDwdp[nsedptr] = (sedW[sedptr + closestw + 1] - sedW[sedptr + closestw])/sedC[sedptr + closestw];
    } else if (closestw == (lsedptr-1)) {
        sedDwdp[nsedptr] = (sedW[sedptr + closestw] - sedW[sedptr + closestw - 1])/sedC[sedptr + closestw];
    } else {
        sedDwdp[nsedptr] = (sedW[sedptr + closestw + 1] - sedW[sedptr + closestw - 1])/sedC[sedptr + closestw]/2.0;
    }

    //    printf("DW %e %e %e\n",sedDwdp[nsedptr],sedW[sedptr + closestw + 1] - sedW[sedptr + closestw - 1],sedC[sedptr + closestw]);
    if (sedC[sedptr + closestw] <= 0.0) {
        printf("Error in SED file; 0 value at normalizing wavelength in %s\n", filename.c_str());
        sedDwdp[nsedptr] = 0.0;
    }
    sedDwdp[nsedptr] = sedDwdp[nsedptr]*monochromaticFactor;

    for (j = 1; j < lsedptr; j++) {
        sedC[sedptr + j] += sedC[sedptr + j - 1];
    }

    sedN[nsedptr] = lsedptr;
    sedPtr[nsedptr] = sedptr;
    sedMode[nsedptr] = mode;
    sedptr += lsedptr;
    nsedptr++;
};
