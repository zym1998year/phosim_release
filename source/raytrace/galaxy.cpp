///
/// @package phosim
/// @file galaxy.cpp
/// @brief galaxy photon sampler
///
/// @brief Created by:
/// @author Suzanne Lorenz (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"
#include "galaxy.h"
#include "helpers.h"
#include "constants.h"
#include "ancillary/random.h"

void Galaxy::sampleSersic(char *sersicdata) {
    char temp[20];
    float bN[99];
    double t = 0;
    double sum = 0;
    FILE *fp;

    fAll = static_cast<double**>(calloc(99, sizeof(double*)));
    for (int n = 0; n < 99; n++) {
        fAll[n] = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));
    }
    cAll = static_cast<double**>(calloc(99, sizeof(double*)));
    for (int n = 0; n < 99; n++) {
        cAll[n] = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));
    }
    dAll = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));

    fp = fopen(sersicdata, "r");
    for (int i = 0; i < 99; i++) {
        fgets(temp, 20, fp);
        sscanf(temp, "%f", &bN[i]);
    }
    for (int n = 0; n < 99; n++) {
        t = 0.0;
        sum = 0.0;
        for (int i = 0; i < GALAXY_POINT; i++) {
            fAll[n][i] = t*t*exp(-bN[n]*(pow(t, 1/((n + 1)/GALAXY_SCALE)) - 1));
            cAll[n][i] = fAll[n][i] + sum;
            sum = cAll[n][i];
            dAll[i] = t;
            t = t + 0.01;
        }
        for (int i = 0; i < GALAXY_POINT; i++) {
            cAll[n][i] = cAll[n][i]/sum;
        }
    }
    fclose(fp);

}

void Galaxy::sersic(double a, double b, double c, double alpha, double beta, double gamma,
                    double n, double *x_out, double *y_out, int thread) {

    double theta = PI*2*random[thread].uniform();
    double phi = acos(2*random[thread].uniform() - 1.0);

    double dummy = n*GALAXY_SCALE;
    int nn = static_cast<int>(dummy);
    int q = 0;
    double randomValue = random[thread].uniform();
    find(*(cAll + nn), GALAXY_POINT, randomValue, &q);
    int down = q;
    int up = q + 1;
    double frak = (randomValue - cAll[nn][down])/(cAll[nn][up] - cAll[nn][down]);
    double r_a = a*(dAll[down]*(1 - frak) + dAll[up]*frak);
    double r_b = b*(dAll[down]*(1 - frak) + dAll[up]*frak);
    double r_c = c*(dAll[down]*(1 - frak) + dAll[up]*frak);

    double x, y, z;
    x = sin(phi)*cos(theta)*r_a;
    y = sin(phi)*sin(theta)*r_b;
    z = cos(phi)*r_c;

    // // *z_out = sin(beta)*sin(alpha)*x - sin(beta)*cos(alpha)*y + cos(beta)*z;
    // *x_out = cos(alpha)*x + sin(alpha)*y;
    // *y_out = -cos(beta)*sin(alpha)*x + cos(beta)*cos(alpha)*y + sin(beta)*z;

    *x_out = (cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma))*x +
        (cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma))*y +
        (sin(gamma)*sin(beta))*z;
    *y_out = (-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma))*x +
        (-sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma))*y +
        (cos(gamma)*sin(beta))*z;

}

void Galaxy::sampleSersic2d(char *sersicdata, int numthread) {
    int i;
    int n;
    char temp[20];
    float bN[99];
    double t = 0;
    double sum = 0;
    FILE *fp;

    fAll2d = static_cast<double**>(calloc(99, sizeof(double*)));
    for (n = 0; n < 99; n++) {
        fAll2d[n] = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));
    }
    cAll2d = static_cast<double**>(calloc(99, sizeof(double*)));
    for (n = 0; n < 99; n++) {
        cAll2d[n] = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));
    }
    dAll2d = static_cast<double*>(calloc(GALAXY_POINT, sizeof(double)));
    clumpX = static_cast<double*>(calloc(1000*numthread, sizeof(double)));
    clumpY = static_cast<double*>(calloc(1000*numthread, sizeof(double)));
    clumpZ = static_cast<double*>(calloc(1000*numthread, sizeof(double)));

    fp = fopen(sersicdata, "r");
    for (i = 0; i < 99; i++) {
        fgets(temp, 20, fp);
        sscanf(temp, "%f", &bN[i]);
    }
    for (n = 0; n < 99; n++) {
        t = 0.0;
        sum = 0.0;
        for (i = 0; i < GALAXY_POINT; i++) {
            fAll2d[n][i] = t*exp(-bN[n]*(pow(t, 1/((n + 1)/GALAXY_SCALE)) - 1));
            cAll2d[n][i] = fAll2d[n][i] + sum;
            sum = cAll2d[n][i];
            dAll2d[i] = t;
            t = t + 0.01;
        }
        for (i = 0; i < GALAXY_POINT; i++) {
            cAll2d[n][i] = cAll2d[n][i]/sum;
        }
    }
    fclose(fp);

}

void Galaxy::sersic2d(double a, double b, double beta, double n, double *x_out, double *y_out, int thread) {

    double theta = PI*2*(random[thread].uniform());
    double randomValue = random[thread].uniform();
    int nn = (int)(n*GALAXY_SCALE);
    int q = 0;
    find(*(cAll2d + nn), GALAXY_POINT, randomValue, &q);
    double r = interpolate(dAll2d, *(cAll2d + nn), randomValue, q);
    double x = cos(theta)*a*r;
    double y = sin(theta)*b*r;

    *x_out = cos(beta)*x + sin(beta)*y;
    *y_out = -sin(beta)*x + cos(beta)*y;
}


void Galaxy::sersicComplex(double a, double b, double c, double alpha, double beta, double gamma,
                          double n, double clumpFrac, double clump, double clumpWidth, double spiralFrac,
                           double alphaSpiral, double bar, double spiralWidth, double phi0, int thread, int *init,
                          double *x_out, double *y_out) {

    if (*init == 0) {
        int intClump = static_cast<int>(clump);
        if (intClump > 1000) intClump = 1000;
        if (intClump > 0) {
            for (int i = 0; i < intClump; i++) {
                double theta = PI*2*random[thread].uniform();
                double phi = acos(2*random[thread].uniform() - 1.0);

                double dummy = n*GALAXY_SCALE;
                int nn = static_cast<int>(dummy);
                int q = 0;
                double randomValue = random[thread].uniform();
                find(*(cAll + nn), GALAXY_POINT, randomValue, &q);
                int down = q;
                int up = q + 1;
                double frak = (randomValue - cAll[nn][down])/(cAll[nn][up] - cAll[nn][down]);
                double r_a = a*(dAll[down]*(1 - frak) + dAll[up]*frak);
                double r_b = b*(dAll[down]*(1 - frak) + dAll[up]*frak);
                double r_c = c*(dAll[down]*(1 - frak) + dAll[up]*frak);

                double x, y, z;
                if (random[thread].uniform() < spiralFrac) {
                    double rPol = sqrt(r_a*r_a + r_b*r_b);
                    if (rPol < bar) rPol = bar;
                    if (random[thread].uniform() < 0.5) {
                        theta = PI + log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
                    } else {
                        theta = log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
                    }
                    x = sin(phi)*cos(theta)*r_a + random[thread].normal()*spiralWidth;
                    y = sin(phi)*sin(theta)*r_b + random[thread].normal()*spiralWidth;
                    z = cos(phi)*r_c;
                } else {
                    x = sin(phi)*cos(theta)*r_a;
                    y = sin(phi)*sin(theta)*r_b;
                    z = cos(phi)*r_c;
                }
                clumpX[intClump*thread + i] = x;
                clumpY[intClump*thread + i] = y;
                clumpZ[intClump*thread + i] = z;
            }
        }
        *init = 1;
    }

    double theta = PI*2*random[thread].uniform();
    double phi = acos(2*random[thread].uniform() - 1.0);

    double dummy = n*GALAXY_SCALE;
    int nn = static_cast<int>(dummy);
    int q = 0;
    double randomValue = random[thread].uniform();
    find(*(cAll + nn), GALAXY_POINT, randomValue, &q);
    int down = q;
    int up = q + 1;
    double frak = (randomValue - cAll[nn][down])/(cAll[nn][up] - cAll[nn][down]);
    double r_a = a*(dAll[down]*(1 - frak) + dAll[up]*frak);
    double r_b = b*(dAll[down]*(1 - frak) + dAll[up]*frak);
    double r_c = c*(dAll[down]*(1 - frak) + dAll[up]*frak);

    double x, y, z;
    if (random[thread].uniform() < clumpFrac) {
        int intClump = static_cast<int>(clump);
        if (intClump > 1000) intClump = 1000;
        int index = floor(random[thread].uniform()*intClump);
        x = clumpX[intClump*thread + index] + random[thread].normal()*clumpWidth;
        y = clumpY[intClump*thread + index] + random[thread].normal()*clumpWidth;
        z = clumpZ[intClump*thread + index] + random[thread].normal()*clumpWidth;
    } else {
        if (random[thread].uniform() < spiralFrac) {
            double rPol = sqrt(r_a*r_a + r_b*r_b);
            if (rPol < bar) rPol = bar;
            if (random[thread].uniform() < 0.5) {
                theta = PI + log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
            } else {
                theta = log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
            }
            x = sin(phi)*cos(theta)*r_a + random[thread].normal()*spiralWidth;
            y = sin(phi)*sin(theta)*r_b + random[thread].normal()*spiralWidth;
            z = cos(phi)*r_c;
        } else {
            x = sin(phi)*cos(theta)*r_a;
            y = sin(phi)*sin(theta)*r_b;
            z = cos(phi)*r_c;
        }
    }

    // // *z_out = sin(beta)*sin(alpha)*x - sin(beta)*cos(alpha)*y + cos(beta)*z;
    // *x_out = cos(alpha)*x + sin(alpha)*y;
    // *y_out = -cos(beta)*sin(alpha)*x + cos(beta)*cos(alpha)*y + sin(beta)*z;

    *x_out = (cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma))*x +
        (cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma))*y +
        (sin(gamma)*sin(beta))*z;
    *y_out = (-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma))*x +
        (-sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma))*y +
        (cos(gamma)*sin(beta))*z;
}


void Galaxy::sersicDisk(double a, double b, double c, double alpha, double beta, double gamma,
                        double n, double *x_out, double *y_out, int thread) {

    double theta = PI*2*random[thread].uniform();

    double dummy = n*GALAXY_SCALE;
    int nn = static_cast<int>(dummy);
    int q = 0;
    double randomValue = random[thread].uniform();
    find(*(cAll2d + nn), GALAXY_POINT, randomValue, &q);
    double r = interpolate(dAll2d, *(cAll2d + nn), randomValue, q);
    double r_a = cos(theta)*a*r;
    double r_b = sin(theta)*b*r;
    double r_c = -c * log(random[thread].uniform());

    double x = r_a;
    double y = r_b;
    double z = r_c;

    // // *z_out = sin(beta)*sin(alpha)*x - sin(beta)*cos(alpha)*y + cos(beta)*z;
    // *x_out = cos(alpha)*x + sin(alpha)*y;
    // *y_out = -cos(beta)*sin(alpha)*x + cos(beta)*cos(alpha)*y + sin(beta)*z;

    *x_out = (cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma))*x +
        (cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma))*y +
        (sin(gamma)*sin(beta))*z;
    *y_out = (-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma))*x +
        (-sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma))*y +
        (cos(gamma)*sin(beta))*z;
}

void Galaxy::sersicDiskComplex(double a, double b, double c, double alpha, double beta, double gamma,
                          double n, double clumpFrac, double clump, double clumpWidth, double spiralFrac,
                               double alphaSpiral, double bar, double spiralWidth, double phi0, int thread, int *init,
                          double *x_out, double *y_out) {

    if (*init == 0) {
        int intClump = static_cast<int>(clump);
        if (intClump > 1000) intClump = 1000;
        if (intClump > 0) {
            for (int i = 0; i < intClump; i++) {
                double theta = PI*2*random[thread].uniform();

                double dummy = n*GALAXY_SCALE;
                int nn = static_cast<int>(dummy);
                int q = 0;
                double randomValue = random[thread].uniform();
                find(*(cAll2d + nn), GALAXY_POINT, randomValue, &q);
                double r = interpolate(dAll2d, *(cAll2d + nn), randomValue, q);
                double r_a = cos(theta)*a*r;
                double r_b = sin(theta)*b*r;
                double r_c = -c * log(random[thread].uniform());

                double x, y, z;
                if (random[thread].uniform() < spiralFrac) {
                    double rPol = sqrt(r_a*r_a + r_b*r_b);
                    if (rPol < bar) rPol = bar;
                    if (random[thread].uniform() < 0.5) {
                        theta = PI + log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
                    } else {
                        theta = log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
                    }
                    x = r_a + random[thread].normal()*spiralWidth;
                    y = r_b + random[thread].normal()*spiralWidth;
                    z = r_c;
                } else {
                    x = r_a;
                    y = r_b;
                    z = r_c;
                }
                clumpX[intClump*thread+i] = x;
                clumpY[intClump*thread+i] = y;
                clumpZ[intClump*thread+i] = z;
            }
        }
        *init = 1;
    }

    double theta = PI*2*random[thread].uniform();

    double dummy = n*GALAXY_SCALE;
    int nn = static_cast<int>(dummy);
    int q = 0;
    double randomValue = random[thread].uniform();
    find(*(cAll2d + nn), GALAXY_POINT, randomValue, &q);
    double r = interpolate(dAll2d, *(cAll2d + nn), randomValue, q);
    double r_a = cos(theta)*a*r;
    double r_b = sin(theta)*b*r;
    double r_c = -c * log(random[thread].uniform());

    double x, y, z;
    if (random[thread].uniform() < clumpFrac) {
        int intClump = static_cast<int>(clump);
        if (intClump > 1000) intClump = 1000;
        int index = floor(random[thread].uniform()*intClump);
        x = clumpX[intClump*thread+index] + random[thread].normal()*clumpWidth;
        y = clumpY[intClump*thread+index] + random[thread].normal()*clumpWidth;
        z = clumpZ[intClump*thread+index] + random[thread].normal()*clumpWidth;
    } else {
        if (random[thread].uniform() < spiralFrac) {
            double rPol = sqrt(r_a*r_a + r_b*r_b);
            if (rPol < bar) rPol = bar;
            if (random[thread].uniform() < 0.5) {
                theta = PI + log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
            } else {
                theta = log(rPol+1e-3)/tan(alphaSpiral) - log(bar+1e-3)/tan(alphaSpiral) + phi0;
            }
            x = r_a + random[thread].normal()*spiralWidth;
            y = r_b + random[thread].normal()*spiralWidth;
            z = r_c;
        } else {
            x = r_a;
            y = r_b;
            z = r_c;
        }
    }

    // // *z_out = sin(beta)*sin(alpha)*x - sin(beta)*cos(alpha)*y + cos(beta)*z;
    // *x_out = cos(alpha)*x + sin(alpha)*y;
    // *y_out = -cos(beta)*sin(alpha)*x + cos(beta)*cos(alpha)*y + sin(beta)*z;

    *x_out = (cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma))*x +
        (cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma))*y +
        (sin(gamma)*sin(beta))*z;
    *y_out = (-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma))*x +
        (-sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma))*y +
        (cos(gamma)*sin(beta))*z;
}

void Galaxy::distortedSphere(double alpha, double coeff0, double coeff1, double coeff2, double coeff3,
                             double coeff4, double coeff5, double coeff6, double coeff7, double coeff8,
                             double *x_out, double *y_out, int thread) {



    double r0;
    double phi = PI*2*random[thread].uniform();
    double theta = acos(2*random[thread].uniform() - 1.0);

    r0 = pow(random[thread].uniform(), 1.0/(alpha + 1.0));

    double r = coeff0*r0;
    r += coeff1*r0*sqrt(3.0/2.0)*sin(theta)*sin(phi);
    r += coeff2*r0*sqrt(3.0)*cos(theta);
    r -= coeff3*r0*sqrt(3.0/2.0)*sin(theta)*cos(phi);
    r += coeff4*r0*sqrt(15/2)/2.0*sin(theta)*sin(theta)*sin(2*phi);
    r += coeff5*r0*sqrt(15/2)*sin(theta)*cos(theta)*sin(phi);
    r += coeff6*r0*sqrt(5)/2.0*(3.0*cos(theta)*cos(theta)-1.0);
    r -= coeff7*r0*sqrt(15/2)*sin(theta)*cos(theta)*cos(phi);
    r += coeff8*r0*sqrt(15/2)/2.0*sin(theta)*sin(theta)*cos(2*phi);

    *x_out = r*sin(theta)*cos(phi);
    *y_out = r*sin(theta)*sin(phi);

}
