///
/// @package phosim
/// @file contamination.cpp
/// @brief contamination class
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

#include <math.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "contamination.h"
#include "constants.h"
#include "ancillary/random.h"

void Contamination::setup (double *innerRadius, double *outerRadius, long totalSurface, long points,
                           double pixsize, long nx, long ny, float qevariation, long seedchip) {

    double mmd = 0.3;
    double st = 0.0025;
    double stw = 0.0025;
    double std = 0.00025;
    double sd, sdw, sdd;
    double sigma = 1.0;
    double mmdw = 1.0;
    double sigmaw = 1.0;
    double sizetol = 0.3;
    double gridsize;
    double smallgridsize;
    double area;
    long particles, particlesd;
    double size;
    double fraction;
    // int r0, phi0, x0, y0;
    // double deltar, deltaphi, deltax, deltay;
    // int rd, pd, xd, yd;
    // int rc, pc, xc, yc;
    // int dx, dy, dr, dphi;
    // double x, y, r, phi;
    long number;
    double oversim;
    double totalloss;

    // double x1 = 0.7 + 2.6*random.normal();
    // double y1 = (x1 - 0.7)*(-0.08) + 0.93 + 0.40*random.normalFixed();
    // x1 = exp(x1);
    // mmdw = x1;
    // sigmaw = y1;

    // double x3 = -3.6  +  1.2*random.normal();
    // double y3 = (x3 + 3.6)*(-0.22) + 1.16 + 0.10*random.normalFixed();
    // x3= exp(x3);
    // mmd=x3;
    // sigma=y3;

    //much larger coarse size distributions
    //20 microns in clouds or fog
    mmd=10.0;
    sigma=1.0;
    mmdw = 20.0;
    sigmaw=1.0;

    random.setSeed32Fixed(1000);
    oversim = 1.0;

    //optics
    surfacenumber = new long[totalSurface]();
    for (long i = 0; i < totalSurface; i++) {
        surfacenumber[i] = 0;
    }

    transmission = new float[totalSurface*points*points]();
    for (long i = 0; i < totalSurface; i++) {
        for (long j = 0; j < points; j++) {
            for (long k = 0; k < points; k++) {
                transmission[i*points*points + j*points + k] = 1.0;
            }
        }
    }

    chiptransmission = new float[nx*ny]();
    for (long j = 0; j < nx; j++) {
        for (long k = 0; k < ny; k++) {
            chiptransmission[k*nx + j] = 1.0 - qevariation*fabs(random.normalFixed());
        }
    }

    surfacelistmap = new long[totalSurface*points*points]();
    surfacelist2ndmap = new long[totalSurface*points*points]();
    for (long i = 0; i < totalSurface; i++) {
        for (long j = 0; j < points; j++) {
            for (long k = 0; k < points; k++) {
                surfacelistmap[i*points*points + j*points + k] = -1;
                surfacelist2ndmap[i*points*points + j*points + k] = -1;
            }
        }
    }

    surfacelistx.resize(totalSurface);
    surfacelisty.resize(totalSurface);
    surfacelists.resize(totalSurface);
    surfacelistt.resize(totalSurface);

#ifdef PERTURBATION_COORDINATE
    int dphi, dr, pd, rd, rc, r0, pc, phi0;
    double deltaphi, deltar, r, phi;
    for (long i = 0; i<totalSurface; i++) {

        st = surfaceDustTransmission[i];
        stw = surfaceCondensationTransmission[i];
        sd = st/(PI*exp(2.0*log(mmd) + 0.5*2.0*2.0*sigma*sigma)*1e-3*1e-3);
        sdw = stw/(PI*exp(2.0*log(mmdw) + 0.5*2.0*2.0*sigmaw*sigmaw)*1e-3*1e-3);

        
        deltar = outerRadius[i] - sqrt((outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i])*(points - 2)/(points - 1)+
                                       innerRadius[i]*innerRadius[i]);
        deltaphi = 2*PI/(points - 1);
        totalloss = 0.0;
        area = PI*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]);
        particles = static_cast<long>(ceil(sd*area + sdw*area));
        particlesd = static_cast<long>(ceil(sd*area));
        gridsize = sizetol*sqrt(area/PI)/(static_cast<double>(points));
        smallgridsize = sizetol*gridsize;
        long j = 0;
        //        long c1 = 0;
        //        long c2 = 0;
        //        long c3 = 0;
        while (j < particles) {
            if (j < particlesd) {
                size = exp(log(mmd) + sigma*random.normalFixed())*1e-3;
            } else {
                size = exp(log(mmdw) + sigmaw*random.normalFixed())*1e-3;
            }
            number = oversim;
            if (number < 1) number = 1;
            if (size > gridsize) {
                // large contaminant
                for (int k = 0; k < number; k++) {
                    //                    c1++;
                    r = sqrt(random.uniformFixed()*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]) +
                             innerRadius[i]*innerRadius[i]);
                    phi = random.uniformFixed()*2*PI;
                    r0 = static_cast<long>((r*r - innerRadius[i]*innerRadius[i])/(outerRadius[i]*outerRadius[i]
                                                                                  - innerRadius[i]*innerRadius[i])*(static_cast<double>(points)-1));
                    phi0 = static_cast<long>((phi)/(2*PI)*(static_cast<double>(points) - 1));
                    x = r*cos(phi);
                    y = r*sin(phi);
                    surfacelistx[i].push_back(x);
                    surfacelisty[i].push_back(y);
                    surfacelists[i].push_back(size);
                    if (j < particlesd) {
                        surfacelistt[i].push_back(0);
                    } else {
                        surfacelistt[i].push_back(1);
                    }
                    dr = ceil((size/deltar) + 1);
                    dphi = ceil((size/(deltaphi*sqrt((outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i])*
                                                     (static_cast<double>(r0))/(static_cast<double>(points) - 1) +
                                                     innerRadius[i]*innerRadius[i]))) + 1);
                    if ((dphi < 0) || (dphi > points/2)) dphi = points/2;
                    for (rd = r0 - dr; rd <= r0 + dr; rd++) {
                        for (pd = phi0 - dphi; pd <= phi0 + dphi; pd++) {
                            rc = rd;
                            pc = pd;
                            if (rc < 0) rc = 0;
                            if (rc > points - 1) rc = points - 1;
                            if (pc < 0) pc = pc + points;
                            if (pc > points - 1) pc = pc - points;
                            long index=i*points*points + rc*points + pc;
                            if (surfacelistmap[index] != -1) {
                                if (surfacelists[i][surfacenumber[i]] > surfacelists[i][surfacelistmap[index]]) {
                                    surfacelist2ndmap[index] = surfacelistmap[index];
                                    surfacelistmap[index] = surfacenumber[i];
                                } else {
                                    if (surfacelists[i][surfacenumber[i]] > surfacelists[i][surfacelist2ndmap[index]]) {
                                        surfacelist2ndmap[index] = surfacenumber[i];
                                    }
                                }
                            } else {
                                    surfacelistmap[index] = surfacenumber[i];
                            }
                        }
                    }
                    surfacenumber[i]++;
                }
            } else if (size > smallgridsize) {

                // small contaminant
                for (int k = 0; k < number; k++) {
                    //                    c2++;
                    r = sqrt(random.uniformFixed()*(outerRadius[i]*outerRadius[i] - innerRadius[i]*innerRadius[i]) +
                             innerRadius[i]*innerRadius[i]);
                    phi = random.uniformFixed()*2*PI;
                    r0 = static_cast<long>((r*r - innerRadius[i]*innerRadius[i])/(outerRadius[i]*outerRadius[i]
                                                                                  - innerRadius[i]*innerRadius[i])*(static_cast<double>(points)-1));
                    phi0 = static_cast<long>((phi)/(2*PI)*(static_cast<double>(points) - 1));
                    if (r0 < 0) r0 = 0;
                    if (r0 > points - 1) r0 = points - 1;
                    if (phi0 < 0) phi0 = phi0 + points;
                    if (phi0 > points - 1) phi0 = phi0 - points;
                    fraction = PI*size*size/(area/(static_cast<double>(points))/(static_cast<double>(points)));
                    long index = i*points*points + r0*points + phi0;
                    transmission[index] -= fraction;
                    if (transmission[index] < 0) {
                        transmission[index] = 0.0;
                    }
                }
            } else {
                // super small contaminant
                totalloss += number*PI*size*size;
                //                c3 += number;
            }

            j += number;
        }
        double avg = 0.0;
        for (long k = 0; k < points; k++) {
            for (long l = 0; l < points; l++) {
                long index = i*points*points + k*points + l;
                transmission[index] -= totalloss/area;
                if (transmission[index] < 0) {
                    transmission[index] = 0.0;
                }
                avg += transmission[index];
            }
        }

    }
#else
    int x0, y0, xd, yd, xc, yc, dx, dy;
    double x, y, deltax, deltay;
    for (long i = 0; i<totalSurface; i++) {

        st = surfaceDustTransmission[i];
        stw = surfaceCondensationTransmission[i];
       sd = st/(PI*exp(2.0*log(mmd) + 0.5*2.0*2.0*sigma*sigma)*1e-3*1e-3);
        sdw = stw/(PI*exp(2.0*log(mmdw) + 0.5*2.0*2.0*sigmaw*sigmaw)*1e-3*1e-3);
        deltax = outerRadius[i]*2.0/(points - 1);
        deltay = outerRadius[i]*2.0/(points - 1);
        totalloss = 0.0;
        area = 4.0*(outerRadius[i]*outerRadius[i]);
        particles = static_cast<long>(ceil(sd*area + sdw*area));
        particlesd = static_cast<long>(ceil(sd*area));
        gridsize = sizetol*sqrt(area)/(static_cast<double>(points));
        smallgridsize = sizetol*gridsize;
        long j = 0;
        // long c1 = 0;
        // long c2 = 0;
        // long c3 = 0;
        while (j < particles) {
            if (j < particlesd) {
                size = exp(log(mmd) + sigma*random.normalFixed())*1e-3;
            } else {
                size = exp(log(mmdw) + sigmaw*random.normalFixed())*1e-3;
            }
            number = oversim;
            if (number < 1) number = 1;
            if (size > gridsize) {
                // large contaminant
                for (int k = 0; k < number; k++) {
                                        // c1++;
                    x = random.uniformFixed()*(outerRadius[i]*2.0)-outerRadius[i];
                    y = random.uniformFixed()*(outerRadius[i]*2.0)-outerRadius[i];
                    x0 = static_cast<long>((x+outerRadius[i])/(outerRadius[i]*2.0)*(static_cast<double>(points)-1));
                    y0 = static_cast<long>((y+outerRadius[i])/(outerRadius[i]*2.0)*(static_cast<double>(points)-1));
                    surfacelistx[i].push_back(x);
                    surfacelisty[i].push_back(y);
                    surfacelists[i].push_back(size);
                    if (j < particlesd) {
                        surfacelistt[i].push_back(0);
                    } else {
                        surfacelistt[i].push_back(1);
                    }
                    dx = ceil((size/deltax) + 1);
                    dy = ceil((size/deltay) + 1);
                    for (xd = x0 - dx; xd <= x0 + dx; xd++) {
                        for (yd = y0 - dy; yd <= y0 + dy; yd++) {
                            if (sqrt((xd-x0)*(xd-x0)+(yd-y0)*(yd-y0)) <= dx) {
                            xc = xd;
                            yc = yd;
                            if (xc < 0) xc = 0;
                            if (xc >= points - 1) xc = points - 1;
                            if (yc < 0) yc = 0;
                            if (yc >= points - 1) yc = points - 1;
                            long index=i*points*points + yc*points + xc;
                            if (surfacelistmap[index] != -1) {
                                if (surfacelists[i][surfacenumber[i]] > surfacelists[i][surfacelistmap[index]]) {
                                    surfacelist2ndmap[index] = surfacelistmap[index];
                                    surfacelistmap[index] = surfacenumber[i];
                                } else {
                                    if (surfacelists[i][surfacenumber[i]] > surfacelists[i][surfacelist2ndmap[index]]) {
                                        surfacelist2ndmap[index] = surfacenumber[i];
                                    }
                                }
                            } else {
                                    surfacelistmap[index] = surfacenumber[i];
                            }
                            }
                        }
                    }
                    surfacenumber[i]++;
                }
            } else if (size > smallgridsize) {

                // small contaminant
                for (int k = 0; k < number; k++) {
                                        // c2++;
                    x = random.uniformFixed()*(outerRadius[i]*2.0)-outerRadius[i];
                    y = random.uniformFixed()*(outerRadius[i]*2.0)-outerRadius[i];
                    x0 = static_cast<long>((x+outerRadius[i])/(outerRadius[i]*2.0)*(static_cast<double>(points)-1));
                    y0 = static_cast<long>((y+outerRadius[i])/(outerRadius[i]*2.0)*(static_cast<double>(points)-1));
                    if (x0 < 0) x0 = 0;
                    if (x0 > points - 1) x0 = points - 1;
                    if (y0 < 0) y0 = 0;
                    if (y0 > points - 1) y0 = points - 1;
                    fraction = PI*size*size/(area/(static_cast<double>(points))/(static_cast<double>(points)));
                    long index = i*points*points + y0*points + x0;
                    transmission[index] -= fraction;
                    if (transmission[index] < 0) {
                        transmission[index] = 0.0;
                    }
                }
            } else {
                // super small contaminant
                totalloss += number*PI*size*size;
                // c3 += number;
            }

            j += number;
        }

        //        double avg = 0.0;
        for (long k = 0; k < points; k++) {
            for (long l = 0; l < points; l++) {
                long index = i*points*points + k*points + l;
                transmission[index] -= totalloss/area;
                if (transmission[index] < 0) {
                    transmission[index] = 0.0;
                }
                //    avg += transmission[index];
            }
        }

    }
#endif

    // int status;
    // long naxesa[2];
    // fitsfile *faptr;
    // status=0;
    // fits_create_file(&faptr,"!cont_m.fits",&status);
    // naxesa[0]=points; naxesa[1]=points;
    // fits_create_img(faptr,FLOAT_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TFLOAT,1,points*points,transmission,&status);
    // fits_close_file(faptr,&status);

    // status=0;
    // fits_create_file(&faptr,"!cont_p.fits",&status);
    // naxesa[0]=points; naxesa[1]=points;
    // fits_create_img(faptr,LONG_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TLONG,1,points*points,surfacelistmap,&status);
    // fits_close_file(faptr,&status);

    //chip
    std = surfaceDustTransmission[totalSurface];
    sdd = std/(PI*exp(2.0*log(mmd) + 0.5*2.0*2.0*sigma*sigma)*1e-3*1e-3);

    random.setSeed32Fixed(1000 + seedchip);
    chiplistmap = new int[(nx)*(ny)]();
    chiplist2ndmap = new int[(nx)*(ny)]();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            chiplistmap[j*nx + i] = -1;
            chiplist2ndmap[j*nx + i] = -1;
        }
    }
    chipnumber = 0;
    area = (nx*pixsize*ny*pixsize);
    particles = static_cast<int>(ceil(sdd*area));
    gridsize = sizetol*pixsize/sqrt(PI);
    smallgridsize = sizetol*gridsize;
    deltax = pixsize;
    deltay = pixsize;
    for (long j = 0; j < particles; j++) {
        size = exp(log(mmd) + sigma*random.normalFixed())*1e-3;
        x = random.uniformFixed()*nx*pixsize;
        y = random.uniformFixed()*ny*pixsize;
        x0 = round((int)(x/pixsize));
        y0 = round((int)(y/pixsize));
        x = x - nx/2*pixsize;
        y = y - ny/2*pixsize;

        if (size > gridsize) {
            chiplistx.push_back(x);
            chiplisty.push_back(y);
            chiplists.push_back(size);
            dx = ceil((size/deltax) + 1);
            dy = ceil((size/deltay) + 1);
            for (xd = x0 - dx; xd <= x0 + dx; xd++) {
                for (yd = y0 - dy; yd <= y0 + dy; yd++) {
                    if (sqrt((xd-x0)*(xd-x0)+(yd-y0)*(yd-y0)) <= dx) {
                    xc = xd;
                    yc = yd;
                    if (xc < 0) xc = 0;
                    if (xc >= nx) xc = nx - 1;
                    if (yc < 0) yc = 0;
                    if (yc >= ny) yc = ny - 1;
                    if (chiplistmap[nx*yc + xc] != -1) {
                        if (chiplists[chipnumber] > chiplists[chiplistmap[nx*yc + xc]]) {
                                chiplist2ndmap[nx*yc + xc] = chiplistmap[nx*yc + xc];
                                chiplistmap[nx*yc + xc] = chipnumber;
                        } else {
                            if (chiplists[chipnumber] > chiplists[chiplist2ndmap[nx*yc + xc]]) {
                                chiplist2ndmap[nx*yc + xc] = chipnumber;
                            }
                        }
                    } else {
                        chiplistmap[nx*yc + xc] = chipnumber;
                    }
                }
                }
            }
            chipnumber++;
        } else {
            fraction = PI*size*size/(pixsize*pixsize);
            if (x0 < 0) x0 = 0;
            if (y0 < 0) y0 = 0;
            if (x0 > nx - 1) x0 = nx - 1;
            if (y0 > ny - 1) y0 = ny - 1;
            chiptransmission[y0*nx + x0] -= fraction;
            if (chiptransmission[y0*nx + x0] < 0) {
                chiptransmission[y0*nx + x0] = 0.0;
            }
        }

    }

    elements = 1024;
    henyey_greenstein = new double[elements]();
    henyey_greenstein_mu = new double[elements]();
    double g = 0.6;
    for (int i = 0; i < elements; i++) {
        henyey_greenstein_mu[i] = -1 + 2.0*(static_cast<double>(i))/(elements - 1);
        henyey_greenstein[i] = (1 - g*g)/pow(1 + g*g - 2*g*henyey_greenstein_mu[i], 1.5);
    }
    double total = 0.0;
    for (int i = 0; i < elements; i++) {
        total += henyey_greenstein[i];
    }
    for (int i = 0; i < elements; i++) {
        henyey_greenstein[i] = henyey_greenstein[i]/total;
    }
    for (int i = 1; i < elements; i++) {
        henyey_greenstein[i] = henyey_greenstein[i] + henyey_greenstein[i - 1];
    }
    absorptionLength = 10.0;

    henyey_greenstein_w = new double[elements]();
    henyey_greenstein_mu_w = new double[elements]();
    g = 0.9;
    for (int i = 0; i < elements; i++) {
        henyey_greenstein_mu_w[i] = -1 + 2.0*(static_cast<double>(i))/(elements - 1);
        henyey_greenstein_w[i] = (1 - g*g)/pow(1 + g*g - 2*g*henyey_greenstein_mu_w[i], 1.5);
    }
    total = 0.0;
    for (int i = 0; i < elements; i++) {
        total += henyey_greenstein_w[i];
    }
    for (int i = 0; i < elements; i++) {
        henyey_greenstein_w[i] = henyey_greenstein_w[i]/total;
    }
    for (int i = 1; i < elements; i++) {
        henyey_greenstein_w[i] = henyey_greenstein_w[i] + henyey_greenstein_w[i - 1];
    }

    // status = 0;
    // fits_create_file(&faptr, "!cont.fits", &status);
    // naxesa[0] = nx; naxesa[1] = ny;
    // fits_create_img(faptr, LONG_IMG, 2, naxesa, &status);
    // fits_write_img(faptr, TLONG, 1, nx*ny, chiplistmap, &status);
    // fits_close_file(faptr, &status);

    // status = 0;
    // fits_create_file(&faptr, "!cont_t.fits", &status);
    // naxesa[0] = nx; naxesa[1] = ny;
    // fits_create_img(faptr, FLOAT_IMG, 2, naxesa, &status);
    // fits_write_img(faptr, TFLOAT, 1, nx*ny, chiptransmission, &status);
    // fits_close_file(faptr, &status);

}
