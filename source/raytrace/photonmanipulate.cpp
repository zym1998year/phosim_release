///
/// @package phosim
/// @file photonmanipulate.cpp
/// @brief photon manipulation routines (part of Image class)
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include "phosim_compat.h"

int Image::domeSeeing(Vector *angle, Photon *aph) {

    double phi, r;

    phi = 2*PI*random[aph->thread].uniform();
    r = sqrt(domeseeing*domeseeing + toypsf*toypsf)*ARCSEC/2.35482*random[aph->thread].normal();

    double dx = r*cos(phi);
    double dy = r*sin(phi);


    if (toyg1 != 0.0 || toyg2 != 0) {
        double dxp = dx*(1 + toyg1) - dy*(toyg2);
        double dyp = dy*(1 - toyg1) - dx*(toyg2);
        dx = - dxp*cos(rotatez) - dyp*sin(rotatez);
        dy = - dxp*sin(rotatez) + dyp*cos(rotatez);
    }

    angle->x = angle->x + dx;
    angle->y = angle->y + dy;
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);

}

int Image::tracking(Vector *angle, double time) {

    double rindex, vxp, vyp;
    int index;

    index = find_linear(perturbation.jittertime, trackinglines, time, &rindex);
    rotation(&vxp, &vyp, spiderangle, perturbation.jitterele[index], perturbation.jitterazi[index]);
    vxp = angle->x - vxp;
    vyp = angle->y + vyp;
    rotation(&(angle->x), &(angle->y), perturbation.jitterrot[index], vxp, vyp);
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    //wind shake
    double wx = perturbation.windshake[index]*cos(winddir[7]*DEGREE - azimuth);
    double wy = perturbation.windshake[index]*sin(winddir[7]*DEGREE - azimuth);

    angle->x += wx - wy*screen.jitterwind[index];
    angle->y += wy + wx*screen.jitterwind[index];
    angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    return(0);
}

int Image::largeAngleScattering(Vector *largeAngle, Photon *aph, int surf) {

    int index;
    double r, phi;

    double cosThetaI = surface.costotal[surf]/surface.coscount[surf];
    double scatteringFraction = 1 - exp(-pow((4*PI*perturbation.scatteringsigma[surf]/(aph->wavelength*1e-6)*cosThetaI),2.0));

    r = random[aph->thread].uniform();
    if (r < scatteringFraction) {
    redoScatter:;
        r = random[aph->thread].uniform();
        find(perturbation.scatteringprofile+surf*SCATTERING_POINTS, SCATTERING_POINTS, r, &index);
        r = ((double)(index))/(static_cast<double>(SCATTERING_POINTS))*SCATTERING_LIMIT*(aph->wavelength*1e-6/(perturbation.scatteringb[surf]));

        double phiInternal = 2.0*PI*random[aph->thread].uniform();
        double cosTheta = sqrt( cosThetaI*cosThetaI - r*r - 2*r*sqrt(1-cosThetaI*cosThetaI)*cos(phiInternal) );
        double factor = pow((cosThetaI+cosTheta),2.0)/pow((1+cosThetaI),2.0)*0.9*cosThetaI/cosTheta;
        if (random[aph->thread].uniform() > factor) goto redoScatter;

        phi = 2.0*PI*random[aph->thread].uniform();
        largeAngle->x += r*cos(phi);
        largeAngle->y += r*sin(phi);
    }

    return(0);

}

int Image::secondKick(Vector *largeAngle, Photon *aph) {

    int index;
    double r, phi;
    if (pupilscreenMode == 1 && atmospheremode < 2) {
        int i;
        int j;
        r = random[aph->thread].uniform();
        // 2D sample from normalized cumulative FFT distribution
        for (i = 0; i < SCREEN_SIZE; i++) {
            if (screen.focalscreencum[i*SCREEN_SIZE] > r) {
                i--;
                goto jumpi;
            }
        }
        i--;
    jumpi:;
        for (j = 0; j < SCREEN_SIZE; j++) {
            if (screen.focalscreencum[i*SCREEN_SIZE + j] > r) {
                j--;
                goto jumpj;
            }
        }
        j--;
    jumpj:;
        i -= SCREEN_SIZE/2;
        j -= SCREEN_SIZE/2;
        r = (static_cast<double>(sqrt(i*i + j*j)))*screen.pupilscreenscale*1e-3/
            (SCREEN_SIZE*screen.fine_sizeperpixel*screen.paddingfactor)*aph->wavelength;
        phi = static_cast<double>(atan2(j, i));
        largeAngle->x += r*cos(phi);
        largeAngle->y += r*sin(phi);
    } else {
        r = random[aph->thread].uniform();
        find(screen.hffunc, SCREEN_SIZE*2, r, &index);
        // If turbulence off, do pupil diffraction. Otherwise, do Kolmogorov diffraction.
        if (atmospheremode < 2) {
            r = ((double)(index))*0.5*1e-3/(surface.outerRadius[0]*screen.paddingfactor)*aph->wavelength;
        }
        else
            r = ((double)(index))*0.5*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel)*aph->wavelengthFactor;
        phi = 2*PI*random[aph->thread].uniform();
        largeAngle->x += r*cos(phi);
        largeAngle->y += r*sin(phi);
    }
    return(0);

}


int Image::diffraction(Vector *position, Vector angle, Vector *largeAngle, Photon *aph) {

    double distance;
    double xl, yl, zorig, r;
    int i;
    double diffdist, adiffdist;
    int signv = 0;
    double diffx = 0.0, diffy = 0.0, mindiffdist;

    mindiffdist = 1e30;
    zorig = position->z;
    double cosShiftedAngle = cos(aph->shiftedAngle);
    double sinShiftedAngle = sin(aph->shiftedAngle);
    double mindist = 1e-3;

    for (i = 0; i < obstruction.nspid; i++) {


        distance = (obstruction.height[i] - (position->z))/angle.z;
        propagate(position, angle, distance);
        xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
        yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

        if (obstruction.type[i] == 1) {

            if (obstruction.angle[i] != 0) {
                distance = (obstruction.height[i] +
                            (fabs(yl) - obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = xl - obstruction.center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = sinShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.height[i] - obstruction.depth[i] +
                        (fabs(yl) - obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = xl - obstruction.center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);

            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = cosShiftedAngle*diffdist/adiffdist;
                diffy = sinShiftedAngle*diffdist/adiffdist;
            }

        }

        if (obstruction.type[i] == 2) {

            if (obstruction.angle[i] != 0) {
                distance = (obstruction.height[i] +
                            (fabs(xl)-obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                            (position->z))/angle.z;
                propagate(position, angle, distance);
                xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
                yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;
            }

            diffdist = yl - obstruction.center[i];
            if (diffdist > 0) {
                signv = 1;
            } else {
                signv = -1;
            }
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = -sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }

            distance = (obstruction.height[i] - obstruction.depth[i] +
                        (fabs(xl)-obstruction.reference[i])*sin(DEGREE*obstruction.angle[i]) -
                        (position->z))/angle.z;
            propagate(position, angle, distance);
            xl = (position->x)*cosShiftedAngle + (position->y)*sinShiftedAngle;
            yl = (position->x)*(-sinShiftedAngle) + (position->y)*cosShiftedAngle;

            diffdist = yl - obstruction.center[i];
            if (diffdist > 0 && signv == -1) return(1);
            if (diffdist < 0 && signv == 1) return(1);
            adiffdist = fabs(diffdist);
            if (adiffdist < obstruction.width[i] + mindist) return(1);
            if (adiffdist - obstruction.width[i] < (mindiffdist)) {
                mindiffdist = fabs(diffdist) - obstruction.width[i];
                diffx = -sinShiftedAngle*diffdist/adiffdist;
                diffy = cosShiftedAngle*diffdist/adiffdist;
            }
        }
    }

    distance = (zorig - (position->z))/angle.z;
    propagate(position, angle, distance);

    if (obstruction.pupil == 1) {
        r = sqrt((position->x)*(position->x) + (position->y)*(position->y));
        diffdist = surface.outerRadius[0] - r;
        if (diffdist < mindist) return(1);
        if (diffdist < mindiffdist) {
            mindiffdist = diffdist;
            diffx = (position->x)/r;
            diffy = (position->y)/r;
        }
        diffdist = r - surface.innerRadius[0];
        if (diffdist < mindist) return(1);
        if (diffdist < mindiffdist) {
            mindiffdist = diffdist;
            diffx = (position->x)/r;
            diffy = (position->y)/r;
        }
    }

    if (mindiffdist < surface.outerRadius[0] - surface.innerRadius[0]) {
        r = (aph->wavelength/1000.0)/(4*PI*mindiffdist);
        if (diffractionMode != 5) {
            largeAngle->x += r*diffx;
            largeAngle->y += r*diffy;
        }
    }


    return(0);

}

double Image::airIndexRefraction(Photon *aph, int layer, int skipADC, double *ADC) {

    if (airrefraction) {
        if (layer == -1) {
            *ADC = 0.0;
            return(0.0);
        } else {


            double currentheight=height[layer]+groundlevel/1e3;
            int index;
            find_v(air.altitudeProfile, currentheight, &index);
            double temp = interpolate_v(air.temperatureProfile,air.altitudeProfile,currentheight,index);
            double press =  interpolate_v(air.pressureProfile,air.altitudeProfile,currentheight,index)/PRESS_ATMOSPHERE*PRESS_ATMOSPHERE_MMHG;
            double water =  interpolate_v(air.h2oProfile,air.altitudeProfile,currentheight,index)*press;

 
            double sigma=1.0/aph->wavelength;
            double ps=press/760.00*1013.25;
            double pw=water/760.00*1013.25;
            double dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
            double ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
            double n=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
            n=n+(6478.31+58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
            n = 1e-8*n;

            if (skipADC == 0) {
                sigma=1.0/centralwavelength;
                ps=press/760.00*1013.25;
                pw=water/760.00*1013.25;
                dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp;
                ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp;
                *ADC=(2371.34+683939.7/(130.0-pow(sigma,2))+4547.3/(38.9-pow(sigma,2)))*ds;
                *ADC=*ADC+(6478.31+58.058*pow(sigma,2)-0.71150*pow(sigma,4)+0.08851*pow(sigma,6))*dw;
                *ADC = 1e-8*(*ADC);
            }

            return(n);
        }

    } else {
        return(0.0);
    }

}

int Image::atmosphericDispersion (Vector *angle, Photon *aph, int layer) {

    double dx, dy, phi, theta, psi;
    double signa=1.0;
    double signb=1.0;

    if (atmospheric_dispersion) {
    if (layer == -1) {
        double nn, nn2;
        if (atmosphericdispcenter) {
            nn=airIndexRefraction(aph,natmospherefile-1,0,&nn2);
        } else {
            nn=airIndexRefraction(aph,natmospherefile-1,1,&nn2);
        }
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(-angle->x) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(-angle->y) +
            (sin(psi)*sin(theta))*fabs(angle->z);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(-angle->x) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(-angle->y) +
            (cos(psi)*sin(theta))*fabs(angle->z);
        aph->refractionInvariant=(1+nn)*(height[natmospherefile-1]+groundlevel/1e3+RADIUS_EARTH)*(sqrt(dx*dx + dy*dy));
        if (atmosphericdispcenter) {
            phi = aph->shiftedAngle*signa;
            theta = zenith;
            psi = -aph->shiftedAngle*signb;
            dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(0.0) +
                (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(0.0) +
                (sin(psi)*sin(theta))*fabs(1.0);
            dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(0.0) +
                (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(0.0) +
                (cos(psi)*sin(theta))*fabs(1.0);
            aph->refractionInvariantCenter=(1+nn2)*(height[natmospherefile-1]+groundlevel/1e3+RADIUS_EARTH)*(sqrt(dx*dx + dy*dy));
            aph->adcx = 0.0;
            aph->adcy = 0.0;
            aph->adcz = smallAnglePupilNormalize(aph->adcx, aph->adcy);
        }
    } else {

        double currentheight=height[layer]+groundlevel/1e3+RADIUS_EARTH;
        double correctionCenter=0.0;
        if (atmosphericdispcenter) {
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(0.0) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(0.0) +
            (sin(psi)*sin(theta))*fabs(1.0);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(0.0) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(0.0) +
            (cos(psi)*sin(theta))*fabs(1.0);
        correctionCenter=aph->refractionInvariantCenter/sqrt(pow(currentheight*(1+aph->airRefractionADC),2.0) -
                                                             pow(aph->refractionInvariantCenter,2.0))*
            (aph->airRefractionADC - aph->airRefractionPreviousADC)/(1 + aph->airRefractionADC);
            if ((dx == 0.0) && (dy == 0.0)) {
                //                aph->adcx = correctionCenter;
                aph->adcx = tan(correctionCenter);
                aph->adcy = 0.0;
            } else {
                aph->adcx = tan(correctionCenter)*dx/sqrt(dx*dx + dy*dy);
                aph->adcy = tan(correctionCenter)*dy/sqrt(dx*dx + dy*dy);
                //                aph->adcx = correctionCenter*dx/sqrt(dx*dx + dy*dy);
                //aph->adcy = correctionCenter*dy/sqrt(dx*dx + dy*dy);
            }
            aph->adcz = smallAnglePupilNormalize(aph->adcx, aph->adcy);
        }
        double correction=aph->refractionInvariant/sqrt(pow(currentheight*(1+aph->airRefraction),2.0) -
                                                        pow(aph->refractionInvariant,2.0))*
            (aph->airRefraction - aph->airRefractionPrevious)/(1 + aph->airRefraction);
        phi = aph->shiftedAngle*signa;
        theta = zenith;
        psi = -aph->shiftedAngle*signb;
        dx = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*(-angle->x) +
            (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*(-angle->y) +
            (sin(psi)*sin(theta))*fabs(angle->z);
        dy = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*(-angle->x) +
            (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*(-angle->y) +
            (cos(psi)*sin(theta))*fabs(angle->z);
        if ((dx == 0.0) && (dy == 0.0)) {
            //            angle->x += correction;
            angle->x += tan(correction);
        } else {
            //            angle->x += correction*dx/sqrt(dx*dx + dy*dy);
            //angle->y += correction*dy/sqrt(dx*dx + dy*dy);
            angle->x += tan(correction)*dx/sqrt(dx*dx + dy*dy);
            angle->y += tan(correction)*dy/sqrt(dx*dx + dy*dy);
        }
        // printf("a %ld %e\n",layer,(correction-correctionCenter)/ARCSEC);
        if (atmosphericdispcenter) {
            angle->x -= aph->adcx;
            angle->y -= aph->adcy;
        }
        // printf("b %ld %e %e %e\n",layer,correctionCenter,aph->refractionInvariantCenter,aph->airRefraction);
        // printf("%e %e %e %e\n",correction,correctionCenter,aph->airRefraction,aph->airRefractionADC);
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }
    }
    // if (layer == 0) {

    //     if (atmosphericdispcenter) {
    //         dx = zenith*sin(aph->shiftedAngle);
    //         dy = zenith*cos(aph->shiftedAngle);
    //         adcx = tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*air.air_refraction_adc;
    //         adcy = tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*air.air_refraction_adc;
    //         if (zenith == 0.0) {
    //             adcx = 0.0;
    //             adcy = 0.0;
    //         }
    //     }
    //     if (atmospheric_dispersion) {
    //         dx = -angle->x + zenith*sin(aph->shiftedAngle);
    //         dy = -angle->y + zenith*cos(aph->shiftedAngle);
    //         angle->x = angle->x + tan(sqrt(dx*dx + dy*dy))*dx/sqrt(dx*dx + dy*dy)*aph->airRefraction;
    //         angle->y = angle->y + tan(sqrt(dx*dx + dy*dy))*dy/sqrt(dx*dx + dy*dy)*aph->airRefraction;
    //         if (atmosphericdispcenter) {
    //             angle->x = angle->x - adcx;
    //             angle->y = angle->y - adcy;
    //         }
    //         angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    //     }

    // }
    return(0);

}


int Image::samplePupil (Vector *position, long ray, Photon *aph) {

    double r, phi; // mm, radians

    if (aperturemode == 1) {
        double x, y;
        x = ((static_cast<double>((ray % 40000) % 200))/200.0)*2.0*maxr - maxr;
        y = ((static_cast<double>((ray % 40000) / 200))/200.0)*2.0*maxr - maxr;
        r = sqrt(x*x + y*y);
        phi = atan2(y, x);
        if (r < minr || r > maxr) return(1);
    } else if (aperturemode == 2) {
        if (ray < 2) {
            r = 1e-14;
            phi = 0.0;
        } else {
            double x, y;
            long opdsize2 = opdsize*opdsize*opdsampling*opdsampling;
            x = ((static_cast<double>((ray % opdsize2) % (opdsampling*opdsize)))/(opdsampling*opdsize - 1))*2.0*maxr - maxr;
            y = ((static_cast<double>((ray % opdsize2) / (opdsampling*opdsize)))/(opdsampling*opdsize - 1))*2.0*maxr - maxr;
            r = sqrt(x*x + y*y);
            phi = atan2(y, x);
            if (r < minr || r > maxr) return(1);
        }
    } else {
        r = sqrt(random[aph->thread].uniform()*(maxr*maxr - minr*minr) + minr*minr);
        phi = random[aph->thread].uniform()*2*PI;
    }

    int surfaceIndex = 0; // start at first surface in optics file

    double R = surface.radiusCurvature[surfaceIndex]; // radius of curvature (mm)
    double asphere; // sum of asphere coefficients (mm)
    double third = -surface.three[surfaceIndex]*1e3; // 3rd order asphere coefficient (mm)
    double fourth = -surface.four[surfaceIndex]*1e3; // 4th order asphere coefficient (mm)
    double fifth = -surface.five[surfaceIndex]*1e3; // 5th order asphere coefficient (mm)
    double sixth = -surface.six[surfaceIndex]*1e3; // 6th order asphere coefficient (mm)
    double seventh = -surface.seven[surfaceIndex]*1e3; // 7th order asphere coefficient (mm)
    double eighth = -surface.eight[surfaceIndex]*1e3; // 8th order asphere coefficient (mm)
    double ninth = -surface.nine[surfaceIndex]*1e3; // 9th order asphere coefficient (mm)
    double tenth = -surface.ten[surfaceIndex]*1e3; // 10th order asphere coefficient (mm)

    asphere = third*pow(r, 3.0) + fourth*pow(r, 4.0) + fifth*pow(r, 5.0) + sixth*pow(r, 6.0) + seventh*pow(r, 7.0) + eighth*pow(r, 8.0) + ninth*pow(r,9.0) + tenth*pow(r, 10.0);

    position->x = r*cos(phi);
    position->y = r*sin(phi);
    if (R != 0.0) {
        double k = surface.conic[surfaceIndex]; // conic constant (unitless)
        // sagitta equation
        position->z = surface.height[surfaceIndex] + r*r/(R*(1.0 + sqrt(1.0 - (1.0 + k)*r*r/(R*R)))) + asphere;
    } else {
        position->z = surface.height[surfaceIndex] + asphere;
    }

    return(0);

}


inline int Image::transmissionCheck(float transmission, int surfaceIndex, int waveSurfaceIndex, Photon *aph) {

    (aph->counter)++;

    if (aph->ghostFlag > 0) {
        int flag = aph->ghostFlag;
        if (flag > 3) flag = 3;
        waveSurfaceIndex += flag*state.waveLimit*state.surfaceLimit;
    }
    if (transmission > state.dynamicTransmission[waveSurfaceIndex + surfaceIndex]) {
        if ((!std::isinf(transmission)) && (!std::isnan(transmission))) {
            state.dynamicTransmission[waveSurfaceIndex + surfaceIndex] = transmission;
            state.dynamicTransmissionInit[waveSurfaceIndex + surfaceIndex] += 1;
            // int maxw=static_cast<int>((maxwavelength - minwavelength));
            // int d=0;
            //  for (int j=0; j< 65;j++) {
            //      for (int i =0; i < maxw; i++) {
            //          int q=state.dynamicTransmissionInit[i*65+j];
            //          if (q==0) d++;
            //      }
            //  }
            //  if ((d%100)==0) printf("R %d\n",d);

        }
    }
    if (transmission < state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex]) {
        if ((!std::isinf(transmission)) && (!std::isnan(transmission))) {
            state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex] = transmission;
            state.dynamicTransmissionInit[waveSurfaceIndex + surfaceIndex] += 1;
        }
    }
    if (aph->counter <= aph->maxcounter) {
        aph->transmissionCheckRand = aph->saveRand[aph->counter];
    } else {
        aph->transmissionCheckRand = random[aph->thread].uniform();
    }
    if (aph->transmissionCheckRand > transmission) {
        return(1);
    } else {
        return(0);
    }

}

inline int Image::transmissionPreCheck(int index, Photon *aph, int rewrite) {

    (aph->counter)++;
    if (rewrite==0) aph->saveRand[aph->counter] = random[aph->thread].uniform();

    // float f= state.dynamicTransmission[index] ;
    // int d = state.dynamicTransmissionInit[index];
    // if (d==0) {
    //     int waveIndex = static_cast<int>((aph->wavelength*1000 - minwavelength)*10);
    //     int waveSurfaceIndex = waveIndex*65;
    //     printf("Y %d %f %d %lf\n",index-waveSurfaceIndex,f,d,aph->wavelength);
    // }

    if (aph->preGhostFlag > 0) {
        int flag = aph->preGhostFlag;
        if (flag > 3) flag = 3;
        index += flag*state.waveLimit*state.surfaceLimit;
    }
    if (state.dynamicTransmissionInit[index]==0) return(2);
    if (aph->saveRand[aph->counter] > (state.dynamicTransmission[index] + transtol)) {
        return(1);
    } else {
        if (aph->saveRand[aph->counter] < (state.dynamicTransmissionLow[index] - transtol)) {
            return(0);
        } else {
            return(2);
        }
    }

}

inline int Image::transmissionPreCheckZero(int index, Photon *aph) {

    (aph->counter)++;
    aph->saveRand[aph->counter] = random[aph->thread].uniform();

    // float f1=state.dynamicTransmission[index];
    // float f2=state.dynamicTransmissionLow[index];
    // int d1=state.dynamicTransmissionInit[index];
    // printf("%d %f %f %d\n",index,f1,f2,d1);
    if (state.dynamicTransmissionInit[index]==0) return(2);
    if (aph->saveRand[aph->counter] > (state.dynamicTransmission[index] + transtol)) {
        return(1);
    } else {
        if (aph->saveRand[aph->counter] < (state.dynamicTransmissionLow[index] - transtol)) {
            return(0);
        } else {
            return(2);
        }
    }

}

inline int Image::transmissionPreCheckSuppression(int index, Photon *aph, int rewrite) {

    (aph->counter)++;
    if (rewrite==0) aph->saveRand[aph->counter] = random[aph->thread].uniform();

    if (aph->preGhostFlag > 0) {
        int flag = aph->preGhostFlag;
        if (flag > 3) flag = 3;
        index += flag*state.waveLimit*state.surfaceLimit;
    }
    if (state.dynamicTransmissionInit[index]==0) return(2);
    if (aph->saveRand[aph->counter] > tauSuppression(state.dynamicTransmission[index] + transtol)) {
        return(1);
    } else {
        if (aph->saveRand[aph->counter] < tauSuppression(state.dynamicTransmissionLow[index] - transtol)) {
            return(0);
        } else {
            return(2);
        }
    }

}

inline int Image::transmissionCheckNoRoll(float transmission, int surfaceIndex, int waveSurfaceIndex, Photon *aph) {

    float randNum;

    if (aph->ghostFlag > 0) {
        int flag = aph->ghostFlag;
        if (flag > 3) flag = 3;
        waveSurfaceIndex += flag*state.waveLimit*state.surfaceLimit;
    }
    if (transmission > state.dynamicTransmission[waveSurfaceIndex + surfaceIndex]) {
        if ((!std::isinf(transmission)) && (!std::isnan(transmission))) {
            state.dynamicTransmission[waveSurfaceIndex + surfaceIndex] = transmission;
            state.dynamicTransmissionInit[waveSurfaceIndex + surfaceIndex] += 1;
        }
    }
    if (transmission < state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex]) {
        if ((!std::isinf(transmission)) && (!std::isnan(transmission))) {
            state.dynamicTransmissionLow[waveSurfaceIndex + surfaceIndex] = transmission;
            state.dynamicTransmissionInit[waveSurfaceIndex + surfaceIndex] += 1;
        }
    }
    if (aph->counter <= aph->maxcounter) {
        randNum = aph->saveRand[aph->counter];
    } else {
        randNum = random[aph->thread].uniform();
    }
    if (randNum > transmission) {
        return(1);
    } else {
        return(0);
    }

}

inline int Image::transmissionPreCheckNoRoll(int index, Photon *aph) {
    
    if (aph->preGhostFlag > 0) {
        int flag = aph->preGhostFlag;
        if (flag > 3) flag = 3;
        index += flag*state.waveLimit*state.surfaceLimit;
    }
    if (state.dynamicTransmissionInit[index]==0) return(2);
    if (aph->saveRand[aph->counter] > (state.dynamicTransmission[index] + transtol)) {
        return(1);
    } else {
        if (aph->saveRand[aph->counter] < (state.dynamicTransmissionLow[index] - transtol)) {
            return(0);
        } else {
            return(2);
        }
    }

}


double Image::surfaceCoating (double wavelength, Vector angle, Vector normal, int newSurf, float *reflection, Photon *aph) {

    double filterAngle, rindex, crindex;
    int index, cindex;

    if (surface.surfacecoating[newSurf] != 0 && coatingmode == 1) {
        if (coating.angleNumber[newSurf] > 1) {
            if (coating.wavelengthNumber[newSurf] > 1) {
                Vector tempangle;
                vectorCopy(angle, &tempangle);
                refract(&tempangle, normal, aph->ncurr, 1.0);
                double arg = fabs(normal.x*tempangle.x + normal.y*tempangle.y + normal.z*tempangle.z);
                if (arg > 1) arg = 1.0;
                filterAngle = acos(arg)/DEGREE;
                index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
                cindex = find_linear(coating.angle[newSurf], coating.angleNumber[newSurf], filterAngle, &crindex);
                if (crindex - cindex < 0.0) crindex = static_cast<double>(cindex);
                if (crindex - cindex > 1.0) crindex = static_cast<double>(cindex + 1);
                if (rindex - index < 0.0) rindex = static_cast<double>(index);
                if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
                //                if (newSurf==1) printf("%lf %d %lf\n",filterAngle,cindex,interpolate_bilinear(coating.transmission[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex));
                *reflection = interpolate_bilinear(coating.reflection[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex);
                return(interpolate_bilinear(coating.transmission[newSurf], coating.wavelengthNumber[newSurf], cindex, crindex, index, rindex));
            } else {
                Vector tempangle;
                vectorCopy(angle, &tempangle);
                refract(&tempangle, normal, aph->ncurr, 1.0);
                double arg = fabs(normal.x*tempangle.x + normal.y*tempangle.y + normal.z*tempangle.z);
                if (arg > 1) arg = 1.0;
                filterAngle = acos(arg)/DEGREE;
                cindex = find_linear(coating.angle[newSurf], coating.angleNumber[newSurf], filterAngle, &crindex);
                if (crindex - cindex < 0.0) crindex = static_cast<double>(cindex);
                if (crindex - cindex > 1.0) crindex = static_cast<double>(cindex + 1);
                *reflection = interpolate_linear(coating.reflection[newSurf], cindex, crindex);
                return(interpolate_linear(coating.transmission[newSurf], cindex, crindex));
            }
        } else {
            if (coating.wavelengthNumber[newSurf]==1) {
                *reflection = *(coating.reflection[newSurf]);
                return(*(coating.transmission[newSurf]));
            } else {
                index = find_linear(coating.wavelength[newSurf], coating.wavelengthNumber[newSurf], wavelength, &rindex);
                if (rindex - index < 0.0) rindex = static_cast<double>(index);
                if (rindex - index > 1.0) rindex = static_cast<double>(index + 1);
                *reflection = interpolate_linear(coating.reflection[newSurf], index, rindex);
                return(interpolate_linear(coating.transmission[newSurf], index, rindex));
            }
        }
    } else {
        return(1.0);
    }

}

void Image::newRefractionIndex(int surfaceIndex, Photon *aph) {

    int index, newMedium;
    double rindex;

    aph->nprev = aph->ncurr;
    if (aph->direction == 1) {
        newMedium = surfaceIndex;
    } else {
        newMedium = surfaceIndex - 1;
    }
    if (newMedium >= 0) {
        if (surface.surfacemed[newMedium] == 0) {
            aph->ncurr = 1.0;
        } else if (surface.surfacemed[newMedium] == 2) {
            aph->ncurr = 1.0 + aph->airRefraction;
        } else {
            index = find_linear(medium.indexRefractionWavelength[newMedium], medium.indexRefractionNumber[newMedium], aph->wavelength, &rindex);
            aph->ncurr = interpolate_linear(medium.indexRefraction[newMedium], index, rindex);
        }
    } else {
        aph->ncurr = 1.0 + aph->airRefraction;
    }

}


void Image::atmospherePropagate(Vector *position, Vector angle, int layer, int mode, Photon *aph) {

    if (layer == -1) {

        propagate(position, angle, (1e6*air.topAtmosphere - position->z)/(angle.z));
        aph->xporig = position->x;
        aph->yporig = position->y;
        aph->zporig = position->z;

        if (mode == 2) {

            if (fabs(aph->time - aph->prtime) > screentol) {
                for (long k = 0; k < SCREEN_SIZE; k++) {
                    for (long l = 0; l < SCREEN_SIZE; l++) {
                        long index = aph->thread*SCREEN_SIZE*SCREEN_SIZE;
                        *(screen.phasescreen + index + k*SCREEN_SIZE + l) = 0;
                    }
                }
            }

        }

    } else {

        double distance = (1e6*height[layer] - position->z)/(angle.z);
        propagate(position, angle, distance);

    }

}

void Image::atmosphereIntercept(Vector *position, int layer, Photon *aph) {

    double rindex;
    int index;
    double wx = wind[layer]*1.0e3*(aph->absoluteTime)*cos(winddir[layer]*DEGREE - azimuth);
    double wy = wind[layer]*1.0e3*(aph->absoluteTime)*sin(winddir[layer]*DEGREE - azimuth);


    // wind blur
    index = find_linear(perturbation.jittertime, trackinglines, aph->absoluteTime, &rindex);

    aph->windx = wx - wy*screen.jitterwind[index]*DEGREE;
    aph->windy = wy + wx*screen.jitterwind[index]*DEGREE;

    aph->lindex = SCREEN_SIZE*SCREEN_SIZE*layer;

    aph->xpos = position->x + xtelloc + aph->windx;
    aph->ypos = position->y + ytelloc + aph->windy;

    find_linear_wrapFloat(aph->xpos, screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexlx0), &(aph->indexlx1), &(aph->dlx));
    find_linear_wrapFloat(aph->ypos, screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexly0), &(aph->indexly1), &(aph->dly));
    find_linear_wrapFloat(aph->xpos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcx0), &(aph->indexcx1), &(aph->dcx));
    find_linear_wrapFloat(aph->ypos, screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcy0), &(aph->indexcy1), &(aph->dcy));
    find_linear_wrapFloat(aph->xpos, screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmx0), &(aph->indexmx1), &(aph->dmx));
    find_linear_wrapFloat(aph->ypos, screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmy0), &(aph->indexmy1), &(aph->dmy));
    find_linear_wrapFloat(aph->xpos, screen.fine_sizeperpixel, SCREEN_SIZE, &(aph->indexfx0), &(aph->indexfx1), &(aph->dfx));
    find_linear_wrapFloat(aph->ypos, screen.fine_sizeperpixel, SCREEN_SIZE, &(aph->indexfy0), &(aph->indexfy1), &(aph->dfy));


}

void Image::atmosphereRefraction (Vector *angle, int layer, int mode, Photon *aph) {

    double scaleOuter;

    if (mode == 1 || mode == 5) {
        scaleOuter = aph->wavelengthFactor*ARCSEC*screen.secondKickSize;
    } else {
        scaleOuter = aph->wavelengthFactor*ARCSEC;
    }

    if (mode <= 1 || mode >= 4) {

        if (atmdebug == 0) {

            (angle->x) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeX + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                           aph->dlx, aph->indexly0, aph->indexly1, aph->dly) +
                           interpolate_bilinear_float_wrap(screen.turbulenceCoarseX + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                           aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy) +
                           interpolate_bilinear_float_wrap(screen.turbulenceMediumX + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                           aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy))*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeY + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                           aph->dlx, aph->indexly0, aph->indexly1, aph->dly) +
                           interpolate_bilinear_float_wrap(screen.turbulenceCoarseY + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                           aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy) +
                           interpolate_bilinear_float_wrap(screen.turbulenceMediumY + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                           aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy))*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);
            
        } else {

            aph->indexlx0 = static_cast<int>(floor((aph->indexlx0)/largeGrid)*largeGrid);
            aph->indexly0 = static_cast<int>(floor((aph->indexly0)/largeGrid)*largeGrid);
            aph->indexcx0 = static_cast<int>(floor((aph->indexcx0)/coarseGrid)*coarseGrid);
            aph->indexcy0 = static_cast<int>(floor((aph->indexcy0)/coarseGrid)*coarseGrid);
            aph->indexmx0 = static_cast<int>(floor((aph->indexmx0)/mediumGrid)*mediumGrid);
            aph->indexmy0 = static_cast<int>(floor((aph->indexmy0)/mediumGrid)*mediumGrid);
            aph->indexfx0 = static_cast<int>(floor((aph->indexfx0)/fineGrid)*fineGrid);
            aph->indexfy0 = static_cast<int>(floor((aph->indexfy0)/fineGrid)*fineGrid);
            aph->indexlx1 = static_cast<int>(floor((aph->indexlx1)/largeGrid)*largeGrid);
            aph->indexly1 = static_cast<int>(floor((aph->indexly1)/largeGrid)*largeGrid);
            aph->indexcx1 = static_cast<int>(floor((aph->indexcx1)/coarseGrid)*coarseGrid);
            aph->indexcy1 = static_cast<int>(floor((aph->indexcy1)/coarseGrid)*coarseGrid);
            aph->indexmx1 = static_cast<int>(floor((aph->indexmx1)/mediumGrid)*mediumGrid);
            aph->indexmy1 = static_cast<int>(floor((aph->indexmy1)/mediumGrid)*mediumGrid);
            aph->indexfx1 = static_cast<int>(floor((aph->indexfx1)/fineGrid)*fineGrid);
            aph->indexfy1 = static_cast<int>(floor((aph->indexfy1)/fineGrid)*fineGrid);

            (angle->x) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeX + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseX + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumX + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale)*scaleOuter;

            (angle->y) += (interpolate_bilinear_float_wrap(screen.turbulenceLargeY + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                         aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceCoarseY + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                         aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                         interpolate_bilinear_float_wrap(screen.turbulenceMediumY + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                         aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale)*scaleOuter;

            angle->z = smallAnglePupilNormalize(angle->x, angle->y);


        }

    } else {

        if (fabs(aph->time - aph->prtime) > screentol) {
            double randomi = random[aph->thread].uniform();
            double randomj = random[aph->thread].uniform();
            for (int i = 0; i < SCREEN_SIZE; i++) {
                for (int j = 0; j < SCREEN_SIZE; j++) {

                    find_linear_wrapFloat(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexlx0), &(aph->indexlx1), &(aph->dlx));
                    find_linear_wrapFloat(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcx0), &(aph->indexcx1), &(aph->dcx));
                    find_linear_wrapFloat(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmx0), &(aph->indexmx1), &(aph->dmx));
                    find_linear_wrapFloat(aph->xpos - aph->xp + (i + randomi - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE, &(aph->indexfx0), &(aph->indexfx1), &(aph->dfx));
                    find_linear_wrapFloat(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.large_sizeperpixel, SCREEN_SIZE, &(aph->indexly0), &(aph->indexly1), &(aph->dly));
                    find_linear_wrapFloat(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.coarse_sizeperpixel, SCREEN_SIZE, &(aph->indexcy0), &(aph->indexcy1), &(aph->dcy));
                    find_linear_wrapFloat(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.medium_sizeperpixel, SCREEN_SIZE, &(aph->indexmy0), &(aph->indexmy1), &(aph->dmy));
                    find_linear_wrapFloat(aph->ypos - aph->yp + (j + randomj - 0.5 - ((double)(SCREEN_SIZE/2) - 0.5))*screen.fine_sizeperpixel,
                                     screen.fine_sizeperpixel,  SCREEN_SIZE,  &(aph->indexfy0),  &(aph->indexfy1), &(aph->dfy));

                    if (atmdebug == 0) {

                        if (mode == 2) {

                            long index = aph->thread*SCREEN_SIZE*SCREEN_SIZE;
                            *(screen.phasescreen + index + i*SCREEN_SIZE + j) +=

                                (interpolate_bilinear_float_wrap(screen.phaseLarge + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                                 aph->dlx, aph->indexly0, aph->indexly1, aph->dly)+
                                 interpolate_bilinear_float_wrap(screen.phaseCoarse + aph->lindex,SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                                 aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)+
                                 interpolate_bilinear_float_wrap(screen.phaseMedium + aph->lindex,SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)+
                                 interpolate_bilinear_float_wrap(screen.phaseFine + aph->lindex  ,SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {

                            long index = aph->thread*SCREEN_SIZE*SCREEN_SIZE;
                            *(screen.phasescreen + index + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseMediumH + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)+
                                 interpolate_bilinear_float_wrap(screen.phaseFineH + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx,  aph->indexfy0, aph->indexfy1, aph->dfy))*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }

                    } else {

                        aph->indexlx0 = static_cast<int>(floor((aph->indexlx0)/largeGrid)*largeGrid);
                        aph->indexly0 = static_cast<int>(floor((aph->indexly0)/largeGrid)*largeGrid);
                        aph->indexcx0 = static_cast<int>(floor((aph->indexcx0)/coarseGrid)*coarseGrid);
                        aph->indexcy0 = static_cast<int>(floor((aph->indexcy0)/coarseGrid)*coarseGrid);
                        aph->indexmx0 = static_cast<int>(floor((aph->indexmx0)/mediumGrid)*mediumGrid);
                        aph->indexmy0 = static_cast<int>(floor((aph->indexmy0)/mediumGrid)*mediumGrid);
                        aph->indexfx0 = static_cast<int>(floor((aph->indexfx0)/fineGrid)*fineGrid);
                        aph->indexfy0 = static_cast<int>(floor((aph->indexfy0)/fineGrid)*fineGrid);
                        aph->indexlx1 = static_cast<int>(floor((aph->indexlx1)/largeGrid)*largeGrid);
                        aph->indexly1 = static_cast<int>(floor((aph->indexly1)/largeGrid)*largeGrid);
                        aph->indexcx1 = static_cast<int>(floor((aph->indexcx1)/coarseGrid)*coarseGrid);
                        aph->indexcy1 = static_cast<int>(floor((aph->indexcy1)/coarseGrid)*coarseGrid);
                        aph->indexmx1 = static_cast<int>(floor((aph->indexmx1)/mediumGrid)*mediumGrid);
                        aph->indexmy1 = static_cast<int>(floor((aph->indexmy1)/mediumGrid)*mediumGrid);
                        aph->indexfx1 = static_cast<int>(floor((aph->indexfx1)/fineGrid)*fineGrid);
                        aph->indexfy1 = static_cast<int>(floor((aph->indexfy1)/fineGrid)*fineGrid);

                        if (mode == 2) {

                            long index = aph->thread*SCREEN_SIZE*SCREEN_SIZE;
                            *(screen.phasescreen + index + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseLarge + aph->lindex, SCREEN_SIZE, aph->indexlx0, aph->indexlx1,
                                                                 aph->dlx, aph->indexly0, aph->indexly1, aph->dly)*largeScale +
                                 interpolate_bilinear_float_wrap(screen.phaseCoarse + aph->lindex, SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                                 aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)*coarseScale +
                                 interpolate_bilinear_float_wrap(screen.phaseMedium + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale +
                                 interpolate_bilinear_float_wrap(screen.phaseFine + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy)*fineScale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        } else {
                            long index = aph->thread*SCREEN_SIZE*SCREEN_SIZE;
                            *(screen.phasescreen + index + i*SCREEN_SIZE + j) +=
                                (interpolate_bilinear_float_wrap(screen.phaseMediumH + aph->lindex, SCREEN_SIZE, aph->indexmx0, aph->indexmx1,
                                                                 aph->dmx, aph->indexmy0, aph->indexmy1, aph->dmy)*mediumScale +
                                 interpolate_bilinear_float_wrap(screen.phaseFineH + aph->lindex, SCREEN_SIZE, aph->indexfx0, aph->indexfx1,
                                                                 aph->dfx, aph->indexfy0, aph->indexfy1, aph->dfy)*fineScale)*
                                screen.phase_norm[layer]*(seefactor[layer]/(totalseeing*pow(1/cos(zenith), 0.6)/2.35));

                        }


                    }
                }
            }
        }

    }

}


double Image::cloudOpacity (int layer, Photon *aph) {

    double transmission, transmissionC, airmassl, densityScale5, densityScale6, densityScale7, densityScale8;

    transmission = 1.0;
    transmissionC = 1.0;
    if (aph->dvr == 0.0) {
        airmassl = 1.0;
    } else {
        airmassl = RADIUS_EARTH*(sqrt(pow(1 + height[layer]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                 sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
            sin(aph->dvr)/(height[layer] - height[layer + 1]);
    }
    int index;
    find(mie.muWavelengthWater, MIE_IN_POINT, aph->wavelength, &index);
    double wavelengthF4 = interpolate(mie.totalExtinctionWater, mie.muWavelengthWater, aph->wavelength, index);
    find(mie.muWavelengthSeasalt, MIE_IN_POINT, aph->wavelength, &index);
    double wavelengthF0 = interpolate(mie.totalExtinctionSeasalt, mie.muWavelengthSeasalt, aph->wavelength, index);
    find(mie.muWavelengthDust, MIE_IN_POINT, aph->wavelength, &index);
    double wavelengthF1 = interpolate(mie.totalExtinctionDust, mie.muWavelengthDust, aph->wavelength, index);
    find(mie.muWavelengthSmoke, MIE_IN_POINT, aph->wavelength, &index);
    double wavelengthF2 = interpolate(mie.totalExtinctionSmoke, mie.muWavelengthSmoke, aph->wavelength, index);
    find(mie.muWavelengthPollution, MIE_IN_POINT, aph->wavelength, &index);
    double wavelengthF3 = interpolate(mie.totalExtinctionPollution, mie.muWavelengthPollution, aph->wavelength, index);
    densityScale5 = 1.0 + aerosolSeasaltgradient*(aph->xpos*cos(aerosolSeasaltAngle) + aph->ypos*sin(aerosolSeasaltAngle))/(1000000.0);
    densityScale6 = 1.0 + aerosolDustgradient*(aph->xpos*cos(aerosolDustAngle) + aph->ypos*sin(aerosolDustAngle))/(1000000.0);
    densityScale7 = 1.0 + aerosolSmokegradient*(aph->xpos*cos(aerosolSmokeAngle) + aph->ypos*sin(aerosolSmokeAngle))/(1000000.0);
    densityScale8 = 1.0 + aerosolPollutiongradient*(aph->xpos*cos(aerosolPollutionAngle) + aph->ypos*sin(aerosolPollutionAngle))/(1000000.0);

    transmission *= exp(-(*(air.tau[12*(layer + 1) + 7] + aph->oindex)*densityScale5*airmassl/air.airmassLayer[layer + 1]*wavelengthF0));
    aph->transmissionSave[0] = transmission;
    transmission *= exp(-(*(air.tau[12*(layer + 1) + 8] + aph->oindex)*densityScale6*airmassl/air.airmassLayer[layer + 1]*wavelengthF1));
    aph->transmissionSave[1] = transmission;
    transmission *= exp(-(*(air.tau[12*(layer + 1) + 9] + aph->oindex)*densityScale7*airmassl/air.airmassLayer[layer + 1]*wavelengthF2));
    aph->transmissionSave[2] = transmission;
    transmission*= exp(-(*(air.tau[12*(layer + 1) + 10] + aph->oindex)*densityScale8*airmassl/air.airmassLayer[layer + 1]*wavelengthF3));
    aph->transmissionSave[3] = transmission;

    if (cloudmean[layer] != 0 || cloudvary[layer] != 0) {

        transmissionC = pow(10.0, -0.4*(cloudmean[layer] +
                                        cloudvary[layer]*(interpolate_bilinear_float_wrap(screen.cloud[layer], SCREEN_SIZE, aph->indexcx0, aph->indexcx1,
                                                                                          aph->dcx, aph->indexcy0, aph->indexcy1, aph->dcy)))*wavelengthF4*airmassl/air.airmassLayer[layer + 1]);
        if (transmissionC > 1.0) transmissionC = 1.0;
    }
    transmission *= transmissionC;
    aph->transmissionSave[4] = transmission;
    return(transmission);
}

double Image::cloudOpacityMoon (int layer, Photon *aph) {

    double transmission;

    transmission = 1.0;
    transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 7] + aph->oindex)));
    aph->transmissionSave[0] = transmission;
    transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 8] + aph->oindex)));
    aph->transmissionSave[1] = transmission;
    transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 9] + aph->oindex)));
    aph->transmissionSave[2] = transmission;
    transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 10] + aph->oindex)));
    aph->transmissionSave[3] = transmission;
    transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 11] + aph->oindex)));
    aph->transmissionSave[4] = transmission;

    return(transmission);
}

double Image::cloudOpacitySun (int layer, Photon *aph) {

    double transmission;

    transmission = 1.0;
    transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 7] + aph->oindex)));
    aph->transmissionSave[0] = transmission;
    transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 8] + aph->oindex)));
    aph->transmissionSave[1] = transmission;
    transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 9] + aph->oindex)));
    aph->transmissionSave[2] = transmission;
    transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 10] + aph->oindex)));
    aph->transmissionSave[3] = transmission;
    transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 11] + aph->oindex)));
    aph->transmissionSave[4] = transmission;

    return(transmission);
}


double Image::atmosphereOpacity (Vector angle, int layer, Photon *aph, int waveSurfaceIndex, int surfaceIndex) {

    double dvx, dvy, airmassl;
    double rindex;
    double densityScale, densityScale2, densityScale3, densityScale4;
    double densityScale5, densityScale6, densityScale7;

    if (layer == -1) {

        dvx = -angle.x + zenith*sin(aph->shiftedAngle);
        dvy = -angle.y + zenith*cos(aph->shiftedAngle);
        aph->dvr = sqrt(dvx*dvx + dvy*dvy);
        if (aph->dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + air.topAtmosphere/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
                sin(aph->dvr)/(air.topAtmosphere - height[layer + 1]);
        }

        aph->oindex = find_linear(air.tauWavelength, AIR_OPACITY_POINT, aph->wavelength, &rindex);
        return(exp(-(*(air.tau[0] + aph->oindex))*airmassl/air.airmassLayer[layer + 1]));


    } else {

        if (aph->dvr == 0.0) {
            airmassl = 1.0;
        } else {
            airmassl = RADIUS_EARTH*(sqrt(pow(1 + height[layer]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1) -
                                     sqrt(pow(1 + height[layer + 1]/RADIUS_EARTH, 2.0)/sin(aph->dvr)/sin(aph->dvr) - 1))*
                sin(aph->dvr)/(height[layer] - height[layer + 1]);
        }

        densityScale = 1.0 + raygradient*(aph->xpos*cos(rayAngle) + aph->ypos*sin(rayAngle))/(1000000.0);
        densityScale2 = 1.0 + o2gradient*(aph->xpos*cos(o2Angle) + aph->ypos*sin(o2Angle))/(1000000.0);
        densityScale3 = 1.0 + h2ogradient*(aph->xpos*cos(h2oAngle) + aph->ypos*sin(h2oAngle))/(1000000.0);
        densityScale4 = 1.0 + o3gradient*(aph->xpos*cos(o3Angle) + aph->ypos*sin(o3Angle))/(1000000.0);
        densityScale5 = 1.0 + n2ogradient*(aph->xpos*cos(n2oAngle) + aph->ypos*sin(n2oAngle))/(1000000.0);
        densityScale6 = 1.0 + cogradient*(aph->xpos*cos(coAngle) + aph->ypos*sin(coAngle))/(1000000.0);
        densityScale7 = 1.0 + ch4gradient*(aph->xpos*cos(ch4Angle) + aph->ypos*sin(ch4Angle))/(1000000.0);
 
        return(exp(-(*(air.tau[12*(layer + 1) + 0] + aph->oindex)*densityScale +
                     *(air.tau[12*(layer + 1) + 1] + aph->oindex)*densityScale2 +
                     *(air.tau[12*(layer + 1) + 2] + aph->oindex)*densityScale3 +
                     *(air.tau[12*(layer + 1) + 3] + aph->oindex)*densityScale4 +
                     *(air.tau[12*(layer + 1) + 4] + aph->oindex)*densityScale5 +
                     *(air.tau[12*(layer + 1) + 5] + aph->oindex)*densityScale6 +
                     *(air.tau[12*(layer + 1) + 6] + aph->oindex)*densityScale7)*airmassl/air.airmassLayer[layer + 1]));

    }

}

double Image::atmosphereOpacityMoon (Vector angle, int layer, Photon *aph) {

    double rindex;

    if (layer == -1) {

        aph->oindex = find_linear(air.tauWavelength, AIR_OPACITY_POINT, aph->wavelength, &rindex);
        return(exp(-(*(air.tau[0] + aph->oindex))));

    } else {

        double transmission = 1.0;
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 0] + aph->oindex)));
        aph->transmissionSave[0] = transmission;
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 1] + aph->oindex)));
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 2] + aph->oindex)));
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 3] + aph->oindex)));
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 4] + aph->oindex)));
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 5] + aph->oindex)));
        transmission *= exp(-(*(air.tauMoonScatter[12*(layer + 1) + 6] + aph->oindex)));
        aph->transmissionSave[1] = transmission;

        return(transmission);

    }

}

double Image::atmosphereOpacitySun (Vector angle, int layer, Photon *aph) {

    double rindex;

    if (layer == -1) {

        aph->oindex = find_linear(air.tauWavelength, AIR_OPACITY_POINT, aph->wavelength, &rindex);
        return(exp(-(*(air.tau[0] + aph->oindex))));

    } else {

         double transmission = 1.0;
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 0] + aph->oindex)));
        aph->transmissionSave[0] = transmission;
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 1] + aph->oindex)));
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 2] + aph->oindex)));
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 3] + aph->oindex)));
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 4] + aph->oindex)));
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 5] + aph->oindex)));
        transmission *= exp(-(*(air.tauSunScatter[12*(layer + 1) + 6] + aph->oindex)));
        aph->transmissionSave[1] = transmission;

        return(transmission);

    }

}


// transform to optic frame
void Image::transform( Vector *position, Vector *angle, int surfaceIndex, int focusFlag, Photon *aph) {

    double vprime1, vprime2, vprime3;

    // defocus
    position->z  =  position->z - perturbation.defocus[surfaceIndex] - surface.height[surfaceIndex];
    // decenter
    position->x = position->x - perturbation.decenterX[surfaceIndex] - surface.centerx[surfaceIndex];
    position->y = position->y - perturbation.decenterY[surfaceIndex] - surface.centery[surfaceIndex];

    if (focusFlag == 1) {
        int step;
        double shiftedValue;
        step = floor(aph->time/exptime*focusSteps);
        shiftedValue = step - (focusSteps - 1)/2.0;
        position->z += shiftedValue*focusStepZ;
        if (step == 0) shiftedValue -= 1.0;
        position->x += shiftedValue*focusStepX*cos(aph->shiftedAngle);
        position->y += shiftedValue*focusStepX*sin(aph->shiftedAngle);
    }

    // angles
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex+0*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.rotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;

    position->z = position->z + surface.height[surfaceIndex];
    position->x = position->x + surface.centerx[surfaceIndex];
    position->y = position->y + surface.centery[surfaceIndex];

}

// transform back to lab frame
void Image::transformInverse( Vector *position, Vector *angle, int surfaceIndex) {

    double vprime1, vprime2, vprime3;


    position->z = position->z - surface.height[surfaceIndex];
    position->x = position->x - surface.centerx[surfaceIndex];
    position->y = position->y - surface.centery[surfaceIndex];

    // angles
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(angle->z);
    vprime2 = (*(perturbation.inverserotationmatrix+9*surfaceIndex + 1*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 2))*(angle->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(angle->x)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(angle->y)+
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(angle->z);

    angle->x = vprime1;
    angle->y = vprime2;
    angle->z = vprime3;

    // position
    vprime1 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 0*3 + 2))*(position->z);
    vprime2 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 1*3 + 2))*(position->z);
    vprime3 = (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 0))*(position->x) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 1))*(position->y) +
        (*(perturbation.inverserotationmatrix + 9*surfaceIndex + 2*3 + 2))*(position->z);

    position->x = vprime1;
    position->y = vprime2;
    position->z = vprime3;


    // defocus
    position->z = position->z + perturbation.defocus[surfaceIndex] + surface.height[surfaceIndex];

    // decenter
    position->x = position->x + perturbation.decenterX[surfaceIndex] + surface.centerx[surfaceIndex];
    position->y = position->y + perturbation.decenterY[surfaceIndex] + surface.centery[surfaceIndex];

}

void Image::interceptDerivatives(Vector *normal, Vector position, int surfaceIndex, Photon *aph) {

    double rxd, udmin = 0.0, wdmin = 0.0;
    int umin = 0, rx, wmin = 0;
    double normal3, normal2;
    double r;
    double dx, dy, dxp, dyp;

    dx = position.x - surface.centerx[surfaceIndex];
    dy = position.y - surface.centery[surfaceIndex];
    if (surface.surfacetype[surfaceIndex]==MIRROR) {
        dxp = cos(aph->shiftedAngle)*dx + sin(aph->shiftedAngle)*dy;
        dyp = -sin(aph->shiftedAngle)*dx + cos(aph->shiftedAngle)*dy;
    } else {
        dxp=dx;
        dyp=dy;
    }

    r = sqrt(dxp*dxp + dyp*dyp);
    rx = find_linear(&surface.radius[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &rxd);
    if (perturbation.zernikeflag[surfaceIndex] == 1) {
#ifdef  PERTURBATION_COORDINATE
        double phi;
        if (r < RAYTRACE_TOLERANCE) {
            phi = atan2(dyp, dxp);
            if (phi < 0) phi += 2*PI;
            dxp += RAYTRACE_TOLERANCE*cos(phi);
            dyp += RAYTRACE_TOLERANCE*sin(phi);
            r = sqrt(dxp*dxp + dyp*dyp);
            rx = find_linear(&surface.radius[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &rxd);
        }
        wmin = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &wdmin);
        phi = atan2(dyp, dxp);
        if (phi < 0) phi += 2*PI;
        umin = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, phi, &udmin);

#else
        umin = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, dxp/perturbation.rmax[surfaceIndex], &udmin);
        wmin = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, dyp/perturbation.rmax[surfaceIndex], &wdmin);
#endif
    }

    normal3 = interpolate_linear(&surface.normal[SURFACE_POINTS*surfaceIndex], rx, rxd);
    if (perturbation.zernikeflag[surfaceIndex] == 0) {
        normal->x = -normal3*dx/r;
        normal->y = -normal3*dy/r;
        normal->z = 1.0;
    } else {

#ifdef  PERTURBATION_COORDINATE
        normal3 += (interpolate_bilinear(perturbation.zernike_summed_nr_p +
                                         surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin))
            /perturbation.rmax[surfaceIndex];
        normal2 = interpolate_bilinear(perturbation.zernike_summed_np_r +
                                       surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin);
        normal->x = -normal3*dx/r + dy*normal2/(r*r);
        normal->y = -normal3*dy/r - dx*normal2/(r*r);
        normal->z = 1.0;
#else
        double normal1;
        normal1 = interpolate_bilinear(perturbation.zernike_summed_nr_p +
                                       surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin)
            /perturbation.rmax[surfaceIndex];
        normal2 = interpolate_bilinear(perturbation.zernike_summed_np_r +
                                       surfaceIndex*PERTURBATION_POINTS*PERTURBATION_POINTS, PERTURBATION_POINTS, umin, udmin, wmin, wdmin)
            /perturbation.rmax[surfaceIndex];
        if (surface.surfacetype[surfaceIndex]==MIRROR) {
            normal->x = -normal3*dx/r - (normal1*cos(aph->shiftedAngle) - normal2*sin(aph->shiftedAngle));
            normal->y = -normal3*dy/r - (normal1*(sin(aph->shiftedAngle)) + normal2*cos(aph->shiftedAngle));
        } else {
            normal->x = -normal3*dx/r - normal1;
            normal->y = -normal3*dy/r - normal2;
        }
        normal->z = 1.0;
#endif
    }
    normalize(normal);

}

int Image::chooseSurface (int *newSurf, int *oldSurf, Photon *aph) {

    if (aph->direction == 1) {
    trynextsurface:;
        if ((*newSurf) <= (*oldSurf) + 1) {
            (*newSurf)--;
            if ((*newSurf) == (*oldSurf)) (*newSurf) = (*oldSurf) - 1;
            if ((*newSurf) <= -1) (*newSurf) = (*oldSurf) + 2;
        } else {
            (*newSurf)++;
        }
        if ((*newSurf) > nsurf - 1) {
            return(1);
        } else {
            if (ghost[(*newSurf)] == 1) goto trynextsurface;
            return(0);
        }
    } else {
    trynextsurfaceback:;
        if ((*newSurf) >= (*oldSurf) - 1) {
            (*newSurf)++;
            if ((*newSurf) == (*oldSurf)) (*newSurf) = (*oldSurf) + 1;
            if ((*newSurf) > nsurf - 1) (*newSurf) = (*oldSurf) - 2;
        } else {
            (*newSurf)--;
        }
        if ((*newSurf) <= -1) {
            return(1);
        } else {
            if (ghost[(*newSurf)] == 1) goto trynextsurfaceback;
            return(0);
        }
    }

}


void Image::atmosphereDiffraction (Vector *angle, Photon *aph) {

    double radius, cr, cc, tf;
    fftw_plan pb;
    long i, j = 0, ix, jx;
    //    double dx, dy;
    double norm;
    double seeing;
    double rinner = 0;
    double router = SCREEN_SIZE;
    double datamax = 255.0;

    long index;
    if (atmospheremode >=2) index = aph->thread*SCREEN_SIZE*SCREEN_SIZE; else index = 0;


    if (atmospheremode < 2) { //turbulence off or atmosphere completely off
        seeing = 0.0;
        norm = 0.0;
        if (pupilscreenMode == 1) {
            char tempstring[4096];
            snprintf(tempstring, sizeof(tempstring), "%s/pupilscreen.fits",instrdir.c_str());
            datamax = 255.0;
            datamax = screen.readScreen(9, screen.pupil_temp, tempstring);
            screen.paddingfactor = 8;
            screen.paddingfactor = screen.readScreen(10, screen.pupil_temp, tempstring);
            screen.pupilscreenscale = 1.0;
            screen.pupilscreenscale = screen.readScreen(11, screen.pupil_temp, tempstring);
            if (ATM_RESAMPLE != 1) screen.downsample(screen.pupil_values, screen.pupil_temp);
        } else {
            screen.paddingfactor = 32;
            rinner = screen.fine_sizeperpixel*SCREEN_SIZE*surface.innerRadius[0]/(2*surface.outerRadius[0]*screen.paddingfactor);
            router = screen.fine_sizeperpixel*SCREEN_SIZE/(2*screen.paddingfactor);
        }
    } else {
        seeing = (totalseeing + 1e-6)*pow(1/cos(zenith), 0.6)*aph->wavelengthFactor;
        norm = sqrt(0.0229*pow(0.98*(1e-4*0.5)/(seeing*ARCSEC), -5.0/3.0))*3.25e8;
        screen.paddingfactor = 1;
        rinner = surface.innerRadius[0]/screen.paddingfactor;
        router = surface.outerRadius[0]/screen.paddingfactor;
    }

    for (i = 0; i < SCREEN_SIZE; i++) {
        for (j = 0; j < SCREEN_SIZE; j++) {
            radius = sqrt((i - SCREEN_SIZE/2 + 0.5)*(i - SCREEN_SIZE/2 + 0.5)+
                          (j - SCREEN_SIZE/2 + 0.5)*(j - SCREEN_SIZE/2 + 0.5))*screen.fine_sizeperpixel;
             if (atmospheremode < 2 && pupilscreenMode == 1) {
                screen.inscreen[index + SCREEN_SIZE*i + j][0] = screen.pupil_values[index + j*SCREEN_SIZE + i]/datamax;
                screen.inscreen[index + SCREEN_SIZE*i + j][1] = 0.0;
            } else {
                if (radius > rinner && radius < router) {
                    screen.inscreen[index + SCREEN_SIZE*i + j][0] = cos(screen.phasescreen[index + i*SCREEN_SIZE + j]*norm);
                    screen.inscreen[index + SCREEN_SIZE*i + j][1] = sin(screen.phasescreen[index + i*SCREEN_SIZE + j]*norm);
                } else {
                    screen.inscreen[index + SCREEN_SIZE*i + j][0] = 0.0;
                    screen.inscreen[index + SCREEN_SIZE*i + j][1] = 0.0;
                }
            }
        }
    }

    pthread_mutex_lock(&lock.lock1);
    pb = fftw_plan_dft_2d(SCREEN_SIZE, SCREEN_SIZE, screen.inscreen + index, screen.outscreen + index, FFTW_FORWARD, FFTW_ESTIMATE);
    pthread_mutex_unlock(&lock.lock1);
    fftw_execute(pb);
    pthread_mutex_lock(&lock.lock3);
    fftw_destroy_plan(pb);
    pthread_mutex_unlock(&lock.lock3);

    double total = 0.0;
    double value;
    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            ix = i + SCREEN_SIZE/2;
            jx = j + SCREEN_SIZE/2;
            ix = ix % (SCREEN_SIZE);
            jx = jx % (SCREEN_SIZE);
            if (atmospheremode < 2 && pupilscreenMode == 1) {
                value = pow(screen.outscreen[index + i*SCREEN_SIZE + j][0], 2) + pow(screen.outscreen[index + i*SCREEN_SIZE + j][1], 2);
                total += value;
                screen.focalscreen[index + ix*SCREEN_SIZE + jx] = value;
            } else {
                screen.focalscreen[index + ix*SCREEN_SIZE + jx] = pow(screen.outscreen[index + i*SCREEN_SIZE + j][0], 2)+
                    pow(screen.outscreen[index + i*SCREEN_SIZE + j][1], 2);
            }
        }
    }


    double cump = 0.0;
    tf = 0.0;
    for (i = 0; i < SCREEN_SIZE; i++) {
        for (j = 0; j < SCREEN_SIZE; j++) {
            tf += screen.focalscreen[index + i*SCREEN_SIZE + j];
            if (atmospheremode < 2 && pupilscreenMode == 1) { // normalize
                screen.focalscreencum[index + i*SCREEN_SIZE + j] = cump;
                cump += screen.focalscreen[index + i*SCREEN_SIZE + j]/total;
            }
        }
    }
    cr = random[aph->thread].uniform()*tf;
    cc = 0.0;
    for (i = 0; i<SCREEN_SIZE; i++) {
        for (j = 0; j<SCREEN_SIZE; j++) {
            cc += screen.focalscreen[index + i*SCREEN_SIZE + j];
            if (cr < cc) goto breakphoton;
        }
    }

breakphoton:;



    if (atmospheremode > 1) { // else, do properly in secondKick instead
        angle->x = angle->x + (i-SCREEN_SIZE/2 + random[aph->thread].uniform()-0.5)*aph->wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
        angle->y = angle->y + (j-SCREEN_SIZE/2 + random[aph->thread].uniform()-0.5)*aph->wavelength*1e-3/(SCREEN_SIZE*screen.fine_sizeperpixel);
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }

    aph->prtime = aph->time;

}


int Image::bloom (int saturatedFlag, Photon *aph) {

    int newxpos, oldxpos;
    int oldypos;
    int startstep;
    int location;
    int side, origside;
    double rd = random[aph->thread].uniform();

    location = chip.extraX*(aph->pos.y - minyDrift) + aph->pos.x - minxDrift;
    if (*(state.satmap + location) == -1) return(1);

    // move in x
    if ((silicon.xTransferFlag==1) && (silicon.yTransferFlag==0)) {
        oldxpos = aph->pos.x;
        startstep = *(state.satmap + location);
        if (startstep < 1) startstep = 1;
        if (aph->pos.x >= chip.midpoint) {
            origside = 1;
        } else {
            origside = 0;
        }
        for (int iii = startstep ; iii < chip.midpoint + 1; iii++) {
            if (rd > 0.5) {
                newxpos = oldxpos + iii;
            } else {
                newxpos = oldxpos - iii;
            }
            if (newxpos >= chip.midpoint) side=1; else side=0;
            if (((side == origside) || (silicon.midlineFlag==0)) && (newxpos >= minxDrift) && (newxpos <= maxxDrift)) {
                if (*(state.focalPlane + chip.extraX*(aph->pos.y - minyDrift) + newxpos - minxDrift) < well_depth) {
                    aph->pos.x = newxpos;
                    *(state.satmap + location) = iii;
                    return(0);
                }
            }
            if (rd > 0.5) {
                newxpos = oldxpos - iii;
            } else {
                newxpos = oldxpos + iii;
            }
            if (newxpos >= chip.midpoint) side=1; else side=0;
            if (((side == origside) || (silicon.midlineFlag==0)) && (newxpos >= minxDrift) && (newxpos <= maxxDrift)) {
                if (*(state.focalPlane + chip.extraX*(aph->pos.y - minyDrift) + newxpos - minxDrift) < well_depth) {
                    aph->pos.x = newxpos;
                    *(state.satmap + location) = iii;
                    return(0);
                }
            }
        }
        *(state.satmap + location) = -1;
        return(1);
    // move in y
    }
    if ((silicon.xTransferFlag==0) && (silicon.yTransferFlag==1)) {
        oldxpos = aph->pos.y;
        startstep = *(state.satmap + location);
        if (startstep < 1) startstep = 1;
        if (aph->pos.y >= chip.midpoint) {
            origside = 1;
        } else {
            origside = 0;
        }
        for (int iii = startstep ; iii < chip.midpoint + 1; iii++) {
            if (rd > 0.5) {
                newxpos = oldxpos + iii;
            } else {
                newxpos = oldxpos - iii;
            }
            if (newxpos >= chip.midpoint) side=1; else side=0;
            if (((side == origside) || (silicon.midlineFlag==0)) && (newxpos >= minyDrift) && (newxpos <= maxyDrift)) {
                if (*(state.focalPlane + chip.extraX*(newxpos - minyDrift) + aph->pos.x - minxDrift) < well_depth) {
                    aph->pos.y = newxpos;
                    *(state.satmap + location) = iii;
                    return(0);
                }
            }
            if (rd > 0.5) {
                newxpos = oldxpos - iii;
            } else {
                newxpos = oldxpos + iii;
            }
            if (newxpos >= chip.midpoint) side=1; else side=0;
            if (((side == origside) || (silicon.midlineFlag==0)) && (newxpos >= minyDrift) && (newxpos <= maxyDrift)) {
                if (*(state.focalPlane + chip.extraX*(newxpos - minyDrift) + aph->pos.x - minxDrift) < well_depth) {
                    aph->pos.y = newxpos;
                    *(state.satmap + location) = iii;
                    return(0);
                }
            }
        }
        *(state.satmap + location) = -1;
        return(1);
    // move in both
    }
    if ((silicon.xTransferFlag==1) && (silicon.yTransferFlag==1)) {
        oldxpos = aph->pos.x;
        oldypos = aph->pos.y;
        startstep = *(state.satmap + location);
        if (startstep < 1) startstep = 1;
        for (int iii = startstep ; iii < chip.midpoint + 1; iii++) {
            double innerSquare = iii/sqrt(2.0);
            double outerSquare = iii + 1;
            if (rd < 0.25) {
                for (int xx =-outerSquare; xx <= outerSquare; xx++) {
                    for (int yy =-outerSquare; yy<=outerSquare; yy++) {
                        if ((abs(xx) >= innerSquare) || (abs(yy) >= innerSquare)) {
                            float rr = sqrt(xx*xx+yy*yy);
                            if ((rr >= iii) && (rr <= iii+1)) {
                                int xxx = xx + oldxpos;
                                int yyy = yy + oldypos;
                                if ((xxx >= minxDrift) && (xxx <= maxxDrift) &&
                                    (yyy >= minyDrift) && (yyy <= maxyDrift)) {
                                    if ((*(state.focalPlane + chip.extraX*(yyy - minyDrift) + xxx - minxDrift) < well_depth)) {
                                        aph->pos.y = yyy;
                                        aph->pos.x = xxx;
                                        *(state.satmap + location) = iii;
                                        return(0);
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (rd < 0.5) {
                for (int xx =-outerSquare; xx <= outerSquare; xx++) {
                    for (int yy =outerSquare; yy>=-outerSquare; yy--) {
                        if ((abs(xx) >= innerSquare) || (abs(yy) >= innerSquare)) {
                            float rr = sqrt(xx*xx+yy*yy);
                            if ((rr >= iii) && (rr <= iii+1)) {
                                int xxx = xx + oldxpos;
                                int yyy = yy + oldypos;
                                if ((xxx >= minxDrift) && (xxx <= maxxDrift) &&
                                    (yyy >= minyDrift) && (yyy <= maxyDrift)) {
                                    if ((*(state.focalPlane + chip.extraX*(yyy - minyDrift) + xxx - minxDrift) < well_depth)) {
                                        aph->pos.y = yyy;
                                        aph->pos.x = xxx;
                                        *(state.satmap + location) = iii ;
                                        return(0);
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (rd < 0.75) {
                for (int xx =outerSquare; xx >= -outerSquare; xx--) {
                    for (int yy =-outerSquare; yy<=outerSquare; yy++) {
                        if ((abs(xx) >= innerSquare) || (abs(yy) >= innerSquare)) {
                            float rr = sqrt(xx*xx+yy*yy);
                            if ((rr >= iii) && (rr <= iii+1)) {
                                int xxx = xx + oldxpos;
                                int yyy = yy + oldypos;
                                if ((xxx >= minxDrift) && (xxx <= maxxDrift) &&
                                    (yyy >= minyDrift) && (yyy <= maxyDrift)) {
                                    if ((*(state.focalPlane + chip.extraX*(yyy - minyDrift) + xxx - minxDrift) < well_depth)) {
                                        aph->pos.y = yyy;
                                        aph->pos.x = xxx;
                                        *(state.satmap + location) = iii;
                                        return(0);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for (int xx =outerSquare; xx>=-outerSquare; xx--) {
                    for (int yy =outerSquare; yy>=-outerSquare; yy--) {
                        if ((abs(xx) >= innerSquare) || (abs(yy) >= innerSquare)) {
                            float rr = sqrt(xx*xx+yy*yy);
                            if ((rr >= iii) && (rr <= iii+1)) {
                                int xxx = xx + oldxpos;
                                int yyy = yy + oldypos;
                                if ((xxx >= minxDrift) && (xxx <= maxxDrift) &&
                                    (yyy >= minyDrift) && (yyy <= maxyDrift)) {
                                    if ((*(state.focalPlane + chip.extraX*(yyy - minyDrift) + xxx - minxDrift) < well_depth)) {
                                        aph->pos.y = yyy;
                                        aph->pos.x = xxx;
                                        *(state.satmap + location) = iii;
                                        return(0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        *(state.satmap + location) = -1;
        return(1);
    }

    return(1);


}



void Image::saturateIn (int source, Vector* largeAngle, Photon *aph, double shiftX, double shiftY, double extraLargeAngle, Vector *largeAngle2)

{
    int location, origx, origy;
    unsigned int leftover;
    int minrad;
    int saturated = 0;

    leftover = aph->sourceOver_m;
    location = chip.extraX*(aph->pos.y - minyDrift) + aph->pos.x - minxDrift;

    rebloom:;
        if (aph->ghostFlag == 1) {
            if ((leftover > sqrt((*(state.focalPlane + location))+1)) &&
                ((*(state.focalPlane + location) + leftover) < well_depth)) {
                leftover=sqrt((*(state.focalPlane + location))+1);
            }
            *(state.focalPlane + location) += leftover;
        } else {
            *(state.focalPlane + location) += leftover;
        }
        if (*(state.focalPlane + location) > well_depth) {
            if (saturated==0) {
                origx = aph->pos.x;
                origy = aph->pos.y;
                saturated = 1;
            }
            leftover = *(state.focalPlane + location) - well_depth;
            *(state.focalPlane + location) = well_depth;
            if (blooming == 1) {
                if (bloom(1, aph)) goto fullysat;
                location = chip.extraX*(aph->pos.y - minyDrift) + aph->pos.x - minxDrift;
                goto rebloom;
            }
        }

    fullysat:;

        if (saturated) {
            aph->saturationFlag = 1;
            if (aph->ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1  && sources.type[source] >= 9) {
                minrad = 0;
                int xv = 0, yv = 0;
                //int rad = static_cast<int>(sqrt(pow(largeAngle->x + shiftX, 2.0) + pow(largeAngle->y + shiftY, 2.0) + pow(extraLargeAngle, 2.0))/DEGREE*platescale/pixsize);
                               int rad = static_cast<int>(sqrt(pow(largeAngle->x + shiftX + largeAngle2->x, 2.0) + pow(largeAngle->y + shiftY + largeAngle2->y, 2.0))/DEGREE*platescale/pixsize);
                // if (silicon.xTransferFlag != silicon.yTransferFlag) {
                    int xn, yn;
                    xn = origx + rad;
                    if (xn <= maxxDrift) {
                        if  (*(state.focalPlane+chip.extraX*(origy - minyDrift) + xn - minxDrift) >= well_depth) xv++;
                    }
                    xn = origx - rad;
                    if (xn >= minxDrift) {
                        if  (*(state.focalPlane+chip.extraX*(origy - minyDrift) + xn - minxDrift) >= well_depth) xv++;
                    }
                    yn = origy + rad;
                    if (yn <= maxyDrift) {
                        if  (*(state.focalPlane+chip.extraX*(yn - minyDrift) + origx - minxDrift) >= well_depth) yv++;
                    }
                    yn = origy - rad;
                    if (yn >= minyDrift) {
                        if  (*(state.focalPlane+chip.extraX*(yn - minyDrift) + origx - minxDrift) >= well_depth) yv++;
                    }
                // } else {
                //     xv = 2;
                //     yv = 2;
                // }
                 if ((xv>=2) && (yv>=2)) {
                     minrad = rad - satbuffer;
                     aph->sat_xpos += origx;
                        aph->sat_ypos += origy;
                        aph->sat_n += 1;
                        aph->sat_xx += origx*origx;
                        aph->sat_yy += origy*origy;
                        aph->sat_r += rad;
                        aph->sat_sigma = aph->sat_xx /aph->sat_n - pow(aph->sat_xpos/aph->sat_n, 2.0) + aph->sat_yy /aph->sat_n - pow(aph->sat_ypos/aph->sat_n, 2.0) - pow(aph->sat_r/aph->sat_n,2.0);
                        if (aph->sat_sigma < 0) aph->sat_sigma=0.0;
                        aph->sat_sigma = sqrt(aph->sat_sigma);
                        if (aph->sat_n > 10000) {
                            double r= sqrt(pow(origx-aph->sat_xpos/aph->sat_n, 2.0) + pow(origy-aph->sat_ypos/aph->sat_n, 2.0) );
                            if (r > aph->sat_maxr) {
                                aph->sat_maxr = r;
                            }
                        }
                        if (aph->sat_n > 100) {
                          minrad = rad - static_cast<int>(aph->sat_sigma*2.0);
                        }
                        if (minrad == (aph->sourceSaturationRadius + 1)) {
                            aph->sourceSaturationRadius = minrad;
                            //printf("%d %lf %lf %d %lf\n",static_cast<int>(minrad),aph->sat_sigma,aph->sat_n,satbuffer,aph->sat_r/aph->sat_n);
                        }
                 }
            }
        }
}

void Image::saturateOut (int source, Vector* largeAngle, Photon *aph,double shiftX, double shiftY, double extraLargeAngle) {

        int minrad;

        if (aph->ghostFlag == 0 && sources.spatialtype[source] != 4 && sources.spatialtype[source] != 1
            && sources.type[source] >= 9 && aph->saturationFlag == 0) {
            if (aph->pos.x > 1e6 || aph->pos.x < -1e6 || aph->pos.y > 1e6 || aph->pos.y < -1e6) {

            } else {
            minrad = 0;
            double deltaX, deltaY;
            double extra = extraLargeAngle/DEGREE*platescale/pixsize;
            deltaX = (largeAngle->x + shiftX)/DEGREE*platescale/pixsize;
            deltaY = (largeAngle->y + shiftY)/DEGREE*platescale/pixsize;
            if (aph->pos.y - deltaY - extra > maxyDrift && aph->pos.y > maxyDrift) {
                if (-deltaY - satbuffer - extra > minrad) minrad = static_cast<int>(-deltaY - satbuffer - extra);
            }
            if (aph->pos.y - deltaY - extra < minyDrift && aph->pos.y < minyDrift) {
                if (deltaY - satbuffer - extra > minrad) minrad = static_cast<int>(deltaY - satbuffer - extra);
            }
            if (aph->pos.x - deltaX - extra > maxxDrift && aph->pos.x > maxxDrift) {
                if (-deltaX - satbuffer - extra > minrad) minrad = static_cast<int>(-deltaX - satbuffer - extra);
            }
            if (aph->pos.x - deltaX - extra < minxDrift && aph->pos.x < minxDrift) {
                if (deltaX - satbuffer - extra > minrad) minrad = static_cast<int>(deltaX - satbuffer - extra);
            }
            if (minrad < 0) minrad = 0;
            //            printf("%d %d %lf %lf %lf %e\n",minrad,aph->sourceSaturationRadius,extra,deltaX,deltaY,extraLargeAngle);
            if (minrad == (aph->sourceSaturationRadius + 1)) {
                aph->sourceSaturationRadius = minrad;
            }
            }
        }


}


int Image::findSurface (Vector angle, Vector position, double *distance, int surfaceIndex, Photon *aph) {


    double zmin;
    int miss;
    double a, b, c, l;
    double disc, l1, l2;
    double tol;
    double discrepancy;
    double zv = 0.0;
    double radcurv;
    int iter = 0;

    if (surface.radiusCurvature[surfaceIndex] == 0.0) {

        l = (surface.height[surfaceIndex] - position.z)/angle.z;
        miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
        discrepancy = fabs(position.z + angle.z*l - zv);
        while (discrepancy > RAYTRACE_TOLERANCE) {
            l = l + (zv - (position.z + angle.z*l))/angle.z;
            miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
            discrepancy = fabs(position.z + angle.z*l - zv);
            iter++;
            if (iter > RAYTRACE_MAX_ITER) goto iterjump;
        }
        zmin = discrepancy;
        *distance = l;
    } else {
        radcurv = surface.radiusCurvature[surfaceIndex];
        a = angle.x*angle.x + angle.y*angle.y + (1 + surface.conic[surfaceIndex])/2.0*angle.z*angle.z;
        b = (1 + surface.conic[surfaceIndex])*angle.z*(position.z - surface.height[surfaceIndex]) +
            2.0*angle.x*position.x + 2.0*angle.y*position.y - 2.0*radcurv*angle.z;
        c = (1 + surface.conic[surfaceIndex])/2.0*(position.z*position.z - 2.0*position.z*surface.height[surfaceIndex] +
                                                   surface.height[surfaceIndex]*surface.height[surfaceIndex]) +
            position.x*position.x + position.y*position.y - 2.0*radcurv*(position.z - surface.height[surfaceIndex]);
        disc = b*b - 4*a*c;
        if (a == 0.0) {
            l = -c/b;
            if (l <= 0.0 || b == 0.0) {
                tol = RAYTRACE_MAX_LENGTH;
                miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
            } else {
                miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                discrepancy = fabs(position.z + angle.z*l - zv);
                while (discrepancy > RAYTRACE_TOLERANCE) {
                    l = l + (zv - (position.z + angle.z*l))/angle.z;
                    miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                    discrepancy = fabs(position.z + angle.z*l - zv);
                    iter++;
                    if (iter > RAYTRACE_MAX_ITER) goto iterjump;
                }
                zmin = discrepancy;
                *distance = l;
            }
        } else {
            if (disc > 0) {
                l1 = ((-b + sqrt(disc))/2.0/a);
                l2 = ((-b - sqrt(disc))/2.0/a);
                l = 0.0;
                if ((l2 > 0) && (l1 > 0)) {
                    if (l2 < l1) {
                        l = l2;
                    } else {
                        l = l1;
                    }
                }
                if ((l2 > 0) && (l1 < 0)) l = l2;
                if ((l2 < 0) && (l1 > 0)) l = l1;
                if (l == 0.0) {
                    tol = RAYTRACE_MAX_LENGTH;
                    miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
                } else {
                    miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                    discrepancy = fabs(position.z + angle.z*l - zv);
                    while (discrepancy > RAYTRACE_TOLERANCE) {
                        l = l + (zv - (position.z + angle.z*l))/angle.z;
                        miss = getIntercept(position.x + angle.x*l, position.y + angle.y*l, &zv, surfaceIndex, aph);
                        discrepancy = fabs(position.z + angle.z*l - zv);
                        iter++;
                        if (iter > RAYTRACE_MAX_ITER || miss == 1) goto iterjump;
                    }
                    zmin = discrepancy;
                    *distance = l;
                }
            } else {
            iterjump:;
                l = 0.0;
                tol = RAYTRACE_MAX_LENGTH;
                miss = goldenBisectSurface(l - tol, l, l + tol, &zmin, angle, position, distance, surfaceIndex, aph);
            }
        }
    }

    if (zmin > RAYTRACE_ERROR || miss == 1) {
        return(1);
    } else {
        return(0);
    }

}


int Image::goldenBisectSurface(double a, double b, double c, double *z, Vector angle, Vector position,
                               double *distance, int surfaceIndex, Photon *aph) {

    //double f0;
    double f1, f2;
    //double f3;
    double x0, x1, x2, x3;
    double zv=0.0;
    int miss = 0, miss1, miss2;

    x0 = a;
    x3 = c;
    if (fabs(c - b) > fabs(b - a)) {
        x1 = b;
        x2 = (b + HALFSQ5M1*(c - b));
    } else {
        x2 = b;
        x1 = (b - HALFSQ5M1*(b - a));
    }

    miss1 = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex, aph);
    f1 = fabs(position.z + angle.z*x1 - zv);
    miss2 = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex, aph);
    f2 = fabs(position.z + angle.z*x2 - zv);

    while (fabs(x2 - x1) > RAYTRACE_TOLERANCE) {

        if (f2 < f1) {
            x0 = x1;
            x1 = x2;
            x2 = HALF3MSQ5*x1 + HALFSQ5M1*x3;
            miss1 = getIntercept(position.x + angle.x*x2, position.y + angle.y*x2, &zv, surfaceIndex, aph);
            //f0 = f1;
            f1 = f2;
            f2 = fabs(position.z + angle.z*x2 - zv);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = HALF3MSQ5*x2 + HALFSQ5M1*x0;
            miss2 = getIntercept(position.x + angle.x*x1, position.y + angle.y*x1, &zv, surfaceIndex, aph);
            //f3 = f2;
            f2 = f1;
            f1 = fabs(position.z + angle.z*x1 - zv);
        }
    }

    if (f1 < f2) {
        *z = f1;
        *distance = x1;
    } else {
        *z = f2;
        *distance = x2;
    }

    if ((miss1==1) || (miss2==1)) {
        miss=1;
    }
    return(miss);

}

inline int Image::getIntercept(double x, double y, double *z, int surfaceIndex, Photon *aph) {

    double r;
    double uu, ww, tt;
    double dx, dy;
    int ttint;
    int lSurfaceIndex;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    if (surface.surfacetype[surfaceIndex]==MIRROR) {
        double dxp = cos(aph->shiftedAngle)*dx + sin(aph->shiftedAngle)*dy;
        double dyp = -sin(aph->shiftedAngle)*dx + cos(aph->shiftedAngle)*dy;
        dx=dxp;
        dy=dyp;
    }

    r = sqrt(dx*dx + dy*dy);
    lSurfaceIndex = SURFACE_POINTS*surfaceIndex;
    if ((r > surface.radius[lSurfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[lSurfaceIndex + 0])) {
        return(1);
    }

    ttint = find_linear(&surface.radius[lSurfaceIndex], SURFACE_POINTS, r, &tt);
    *z = interpolate_linear(&surface.profile[lSurfaceIndex], ttint, tt);

    find(&surface.radiusArea[lSurfaceIndex], PERTURBATION_POINTS, r, &(aph->vvint));


#ifdef  PERTURBATION_COORDINATE
    double phi;
    aph->wwint = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &ww);
    phi = atan2(dy, dx);
    if (phi < 0) phi += 2*PI;
    aph->uuint = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, phi, &uu);
#else
    aph->uuint = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, dx/perturbation.rmax[surfaceIndex], &uu);
    aph->wwint = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, dy/perturbation.rmax[surfaceIndex], &ww);
    //    printf("%lf %lf %d %d\n",dx,dy,aph->uuint,aph->wwint);
#endif

    if (perturbation.zernikeflag[surfaceIndex] == 1) {
        double dz;
        dz=interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, aph->uuint, uu, aph->wwint, ww);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + aph->uuint*PERTURBATION_POINTS + aph->wwint] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + (aph->uuint+1)*PERTURBATION_POINTS + aph->wwint] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + aph->uuint*PERTURBATION_POINTS + aph->wwint + 1] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + (aph->uuint+1)*PERTURBATION_POINTS + aph->wwint + 1] == 1) return(1);
        *z += dz;
    }
    // jump1:;
    return(0);

}

int Image::getDeltaIntercept(double x, double y, double *z, int surfaceIndex, Photon *aph) {

    double r;
    double uu, ww;
    int uuint, wwint;
    double dx, dy;
    int miss;

    miss = 0;

    dx = x - surface.centerx[surfaceIndex];
    dy = y - surface.centery[surfaceIndex];
    if (surface.surfacetype[surfaceIndex]==MIRROR) {
        double dxp = cos(aph->shiftedAngle)*dx + sin(aph->shiftedAngle)*dy;
        double dyp = -sin(aph->shiftedAngle)*dx + cos(aph->shiftedAngle)*dy;
        dx=dxp;
        dy=dyp;
    }
    r = sqrt(dx*dx + dy*dy);


    find(&surface.radiusArea[SURFACE_POINTS*surfaceIndex], SURFACE_POINTS, r, &(aph->vvint));
    if ((r > surface.radius[SURFACE_POINTS*surfaceIndex + SURFACE_POINTS - 1]) ||
        (r < surface.radius[SURFACE_POINTS*surfaceIndex + 0])) miss = 1;
    *z = 0.0;

    if (perturbation.zernikeflag[surfaceIndex] == 1) {

#ifdef  PERTURBATION_COORDINATE
        double phi;
        wwint = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, r/perturbation.rmax[surfaceIndex], &ww);
        phi = atan2(dy, dx);
        if (phi < 0) phi += 2*PI;
        uuint = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, phi, &uu);
#else
         uuint = find_linear(perturbation.zernike_a_grid, PERTURBATION_POINTS, dx/perturbation.rmax[surfaceIndex], &uu);
         wwint = find_linear(perturbation.zernike_b_grid, PERTURBATION_POINTS, dy/perturbation.rmax[surfaceIndex], &ww);
#endif

        double dz;
        dz = interpolate_bilinear(perturbation.zernike_summed + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex,
                                   PERTURBATION_POINTS, uuint, uu, wwint, ww);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + uuint*PERTURBATION_POINTS + wwint] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + (uuint+1)*PERTURBATION_POINTS + wwint] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + uuint*PERTURBATION_POINTS + wwint + 1] == 1) return(1);
        if (perturbation.zernike_summed_flag[PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex + (uuint+1)*PERTURBATION_POINTS + wwint + 1] == 1) return(1);
        *z += dz;
    }
    return(miss);
}


int Image::getWavelengthTime (Photon *aph, int source) {

    double tempf1;
    int index;
    double dustvalue = 1.0;

    // select wavelength
    if (sources.spatialtype[source] != OPD) {
        tempf1 = random[aph->thread].uniform();
        //         printf("a %d %d\n",sedPtr[sources.sedptr[source]],sedN[sources.sedptr[source]]);
        find_v_n(sedC, sedPtr[sources.sedptr[source]],  sedN[sources.sedptr[source]], tempf1, &index);
         if (sources.type[source] == 0 && domewave != 0.0) {
             aph->wavelength = sedW[index];
         } else {
             //             printf("b %d\n",index);
             aph->wavelength = interpolate_v(sedW, sedC, tempf1, index);
         }
        aph->wavelength = aph->wavelength/1000.0;
    } else {
        aph->wavelength = sources.gamma1[source]/1000.0;
    }

    // dust at source
    if (sources.dusttypez[source] != 0) {
        if (sources.dusttypez[source] == 1) dustvalue = dust.ccm(aph->wavelength, maxwavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (sources.dusttypez[source] == 2) dustvalue = dust.calzetti(aph->wavelength, maxwavelength, sources.dustparz[source][0], sources.dustparz[source][1]);
        if (random[aph->thread].uniform() > dustvalue) return(1);
    }

    // redshift
    aph->wavelength = (aph->wavelength)*(1 + sources.redshift[source]);
    if (aph->wavelength < minwavelength/1000.0 || aph->wavelength >= maxwavelength/1000.0) return(1);

    // dust at z=0
    if (sources.dusttype[source] != 0) {
        if (sources.dusttype[source] == 1) dustvalue = dust.ccm(aph->wavelength, maxwavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (sources.dusttype[source] == 2) dustvalue = dust.calzetti(aph->wavelength, maxwavelength, sources.dustpar[source][0], sources.dustpar[source][1]);
        if (random[aph->thread].uniform() > dustvalue) {return(1);}
    }

    // choose time
    aph->time = (random[aph->thread].uniform())*exptime;

    // choose polarization
    if (random[aph->thread].uniform() > 0.5) aph->polarization=1; else aph->polarization=-1;

    return(0);

}

void Image::getAngle (Vector *angle, double time, int source, Photon *aph) {

    if (rotationTracking) {
        if ((sources.spatialtype[source] == MOVINGPOINT) ||
            (sources.spatialtype[source] == MOVINGSPHERE) ||
            (sources.spatialtype[source] == MOVINGDISK) ||
            (sources.spatialtype[source] == MOVINGREFLECTEDSPHERE) ) {
            angle->x = sources.vx[source] + sources.spatialpar[source][0]*cos(rotatez + sources.spatialpar[source][1]*DEGREE)
                *(time - exptime/2.0 + timeoffset)*ARCSEC;
            angle->y = sources.vy[source] + sources.spatialpar[source][0]*sin(rotatez + sources.spatialpar[source][1]*DEGREE)
                *(time - exptime/2.0 + timeoffset)*ARCSEC;
        } else {
            angle->x = sources.vx[source];
            angle->y = sources.vy[source];
        }
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    } else {
        float rotationIndexFloat = (aph->absoluteTime+vistime/2)/(vistime)*(rotationPoints-1);
        int rotationIndex = static_cast<int>(rotationIndexFloat);
        if (rotationIndex > rotationPoints - 2) rotationIndex = rotationPoints - 2;
        if (rotationIndex < 0) rotationIndex = 0;
        rotationIndexFloat -= rotationIndex;
        double dx= 0.0, dy= 0.0;
        if ((sources.spatialtype[source] == MOVINGPOINT) ||
            (sources.spatialtype[source] == MOVINGSPHERE) ||
            (sources.spatialtype[source] == MOVINGDISK) ||
            (sources.spatialtype[source] == MOVINGREFLECTEDSPHERE) ) {
            dx = sources.spatialpar[source][0]*cos(rotatez + sources.spatialpar[source][1]*DEGREE)
                *(time - exptime/2.0 + timeoffset)*ARCSEC;
            dy = sources.spatialpar[source][0]*cos(rotatez + sources.spatialpar[source][1]*DEGREE)
                *(time - exptime/2.0 + timeoffset)*ARCSEC;
        }
        angle->x = sources.vx[source*rotationPoints + rotationIndex]*(1 - rotationIndexFloat) +
            sources.vx[source*rotationPoints + rotationIndex + 1]*(rotationIndexFloat) + dx;
        angle->y = sources.vy[source*rotationPoints + rotationIndex]*(1 - rotationIndexFloat) +
            sources.vy[source*rotationPoints + rotationIndex + 1]*(rotationIndexFloat) + dy;
        angle->z = smallAnglePupilNormalize(angle->x, angle->y);

    }
}

int Image::getDeltaAngle(Vector *angle, Vector *position, int source, double *shiftX, double *shiftY, int thread, int *initGal, Photon *aph) {

    double dx = 0.0;
    double dy = 0.0;

    if (sources.spatialtype[source] != POINT && sources.spatialtype[source] != MOVINGPOINT
        && sources.spatialtype[source] != OPD) {
        if (sources.spatialtype[source] == SERSIC2D) {
            galaxy.sersic2d(sources.spatialpar[source][0], sources.spatialpar[source][1],
                            sources.spatialpar[source][2]*DEGREE, sources.spatialpar[source][3],
                            &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == GAUSSIAN) {
            dx = random[aph->thread].normal()*ARCSEC*sources.spatialpar[source][0];
            dy = random[aph->thread].normal()*ARCSEC*sources.spatialpar[source][0];
       } else if (sources.spatialtype[source] == PINHOLE) {
            double r = sqrt(random[aph->thread].uniform()*sources.spatialpar[source][3]*sources.spatialpar[source][3]);
            double phi = random[aph->thread].uniform()*2*PI;
            double xt = position->x - r*cos(phi) - sources.spatialpar[source][0];
            double yt = position->y - r*sin(phi) - sources.spatialpar[source][1];
            double zt = position->z - surface.height[0] - sources.spatialpar[source][2];
            double finiteDistanceR = sqrt(xt*xt + yt*yt + zt*zt);
            dx = -xt/finiteDistanceR;
            dy = yt/finiteDistanceR;
        } else if (sources.spatialtype[source] == SERSIC) {
            galaxy.sersic(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5]*DEGREE,
                          sources.spatialpar[source][6],
                          &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == SERSICDISK) {
            galaxy.sersicDisk(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5]*DEGREE,
                              sources.spatialpar[source][6],
                              &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        }else if (sources.spatialtype[source] == DISTORTEDSPHERE) {
            galaxy.distortedSphere(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3],
                          sources.spatialpar[source][4], sources.spatialpar[source][5],
                                   sources.spatialpar[source][6], sources.spatialpar[source][7],
                                   sources.spatialpar[source][8], sources.spatialpar[source][9],
                                   &dx, &dy, thread);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == SERSICCOMPLEX) {
            galaxy.sersicComplex(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5]*DEGREE,
                                 sources.spatialpar[source][6],
                          sources.spatialpar[source][7], sources.spatialpar[source][8],
                          sources.spatialpar[source][9], sources.spatialpar[source][10],
                          sources.spatialpar[source][11]*DEGREE, sources.spatialpar[source][12],
                          sources.spatialpar[source][13], sources.spatialpar[source][14]*DEGREE,
                                 thread, initGal, &dx, &dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        }else if (sources.spatialtype[source] == SERSICDISKCOMPLEX) {
            galaxy.sersicComplex(sources.spatialpar[source][0], sources.spatialpar[source][1],
                          sources.spatialpar[source][2], sources.spatialpar[source][3]*DEGREE,
                          sources.spatialpar[source][4]*DEGREE, sources.spatialpar[source][5]*DEGREE,
                                 sources.spatialpar[source][6],
                                 sources.spatialpar[source][7], sources.spatialpar[source][8],
                          sources.spatialpar[source][9], sources.spatialpar[source][10],
                          sources.spatialpar[source][11]*DEGREE, sources.spatialpar[source][12],
                          sources.spatialpar[source][13], sources.spatialpar[source][14]*DEGREE,
                                 thread, initGal, &dx, &dy);
            dx = -dx*ARCSEC;
            dy = dy*ARCSEC;
        } else if (sources.spatialtype[source] == SPHERE) {
            double phi = random[aph->thread].uniform()*2*PI;
            double theta = acos(random[aph->thread].uniform());
            dx = sources.spatialpar[source][0]*sin(theta)*cos(phi)*ARCSEC;
            dy = sources.spatialpar[source][0]*sin(theta)*sin(phi)*ARCSEC;
        } else if (sources.spatialtype[source] == DISK) {
            double phi = random[aph->thread].uniform()*2*PI;
            double r = sqrt(random[aph->thread].uniform());
            dx = sources.spatialpar[source][0]*r*cos(phi)*ARCSEC;
            dy = sources.spatialpar[source][0]*r*sin(phi)*ARCSEC;
        } else if (sources.spatialtype[source] == REFLECTEDSPHERE) {
            double phi = random[aph->thread].uniform()*2*PI;
            double theta = acos(random[aph->thread].uniform());
            double dxp = sources.spatialpar[source][0]*sin(theta)*cos(phi)*ARCSEC;
            double dyp = sources.spatialpar[source][0]*sin(theta)*sin(phi)*ARCSEC;
            dx = dxp*cos(sources.spatialpar[source][2]) + dyp*sin(sources.spatialpar[source][2]);
            dy = -dxp*sin(sources.spatialpar[source][2]) + dyp*cos(sources.spatialpar[source][2]);
            if (sin(sources.spatialpar[source][1]*DEGREE)*sin(theta)*cos(phi) +
                cos(sources.spatialpar[source][1]*DEGREE)*cos(theta) < 0) return(1);
        } else if (sources.spatialtype[source] == MOVINGSPHERE) {
            double phi = random[aph->thread].uniform()*2*PI;
            double theta = acos(random[aph->thread].uniform());
            dx = sources.spatialpar[source][2]*sin(theta)*cos(phi)*ARCSEC;
            dy = sources.spatialpar[source][2]*sin(theta)*sin(phi)*ARCSEC;
        } else if (sources.spatialtype[source] == MOVINGDISK) {
            double phi = random[aph->thread].uniform()*2*PI;
            double r = sqrt(random[aph->thread].uniform());
            dx = sources.spatialpar[source][2]*r*cos(phi)*ARCSEC;
            dy = sources.spatialpar[source][2]*r*sin(phi)*ARCSEC;
        } else if (sources.spatialtype[source] == MOVINGREFLECTEDSPHERE) {
            double phi = random[aph->thread].uniform()*2*PI;
            double theta = acos(random[aph->thread].uniform());
            double dxp = sources.spatialpar[source][2]*sin(theta)*cos(phi)*ARCSEC;
            double dyp = sources.spatialpar[source][2]*sin(theta)*sin(phi)*ARCSEC;
            dx = dxp*cos(sources.spatialpar[source][4]) + dyp*sin(sources.spatialpar[source][4]);
            dy = -dxp*sin(sources.spatialpar[source][4]) + dyp*cos(sources.spatialpar[source][4]);
            if (sin(sources.spatialpar[source][3]*DEGREE)*sin(theta)*cos(phi) +
                cos(sources.spatialpar[source][3]*DEGREE)*cos(theta) < 0) return(1);
        } else if (sources.spatialtype[source] == IMAGE) {
            double imagepixel, localcumulative;
            long ii, jj, selectedsky;
            selectedsky = static_cast<long>(sources.spatialpar[source][2]);
            imagepixel = (random[aph->thread].uniform()*cumulative[selectedsky]);
            localcumulative = 0;
            for (ii = 0; ii < naxesb[selectedsky][0]; ii++) {
                if (cumulativex[selectedsky][ii] > imagepixel) {
                    ii--;
                    goto jumpa;
                }
            }
            ii--;
        jumpa:;
            localcumulative = cumulativex[selectedsky][ii];
            for (jj = 0; jj < naxesb[selectedsky][1]; jj++) {
                localcumulative += (*(tempptr[selectedsky] + (ii)*naxesb[selectedsky][1] + jj));
                if (localcumulative > imagepixel) {
                    jj--;
                    goto jumpb;
                }
            }
            jj--;
        jumpb:;
            dx = ((jj - (naxesb[selectedsky][1]/2.0))*cos(sources.spatialpar[source][1]*DEGREE)+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*sin(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
            dy = -((jj - (naxesb[selectedsky][1]/2.0))*(-sin(sources.spatialpar[source][1]*DEGREE))+
                (-(ii - (naxesb[selectedsky][0]/2.0)))*cos(sources.spatialpar[source][1]*DEGREE))*
                sources.spatialpar[source][0]*ARCSEC;
        }
        if (sources.gamma1[source] != 0.0 || sources.gamma2[source] != 0.0) {
            double dxp = dx*(1 + sources.gamma1[source] - sources.kappa[source]) - dy*(sources.gamma2[source]);
            double dyp = dy*(1 - sources.gamma1[source] - sources.kappa[source]) - dx*(sources.gamma2[source]);
            dx = dxp;
            dy = dyp;
        }
        *shiftX = - dx*cos(rotatez) - dy*(sin(rotatez));
        *shiftY = - dx*(sin(rotatez)) + dy*cos(rotatez);
        // angle->x += *shiftX;
        // angle->y += *shiftY;
        // angle->z = smallAnglePupilNormalize(angle->x, angle->y);
    }

    return(0);
}



int Image::photonSiliconPropagate(Vector *angle, Vector *position, double lambda, Vector normal, double dh, int waveSurfaceIndex, Photon *aph) {

    double travel, dtravel;
    int yindex=0;
    Vector origAngle;
    double dead;
    double twicevidotvn;

    aph->z0 = position->z;
    if (detectorcollimate == 1) {
        angle->z = angle->z/fabs(angle->z);
        angle->x = 0.0;
        angle->y = 0.0;
    }
    vectorCopy(*angle, &origAngle);

    // silicon refraction
    double nSi = 1.0;
    yindex = find_linear_float_mc(silicon.wavelengthGrid, silicon.numWavelength, aph->wavelength, random[aph->thread].uniform());
    if (sensorrefraction == 1) {
        //        nSi = interpolate_linear(silicon.indexRefraction, yindex, ryindex);
        nSi = silicon.indexRefraction[yindex];
    } else {
        nSi = 1e6;
    }
    refract(angle, normal, 1, nSi);

    // photo-electron conversion
    aph->xindex = find_linear_float_mc(silicon.temperatureGrid, silicon.numTemperature, sensorTempNominal + sensorTempDelta, random[aph->thread].uniform());
    //    double mfp = interpolate_bilinear(silicon.meanFreePath, silicon.numWavelength, aph->xindex, aph->rxindex, yindex, ryindex);
    float mfp = silicon.meanFreePath[aph->xindex*silicon.numWavelength + yindex];
    if (std::isinf(mfp) || std::isnan(mfp)) return(1);
    double randNum;
    double conversion;
    if (photoelectric == 1) {
        conversion = 1.0 - exp(-2*sensorthickness/fabs(angle->z)/mfp);
    } else {
        conversion = 1.0;
        mfp = 1e-6;
    }
    (aph->counter)++;
    //float f1=state.dynamicTransmission[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex];
    //float f2=state.dynamicTransmissionLow[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex];
    if (conversion > state.dynamicTransmission[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex]) {
        //        printf("H:  %e %e %e %e %e\n",aph->wavelength,f1,f2,conversion,mfp);
        if ((!std::isinf(conversion)) && (!std::isnan(conversion))) {
            state.dynamicTransmission[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex] = conversion;
            state.dynamicTransmissionInit[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex] += 1;
        }
    }
    if (conversion < state.dynamicTransmissionLow[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex]) {
        //float f3=state.dynamicTransmissionLow[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex];
        //       printf("L:  %e %e %e %e %e\n",aph->wavelength,f1,f2,conversion,mfp);
        if ((!std::isinf(conversion)) && (!std::isnan(conversion))) {
            state.dynamicTransmissionLow[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex] = conversion;
            state.dynamicTransmissionInit[6*(natmospherefile+1) + 2*nsurf + 2 + waveSurfaceIndex] += 1;
        }
    }
    if (aph->counter <= aph->maxcounter) {
        randNum = aph->saveRand[aph->counter];
    } else {
        randNum = random[aph->thread].uniform();
    }
    if (randNum > conversion) {
        propagate(position, *angle, sensorthickness/fabs(angle->z));
        return(1);
    }
    travel = mfp*(-log(1.0 - randNum*0.9999));
    aph->pos.x += (2.0*random[aph->thread].uniform() - 1.0)*DET_RESAMPLE;
    aph->pos.y += (2.0*random[aph->thread].uniform() - 1.0)*DET_RESAMPLE;
    aph->location = (floor((chip.nampx + 2*activeBuffer-1)/DET_RESAMPLE)+1)*
        floor((aph->pos.y - miny  + activeBuffer)/DET_RESAMPLE) +
        floor((aph->pos.x - minx + activeBuffer)/DET_RESAMPLE);
    if (aph->pos.x < minx - activeBuffer || aph->pos.x > maxx + activeBuffer ||
        aph->pos.y < miny - activeBuffer || aph->pos.y > maxy + activeBuffer) {
        int xL, yL;
        xL = aph->pos.x;
        yL = aph->pos.y;
        if (aph->pos.x <  minx - activeBuffer) xL =  minx - activeBuffer;
        if (aph->pos.x > maxx + activeBuffer) xL = maxx + activeBuffer;
        if (aph->pos.y < miny - activeBuffer) yL = miny - activeBuffer;
        if (aph->pos.y > maxy + activeBuffer) yL = maxy + activeBuffer;
        aph->location = (floor((chip.nampx + 2*activeBuffer-1)/DET_RESAMPLE)+1)*
            floor((yL - miny + activeBuffer)/DET_RESAMPLE) +
            floor((xL - minx + activeBuffer)/DET_RESAMPLE);
    }
    if (deadlayer == 1) {
        dead = silicon.deadLayer[aph->location];
    } else {
        dead = 0.0;
    }
    //    printf("%e %e %e %e %e %d %e %d\n",travel,angle->z,dead/1e6,mfp,conversion,yindex,aph->wavelength,aph->xindex);
    if (fabs(travel*angle->z) < dead/1e6) return(1);
    aph->sensorDirection = angle->z/fabs(angle->z);
    aph->collect_z = aph->z0 + aph->sensorDirection*sensorthickness;
    if (fabs(travel*angle->z) >= sensorthickness) {
        double conversion = fringing(origAngle, normal, aph->wavelength, nSi, static_cast<double>(sensorthickness*1e3) + dh*1000.0*fringeflag, mfp, aph->polarization);
        if (random[aph->thread].uniform() > conversion) return(1);
        aph->ghostFlag = 1;
        while (travel > 0) {
            if (fabs(travel*angle->z) >= sensorthickness) {
                dtravel = fabs(sensorthickness/angle->z);
                propagate(position, *angle, dtravel);
                reflect(angle, normal,&twicevidotvn);
                travel -= dtravel;
            } else {
                propagate(position, *angle, travel);
                travel = 0;
            }
        }
    } else {
        propagate(position, *angle, travel);
    }
    return(0);
}

void Image::pixelShift(double x, double y, double *rhoprime, double *cost, double *sint) {

    double rho = sqrt(x*x + y*y);
    *cost = x/rho;
    *sint = y/rho;
    *rhoprime = sqrt(rho*rho +
                     (*cost)*(*cost)*silicon.spaceChargeSpreadX +
                     (*sint)*(*sint)*silicon.spaceChargeSpreadY);

}


void Image::pixelPosition(Vector *position, Photon *aph) {

    aph->pos.x = static_cast<int>(floor(position->x*inversePixsize + halfpixelsx));
    aph->pos.y = static_cast<int>(floor(position->y*inversePixsize + halfpixelsy));

}

void Image::pixelPositionDrift(Vector *position, Photon *aph) {

    aph->pos.x = static_cast<int>(floor(position->x*inversePixsize + halfpixelsx + driftX*aph->time));
    aph->pos.y = static_cast<int>(floor(position->y*inversePixsize + halfpixelsy  + driftY*aph->time));

}

void Image::pixelFraction(Vector *position, Photon *aph) {

    aph->xPosR = position->x/pixsize - floor(position->x/pixsize) - 0.5;
    aph->yPosR = position->y/pixsize - floor(position->y/pixsize) - 0.5;

}


int Image::electronSiliconPropagate(Vector *angle, Vector *position, Photon *aph, int counter, int nopsf) {

    if (counter == 0) {
        // aph->zindex = find_linear_float_mc(silicon.thicknessGrid, silicon.numThickness, fabsf(position->z - aph->z0), random[aph->thread].uniform());
        float temp;
        aph->zindex = find_linearFloat(silicon.thicknessGrid, silicon.numThickness, std::abs(position->z - aph->z0), &temp);
    } else {
        aph->zindex += siliconSubSteps;
    }

    // impurities
    if (counter == 0) {
        int dopantIndex;
        if (impurityvariation) {
            dopantIndex = find_linear_float_mc(silicon.dopantGrid, silicon.numDopant, silicon.nbulkmap[aph->location]*nbulk, random[aph->thread].uniform());
        } else {
            dopantIndex = find_linear_float_mc(silicon.dopantGrid, silicon.numDopant, nbulk, random[aph->thread].uniform());
        }
        aph->dopantThicknessIndex = dopantIndex*silicon.numThickness;
        aph->dopantTemperatureThicknessIndex = silicon.numTemperature*(aph->dopantThicknessIndex) +silicon.numThickness*aph->xindex;
        aph->dopantRhoThicknessIndex = silicon.numRho*(aph->dopantThicknessIndex);
    }

    // charge diffusion
    if ((chargediffusion) && (nopsf==0)) {
        float sg = silicon.sigma[aph->dopantTemperatureThicknessIndex+aph->zindex];
        double pair;
        position->x += sg*random[aph->thread].normalPair(&pair);
        position->y += sg*pair;
    }

    // lateral fields
    if (fieldanisotropy) {
        //        printf("x %d %d\n",aph->dopantThicknessIndex,aph->zindex);
        int index = aph->dopantThicknessIndex + aph->zindex;
        float fsg = silicon.fsigma[index];
        float gsg = silicon.gsigma[index];
        //        printf("y %f %f  %d\n",fsg,gsg,aph->location);
        position->x += fsg*silicon.sigmaX[aph->location] + gsg*silicon.deltaX[aph->location];
        position->y += fsg*silicon.sigmaY[aph->location] + gsg*silicon.deltaY[aph->location];
    }

    //    printf("ae\n");
    // pixel lithography error
    if (counter == 0) {
        if (pixelerror) {
            position->x += silicon.gammaX[aph->location]*0.0;
            position->y += silicon.gammaY[aph->location]*0.0;
        }
    }

    if ((chargesharing) && (aph->zindex >= (SILICON_STEPS-siliconSubSteps))) {
        aph->pos.x = static_cast<int>(floor(position->x*inversePixsize + halfpixelsx + driftX*aph->time));
        aph->pos.y = static_cast<int>(floor(position->y*inversePixsize + halfpixelsy + driftY*aph->time));
        if (aph->pos.y >= minyDrift && aph->pos.y <= maxyDrift) {
            if (aph->pos.x <= maxxDrift && aph->pos.x >= minxDrift) {
                int location = chip.extraX*(aph->pos.y - minyDrift) + aph->pos.x - minxDrift;
                unsigned int chp = *(state.focalPlane + location);
                if (chp >= 1000) {
                    int uindex;
                    aph->xPosR = position->x/pixsize - floor(position->x/pixsize) - 0.5;
                    aph->yPosR = position->y/pixsize - floor(position->y/pixsize) - 0.5;
                    double rho, rhoprime, cost, sint, chv;
                    pixelShift(aph->xPosR, aph->yPosR, &rhoprime, &cost, &sint);
                    uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-1, random[aph->thread].uniform());
                    chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                    position->x += chp*chv*cost;
                    position->y += chp*chv*sint;

                    if (aph->pos.x + 1 <= maxxDrift && aph->pos.x + 1 >= minxDrift) {
                        chp = *(state.focalPlane + location + 1);
                        pixelShift(aph->xPosR - 1, aph->yPosR, &rhoprime, &cost, &sint);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                        position->x += chp*chv*cost;
                        position->y += chp*chv*sint;
                    }

                    if (aph->pos.x - 1 <= maxxDrift && aph->pos.x - 1 >= minxDrift) {
                        chp = *(state.focalPlane + location - 1);
                        pixelShift(aph->xPosR + 1, aph->yPosR, &rhoprime, &cost, &sint);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                        position->x += chp*chv*cost;
                        position->y += chp*chv*sint;
                    }

                    if (aph->pos.y + 1 >= minyDrift && aph->pos.y + 1 <= maxyDrift) {
                        chp = *(state.focalPlane + location + chip.extraX);
                        pixelShift(aph->xPosR, aph->yPosR - 1, &rhoprime, &cost, &sint);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                        position->x += chp*chv*cost;
                        position->y += chp*chv*sint;
                    }

                    if (aph->pos.y - 1 >= minyDrift && aph->pos.y - 1 <= maxyDrift) {
                        chp = *(state.focalPlane + location - chip.extraX);
                        pixelShift(aph->xPosR, aph->yPosR + 1, &rhoprime, &cost, &sint);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rhoprime*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                        position->x += chp*chv*cost;
                        position->y += chp*chv*sint;
                    }

                    if (aph->pos.y >= minyDrift && aph->pos.y - 1 <= maxyDrift) {
                        chp = silicon.chargeStopCharge;
                        rho = fabs(aph->yPosR + 0.5);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rho*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex  + silicon.numThickness*uindex + aph->zindex];
                        position->y += chp*chv*(aph->yPosR + 0.5)/rho;
                    }

                    if (aph->pos.y + 1 >= minyDrift && aph->pos.y <= maxyDrift) {
                        chp = silicon.chargeStopCharge;
                        rho = fabs(aph->yPosR - 0.5);
                        uindex = find_linear_float_mc(silicon.rho, silicon.numTemperature, rho*pixsize*1e-1, random[aph->thread].uniform());
                        chv = silicon.hsigma[aph->dopantRhoThicknessIndex + silicon.numThickness*uindex + aph->zindex];
                        position->y += chp*chv*(aph->yPosR - 0.5)/rho;
                    }

                }
            }
        }
    }

    //    position->z += sensorthickness/siliconSegments*aph->sensorDirection;

    return(0);

}

double Image::monolayer (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization, double nMid, double nMidimag) {

    double arg = fabs(normal.x*angle.x + normal.y*angle.y + normal.z*angle.z);
    if (arg > 1) arg = 1.0;
    double inputAngle = acos(arg);
    double midAngle = asin(sin(inputAngle)/nMid);
    double outputAngle = asin(sin(inputAngle)/nSi);

    std::complex<double> midIndex (nMid, nMidimag);
    std::complex<double> siIndex (nSi, -wavelength/meanFreePath/1e3/4/PI);

    std::complex<double> rhot1, rhot2;
    rhot1 = ((pow(cos(inputAngle), polarization) - midIndex*pow(cos(midAngle), polarization))/
            (pow(cos(inputAngle), polarization) + midIndex*pow(cos(midAngle), polarization)));
    rhot2 = ((midIndex*pow(cos(midAngle), polarization) - siIndex*pow(cos(outputAngle), polarization))/
             (midIndex*pow(cos(midAngle), polarization) + siIndex*pow(cos(outputAngle), polarization)));

    std::complex<double> delta, phase;
    delta= 2*PI/wavelength*thickness*midIndex*sqrt(1.0 - sin(inputAngle)*sin(inputAngle)/(nMid*nMid));

    std::complex<double> imaginary (0, 1);
    phase = exp(-2.0*imaginary*delta);
    std::complex<double> gamma;
    std::complex<double> one (1, 0);
    gamma = (rhot1 + rhot2*phase)/(one + rhot1*rhot2*phase);
    double reflection = real(gamma*conj(gamma));
    if (std::isnan(reflection)) reflection=1.0;

    return(reflection);
}


double Image::interfaceMonolayer (double inputAngle, double wavelength, double nIn, double nInimag, double nOut, double nOutimag, double thickness, int polarization, double nMid, double nMidimag) {

    double midAngle = asin(nIn*sin(inputAngle)/nMid);
    double outputAngle = asin(nIn*sin(inputAngle)/nOut);

    std::complex<double> inIndex (nIn, nInimag);
    std::complex<double> midIndex (nMid, nMidimag);
    std::complex<double> outIndex (nOut, nOutimag);

    std::complex<double> rhot1, rhot2;
    rhot1 = ((inIndex*pow(cos(inputAngle), polarization) - midIndex*pow(cos(midAngle), polarization))/
            (inIndex*pow(cos(inputAngle), polarization) + midIndex*pow(cos(midAngle), polarization)));
    rhot2 = ((midIndex*pow(cos(midAngle), polarization) - outIndex*pow(cos(outputAngle), polarization))/
             (midIndex*pow(cos(midAngle), polarization) + outIndex*pow(cos(outputAngle), polarization)));

    std::complex<double> delta, phase;
    delta= 2*PI/wavelength*thickness*midIndex*sqrt(1.0 - sin(inputAngle)*sin(inputAngle)/(nMid*nMid));

    std::complex<double> imaginary (0, 1);
    phase = exp(-2.0*imaginary*delta);
    std::complex<double> gamma;
    std::complex<double> one (1, 0);
    gamma = (rhot1 + rhot2*phase)/(one + rhot1*rhot2*phase);
    double reflection = real(gamma*conj(gamma));

    if (std::isnan(reflection)) reflection=1.0;

    return(reflection);
}

double Image::interfaceMonolayerAveragePolarization(double inputAngle, double wavelength, double nIn, double nInimag, double nOut, double nOutimag, double thickness, double nMid, double nMidimag) {

    double reflection1 = interfaceMonolayer(inputAngle, wavelength, nIn, nInimag, nOut, nOutimag, thickness, 1, nMid, nMidimag);
    double reflection2 = interfaceMonolayer(inputAngle, wavelength, nIn, nInimag, nOut, nOutimag, thickness, -1, nMid, nMidimag);
    return(0.5*(reflection1 + reflection2));
}


double Image::interfaceAveragePolarization (double inputAngle, double nIn, double nInImag, double nOut, double nOutImag) {

    double reflection1 = interface(inputAngle, nIn, nInImag, nOut, nOutImag, 1);
    double reflection2 = interface(inputAngle, nIn, nInImag, nOut, nOutImag, -1);
    return(0.5*(reflection1 + reflection2));
}

double Image::interface (double inputAngle, double nIn, double nInImag, double nOut, double nOutImag, int polarization) {

    double outputAngle  = asin(nIn*sin(inputAngle)/nOut);
    if (sin(inputAngle) >= (nIn/nOut)) {
        return(1.0);
    }
    std::complex<double> inIndex (nIn, nInImag);
    std::complex<double> outIndex (nOut, nOutImag);
    std::complex<double> rho;

    rho =  ((inIndex*pow(cos(inputAngle), polarization) - outIndex*pow(cos(outputAngle), polarization))/
            (inIndex*pow(cos(inputAngle), polarization) + outIndex*pow(cos(outputAngle), polarization)));
    double reflection = real(rho*conj(rho));
    if (std::isnan(reflection)) reflection=1.0;
    return(reflection);

}

double Image::fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization) {

    std::complex<double> nSilicon (nSi, -wavelength/meanFreePath/1e3/4/PI);
    std::complex<double> nOx (1.5, 0);
    std::complex<double> one (1, 0);
    std::complex<double> imaginary (0, 1);
    double thicknessOx = 2.0/1e3;
    double thicknessSi = thickness;

    double arg = fabs(normal.x*angle.x + normal.y*angle.y + normal.z*angle.z);
    if (arg > 1) arg = 1.0;
    double airAngle = acos(arg);
    double siliconAngle = asin(sin(airAngle)/real(nSilicon));
    double oxideAngle = asin(sin(airAngle)/real(nOx));
    double airAngleOut = -airAngle;

    std::complex<double> factor1 (sqrt(1.0 - sin(airAngle)*sin(airAngle)/(real(nSilicon)*real(nSilicon))), 0);
    std::complex<double> factor2 (sqrt(1.0 - sin(airAngle)*sin(airAngle)/(real(nOx)*real(nOx))), 0);
    std::complex<double> delta1;
    std::complex<double> delta2;

    delta1 = 2*PI/wavelength*thicknessSi*nSilicon*factor1;
    delta2 = 2*PI/wavelength*thicknessOx*nOx*factor2;

    std::complex<double> rhot1;
    std::complex<double> rhot2;
    std::complex<double> rhot3;
    rhot1 = (pow(cos(airAngle), polarization) - nSilicon*pow(cos(siliconAngle), polarization))/
            (pow(cos(airAngle), polarization) + nSilicon*pow(cos(siliconAngle), polarization));
    rhot2 = (nSilicon*pow(cos(siliconAngle), polarization) - nOx*pow(cos(oxideAngle), polarization))/
            (nSilicon*pow(cos(siliconAngle), polarization) + nOx*pow(cos(oxideAngle), polarization));
    rhot3 = (nOx*pow(cos(oxideAngle), polarization) - pow(cos(airAngleOut), polarization))/
            (nOx*pow(cos(oxideAngle), polarization) + pow(cos(airAngleOut), polarization));

    std::complex<double> phase1;
    std::complex<double> phase2;
    phase1 = exp(-2.0*imaginary*delta1);
    phase2 = exp(-2.0*imaginary*delta2);

    std::complex<double> gamma1;
    std::complex<double> gamma2;
    gamma2 = (rhot2 + rhot3*phase2)/(one + rhot2*rhot3*phase2);
    gamma1 = (rhot1 + gamma2*phase1)/(one + rhot1*gamma2*phase1);

    double convert = real(gamma1*conj(gamma1)) + real(phase1*conj(phase1))*real(phase2*conj(phase2))*(1.0 - real(rhot3*conj(rhot3)));
    return(convert);


}


int Image::contaminationSurfaceCheck(Vector position, Vector *angle, int surfaceIndex, Photon *aph) {

    #ifdef  PERTURBATION_COORDINATE
    int eeint = aph->vvint;
    #else
    int eeint = aph->wwint;
    #endif

    if (random[aph->thread].uniform() > (double)(*(contamination.transmission +
                                      PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                                      eeint*PERTURBATION_POINTS + aph->uuint))) {
        return(1);
    }

    long cc;
    for (int check=0; check<2; check++) {
        if (((check==0) && (*(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                              eeint*PERTURBATION_POINTS + aph->uuint) != -1)) ||
            ((check==1) && (*(contamination.surfacelist2ndmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                              eeint*PERTURBATION_POINTS + aph->uuint) != -1))) {
            if (check==0) {
                cc = *(contamination.surfacelistmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                    eeint*PERTURBATION_POINTS + aph->uuint);
            } else {
                cc = *(contamination.surfacelist2ndmap + PERTURBATION_POINTS*PERTURBATION_POINTS*surfaceIndex +
                    eeint*PERTURBATION_POINTS + aph->uuint);
            }


        if (sqrt(pow(position.x - contamination.surfacelistx[surfaceIndex][cc], 2.0)+
                 pow(position.y - contamination.surfacelisty[surfaceIndex][cc], 2.0)) <
            contamination.surfacelists[surfaceIndex][cc]) {

            double mu, phi;
            int xx = 0, wx = 0;
            double abs = 0.0;
            double sca = 0.0;
            double frac = 0.0;
            int index;
            double xq = 0.0;
            double ss = contamination.surfacelists[surfaceIndex][cc];

            if (contamination.surfacelistt[surfaceIndex][cc] == 0) {
                xq = 2*PI*ss/(aph->wavelength);
                find(mie.muXDust,MIE_IN_POINT, xq, &xx);
                frac = interpolate_fraction(mie.muXDust, xq, xx);
                if (random[aph->thread].uniform() < frac) xx++;
                int cindex = wx*MIE_IN_POINT + xx;
                abs = mie.absorptionDust[cindex];
                sca = mie.scatteringDust[cindex];
                double ext = abs + sca;
                double uu = random[aph->thread].uniform();
                if (uu < abs/ext) {
                    return(1);
                }
                double tempf1 = random[aph->thread].uniform();
                find(mie.cumulativeMuDust + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
                mu =interpolate(mie.muValueDust, mie.cumulativeMuDust + cindex*MIE_OUT_POINT, tempf1, index);
                phi = random[aph->thread].uniform()*2.0*PI;

            } else {

                xq = 2*PI*ss/(aph->wavelength);
                find(mie.muXWater,MIE_IN_POINT, xq, &xx);
                frac = interpolate_fraction(mie.muXWater, xq, xx);
                if (random[aph->thread].uniform() < frac) xx++;
                int cindex = wx*MIE_IN_POINT + xx;
                abs = mie.absorptionWater[cindex];
                sca = mie.scatteringWater[cindex];
                double ext = abs + sca;
                double uu = random[aph->thread].uniform();
                if (uu < abs/ext) {
                    return(1);
                }
                double tempf1 = random[aph->thread].uniform();
                find(mie.cumulativeMuWater + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
                mu =interpolate(mie.muValueWater, mie.cumulativeMuWater + cindex*MIE_OUT_POINT, tempf1, index);
                phi = random[aph->thread].uniform()*2.0*PI;

            }
            
            
            shift_mu(angle, mu, phi);
            if (mu < 0) {
                return(1);
            }
            aph->ghostFlag = 1;
        }
        }
    }
    return(0);


}

int Image::contaminationDetectorCheck (Vector position, Vector *angle, Photon *aph) {

    if (aph->pos.x >= minx && aph->pos.x <= maxx && aph->pos.y >= miny && aph->pos.y <= maxy) {
        aph->location2 =  (floor((chip.nampx-1)/DET_RESAMPLE)+1)*floor((aph->pos.y- miny)/DET_RESAMPLE) + floor((aph->pos.x - minx)/DET_RESAMPLE);
        if (random[aph->thread].uniformFloat()>(*(contamination.chiptransmission +aph->location2))) {
            return(1);
        }

        long cc =0;
        for (int check=0;check<2;check++) {
            if (((check==0) && (*(contamination.chiplistmap + aph->location2) != -1)) ||
                ((check==1) && (*(contamination.chiplist2ndmap + aph->location2) != -1))) {
            if (check==0) {
                cc = *(contamination.chiplistmap + aph->location2);
            } else {
                cc = *(contamination.chiplist2ndmap + aph->location2);
            }
            if (sqrt(pow(position.x - contamination.chiplistx[cc], 2.0) + pow(position.y - contamination.chiplisty[cc], 2.0)) <
                contamination.chiplists[cc]) {

                double mu, phi;
                int xx = 0, wx = 0;
                double abs = 0.0;
                double sca = 0.0;
                double frac = 0.0;
                int index;
                double xq = 0.0;
                double ss = contamination.chiplists[cc];

                xq = 2*PI*ss/(aph->wavelength);
                find(mie.muXDust,MIE_IN_POINT, xq, &xx);
                frac = interpolate_fraction(mie.muXDust, xq, xx);
                if (random[aph->thread].uniform() < frac) xx++;
                int cindex = wx*MIE_IN_POINT + xx;
                abs = mie.absorptionDust[cindex];
                sca = mie.scatteringDust[cindex];
                double ext = abs + sca;
                double uu = random[aph->thread].uniform();
                if (uu < abs/ext) {
                    return(1);
                }
                double tempf1 = random[aph->thread].uniform();
                find(mie.cumulativeMuDust + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
                mu =interpolate(mie.muValueDust, mie.cumulativeMuDust + cindex*MIE_OUT_POINT, tempf1, index);
                phi = random[aph->thread].uniform()*2.0*PI;

                shift_mu(angle, mu, phi);
                if (mu < 0) {
                    return(1);
                }
                aph->ghostFlag = 1;
            }
        }
        }
    }

    return(0);

}



int Image::mieScatter(double tau, Vector *largeAngle, Photon *aph, int actScatter) {

    Vector aa;
    double zz =1.0;
    int index;
    double frac = 0.0;
    int xx = 0, wx = 0;
    int type = 4;

    if (aph->transmissionCheckRand > aph->transmissionSave[0]) {
        type = 0;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[1]) {
        type = 1;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[2]) {
        type = 2;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[3]) {
        type = 3;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[4]) {
        type = 4;
    } else {
        printf("Error in Mie\n");
    }

    // find effs
    double effs = 0.0;
    double sigma =0.0;
    double maxv = 0.0;
    double minv = 0.0;
    double mean = 0.0;
    if (type == 0) {
        find(mie.muWavelengthSeasalt, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthSeasalt, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find(mie.tauArraySeasalt, MIE_TAU_POINT, tau, &index);
        effs = interpolate(mie.effsSeasalt+MIE_TAU_POINT*wx, mie.tauArraySeasalt, tau, index);
        mean = mie.dropMeanSeasalt;
        sigma = mie.dropSigmaSeasalt;
        minv = mie.dropMinSeasalt;
        maxv = mie.dropMaxSeasalt;
    }
    if (type == 1) {
        find(mie.muWavelengthDust, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthDust, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find(mie.tauArrayDust, MIE_TAU_POINT, tau, &index);
        effs = interpolate(mie.effsDust+MIE_TAU_POINT*wx, mie.tauArrayDust, tau, index);
        mean = mie.dropMeanDust;
        sigma = mie.dropSigmaDust;
        minv = mie.dropMinDust;
        maxv = mie.dropMaxDust;
    }
    if (type == 2) {
        find(mie.muWavelengthSmoke, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthSmoke, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find(mie.tauArraySmoke, MIE_TAU_POINT, tau, &index);
        effs = interpolate(mie.effsSmoke+MIE_TAU_POINT*wx, mie.tauArraySmoke, tau, index);
        mean = mie.dropMeanSmoke;
        sigma = mie.dropSigmaSmoke;
        minv = mie.dropMinSmoke;
        maxv = mie.dropMaxSmoke;
    }
    if (type == 3) {
        find(mie.muWavelengthPollution, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthPollution, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find(mie.tauArrayPollution, MIE_TAU_POINT, tau, &index);
        effs = interpolate(mie.effsPollution+MIE_TAU_POINT*wx, mie.tauArrayPollution, tau, index);
        mean = mie.dropMeanPollution;
        sigma = mie.dropSigmaPollution;
        minv = mie.dropMinPollution;
        maxv = mie.dropMaxPollution;
    }
    if (type == 4) {
        find(mie.muWavelengthWater, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthWater, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find(mie.tauArrayWater, MIE_TAU_POINT, tau, &index);
        effs = interpolate(mie.effsWater+MIE_TAU_POINT*wx, mie.tauArrayWater, tau, index);
        mean = mie.dropMeanWater;
        sigma = mie.dropSigmaWater;
        minv = mie.dropMinWater;
        maxv = mie.dropMaxWater;
    }
    //   int count=0;
   double xq = 0.0;
   double abs = 0.0;
   double sca = 0.0;
   double mu = 0.0;
   int count = 0;
reset:;
    // find size
    double ss = exp(log(mean) + random[aph->thread].normal()*sigma);
    if (ss > maxv) ss=maxv;
    if (ss < minv) ss=minv;

    // find cindex
    if (type == 0) {
        xq = 2*PI*ss/(aph->wavelength);
        find(mie.muXSeasalt,MIE_IN_POINT, xq, &xx);
        frac = interpolate_fraction(mie.muXSeasalt, xq, xx);
        if (random[aph->thread].uniform() < frac) xx++;
    }
    if (type == 1) {
        xq = 2*PI*ss/(aph->wavelength);
        find(mie.muXDust,MIE_IN_POINT, xq, &xx);
        frac = interpolate_fraction(mie.muXDust, xq, xx);
        if (random[aph->thread].uniform() < frac) xx++;
    }
    if (type == 2) {
        xq = 2*PI*ss/(aph->wavelength);
        find(mie.muXSmoke,MIE_IN_POINT, xq, &xx);
        frac = interpolate_fraction(mie.muXSmoke, xq, xx);
        if (random[aph->thread].uniform() < frac) xx++;
    }
    if (type == 3) {
        xq = 2*PI*ss/(aph->wavelength);
        find(mie.muXPollution,MIE_IN_POINT, xq, &xx);
        frac = interpolate_fraction(mie.muXPollution, xq, xx);
        if (random[aph->thread].uniform() < frac) xx++;
    }
    if (type == 4 ) {
        xq = 2*PI*ss/(aph->wavelength);
        find(mie.muXWater,MIE_IN_POINT, xq, &xx);
        frac = interpolate_fraction(mie.muXWater, xq, xx);
        if (random[aph->thread].uniform() < frac) xx++;
    }

    int cindex = wx*MIE_IN_POINT + xx;

    if (type == 0) {
        abs = mie.absorptionSeasalt[cindex];
        sca = mie.scatteringSeasalt[cindex];
    }
    if (type == 1) {
        abs = mie.absorptionDust[cindex];
        sca = mie.scatteringDust[cindex];
    }
    if (type == 2) {
        abs = mie.absorptionSmoke[cindex];
        sca = mie.scatteringSmoke[cindex];
    }
    if (type == 3) {
        abs = mie.absorptionPollution[cindex];
        sca = mie.scatteringPollution[cindex];
    }
    if (type == 4) {
        abs = mie.absorptionWater[cindex];
        sca = mie.scatteringWater[cindex];
    }
    double ext = abs + sca;

    zz = 1.0;
    aa.z = -1.0;
    aa.x = 0.0;
    aa.y = 0.0;
    int scattered = 0;
    redo:;
    //printf("%ld %ld %d %lf %lf %lf %lf\n",wx,xx,type,effs,tau,ss,aph->wavelength);
    //          double mfplocal = 1.0/tau*effs*effs/ss/ss*2.0/4.0;
    double mfplocal = 1.0/tau*effs*effs/ss/ss/ext;
    if (mfplocal < 0.01) mfplocal=0.01;
    double step = mfplocal*random[aph->thread].exponential();
    if (step < 0.01) step=0.01;
    if (step > 10.0 || std::isnan(step)) step=10.0;

    zz += aa.z*step;

          if (zz < 0.0) {
              // removed unscattered light because assuming was extincted
              if (scattered==0) {
                  count++;
                  if (count > 100) return(1);
                  goto reset;
              }
              // only scattered (not reflected or absorbed)
              if (actScatter == 1) {
                  largeAngle->x += aa.x;
                  largeAngle->y += aa.y;
              }
              return(0);
          }
          if (zz > 1.0) {
              // reflected
              return(1);
          }
          double uu = random[aph->thread].uniform();
          // if (uu > ext/4.0) {
          //     goto redo;
          // }
          //          if (uu < abs/4.0) {
          if (uu < abs/ext) {
              //absorbed
              return(1);
          }

          // choose phi
          double phi = random[aph->thread].uniform()*2.0*PI;

          // choose mu
          double tempf1 = random[aph->thread].uniform();
          if (type == 0) {
              find(mie.cumulativeMuSeasalt + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
              mu =interpolate(mie.muValueSeasalt, mie.cumulativeMuSeasalt + cindex*MIE_OUT_POINT, tempf1, index);
          }
          if (type == 1) {
              find(mie.cumulativeMuDust + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
              mu =interpolate(mie.muValueDust, mie.cumulativeMuDust + cindex*MIE_OUT_POINT, tempf1, index);
          }
          if (type == 2) {
              find(mie.cumulativeMuSmoke + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
              mu =interpolate(mie.muValueSmoke, mie.cumulativeMuSmoke + cindex*MIE_OUT_POINT, tempf1, index);
          }
          if (type == 3) {
              find(mie.cumulativeMuPollution + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
              mu =interpolate(mie.muValuePollution, mie.cumulativeMuPollution + cindex*MIE_OUT_POINT, tempf1, index);
          }
          if (type == 4) {
              find(mie.cumulativeMuWater + cindex*MIE_OUT_POINT, MIE_OUT_POINT, tempf1, &index);
              mu =interpolate(mie.muValueWater, mie.cumulativeMuWater + cindex*MIE_OUT_POINT, tempf1, index);
          }
          if (mu > 1) mu = 1.0;
          if (mu < -1) mu = -1.0;

          shift_mu(&aa,mu,phi);
          scattered++;

          if (scattered<1000)  goto redo;

          return(1);

}


float Image::tauSuppression(float transmission) {

    if (transmission != 0.0) {
        float tau = -log(transmission);
        float suppression = (0.5 - 0.5*(5.28*tau + 0.62*tau*tau)/(1.0 + 5.28*tau + 0.62*tau*tau));
        return(transmission + (1 - transmission)*(suppression));
    } else {
        return(0.0);
    }

}


int Image::rejectScatterCloud(double bAngle, Photon *aph) {

    int xx = 0, wx = 0;
    double frac, factor = 0.0;

    if (aph->transmissionCheckRand > aph->transmissionSave[0]) {
        find(mie.muWavelengthSeasalt, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthSeasalt, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find_reverse(mie.muValueSeasalt, MIE_OUT_POINT, cos(bAngle), &xx);
        factor = mie.totalProfileSeasalt[MIE_OUT_POINT*wx + xx]/((sin(bAngle)+1e-6)*2/static_cast<double>(MIE_OUT_POINT))/4.0;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[1]) {
        find(mie.muWavelengthDust, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthDust, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find_reverse(mie.muValueDust, MIE_OUT_POINT, cos(bAngle), &xx);
        factor = mie.totalProfileDust[MIE_OUT_POINT*wx + xx]/((sin(bAngle)+1e-6)*2/static_cast<double>(MIE_OUT_POINT))/4.0;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[2]) {
        find(mie.muWavelengthSmoke, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthSmoke, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find_reverse(mie.muValueSmoke, MIE_OUT_POINT, cos(bAngle), &xx);
        factor = mie.totalProfileSmoke[MIE_OUT_POINT*wx + xx]/((sin(bAngle)+1e-6)*2/static_cast<double>(MIE_OUT_POINT))/4.0;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[3]) {
        find(mie.muWavelengthPollution, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthPollution, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find_reverse(mie.muValuePollution, MIE_OUT_POINT, cos(bAngle), &xx);
        factor = mie.totalProfilePollution[MIE_OUT_POINT*wx + xx]/((sin(bAngle)+1e-6)*2/static_cast<double>(MIE_OUT_POINT))/4.0;
    } else if (aph->transmissionCheckRand > aph->transmissionSave[4]) {
        find(mie.muWavelengthWater, MIE_IN_POINT, aph->wavelength, &wx);
        frac = interpolate_fraction(mie.muWavelengthWater, aph->wavelength, wx);
        if (random[aph->thread].uniform() < frac) wx++;
        find_reverse(mie.muValueWater, MIE_OUT_POINT, cos(bAngle), &xx);
        factor = mie.totalProfileWater[MIE_OUT_POINT*wx + xx]/((sin(bAngle)+1e-6)*2/static_cast<double>(MIE_OUT_POINT))/4.0;
    }

    double uu = random[aph->thread].uniform();

    if (uu < factor) {
        return(0);
    } else {
        return(1);
    }
}

int Image::rejectScatterAtmosphere(double bAngle, Photon *aph) {

    double factor = 0.0;

    if (aph->transmissionCheckRand > aph->transmissionSave[0]) {
        factor = 0.75*(1.0 + cos(bAngle)*cos(bAngle))/4.0;
    } else {
        factor = 1.0/4.0;
    }

    double uu = random[aph->thread].uniform();
    if (uu < factor) {
        return(0);
    } else {
        return(1);
    }
}


PHOSIM_NOINLINE void Image::gaussianXY(Vector *newPosition, Vector oldPosition, float scale, Photon *aph) {

    double pair;

    newPosition->x = oldPosition.x + random[aph->thread].normalPair(&pair)*scale;
    newPosition->y = oldPosition.y + pair*scale;
    newPosition->z = oldPosition.z;

}

PHOSIM_NOINLINE int Image::inBounds(TwoVectorInt *position) {

    if ((position->x >= minx) && (position->y >= miny) &&
        (position->x <= maxx) && (position->y <= maxy)) return(1); else return(0);

}

PHOSIM_NOINLINE int Image::inBoundsDrift(TwoVectorInt *position) {

    if ((position->x >= minxDrift) && (position->y >= minyDrift) &&
        (position->x <= maxxDrift) && (position->y <= maxyDrift)) return(1); else return(0);

}

PHOSIM_NOINLINE int Image::inBoundsBuffer(TwoVectorInt *position, int buffer) {

    if ((position->x >= minx-buffer) && (position->y >= miny-buffer) &&
        (position->x <= maxx+buffer) && (position->y <= maxy+buffer)) return(1); else return(0);

}

PHOSIM_NOINLINE void Image::diagnosticCount(long ray, int count) {

     if (ray>0)  if ((ray % 1000)==0) state.diagnosticCounter[count]++;

}

PHOSIM_NOINLINE void Image::backgroundProbCount(int bindex, int backgroundProbOk, int ssource) {

    if ((sources.type[ssource] < 9) && (backgroundProbOk)) {
        state.bpRej[bindex] += 1;
    }

}
