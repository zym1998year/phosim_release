///
/// @package phosim
/// @file raytrace.cpp
/// @brief raytrace functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author J. Garrett Jernigan (UCB)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <math.h>
#include "raytrace.h"
#include "constants.h"
#include "helpers.h"

void vectorCopy (Vector vectorIn, Vector *vectorOut) {

    vectorOut->x = vectorIn.x;
    vectorOut->y = vectorIn.y;
    vectorOut->z = vectorIn.z;

}

void vectorAdd (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x + vectorB.x;
    vectorOut->y = vectorA.y + vectorB.y;
    vectorOut->z = vectorA.z + vectorB.z;

}

void vectorSubtract (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x-vectorB.x;
    vectorOut->y = vectorA.y - vectorB.y;
    vectorOut->z = vectorA.z - vectorB.z;

}

void vectorInit (Vector *vector) {

    vector->x = 0.0;
    vector->y = 0.0;
    vector->z = 0.0;

}

void dcrsph (Vector *vector, double longitude, double latitude) {

    vector->x = cos(longitude)*cos(latitude);
    vector->y = sin(longitude)*cos(latitude);
    vector->z = sin(latitude);

}

void dcarsph (Vector *vector, double *longitude, double *latitude) {

    *latitude = asin(vector->z);
    *longitude = atan2(vector->y, vector->x);
    if (*longitude<0.0) *longitude += 2.0*PI;

}

double vectorDot(Vector *vectorA, Vector *vectorB) {

    return(vectorA->x*vectorB->x + vectorA->y*vectorB->y + vectorA->z*vectorB->z);
}

void vectorCross(Vector *vectorA, Vector *vectorB, Vector *vectorC) {

    vectorC->x = vectorA->y*vectorB->z - vectorA->z*vectorB->y;
    vectorC->y = vectorA->z*vectorB->x - vectorA->x*vectorB->z;
    vectorC->z = vectorA->x*vectorB->y - vectorA->y*vectorB->x;
}


void setup_tangent ( double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz) {

    dcrsph(tpz, ra, dec);
    dcrsph(tpx, ra + (PI/2.0), 0.0);
    vectorCross(tpz, tpx, tpy);
}

void tangent (double ra, double dec, double *x, double *y, Vector *tpx, Vector *tpy, Vector *tpz) {

    double dotz;
    Vector s;

    dcrsph(&s, ra, dec);
    dotz = vectorDot(&s, tpz);
    *x = vectorDot(&s, tpx)/dotz;
    *y = vectorDot(&s, tpy)/dotz;

    // s.x = tpz->x  +   (*x)*tpx->x   +  (*y)*tpy->x;
    // s.y = tpz->y  +   (*x)*tpx->y   +  (*y)*tpy->y;
    // s.z = tpz->z  +   (*x)*tpx->z   +  (*y)*tpy->z;
    // normalize(&S);

}

void shift_mu (Vector *vector, double mu, double phi) {

    double vx = vector->x;
    double vy = vector->y;
    double vz = vector->z;
    if (vx==0.0 && vy==0.0) vy = vy + 1e-9;
    double vr = sqrt(vx*vx + vy*vy);

    double x = sqrt(1 - mu*mu)*cos(phi);
    double y = sqrt(1 - mu*mu)*sin(phi);

    vector->x = (-vy)/vr*x + vz*vx/vr*y + vx*mu;
    vector->y = vx/vr*x + vz*vy/vr*y + vy*mu;
    vector->z = (-vr)*y + vz*mu;

}

void normalize (Vector *vector) {

    double r;

    r = modulus(vector);
    vector->x /= r;
    vector->y /= r;
    vector->z /= r;

}

double rotateInverseX(Vector *vector, double angle) {

    return(vector->x*cos(angle) + vector->y*sin(angle));

}

double rotateInverseY(Vector *vector, double angle) {

    return(vector->y*cos(angle) - vector->x*sin(angle));

}


double modulus (Vector *vector) {

    return(sqrt((vector->x)*(vector->x) + (vector->y)*(vector->y) + (vector->z)*(vector->z)));

}



void propagate(Vector *position, Vector angle, double distance) {

    position->x = position->x + angle.x*distance;
    position->y = position->y + angle.y*distance;
    position->z = position->z + angle.z*distance;

}

void reflect( Vector *vectorIn, Vector normal, double *twicevidotvn) {

    *twicevidotvn = 2*((vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z);

    vectorIn->x -= (*twicevidotvn)*normal.x;
    vectorIn->y -= (*twicevidotvn)*normal.y;
    vectorIn->z -= (*twicevidotvn)*normal.z;

}

void refract (Vector *vectorIn, Vector normal, double n1, double n2) {

    double vidotvn;
    double vidotvn2;
    double twicevidotvn;
    double a, b, b2, d2;

    vidotvn = (vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z;

    if (vidotvn>=0.0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
        vidotvn = -vidotvn;
    }

    vidotvn2 = vidotvn*vidotvn;
    b = n1/n2;
    b2 = b*b;
    d2 = 1.0 - b2*(1.0 - vidotvn2);

    if (d2 >= 0.0) {
        if (vidotvn >= 0.0) {
            a = -b*vidotvn + sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        } else {
            a = -b*vidotvn-sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        }
    } else {
        // Total Internal Reflection
        twicevidotvn = vidotvn + vidotvn;
        vectorIn->x = -twicevidotvn*normal.x + (vectorIn->x);
        vectorIn->y = -twicevidotvn*normal.y + (vectorIn->y);
        vectorIn->z = -twicevidotvn*normal.z + (vectorIn->z);
    }

}

void setupGrids(double *zernikeR, double *zernikePhi, double *zernikeX, double *zernikeY, double *zernikeA, double *zernikeB, long numelements, int coordinate) {

    long index;

    if (coordinate==1) {
        for (long i = 0; i < numelements; i++) {
            for (long j = 0; j < numelements; j++) {
                index= i*numelements + j;
                *(zernikeX + index) = -1.0 + (static_cast<double>(i))/(static_cast<double>(numelements) - 1)*2.0;
                *(zernikeY + index) = -1.0 + (static_cast<double>(j))/(static_cast<double>(numelements) - 1)*2.0;
                *(zernikeR + index) = sqrt(zernikeX[index]*zernikeX[index] + zernikeY[index]*zernikeY[index]);
                double phi=0.0;
                phi = atan2(zernikeY[index], zernikeX[index]);
                if (phi < 0) phi=phi+2*PI;
                *(zernikePhi + index) = phi;
            }
        }
        for (long i = 0; i < numelements; i++) {
            zernikeA[i]=-1.0 + (static_cast<double>(i))/(static_cast<double>(numelements) - 1)*2.0;
            zernikeB[i]=-1.0 + (static_cast<double>(i))/(static_cast<double>(numelements) - 1)*2.0;
        }
    } else {
        for (long i = 0; i < numelements; i++) {
            for (long j = 0; j < numelements; j++) {
                index= i*numelements + j;
                *(zernikeR + index) = (static_cast<double>(j))/(static_cast<double>(numelements) - 1);
                *(zernikePhi + index) = (static_cast<double>(i))/(static_cast<double>(numelements))*2*PI;
                *(zernikeX + index) = zernikeR[index]*cos(zernikePhi[index]);
                *(zernikeY + index) = zernikeR[index]*sin(zernikePhi[index]);
            }
        }
        for (long i = 0; i < numelements; i++) {
            zernikeA[i]=(static_cast<double>(i))/(static_cast<double>(numelements) - 1);
            zernikeB[i]=(static_cast<double>(i))/(static_cast<double>(numelements))*2*PI;
        }
    }
}

void zernikes(double *zernikeR, double *zernikePhi, double *zernikeRGrid, double *zernikePhiGrid,
              double *zernikeNormalR, double *zernikeNormalPhi, long numelements, long nzernikes, int coordinate) {

    double r, phi;

    for (long i = 0; i < numelements; i++) {
        for (long j = 0; j < numelements; j++) {

            long index2 = numelements*numelements;
            long index = i*numelements+j;
            r= *(zernikeRGrid +index);
            phi=*(zernikePhiGrid+index);

        // r = (static_cast<double>(i))/(static_cast<double>(numelements) - 1);
        // phi = (static_cast<double>(i))/(static_cast<double>(numelements))*2*PI;
        // *(zernikeRGrid + i) = r;
        // *(zernikePhiGrid + i) = phi;

        // defocus 2/0
        *(zernikeR + 0*index2 + index) = 1.0;
        *(zernikePhi + 0*index2 + index) = 1.0;
        *(zernikeNormalR + 0*index2 + index) = 0.0;
        *(zernikeNormalPhi + 0*index2 + index) = 0.0;

        // defocus 2/0
        *(zernikeR + 1*index2 + index) = 2*r;
        *(zernikePhi + 1*index2 + index) = cos(phi);
        *(zernikeNormalR + 1*index2 + index) = 2.0;
        *(zernikeNormalPhi + 1*index2 + index) = -sin(phi);

        // defocus 2/0
        *(zernikeR + 2*index2 + index) = 2*r;
        *(zernikePhi + 2*index2 + index) = sin(phi);
        *(zernikeNormalR + 2*index2 + index) = 2.0;
        *(zernikeNormalPhi + 2*index2 + index) = cos(phi);

        // defocus 2/0
        *(zernikeR + 3*index2 + index) = sqrt(3.0)*(2*r*r - 1);
        *(zernikePhi + 3*index2 + index) = 1.0;
        *(zernikeNormalR + 3*index2 + index) = sqrt(3.0)*(4*r);
        *(zernikeNormalPhi + 3*index2 + index) = 0.0;

        // astigmatism 2/2
        *(zernikeR + 4*index2 + index) = sqrt(6.0)*r*r;
        *(zernikePhi + 4*index2 + index) = sin(2*phi);
        *(zernikeNormalR + 4*index2 + index) = sqrt(6.0)*2*r;
        *(zernikeNormalPhi + 4*index2 + index) = cos(2*phi)*2;

        *(zernikeR + 5*index2 + index) = sqrt(6.0)*r*r;
        *(zernikePhi + 5*index2 + index) = cos(2*phi);
        *(zernikeNormalR + 5*index2 + index) = sqrt(6.0)*2*r;
        *(zernikeNormalPhi + 5*index2 + index) = -sin(2*phi)*2;

        // coma 3/1
        *(zernikeR + 6*index2 + index) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernikePhi + 6*index2 + index) = sin(phi);
        *(zernikeNormalR + 6*index2 + index) = sqrt(8.0)*(9*r*r - 2);
        *(zernikeNormalPhi + 6*index2 + index) = cos(phi);

        *(zernikeR + 7*index2 + index) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernikePhi + 7*index2 + index) = cos(phi);
        *(zernikeNormalR + 7*index2 + index) = sqrt(8.0)*(9*r*r - 2);
        *(zernikeNormalPhi + 7*index2 + index) = -sin(phi);

        // trefoil 3/3
        *(zernikeR + 8*index2 + index) = sqrt(8.0)*(r*r*r);
        *(zernikePhi + 8*index2 + index) = sin(3*phi);
        *(zernikeNormalR + 8*index2 + index) = sqrt(8.0)*(3*r*r);
        *(zernikeNormalPhi + 8*index2 + index) = cos(3*phi)*3;

        *(zernikeR + 9*index2 + index) = sqrt(8.0)*(r*r*r);
        *(zernikePhi + 9*index2 + index) = cos(3*phi);
        *(zernikeNormalR + 9*index2 + index) = sqrt(8.0)*(3*r*r);
        *(zernikeNormalPhi + 9*index2 + index) = -sin(3*phi)*3;

        // spherical aberration 4/0
        *(zernikeR + 10*index2 + index) = sqrt(5.0)*(6*r*r*r*r - 6*r*r + 1);
        *(zernikePhi + 10*index2 + index) = 1.0;
        *(zernikeNormalR + 10*index2 + index) = sqrt(5.0)*(24*r*r*r - 12*r);
        *(zernikeNormalPhi + 10*index2 + index) = 0.0;

        // 2nd astigmatism 4/2
        *(zernikeR + 11*index2 + index) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernikePhi + 11*index2 + index) = cos(2*phi);
        *(zernikeNormalR + 11*index2 + index) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernikeNormalPhi + 11*index2 + index) = -sin(2*phi)*2;

        *(zernikeR + 12*index2 + index) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernikePhi + 12*index2 + index) = sin(2*phi);
        *(zernikeNormalR + 12*index2 + index) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernikeNormalPhi + 12*index2 + index) = cos(2*phi)*2;

        // quadfoil 4/4
        *(zernikeR + 13*index2 + index) = sqrt(10.0)*(r*r*r*r);
        *(zernikePhi + 13*index2 + index) = cos(4*phi);
        *(zernikeNormalR + 13*index2 + index) = sqrt(10.0)*(4*r*r*r);
        *(zernikeNormalPhi + 13*index2 + index) = -sin(4*phi)*4;

        *(zernikeR + 14*index2 + index) = sqrt(10.0)*(r*r*r*r);
        *(zernikePhi + 14*index2 + index) = sin(4*phi);
        *(zernikeNormalR + 14*index2 + index) = sqrt(10.0)*(4*r*r*r);
        *(zernikeNormalPhi + 14*index2 + index) = cos(4*phi)*4;

        // 2nd coma 5/1
        *(zernikeR + 15*index2 + index) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernikePhi + 15*index2 + index) = cos(phi);
        *(zernikeNormalR + 15*index2 + index) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernikeNormalPhi + 15*index2 + index) = -sin(phi);

        *(zernikeR + 16*index2 + index) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernikePhi + 16*index2 + index) = sin(phi);
        *(zernikeNormalR + 16*index2 + index) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernikeNormalPhi + 16*index2 + index) = cos(phi);

        // 2nd trefoil 5/3
        *(zernikeR + 17*index2 + index) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernikePhi + 17*index2 + index) = cos(3*phi);
        *(zernikeNormalR + 17*index2 + index) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernikeNormalPhi + 17*index2 + index) = -sin(3*phi)*3;

        *(zernikeR + 18*index2 + index) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernikePhi + 18*index2 + index) = sin(3*phi);
        *(zernikeNormalR + 18*index2 + index) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernikeNormalPhi + 18*index2 + index) = cos(3*phi)*3;

        // pentafoil 5/5
        *(zernikeR + 19*index2 + index) = sqrt(12.0)*(r*r*r*r*r);
        *(zernikePhi + 19*index2 + index) = cos(5*phi);
        *(zernikeNormalR + 19*index2 + index) = sqrt(12.0)*(5*r*r*r*r);
        *(zernikeNormalPhi + 19*index2 + index) = -sin(5*phi)*5;

        *(zernikeR + 20*index2 + index) = sqrt(12.0)*(r*r*r*r*r);
        *(zernikePhi + 20*index2 + index) = sin(5*phi);
        *(zernikeNormalR + 20*index2 + index) = sqrt(12.0)*(5*r*r*r*r);
        *(zernikeNormalPhi + 20*index2 + index) = cos(5*phi)*5;

        *(zernikeR + 21*index2 + index) = sqrt(7.0)*(20.0*r*r*r*r*r*r - 30.0*r*r*r*r + 12.0*r*r - 1.0);
        *(zernikePhi + 21*index2 + index) = 1.0;
        *(zernikeNormalR + 21*index2 + index) = sqrt(7.0)*(120.0*r*r*r*r*r - 120.0*r*r*r + 24.0*r);
        *(zernikeNormalPhi + 21*index2 + index) = 0.0;

        *(zernikeR + 22*index2 + index) = sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        *(zernikePhi + 22*index2 + index) = sin(2*phi);
        *(zernikeNormalR + 22*index2 + index) = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        *(zernikeNormalPhi + 22*index2 + index) = cos(2*phi)*2;

        *(zernikeR + 23*index2 + index) = sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        *(zernikePhi + 23*index2 + index) = cos(2*phi);
        *(zernikeNormalR + 23*index2 + index) = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        *(zernikeNormalPhi + 23*index2 + index) = -sin(2*phi)*2;

        *(zernikeR + 24*index2 + index) = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        *(zernikePhi + 24*index2 + index) = sin(4*phi);
        *(zernikeNormalR + 24*index2 + index) = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        *(zernikeNormalPhi + 24*index2 + index) = cos(4*phi)*4;

        *(zernikeR + 25*index2 + index) = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        *(zernikePhi + 25*index2 + index) = cos(4*phi);
        *(zernikeNormalR + 25*index2 + index) = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        *(zernikeNormalPhi + 25*index2 + index) = -sin(4*phi)*4;

        *(zernikeR + 26*index2 + index) = sqrt(14.0)*(r*r*r*r*r*r);
        *(zernikePhi + 26*index2 + index) = sin(6*phi);
        *(zernikeNormalR + 26*index2 + index) = sqrt(14.0)*(6.0*r*r*r*r*r);
        *(zernikeNormalPhi + 26*index2 + index) = cos(6*phi)*6;

        *(zernikeR + 27*index2 + index) = sqrt(14.0)*(r*r*r*r*r*r);
        *(zernikePhi + 27*index2 + index) = cos(6*phi);
        *(zernikeNormalR + 27*index2 + index) = sqrt(14.0)*(6.0*r*r*r*r*r);
        *(zernikeNormalPhi + 27*index2 + index) = -sin(6*phi)*6;

        if (coordinate!=2) {
            double fr, fphi;
            for (int k=0; k<28; k++) {
                fr= *(zernikeNormalR + k*index2 + index);
                fphi= *(zernikeNormalPhi + k*index2 + index);
                double x = r*cos(phi);
                double y = r*sin(phi);
                *(zernikeNormalR + k*index2 + index) = fr * x/r - fphi *1/y;
                *(zernikeNormalPhi + k*index2 + index) = fr * y/r + fphi *1/x;
            }
        }

        }
    }


}


void zernikeValue(int n, double x, double y, double *t) {

    double r= sqrt(x*x+y*y);
    double phi=atan2(y,x);

    if (n==0) {
        t[0] = 1.0;
        t[1] = 1.0;
        t[2] = 0.0;
        t[3] = 0.0;
    } else if (n==1) {
        t[0] = 2*r;
        t[1] = cos(phi);
        t[2] = 2.0;
        t[3] = -sin(phi);
    } else if (n==2) {
        t[0] = 2*r;
        t[1] = sin(phi);
        t[2] = 2.0;
        t[3] = cos(phi);
    } else if (n==3) {
        t[0] = sqrt(3.0)*(2*r*r - 1);
        t[1] = 1.0;
        t[2] = sqrt(3.0)*(4*r);
        t[3] = 0.0;
    } else if (n==4) {
        t[0] = sqrt(6.0)*r*r;
        t[1] = sin(2*phi);
        t[2] = sqrt(6.0)*2*r;
        t[3] = cos(2*phi)*2;
    } else if (n==5) {
        t[0] = sqrt(6.0)*r*r;
        t[1] = cos(2*phi);
        t[2] = sqrt(6.0)*2*r;
        t[3] = -sin(2*phi)*2;
    } else if (n==6) {
        t[0] = sqrt(8.0)*(3*r*r*r - 2*r);
        t[1] = sin(phi);
        t[2] = sqrt(8.0)*(9*r*r - 2);
        t[3] = cos(phi);
    } else if (n==7) {
        t[0] = sqrt(8.0)*(3*r*r*r - 2*r);
        t[1] = cos(phi);
        t[2] = sqrt(8.0)*(9*r*r - 2);
        t[3] = -sin(phi);
    } else if (n==8) {
        t[0] = sqrt(8.0)*(r*r*r);
        t[1] = sin(3*phi);
        t[2] = sqrt(8.0)*(3*r*r);
        t[3] = cos(3*phi)*3;
    } else if (n==9) {
        t[0] = sqrt(8.0)*(r*r*r);
        t[1] = cos(3*phi);
        t[2] = sqrt(8.0)*(3*r*r);
        t[3] = -sin(3*phi)*3;
    } else if (n==10) {
        t[0] = sqrt(5.0)*(6*r*r*r*r - 6*r*r + 1);
        t[1] = 1.0;
        t[2] = sqrt(5.0)*(24*r*r*r - 12*r);
        t[3] = 0.0;
    } else if (n==11) {
        t[0] = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        t[1] = cos(2*phi);
        t[2] = sqrt(10.0)*(16*r*r*r - 6*r);
        t[3] = -sin(2*phi)*2;
    } else if (n==12) {
        t[0]= sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        t[1] = sin(2*phi);
        t[2]= sqrt(10.0)*(16*r*r*r - 6*r);
        t[3]= cos(2*phi)*2;
    } else if (n==13) {
        t[0] = sqrt(10.0)*(r*r*r*r);
        t[1] = cos(4*phi);
        t[2] = sqrt(10.0)*(4*r*r*r);
        t[3] = -sin(4*phi)*4;
    } else if (n==14) {
        t[0] = sqrt(10.0)*(r*r*r*r);
        t[1] = sin(4*phi);
        t[2] = sqrt(10.0)*(4*r*r*r);
        t[3] = cos(4*phi)*4;
    } else if (n==15) {
        t[0] = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        t[1] = cos(phi);
        t[2] = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        t[3] = -sin(phi);
    } else if (n==16) {
        t[0] = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        t[1] = sin(phi);
        t[2] = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        t[3] = cos(phi);
    } else if (n==17) {
        t[0] = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        t[1] = cos(3*phi);
        t[2] = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        t[3] = -sin(3*phi)*3;
    } else if (n==18) {
        t[0] = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        t[1] = sin(3*phi);
        t[2] = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        t[3] = cos(3*phi)*3;
    } else if (n==19) {
        t[0] = sqrt(12.0)*(r*r*r*r*r);
        t[1] = cos(5*phi);
        t[2] = sqrt(12.0)*(5*r*r*r*r);
        t[3] = -sin(5*phi)*5;
    } else if (n==20) {
        t[0] = sqrt(12.0)*(r*r*r*r*r);
        t[1] = sin(5*phi);
        t[2] = sqrt(12.0)*(5*r*r*r*r);
        t[3] = cos(5*phi)*5;
    } else if (n==21) {
        t[0] = sqrt(7.0)*(20.0*r*r*r*r*r*r - 30.0*r*r*r*r + 12.0*r*r - 1.0);
        t[1] = 1.0;
        t[2] = sqrt(7.0)*(120.0*r*r*r*r*r - 120.0*r*r*r + 24.0*r);
        t[2]= 0.0;
    } else if (n==22) {
        t[0]= sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        t[1] = sin(2*phi);
        t[2] = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        t[3] = cos(2*phi)*2;
    } else if (n==23) {
        t[0] = sqrt(14.0)*(15.0*r*r*r*r*r*r - 20.0*r*r*r*r + 6*r*r);
        t[1] = cos(2*phi);
        t[2] = sqrt(14.0)*(90.0*r*r*r*r*r - 80.0*r*r*r + 12.0*r);
        t[3]= -sin(2*phi)*2;
    } else if (n==24) {
        t[0] = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        t[1] = sin(4*phi);
        t[2] = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        t[3] = cos(4*phi)*4;
    } else if (n==25) {
        t[0] = sqrt(14.0)*(6.0*r*r*r*r*r*r - 5.0*r*r*r*r);
        t[1] = cos(4*phi);
        t[2] = sqrt(14.0)*(36.0*r*r*r*r*r - 20.0*r*r*r);
        t[3]= -sin(4*phi)*4;
    } else if (n==26) {
        t[0] = sqrt(14.0)*(r*r*r*r*r*r);
        t[1] = sin(6*phi);
        t[2] = sqrt(14.0)*(6.0*r*r*r*r*r);
        t[3] = cos(6*phi)*6;
    } else if (n==27) {
        t[0] = sqrt(14.0)*(r*r*r*r*r*r);
        t[1] = cos(6*phi);
        t[2] = sqrt(14.0)*(6.0*r*r*r*r*r);
        t[3] = -sin(6*phi)*6;
    }

}

//Chebyshev polynomial of the first kind
double chebyshevT (int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    return 2*x*chebyshevT(n - 1, x) - chebyshevT(n - 2, x);
}

//differentiation of Chebyshev polynomial of the first kind
double chebyshevT_x (int n, double x) {
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;
    return 2*chebyshevT(n - 1, x) + 2*chebyshevT_x(n - 1, x) - chebyshevT_x(n - 2, x);
}

void chebyshevT2D (int n, double x, double y, double *t) {
    int degree, d(n + 1);
    for (degree = 0; degree <= n; degree++) {
        d -= (degree + 1);
        if (d <= 0) {
            d += degree;
            break;
        }
    }
    double tx = chebyshevT(degree - d, x);
    double ty = chebyshevT(d, y);
    double tx_x = chebyshevT_x(degree - d, x);
    double ty_y = chebyshevT_x(d, y);
    t[0] = tx*ty;
    t[1] = tx_x*ty;
    t[2] = tx*ty_y;
}

void chebyshevs (double *r_grid, double *phi_grid, double *chebyshev, double *chebyshev_r, double *chebyshev_phi, long nPoint, long nTerm, int coordinate) {
    for (long j = 0; j < nPoint; j++) {
        for (long l = 0; l < nPoint; l++) {
            double r = r_grid[nPoint*l+j];
            double x = r*cos(phi_grid[nPoint*l+j]);
            double y = r*sin(phi_grid[nPoint*l+j]);
            for (long i = 0; i < nTerm; i++) {
                double t[3], t_r, t_phi;
                chebyshevT2D(i, x, y, &t[0]);
                t_r = t[1]*x/r + t[2]*y/r;         // dz/dr = dz/dx cos(theta) + dz/dy sin(theta)
                t_phi = r*(-t[1]*y/r + t[2]*x/r);  // dz/dtheta = r [ -dz/dx sin(theta) + dz/dy cos(theta) ]
                chebyshev[i*nPoint*nPoint + l*nPoint + j] = t[0];
                if (coordinate==2) {
                    chebyshev_r[i*nPoint*nPoint + l*nPoint + j] = t_r;
                    chebyshev_phi[i*nPoint*nPoint + l*nPoint + j] = t_phi;
                } else {
                    chebyshev_r[i*nPoint*nPoint + l*nPoint + j] = t[1];
                    chebyshev_phi[i*nPoint*nPoint + l*nPoint + j] = t[2];
                }
            }
        }
    }
}

void polys (double *r_grid, double *phi_grid, double *poly, double *poly_r, double *poly_phi, long nPoint, long nTerm, int coordinate) {
    for (long j = 0; j < nPoint; j++) {
        for (long l = 0; l < nPoint; l++) {
            double r = r_grid[nPoint*l+j];
            double x = r*cos(phi_grid[nPoint*l+j]);
            double y = r*sin(phi_grid[nPoint*l+j]);
            long nn = static_cast<long>(sqrt(nTerm));
            long i = 0;
            for (long nx = 0; nx < nn; nx++) {
                for (long ny = 0; ny < nn; ny++) {
                    poly[i*nPoint*nPoint + l*nPoint + j] = pow(x, nx)*pow(y, ny);
                    if (coordinate==2) { 
                        poly_r[i*nPoint*nPoint + l*nPoint + j] = nx*pow(x, nx - 1)*pow(y, ny)*x/r + ny*pow(x, nx)*pow(y, ny-1)*y/r;
                        poly_phi[i*nPoint*nPoint + l*nPoint + j] = (-nx*pow(x, nx - 1)*pow(y, ny)*y + ny*pow(x, nx)*pow(y, ny-1)*x);
                    } else {
                        poly_r[i*nPoint*nPoint + l*nPoint + j] = nx*pow(x, nx - 1)*pow(y, ny);
                        poly_phi[i*nPoint*nPoint + l*nPoint + j] = ny*pow(x, nx)*pow(y, ny-1);
                    }
                    i++;
                }
            }
        }
    }
}


// df/dr = df/dx dx/dr      dr/dx = d(sqrt(x^2+y^2)/dx = 2x*0.5*(x^2+y^2)^(-0.5) = x/r

// df/dphi = df/dphi dphi/dx       x=r*cos(phi)    dx/dphi = r*d(cos(phi)/dphi = r*sin(phi) = y
