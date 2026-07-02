///
/// @package phosim
/// @file surface.h
/// @brief header file for surface class
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

#ifndef Surface_H
#define Surface_H

#include <math.h>
#include "grating.h"

enum SurfaceTypes {MIRROR=0,LENS=1,FILTER=2,DETECTOR=3,GRATING=4,EXITPUPIL=5};

class Surface {

 public:
    double* radiusCurvature;
    double* conic;
    double* height;
    double* intz;
    double* outerRadius;
    double* innerRadius;
    double* innerRadius0;
    double* two;
    double* three;
    double* four;
    double* five;
    double* six;
    double* seven;
    double* eight;
    double* nine;
    double* ten;
    double* eleven;
    double* twelve;
    double* thirteen;
    double* fourteen;
    double* fifteen;
    double* sixteen;
    double* centerx;
    double* centery;
    double* rmax;
    double* coscount;
    double* costotal;

    double* normal;
    double* radius;
    double* radiusArea;
    double* profile;

    int* surfacecoating;
    int* surfacetype;
    int* surfacemed;
    
    Grating** ppGrating;

    void setup(long surfaceTotal, long points);
    void asphere(long surfaceIndex, long points);
    
};

#endif
