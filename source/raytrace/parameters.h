///
/// @package phosim
/// @file parameters.h
/// @brief parameters
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

#ifndef Par_H
#define Par_H

// general numerical accuracy
int const PERTURBATION_POINTS = 256;     // 1/X of D
int const DET_RESAMPLE = 4;                          // X pixels
int const ATM_RESAMPLE = 2;                         // X cm
double const SED_TOLERANCE = 2.0;            // X nm
int const MIE_OUT_POINT = 16384;                 // split of cos
int const AIR_OPACITY_POINT = 8192;            // split of entire wavelength range
int const SURFACE_POINTS = 65536;                 //  1/X of D
int const GALAXY_POINT = 8192;                     // steps of 0.01 of radii
int const SCATTERING_POINTS = 65536;        // total steps
double const SCATTERING_LIMIT = 40.0;         // limit in units of lambda/b

//should have minimal effect on memory
//#define PERTURBATION_COORDINATE
int const SCREEN_SIZE = 1024/ATM_RESAMPLE;
int const SILICON_STEPS = 1024;
double const GALAXY_SCALE = 10.0;
int const MAX_SURF = 256;
int const MAX_LAYER = 16;
int const MAX_IMAGE = 256;
int const OPD_SCREEN_SIZE = 255;
int const OPD_SAMPLING = 4;
int const NZERN = 28;
int const NCHEB = 21;
int const NPOLY = 49;
int const AIR_HEIGHT = 128;
int const AIR_HEIGHT_SCATTER = 64;
int const MIE_IN_POINT = 8;
int const MIE_TAU_POINT = 256;
int const MAX_CHIP = 512;
int const MAX_CATALOG = 32;
int const MAX_BOUNCE = 128;
int const RAYTRACE_MAX_ITER = 10;
double const RAYTRACE_TOLERANCE = 1e-7;           // mm
double const RAYTRACE_MAX_LENGTH = 40000.0;  // mm
double const RAYTRACE_ERROR = 10.0;                       // mm

#endif
