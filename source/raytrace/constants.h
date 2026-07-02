///
/// @package phosim
/// @file constants.h
/// @brief fundamental constants
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

// MATH

#ifdef M_PI
    #undef M_PI
#endif

#define PI        (3.141592653589793238462643)
#define ARCSEC    (PI/180.0/3600.0)
#define DEGREE    (PI/180.0)
#define HALFSQ5M1 (0.5*(sqrt(5) - 1.0))
#define HALF3MSQ5 (0.5*(3.0 - sqrt(5)))

// BASE PHYSICS

#define G_NEWTON      6.67430e-11             // Newton kg^-2 m^2 (+/- 0.00015e-11)
#define C_LIGHT       2.99792458e8               // m s^-1        (exact)
#define H_PLANCK      6.62607015e-34       // J s                 (exact)
#define E_CHARGE      1.6022176634e-19    // Coloumbs  (exact)
#define K_BOLTZMANN   1.380649e-23     // Joules K^-1 (exact)
#define N_AVOGADRO    6.02214076e23    // mol^-1        (exact)
#define K_CANDELA     683.0                      // lm W^-1       (exact)

// PARTICLES AND NUCLEI

#define M_PROTON      1.672621898e-27  // kg   (+/- 0.000000021e-27)
#define M_ELECTRON    9.10938356e-31   // kg   (+/- 0.00000011e-31)
#define M_MUON        1.883531627e-28  // kg   (+/- 0.000000042e-28)
#define M_ALPHA       6.6446573357e-27 // kg   (+/- 0.0000000020e-27)

// MOLECULES AND ATOMS

#define M_H           1.008            // g mol^-1
#define M_H2          2.01588          // g mol^-1
#define M_HE          4.002602         // g mol^-1
#define M_O           15.999           // g mol^-1
#define M_CH4         16.0425          // g mol^-1
#define M_OH          17.007           // g mol^-1
#define M_NH3         17.013           // g mol^-1
#define M_H2O         18.0153          // g mol^-1
#define M_HF          20.006           // g mol^-1
#define M_NE          20.1797          // g mol^-1
#define M_C2H2        26.04            // g mol^-1
#define M_HCN         27.025           // g mol^-1
#define M_CO          28.01            // g mol^-1
#define M_N2          28.0134          // g mol^-1
#define M_SI          28.085           // g mol^-1
#define M_NO          30.006           // g mol^-1
#define M_H2CO        30.026           // g mol^-1
#define M_C2H6        30.07            // g mol^-1
#define M_PH3         33.998           // g mol^-1
#define M_O2          31.9988          // g mol^-1
#define M_H2O2        34.015           // g mol^-1
#define M_HCL         36.46            // g mol^-1
#define M_AR          39.948           // g mol^-1
#define M_N2O         44.013           // g mol^-1
#define M_CO2         44.0095          // g mol^-1
#define M_NO2         46.006           // g mol^-1
#define M_O3          47.998           // g mol^-1
#define M_CH3CL       50.49            // g mol^-1
#define M_CLO         51.45            // g mol^-1
#define M_HOCL        52.46            // g mol^-1
#define M_OCS         60.08            // g mol^-1
#define M_SO2         64.07            // g mol^-1
#define M_HNO3        63.013           // g mol^-1
#define M_HBR         80.91            // g mol^-1
#define M_KR          83.798           // g mol^-1
#define M_CD          112.41           // g mol^-1
#define M_TE          127.6            // g mol^-1
#define M_HI          127.912          // g mol^-1
#define M_XE          131.293          // g mol^-1
#define M_HG          200.59           // g mol^-1
#define M_CDTE     240.01            // g mol^-1
#define M_HGTE     328.19           // g mol^-1


// DERIVED PHYSICS

#define CELCIUS_KELVIN 273.16  // (exact)
#define EPSILON_0  (1e7/C_LIGHT/C_LIGHT/4.0/PI)  // (exact)
#define ELECTRON_RADIUS (E_CHARGE*E_CHARGE/(4*PI*EPSILON_0*M_ELECTRON*C_LIGHT*C_LIGHT)) // (exact)
#define K_DEDX (4*PI*N_AVOGADRO*ELECTRON_RADIUS*ELECTRON_RADIUS*M_ELECTRON*C_LIGHT*C_LIGHT) // (exact)

// MEASURED PHYSICS
#define EPSILON_SI 11.7                                                        // Relative permittivity for Si
#define PAIRENERGY_SI 3.62
#define PAIRENERGY_MCT 2.62
#define FANO_SI 0.1
#define FANO_MCT 0.1
#define EXCITATION_SI 173.0
#define EXCITATION_HGTE 660.0
#define EXCITATION_CDTE 500.0
#define DENSITY_SI 2.32
#define DENSITY_HGTE 8.1
#define DENSITY_CDTE 5.8
#define DENSITY_H2O 0.997


// CGS UNITS (TO BE REMOVED LATER)
#define H_CGS (H_PLANCK*1e7)                                                            // ergs s  (exact)
#define C_CGS (C_LIGHT*100.0)                                                              // cm s^-1 (exact)
#define EPSILON_0_CM (1e7/C_LIGHT/C_LIGHT/4.0/PI/100.0)    // Farads cm^-1 (exact)

// ASTRONOMY

#define RADIUS_EARTH 6378.1366   // km     (standard value)
#define EARTH_SUN 149597870.700  // km     (standard value)
#define EARTH_MOON 384403.0      // km
#define GRAV_ACCEL_EARTH 9.80665 // m s^-2 (standard value)
#define PRESS_ATMOSPHERE 101325  // Pa     (standard value)
#define PRESS_ATMOSPHERE_MMHG 760 // mmHg (standard value)
#define PRESS_ATMOSPHERE_MB 1000 // mbar (standard value)
