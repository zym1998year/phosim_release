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

#include <cstddef>
#include "parameters.h"

class Coating {

public:
    double** reflection;
    double** transmission;
    double** wavelength;
    double** angle;

    int* wavelengthNumber;
    int* angleNumber;

    float** thickness;
    int* nonUniformity;
    float* nominalThickness;

    
    void setup(int totalSurface);
    void allocate(int surfaceIndex, int lines);

};
