///
/// @package phosim
/// @file coating.cpp
/// @brief coating class
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

#include "coating.h"

void Coating::setup (int totalSurfaces) {

    wavelengthNumber = new int[totalSurfaces]();
    angleNumber = new int[totalSurfaces]();

    reflection = new double*[totalSurfaces]();
    transmission = new double*[totalSurfaces]();
    wavelength = new double*[totalSurfaces]();
    angle = new double*[totalSurfaces]();
    thickness = new float*[totalSurfaces]();
    nominalThickness = new float[totalSurfaces]();
    nonUniformity = new int[totalSurfaces]();
    //Must initally be defined (even if defined as NULL)
    for (int i = 0; i <totalSurfaces; i++) {
        wavelength[i] = NULL;
        angle[i] = NULL;
        nonUniformity[i] = 0;
    }

}

void Coating::allocate (int surfaceIndex, int lines) {

    transmission[surfaceIndex] = new double[lines]();
    reflection[surfaceIndex] = new double[lines]();
    wavelength[surfaceIndex] = new double[lines]();
    angle[surfaceIndex] = new double[lines]();


}
