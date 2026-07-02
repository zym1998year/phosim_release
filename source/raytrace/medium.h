///
/// @package phosim
/// @file medium.h
/// @brief header file for medium class
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

class Medium {

public:
    int indexRefractionNumber[MAX_SURF];
    double *indexRefraction[MAX_SURF];
    double *indexRefractionWavelength[MAX_SURF];
    double *indexRefractionImag[MAX_SURF];

    void setup(int surfaceIndex, std::string mediumFileName);

};
