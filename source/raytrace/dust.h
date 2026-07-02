///
/// @package phosim
/// @file dust.h
/// @brief dust header file
///
/// @brief Created by
/// @author James Pizagno (UW)
///
/// @brief Modified by:
/// @author John Peterson (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

class Dust {

public:

    double *kGrid;
    double *aGrid;
    double *bGrid;
    double *wavelengthGrid;

    void setup (double maxwavelength);
    double calzetti (double wavelength, double maxwavelength, double av, double rv);
    double calzettiSetup(double l);
    double ccm (double wavelength, double maxwavelength, double av, double rv);
    void ccmSetup(double l, double *a, double *b);

};
