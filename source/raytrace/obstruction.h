///
/// @package phosim
/// @file obstruction.h
/// @brief header file for obstruction class
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

class Obstruction {

public:
    double par[MAX_SURF];
    double center[MAX_SURF];
    double width[MAX_SURF];
    double depth[MAX_SURF];
    double angle[MAX_SURF];
    double reference[MAX_SURF];
    double height[MAX_SURF];
    int type[MAX_SURF];
    int nspid;
    int pupil;

};
