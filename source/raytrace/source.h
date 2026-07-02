///
/// @package phosim
/// @file source.h
/// @brief source structure
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///
///

enum SourceTypes {POINT=0, IMAGE=1, GAUSSIAN=2, MOVINGPOINT=4, SERSIC=5, SERSIC2D=6, PINHOLE=7,
                  OPD=8, SERSICCOMPLEX=9, SERSICDISK=10, SERSICDISKCOMPLEX=11, DISTORTEDSPHERE=12,
                  SPHERE=13, REFLECTEDSPHERE=14, DISK=15, MOVINGSPHERE=16, MOVINGDISK=17, MOVINGREFLECTEDSPHERE=18};

struct Source {
    std::vector<double> ra;
    std::vector<double> redshift;
    std::vector<double> gamma1;
    std::vector<double> gamma2;
    std::vector<double> kappa;
    std::vector<double> dec;
    std::vector<int> split;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> norm;
    std::vector<double> mag;
    std::vector<std::vector<double> > spatialpar;
    std::vector<std::vector<double> > dustpar;
    std::vector<std::vector<double> > dustparz;
    std::vector<std::string> id;
    std::vector<int> spatialtype;
    std::vector<int> dusttype;
    std::vector<int> dusttypez;
    std::vector<int> type;
    std::vector<int> skysameas;
    std::vector<int> sedptr;
    std::vector<double> backgroundAngle;
};
