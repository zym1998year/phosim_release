///
/// @package phosim
/// @file contamination.h
/// @brief header file for contamination class
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

#include <vector>

#include "ancillary/random.h"

class Contamination {

 public:

    float* transmission;
    long* surfacelistmap;
    long* surfacelist2ndmap;
    std::vector<std::vector<double> > surfacelistx;
    std::vector<std::vector<double> > surfacelisty;
    std::vector<std::vector<double> > surfacelists;
    std::vector<std::vector<int> > surfacelistt;
    long* surfacenumber;
    float* chiptransmission;
    std::vector<double> chiplistx;
    std::vector<double> chiplisty;
    std::vector<double> chiplists;
    int* chiplistmap;
    int* chiplist2ndmap;
    int chipnumber;
    double *henyey_greenstein;
    double *henyey_greenstein_mu;
    double *henyey_greenstein_w;
    double *henyey_greenstein_mu_w;
    double absorptionLength;
    int elements;
    Random random;
    float *surfaceDustTransmission;
    float *surfaceCondensationTransmission;

    void setup(double *innerRadius, double *outerRadius, long totalSurface, long points, double pixsize, long nx, long ny, float qevariation, long seedchip);

};
