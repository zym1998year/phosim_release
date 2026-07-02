///
/// @package phosim
/// @file state.h
/// @brief header for state structure
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

#include <atomic>
#include "event.h"
#include "counter.h"

struct State {

    std::atomic<int> *satmap;
    double *opd;
    double *opdcount;
    std::atomic<float> *dynamicTransmission;
    std::atomic<float> *dynamicTransmissionLow;
    std::atomic<int> *dynamicTransmissionInit;
    std::atomic<unsigned int> *focalPlane;
    int surfaceLimit;
    int waveLimit;
    float *focalPlaneFloat;
    double *cx;
    double *cy;
    double *cz;
    double *r0;
    double *epR;
    std::atomic<unsigned int> *bpAcc;
    std::atomic<unsigned int> *bpRej;
    std::atomic<long> *diagnosticCounter;
    EventFile* pEventLogging;
    Clog counterLog;
    Tlog throughputLog;
    Clog globalLog;

};
