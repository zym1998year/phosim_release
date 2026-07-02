///
/// @package phosim
/// @file event.h
/// @brief event logger header
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Glenn Sembroski (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <fitsio.h>
#include <iostream>
#include <string>
#include <vector>

struct PhotonLog {

    double* pX;
    double* pY;
    double* pZ;
    int* pSurface;
    long counter;
    long lastRowWritten;
    long maxNumPhotons;
    long bufferSize;
    std::string outputDir;
    std::string eventFitsFileName;
    fitsfile* pFitsFile;
 };

class EventFile {

 private:
    PhotonLog* pfPhotonLog;

 public:
    EventFile(long maxNumPhotons, std::string outputDir);
    EventFile(long maxNumPhotons, std::string outputDir, std::string eventFileName);
    ~EventFile();

    void eventFileInit(long maxNumPhotons, std::string outputDir,std::string eventFileName);
    void logPhoton(double a, double b, double c, int d);
    void eventFileCreate();
    void eventFileWrite(long long firstRow,long long numRowsToWrite);
    void eventFileClose();
};
