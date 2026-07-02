///
/// @package phosim
/// @file trim.h
/// @brief trim header file
///
/// @brief Created by:
/// @author Alan Meert (Purdue)
///
/// @brief Modified by:
/// @author Justin Bankert (Purdue)
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "raytrace/parameters.h"

#include <vector>
#include <string>

/// @addtogroup trim_group trim
/// @{

class Trim {

public:
    std::string instrdir;
    std::vector<std::string> catalog;
    std::vector<std::string> chipid;
    std::string catalogOpd;
    char outputFilename[MAX_CHIP][4096];
    char outputOpdFilename[4096];
    int nChip;
    int nCatalog;
    int buffer;
    int filter;
    int strayLight;
    int minSource;
    int opdMode;
    int flipX;
    int flipY;
    int swap;
    std::string obshistid;
    long flatDirectory;
    double pointingRa;
    double pointingDec;
    double rotatez;
    double focalLength;
    double plateScale;
    double scale;
    double extendedBuffer;
    double maglimit;
    double magconst;
    double xPosition[MAX_CHIP];
    double yPosition[MAX_CHIP];
    double xDimension[MAX_CHIP];
    double yDimension[MAX_CHIP];
    double pixelSize[MAX_CHIP];
    double angle[MAX_CHIP];
    double deltaX[MAX_CHIP];
    double deltaY[MAX_CHIP];

    void xyPosition(double alpha, double delta, double *x, double *y, double *dis);
    void readCatalog();
    void readOpdCatalog();
    void getDetectorProperties(int d);
    void setup();

};


/// @}
