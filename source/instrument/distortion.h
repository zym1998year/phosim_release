#include "ancillary/readtext.h"

int distortion(int surfaceIndex, int secondSurfaceIndex, long points, double tolerance, int control, double wfeError, double actError, long seed, const std::string & obsID, int filter, const std::string & instrdir,
               double zenith, double temperature, double azimuth, double temperatureChange, double angTol, double initStep, int zernikeStart, int actuatorStart, double *moveOffset, int pertDebug, double gravity, long *diagCall, int nthread, int burnin, int maxiter);
