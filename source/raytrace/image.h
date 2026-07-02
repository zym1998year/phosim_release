///
/// @package phosim
/// @file image.h
/// @brief header for image class
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

#include <fitsio.h>
#include <fitsio2.h>
#include <fftw3.h>
#include <vector>

#include "ancillary/random.h"
#include "galaxy.h"
#include "dust.h"
#include "observation.h"
#include "raytrace.h"
#include "surface.h"
#include "coating.h"
#include "silicon.h"
#include "perturbation.h"
#include "screen.h"
#include "air.h"
#include "medium.h"
#include "obstruction.h"
#include "chip.h"
#include "contamination.h"
#include "grating.h"
#include "photon.h"
#include "state.h"
#include "lock.h"
#include "mie.h"

class Image : public Observation {

    public:

    // objects and structures
    Galaxy galaxy;
    Dust dust;
    Surface surface;
    Coating coating;
    Silicon silicon;
    Air air;
    Perturbation perturbation;
    Screen screen;
    Medium medium;
    Obstruction obstruction;
    Contamination contamination;
    Chip chip;
    Grating* pGrating;
    State state;
    Random* random;
    Mie mie;

    // setup and loop methods
    int atmSetup();
    int telSetup();
    int sourceLoop();
    void photonLoop(int ssource, int thread, int finish);
    void particleLoop(double flux, double energyMin, double energyMax, double energyIndex, double angularIndex, int particleType, int downward, long *nparticles);

    // thread and optimization
    int dynamicTransmissionOptimization(int thread, int k, int *lastSurface, int *preGhost, int waveSurfaceIndex, int straylightcurrent, Photon *aph, int sourceType);
    int dynamicTransmissionOptimizationZero(int *lastSurface, int *preGhost, int waveSurfaceIndex, int straylightcurrent, Photon *aph, int sourceType);
    static void* threadFunction(void *voidArgs);
    Lock lock;
    int remain;
    int openthreads;
    int *openthread;

    // physics methods
    int getWavelengthTime(Photon *aph, int source);
    int domeSeeing(Vector *angle, Photon *aph);
    int tracking(Vector *angle, double time);
    int atmosphericDispersion(Vector *angle, Photon *aph, int layer);
    int largeAngleScattering(Vector *largeAngle, Photon *aph, int surf);
    int secondKick(Vector *largeAngle, Photon *aph);
    int diffraction(Vector *position, Vector angle, Vector *largeAngle, Photon *aph);
    int samplePupil(Vector *position, long ray, Photon *aph);
    int transmissionCheck(float transmission, int surfaceIndex, int waveSurfaceIndex, Photon *aph);
    int transmissionPreCheck(int index, Photon *aph, int rewrite);
    int transmissionPreCheckZero(int index, Photon *aph);
    int transmissionPreCheckSuppression(int index, Photon *aph, int rewrite);
    int transmissionCheckNoRoll(float transmission, int surfaceIndex, int waveSurfaceIndex, Photon *aph);
    int transmissionPreCheckNoRoll(int index, Photon *aph);
    int chooseSurface(int *newSurface, int *oldSurface, Photon *aph);
    int findSurface(Vector angle, Vector position, double *distance, int surfaceIndex, Photon *aph);
    int goldenBisectSurface(double a, double b, double c, double *z, Vector angle,
                            Vector position, double *distance, int surfaceIndex, Photon *aph);
    int getIntercept(double x, double y, double *z, int surfaceIndex, Photon *aph);
    int getDeltaIntercept(double x, double y, double *zv, int surfaceIndex, Photon *aph);
    int bloom(int saturatedFlag, Photon *aph);
    void saturateIn(int source, Vector *largeAngle, Photon *aph, double shiftX, double shiftY, double extraLargeAngle, Vector *largeAngle2);
    void saturateOut(int source, Vector *largeAngle, Photon *aph, double shiftX, double shiftY, double extraLargeAngle);
    int photonSiliconPropagate(Vector *angle, Vector *position, double lambda, Vector normal,
                               double dh, int waveSurfaceIndex, Photon *aph);
    void pixelShift(double, double, double*, double*, double*);
    void pixelPosition(Vector *position, Photon *aph);
    void pixelPositionDrift(Vector *position, Photon *aph);
    void pixelFraction(Vector *position, Photon *aph);
    int electronSiliconPropagate(Vector *angle, Vector *position, Photon *aph, int counter, int nopsf);
    int contaminationSurfaceCheck(Vector position, Vector *angle, int surfaceIndex, Photon *aph);
    int contaminationDetectorCheck(Vector position, Vector *angle, Photon *aph);
    double airIndexRefraction(Photon *aph, int layer, int skipADC, double *ADC);
    double surfaceCoating(double wavelength, Vector angle,
                          Vector normal, int surfaceIndex, float *reflection, Photon *aph);
    double cloudOpacity(int layer, Photon *aph);
    double cloudOpacityMoon(int layer, Photon *aph);
    double cloudOpacitySun(int layer, Photon *aph);
    double atmosphereOpacity(Vector angle, int layer, Photon *aph, int t, int r);
    double atmosphereOpacityRayleigh(Vector angle, int layer, Photon *aph);
    double atmosphereOpacityMoon(Vector angle, int layer, Photon *aph);
    double atmosphereOpacitySun(Vector angle, int layer, Photon *aph);
    double fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization);
    double monolayer (Vector angle, Vector normal, double wavelength, double nSi, double thickness, double meanFreePath, int polarization, double nMid, double nMidimag);
    void getAngle(Vector *angle, double time, int source, Photon *aph);
    int getDeltaAngle(Vector *angle, Vector *position, int source, double *shiftX, double *shiftY, int thread, int *initGal, Photon *aph);
    void newRefractionIndex(int surfaceIndex, Photon *aph);
    void atmospherePropagate(Vector *position, Vector angle, int layer, int mode, Photon *aph);
    void atmosphereIntercept(Vector *position, int layer, Photon *aph);
    void atmosphereRefraction(Vector *angle, int layer, int mode, Photon *aph);
    void atmosphereDiffraction(Vector *angle, Photon *aph);
    void transform(Vector *angle, Vector *position, int surfaceIndex, int focusFlag, Photon *aph);
    void transformInverse(Vector *angle, Vector *position, int surfaceIndex);
    void interceptDerivatives(Vector *normal, Vector position, int surfaceIndex, Photon *aph);
    void gaussianXY(Vector *newPosition, Vector oldPosition, float scale, Photon *aph);
    void cosmicRays(long *raynumber);
    int mieScatter(double tau, Vector *largeAngle, Photon *aph, int actScatter);
    float tauSuppression(float transmission);
    int rejectScatterCloud(double bAngle, Photon *aph);
    int rejectScatterAtmosphere(double bAngle, Photon *aph);
    int inBounds(TwoVectorInt *position);
    int inBoundsBuffer(TwoVectorInt *position, int buffer);
    int inBoundsDrift(TwoVectorInt *position);
    void diagnosticCount(long ray, int count);
    void backgroundProbCount(int bindex, int backgroundProbOk, int ssource);
    double interfaceAveragePolarization (double inputAngle, double nIn, double nInImag, double nOut, double nOutImag);
    double interface (double inputAngle, double nIn, double nInImag, double nOut, double nOutImag, int polarization);
    double interfaceMonolayerAveragePolarization(double inputAngle, double wavelength, double nIn, double nInimag, double nOut, double nOutimag, double thickness, double nMid, double nMidimag);
    double interfaceMonolayer (double inputAngle, double wavelength, double nIn, double nInimag, double nOut, double nOutimag, double thickness, int polarization, double nMid, double nMidimag);

    // output
    void writeImageFile();
    void writeOPD();
    void readCheckpoint(int checkpointcount);
    void writeCheckpoint(int checkpointcount);
    void cleanup();

};

struct thread_args {
    Image *instance;
    int runthread;
    int ssource[16384];
    int thread;
};
