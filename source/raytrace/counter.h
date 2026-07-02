///
/// @package phosim
/// @file counter.h
/// @brief header for counter
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

#include <time.h>
#include <string>
#include <vector>
#include <atomic>

struct Clog {
    std::atomic<long> rejected;
    std::atomic<long> removed;
    std::atomic<long> accepted;
    std::atomic<long> removed_dt;
    std::atomic<long> totalPhoton;
    std::atomic<double> previousCPUTime;
    std::atomic<double> previousWallTime;
    std::atomic<double> originalCPUTime;
    std::atomic<double> originalWallTime;
};

struct Tlog {
    double *throughput;
};

void counterClear(Clog *counterLog);
void counterAdd(Clog *localLog, Clog *counterLog);
void counterInit (Clog *counterLog);
void counterCheck (Clog *counterLog, int sourcecounter, char *name);
void counterFinish (Clog *counterFinish);
void countBad (Clog *counterLog, long photons, long *ray);
void countBad_dt (Clog *counterLog, long photons, long *ray);
void countGood (Clog *counterLog, long photons, long *ray);

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename,
                          Tlog *throughpuTlog, double minwavelength, double maxwavelength, int nsurf);
void addThroughput (Tlog *throughpuTlog, double minwavelength, double maxwavelength, int surf, int waveindex, long sourceover);
void initThroughput (Tlog *throughpuTlog, double minwavelength, double maxwavelength, int nsurf);
void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename,
                        std::vector<long> sourceSaturation, std::vector<long> sourceXpos, std::vector<long> sourceYpos,
                        std::vector<std::string> sourceId, int nsource);

// Waits for a new tick. Returns number of ticks since program start.
inline clock_t GetNewTick() {
    clock_t prev = clock();
    clock_t cur = clock();
    while (cur == prev) {
        cur = clock();
    }
    return cur;
}

// Converts ticks to seconds
inline double TicksToSec(clock_t ticks) {
    return static_cast<double>(ticks/CLOCKS_PER_SEC);
}
