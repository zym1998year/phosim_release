///
/// @package phosim
/// @file counter.cpp
/// @brief various logging functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <stdlib.h>

#include "sys/time.h"
#include "constants.h"
#include "counter.h"

//#include "sys/resource.h"
//#include <stdlib.h>
//#include <stddef.h>
//#include <stdio.h>
//#include <malloc/malloc.h>

void counterInit(Clog *counterLog) {


    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Type                 Sources Photons/Electrons  ( Sat , Rem , Rej ,  Acc )%% Time (s)        Events/s\n");
    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");

    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;
    struct timeval tim;
    gettimeofday(&tim, NULL);
    counterLog->previousWallTime = tim.tv_sec + (tim.tv_usec/1000000.0);
    counterLog->previousCPUTime = TicksToSec(GetNewTick());
    double val1, val2;
    val1 = counterLog->previousWallTime;
    val2 = counterLog->previousCPUTime;
    counterLog->originalWallTime = val1;
    counterLog->originalCPUTime = val2;
}

void counterClear(Clog *counterLog) {

    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;

}

void counterAdd(Clog *localLog, Clog *counterLog) {

    counterLog->rejected += localLog->rejected;
    counterLog->accepted += localLog->accepted;
    counterLog->removed += localLog->removed;
    counterLog->removed_dt += localLog->removed_dt;
    counterLog->totalPhoton += localLog->totalPhoton;

}

void counterCheck(Clog *counterLog, int sourcecounter, char *name) {

    double newCpuTime, newWallTime;
    long rate;
    char sourceString[4096];
    char photonString[4096];
    char rateString[4096];
    long tp;

    //struct rlimit limit;
    //getrlimit (RLIMIT_STACK, &limit);
    //printf ("Stack Limit = %llu and %llu max\n", limit.rlim_cur, limit.rlim_max);
    //malloc_stats();
    //getrlimit (RLIMIT_RSS, &limit);
    //printf ("RSS Limit = %llu and %llu max\n", limit.rlim_cur, limit.rlim_max);

    newCpuTime = TicksToSec(clock());
    struct timeval tim;
    gettimeofday(&tim, NULL);
    newWallTime = tim.tv_sec + (tim.tv_usec/1000000.0);
    rate = static_cast<int>(counterLog->totalPhoton/(newWallTime-counterLog->previousWallTime + 1e-2));

    if (sourcecounter < 1000) snprintf(sourceString, sizeof(sourceString), "%7d", sourcecounter);
    if (sourcecounter >= 1000) snprintf(sourceString, sizeof(sourceString), "%3d,%03d", sourcecounter/1000, sourcecounter%1000);

    tp = counterLog->totalPhoton;
    if (tp < 1000) snprintf(photonString, sizeof(photonString), "  %15ld", tp);
    if (tp >= 1000 && tp < 1000000) {
        snprintf(photonString, sizeof(photonString), "  %11ld,%03ld", tp/1000, tp%1000);
    }
    if (tp >= 1000000 && tp < 1000000000) {
        snprintf(photonString, sizeof(photonString), "  %7ld,%03ld,%03ld", tp/1000000,
                (tp/1000)%1000, tp%1000);
    }
    if (tp >= 1000000000) {
        snprintf(photonString, sizeof(photonString), "  %3ld,%03ld,%03ld,%03ld", (tp/1000000000),
                (tp/1000000)%1000, (tp/1000)%1000,
                tp%1000);
    }

    if (rate < 1000) snprintf(rateString, sizeof(rateString), "%15ld", rate);
    if (rate >= 1000 && rate < 1000000) snprintf(rateString, sizeof(rateString), "%11ld,%03ld", rate/1000, rate%1000);
    if (rate >= 1000000 && rate < 1000000000) snprintf(rateString, sizeof(rateString), "%7ld,%03ld,%03ld", rate/1000000,
                                                      (rate/1000)%1000, rate%1000);
    if (rate >= 1000000000) snprintf(rateString, sizeof(rateString), "%3ld,%03ld,%03ld,%03ld", (rate/1000000000),
                                    (rate/1000000)%1000, (rate/1000)%1000, rate%1000);

    if (counterLog->totalPhoton >= 1) {
        printf("%s %s %s  (%5.1f,%5.1f,%5.1f,%6.2f)  %9.2f %s\n",name, sourceString, photonString,
                static_cast<double>(counterLog->rejected)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->removed_dt)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->removed)/static_cast<double>(counterLog->totalPhoton)*100,
                static_cast<double>(counterLog->accepted)/static_cast<double>(counterLog->totalPhoton)*100,
                newWallTime - counterLog->previousWallTime, rateString);
    } else {
        printf("%s %s\n",name, sourceString);
    }

    counterLog->previousWallTime = newWallTime;
    counterLog->previousCPUTime = newCpuTime;
    counterLog->rejected = 0;
    counterLog->accepted = 0;
    counterLog->removed = 0;
    counterLog->removed_dt = 0;
    counterLog->totalPhoton = 0;

}


void counterFinish(Clog *counterLog) {
    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");
    printf("Wall Time = %20.2f CPU Time = %20.2f\n",counterLog->previousWallTime - counterLog->originalWallTime,counterLog->previousCPUTime - counterLog->originalCPUTime);
    fprintf(stdout, "----------------------------------------------------------------------------------------------------\n");

}

void countGood(Clog *counterLog, long photons, long *ray) {

    counterLog->accepted += 1;
    counterLog->rejected += photons - 1;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void countBad(Clog *counterLog, long photons, long *ray) {

    counterLog->removed += photons;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void countBad_dt(Clog *counterLog, long photons, long *ray) {

    counterLog->removed_dt += photons;
    counterLog->totalPhoton += photons;
    *ray += photons;

}

void addThroughput (Tlog *throughputlog, double minwavelength, double maxwavelength, int surf, int waveindex, long sourceover) {

    int windex = static_cast<long>(waveindex/10);
    throughputlog->throughput[(surf + 1)*(static_cast<int>((maxwavelength - minwavelength) + 1)) + windex] += static_cast<double>(sourceover);
}

void initThroughput (Tlog *throughputlog, double minwavelength, double maxwavelength, int nsurf) {

    throughputlog->throughput = static_cast<double*>(calloc((nsurf + 2)*((maxwavelength - minwavelength) + 1), sizeof(double)));

}

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename, Tlog *throughputlog,
                          double minwavelength, double maxwavelength, int nsurf) {

    FILE *outdafile;
    char tempstring[4096];

    snprintf(tempstring, sizeof(tempstring), "%s/throughput_%s.txt", outputdir.c_str(), outputfilename.c_str());
    outdafile = fopen(tempstring, "w");
    for (int k = 0; k < (maxwavelength - minwavelength) + 1; k++) {
        fprintf(outdafile, "%d ", k);
        for (int i = 0; i < nsurf + 2; i++) {
            fprintf(outdafile, "%lf ", throughputlog->throughput[i*(static_cast<long>((maxwavelength - minwavelength) + 1)) + k]);
        }
        fprintf(outdafile, "\n");
    }
    fclose(outdafile); 

}

void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename,
                        std::vector<long> sourceSaturation, std::vector<long> sourceXpos, std::vector<long> sourceYpos,
                        std::vector<std::string> sourceId, int nsource) {

    FILE *outdafile;
    char tempstring[4096];


    snprintf(tempstring, sizeof(tempstring), "%s/centroid_%s.txt", outputdir.c_str(), outputfilename.c_str());
        outdafile = fopen(tempstring, "w");
        fprintf(outdafile, "SourceID Photons AvgX AvgY\n");
        for (int k = 0; k < nsource; k++) {
            fprintf(outdafile, "%s %ld %lf %lf\n", sourceId[k].c_str(), sourceSaturation[k],
                    (static_cast<double>(sourceXpos[k]))/(static_cast<double>(sourceSaturation[k])),
                    (static_cast<double>(sourceYpos[k]))/(static_cast<double>(sourceSaturation[k])));
        }
        fclose(outdafile);

}
