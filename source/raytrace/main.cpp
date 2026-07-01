///
/// @package phosim
/// @file main.cpp
/// @brief main for raytrace code
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

#include "raytrace/image.h"


// #include <filesystem>
// #include <sstream>
// #include <fstream>
// #include <signal.h>     // ::signal, ::raise
// #include <boost/stacktrace.hpp>

// const char* backtraceFileName = "./backtraceFile.dump";

// void signalHandler(int)
// {
//     ::signal(SIGSEGV, SIG_DFL);
//     ::signal(SIGABRT, SIG_DFL);
//     boost::stacktrace::safe_dump_to(backtraceFileName);
//     ::raise(SIGABRT);
// }
// void sendReport()
// {
//     if (std::filesystem::exists(backtraceFileName))
//     {
//         std::ifstream file(backtraceFileName);

//         auto st = boost::stacktrace::stacktrace::from_dump(file);
//         std::ostringstream backtraceStream;
//         backtraceStream << st << std::endl;

//         // sending the code from st

//         file.close();
//         std::filesystem::remove(backtraceFileName);
//     }
// }

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#ifndef _WIN32
#include <execinfo.h>
#include <unistd.h>
#endif

void handler(int sig) {
#ifndef _WIN32
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
#else
  // execinfo/backtrace is glibc-specific; on Windows just report and exit.
  fprintf(stderr, "Error: signal %d\n", sig);
#endif
  exit(1);
}

int main() {
    Image im;

    //     ::signal(SIGSEGV, signalHandler);
    // ::signal(SIGABRT, signalHandler);

    signal(SIGSEGV, handler);
    im.parser();
    im.background();
    im.atmSetup();
    im.telSetup();
    im.sourceLoop();
    im.cleanup();

    return(0);
}
