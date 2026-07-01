///
/// @package phosim
/// @file lock.h
/// @brief header for lock
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

#include <mutex>
#include <condition_variable>

struct Lock {
    std::mutex lock1;
    std::mutex lock2;
    std::mutex lock3;
    std::mutex lock4;
    std::mutex lock5;
    std::mutex lock6;
    std::mutex lock7;
    std::mutex lock8;
    std::mutex lock9;
    std::mutex lock10;
    std::mutex lock11;
    std::condition_variable cond;
};
