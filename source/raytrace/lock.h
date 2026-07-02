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

struct Lock {
    pthread_mutex_t lock1;
    pthread_mutex_t lock2;
    pthread_mutex_t lock3;
    pthread_mutex_t lock4;
    pthread_mutex_t lock5;
    pthread_mutex_t lock6;
    pthread_mutex_t lock7;
    pthread_mutex_t lock8;
    pthread_mutex_t lock9;
    pthread_mutex_t lock10;
    pthread_mutex_t lock11;
    pthread_cond_t cond;
};
