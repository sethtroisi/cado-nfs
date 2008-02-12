//
// Copyright (C) 2008 INRIA (French National Institute for Research in
// Computer Science and Control)
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//

/**
 * \file    stopwatch.c
 * \author  Jerome Milan
 * \date    Tue Feb 12 2008
 * \version 1.0
 *
 * \brief A very basic stopwatch-like timer.
 *
 * This file implements a very basic stopwatch-like timer (with microsecond
 * precision) based on the \c rusage structure and using the \c getrusage
 * function.
 */

#include <stdlib.h>
#include "stopwatch.h"

//----------------------------------------------------------------------------
inline void init_stopwatch(stopwatch_t* const watch) {    
    watch->elapsed_usec = 0;
    watch->is_running   = false;
}
//----------------------------------------------------------------------------
inline void start_stopwatch(stopwatch_t* const watch) {    
    if (! watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->started_usec  = (uint64_t) watch->rsg->ru_utime.tv_sec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_sec;
        watch->started_usec *= 1000000;        
        watch->started_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
        
        watch->is_running = true;
    }
}
//----------------------------------------------------------------------------
inline void stop_stopwatch(stopwatch_t* const watch) {
    if (watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->elapsed_usec += (   (uint64_t) watch->rsg->ru_utime.tv_sec
                                 + (uint64_t) watch->rsg->ru_stime.tv_sec
                               ) * 1000000;                                         
        watch->elapsed_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->elapsed_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
        watch->elapsed_usec -= watch->started_usec;
        
        watch->is_running  = false;
    }
}
//----------------------------------------------------------------------------
inline void reset_stopwatch(stopwatch_t* const watch) {
    watch->elapsed_usec = 0;
    if (watch->is_running) {
        getrusage(RUSAGE_SELF, watch->rsg);
        
        watch->started_usec  = (uint64_t) watch->rsg->ru_utime.tv_sec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_sec;
        watch->started_usec *= 1000000;        
        watch->started_usec += (uint64_t) watch->rsg->ru_utime.tv_usec;
        watch->started_usec += (uint64_t) watch->rsg->ru_stime.tv_usec;
    }
}
//----------------------------------------------------------------------------
inline double get_stopwatch_elapsed(stopwatch_t* const watch) {
    return  (watch->elapsed_usec / 1000000.0);
}
//----------------------------------------------------------------------------
