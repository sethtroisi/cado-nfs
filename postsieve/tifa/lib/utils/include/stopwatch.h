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
 * \file    stopwatch.h
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

#if !defined(_TIFA_STOPWATCH_H_)
#define _TIFA_STOPWATCH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/resource.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdbool.h>

   /**
    * \struct struct_stopwatch_t stopwatch.h lib/utils/include/stopwatch.h
    * \brief  Defines a very basic stopwatch-like timer.
    *
    * This structure defines very basic stopwatch-like timer based on the
    * \c rusage structure.
    */
struct struct_stopwatch_t {
       /**
        * An \c rusage structure.
        */
    struct   rusage rsg[1];
       /**
        * The time (in microseconds) when the stopwatch started.
        */
    uint64_t started_usec;
       /**
        * The elapsed time accumulated (in microseconds).
        */
    uint64_t elapsed_usec;
       /**
        * \c true iif the stopwatch is currently running.
        */
    bool     is_running;
};
   /**
    * \typedef stopwatch_t
    * \brief Equivalent to <tt>struct struct_stopwatch_t</tt>.
    */
typedef struct struct_stopwatch_t stopwatch_t;

   /**
    * \brief Inits a \c stopwatch_t.
    *
    * Initializes the \c stopwatch_t pointed to by \c watch.
    *
    * \param watch The \c stopwatch_t to init.
    */
inline void init_stopwatch(stopwatch_t* const watch);

   /**
    * \brief Starts a \c stopwatch_t.
    *
    * Starts the \c stopwatch_t pointed to by \c watch.
    *
    * \note Consecutive calls to \c start_stopwatch are without effect.
    *
    * \param watch The \c stopwatch_t to start.
    */
inline void start_stopwatch(stopwatch_t* const watch);

   /**
    * \brief Stop a \c stopwatch_t.
    *
    * Stop the \c stopwatch_t pointed to by \c watch. Successive calls to 
    * \c start_stopwatch and \c stop_stopwatch are cumulative. In other words,
    * the stopwatch's state holds the time elapsed during all previous time
    * intervals defined by a call to \c start_stopwatch followed by a call to
    * \c stop_stopwatch (provided that the stopwatch was not reset with
    * \c reset_stopwatch).
    *
    * \note Consecutive calls to \c stop_stopwatch are without effect.
    *
    * \param watch The \c stopwatch_t to stop.
    */
inline void stop_stopwatch(stopwatch_t* const watch);

   /**
    * \brief Reset a \c stopwatch_t.
    *
    * Reset the \c stopwatch_t pointed to by \c watch. The stopwatch is 
    * \e not stopped and will continue to run unless it was already stopped.
    *
    * \param watch The \c stopwatch_t to reset.
    */
inline void reset_stopwatch(stopwatch_t* const watch);

   /**
    * \brief Returns the elapsed time measured.
    *
    * Returns the elapsed time measured by \c watch in seconds.
    *
    * \warning The returned result is only meaningful if the stopwatch is
    * not running (i.e. it has been stopped with the \c stop_stopwatch 
    * function).
    *
    * \param[in] watch The \c stopwatch_t used for timing.
    */
inline double get_stopwatch_elapsed(stopwatch_t* const watch);

#ifdef __cplusplus
}
#endif

#endif
