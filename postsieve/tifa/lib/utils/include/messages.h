//
// Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
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
 * \file    messages.h
 * \author  Jerome Milan
 * \date    Wed Aug 29 2007
 * \version 1.0.1
 *
 * \brief Status / error messages output macros
 *
 * This file defines some macros used to output status or error messages if
 * some of the \c*_VERBOSE and/or \c*_TIMING symbols are set to non-zero.
 *
 * \warning The <tt>__VERBOSE__</tt>, <tt>__TIMING__</tt> and
 * <tt>__PREFIX__</tt> symbols should be defined in the file including this
 * header. It is under the responsability of the including file to check for
 * these symbol definitions or to define them, if needed.
 */

/*
 * History:
 *     1.0.1: Wed Aug 29 2007 by JM:
 *            - Added message for Fermat's algorithm. 
 *     1.0  : Fri Mar 2 2007 by JM:
 *            - Initial version.
 */


#if !defined(_TIFA_MESSAGES_H_)
   /**
    * \def _TIFA_MESSAGES_H_
    * Standard include guard.
    */
#define _TIFA_MESSAGES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "timer.h"
#include "tifa_config.h"

//
// __VERBOSE__, __TIMING__ and __PREFIX__ should be defined in the file
// including this header
//
#if !defined(__VERBOSE__)
#error "__VERBOSE__ must be defined before including this header."
#endif
#if !defined(__TIMING__)
#error "__TIMING__ must be defined before including this header."
#endif
#if !defined(__PREFIX__)
#error "__PREFIX__ must be defined before including this header."
#endif

#if __VERBOSE__ || __TIMING__
    #include <stdio.h>
    #include "exit_codes.h"
#endif

//
// Strings pertaining to all factorization algorithms
//
#define INIT_STRING             "computing start up data..."
#define UPDATE_STRING           "updating context... "
#define UPDATE_GIVEUP_STRING    "updating context... give up"
#define CLEAN_STRING            "cleaning..."
#define FUNC_SUCCEEDED_STRING   "completed"
#define FUNC_FAILED_STRING      "failed"
//
// Strings pertaining to "congruences of squares" methods only
//
#define COLLECT_RELS_STRING      "collecting (X^2 = Y mod N) relations..."
#define COLLECT_RELS_DONE_STRING "collection of (X^2 = Y mod N) relations"
#define FACTOR_RES_STRING        "factoring residues on the factor base..."
#define LIN_ALG_STRING           "resolving linear algebra system..."
#define DED_FACTORS_STRING       "deducing factors..."
#define UPDATE_MORE_RELS_STRING  "updating context to find more relations... "
#define UPDATE_NEW_MULT_STRING   "updating context to change multiplier... "

//
// Strings pertaining to Fermat's factorization only
//
#define FERMAT_FACT_STRING       "performing Fermat's factorization..."
#define FERMAT_FACT_DONE_STRING  "Fermat's factorization..."
#define UPDATE_MULTIPREC_STRING  "updating context to use multi-precision... "
//
// Strings pertaining to SQUFOF only
//
#define FWD_CYCL_STRING    "forward cycling to find a proper form..."
#define INV_SQRT_STRING    "computing inverse square root of form..."
#define REV_CYCL_STRING    "reverse cycling to find a factor..."
#define UPDATE_RACE_STRING "updating context to perform a race... "
//
// Strings pertaining to trial division only
//
#define DIVIDING_STRING_FORMAT "trial dividing by %lu primes..."

//
// Misc. strings
//
#define TAB    "    "
#define SEPARATION_LINE "-----------------------------------" \
                        "------------------------------------"
//
// No messages are printed if __VERBOSE__ and __TIMING__ are set to 0.
//
#if __VERBOSE__ || __TIMING__
    #define PRINTF(...)      printf(__VA_ARGS__); fflush(stdout)
    #define PRFX_PRINTF(...) printf(__PREFIX__ __VA_ARGS__); fflush(stdout)
#else
    #define PRINTF(...)          /* intentionally left empty */
    #define PRFX_PRINTF(...)     /* intentionally left empty */
#endif

//
// If needed, align the right margin by padding with white spaces
//
#define PRINT_LINE(STRING)    PRFX_PRINTF("%s\n",  STRING)
#define PADDED_PRINTF(STRING) PRFX_PRINTF("%-47s", STRING)

//
// No newline after message if __TIMING__ set to 1 to be able to write the
// timing results on the same line
//
#if !__TIMING__
    #define PRINT_MSG(STRING) PRINT_LINE(STRING)
#else
    #define PRINT_MSG(STRING) PADDED_PRINTF(STRING)
#endif

#if __TIMING__
    #define PRINT_SEPARATION_LINE PRINT_LINE(SEPARATION_LINE)
#else
    #define PRINT_SEPARATION_LINE /* intentionally left empty */
#endif

//
// Messages pertaining to all factorization algorithms
//
#define PRINT_FAILURE_NL           PRINTF("failure!\n")
#define PRINT_INIT_MSG             PRINT_MSG(INIT_STRING)
#define PRINT_UPDATE_MSG           PRINT_MSG(UPDATE_STRING)
#define PRINT_UPDATE_GIVEUP_MSG    PRINT_MSG(UPDATE_GIVEUP_STRING)
#define PRINT_CLEAN_MSG            PRINT_MSG(CLEAN_STRING)
   /**
    * \def PRINT_SUCCESS_MSG(ECODE)
    */
#define PRINT_SUCCESS_MSG(ECODE) do {                         \
    char __tmp_str__[64];                                     \
    (void)sprintf(__tmp_str__, FUNC_SUCCEEDED_STRING " (%s)", \
                  ecode_to_str[ECODE]);                       \
    PRINT_MSG(__tmp_str__);                                   \
} while (0)
   /**
    * \def PRINT_FAILURE_MSG(ECODE)
    */
#define PRINT_FAILURE_MSG(ECODE) do {                         \
    char __tmp_str__[64];                                     \
    (void)sprintf(__tmp_str__, FUNC_FAILED_STRING " (%s)",    \
                  ecode_to_str[ECODE]);                       \
    PRINT_MSG(__tmp_str__);                                   \
} while (0)

#if __VERBOSE__ || __TIMING__
    //
    // Introductive messages for all algorithms
    //

    #define PRINT_INTRO_MSG(NAME)                                     \
        do {                                                          \
            PRINTF(__PREFIX__ "attempting %s factorization\n", NAME); \
            PRINT_SEPARATION_LINE;                                    \
        } while (0)
    //
    // Status message giving returned exit code
    //
    #define PRINT_STATUS(SUCCESS, ECODE)      \
        do {                                  \
            PRINT_SEPARATION_LINE;            \
            if (SUCCESS) {                    \
                PRINT_SUCCESS_MSG(ECODE);     \
            } else {                          \
                PRINT_FAILURE_MSG(ECODE);     \
            }                                 \
            PRINT_TOTAL;                      \
        } while (0)
#else
    #define PRINT_STATUS(SUCCESS, ECODE_STR) /* intentionally left empty */
    #define PRINT_INTRO_MSG(NAME)            /* intentionally left empty */
#endif
//
// Messages pertaining to "congruences of squares" methods only
//
#define PRINT_UPDATE_MORE_RELS_MSG PRINT_MSG(UPDATE_MORE_RELS_STRING)
#define PRINT_UPDATE_NEW_MULT_MSG  PRINT_MSG(UPDATE_NEW_MULT_STRING)

#define PRINT_FACTOR_RES_MSG     PRINT_MSG(FACTOR_RES_STRING)
#define PRINT_LIN_ALG_MSG        PRINT_MSG(LIN_ALG_STRING)
#define PRINT_DED_FACTORS_MSG    PRINT_MSG(DED_FACTORS_STRING)

#if __VERBOSE__
    #define PRINT_COLLECT_RELS_MSG PRINT_LINE(COLLECT_RELS_STRING)
#else
    #define PRINT_COLLECT_RELS_MSG PADDED_PRINTF(COLLECT_RELS_STRING)
#endif

#if __VERBOSE__ && __TIMING__
    #define PRINT_COLLECT_RELS_DONE_MSG PRINT_MSG(COLLECT_RELS_DONE_STRING)
#else
    #define PRINT_COLLECT_RELS_DONE_MSG /* intentionally left empty */
#endif

#if __VERBOSE__
    #define PRINT_NRELS_FOUND(FOUND, TO_FIND)                          \
        PRFX_PRINTF(TAB TAB "(%"PRIu32"/%"PRIu32") relations found\n", \
                    FOUND, TO_FIND);
#else
    #define PRINT_NRELS_FOUND(FOUND, TO_FIND) /* intentionally left empty */
#endif

#if __VERBOSE__ && __TIMING__
    #define PRINT_RES_GENERATED_MSG(TIMING)                                   \
        PRFX_PRINTF(TAB "(relations generated in "TIMING_FORMAT" seconds)\n", \
                    TIMING);
    #define PRINT_RES_SELECTED_MSG(TIMING)                                    \
        PRFX_PRINTF(TAB "(relations selected in  "TIMING_FORMAT" seconds)\n", \
                    TIMING);
#else
    #define PRINT_RES_GENERATED_MSG(TIMING) /* intentionally left empty */
    #define PRINT_RES_SELECTED_MSG(TIMING)  /* intentionally left empty */
#endif
//
// Messages pertaining to Fermat's algorithm only
//
#define PRINT_UPDATE_MULTIPREC_MSG  PRINT_MSG(UPDATE_MULTIPREC_STRING)

#if __VERBOSE__
    #define PRINT_FERMAT_FACT_MSG PRINT_LINE(FERMAT_FACT_STRING)
#else
    #define PRINT_FERMAT_FACT_MSG PADDED_PRINTF(FERMAT_FACT_STRING)
#endif

#if __VERBOSE__ && __TIMING__
    #define PRINT_FERMAT_FACT_DONE_MSG PRINT_MSG(FERMAT_FACT_DONE_STRING)
    #define PRINT_NEXTPRIME_MSG(TIMING)                                     \
        PRFX_PRINTF(TAB "primes computed in   "TIMING_FORMAT" seconds\n",   \
                    TIMING);
    #define PRINT_SQRTM_MSG(TIMING)                                         \
        PRFX_PRINTF(TAB "sqrtm computed in    "TIMING_FORMAT" seconds\n",   \
                    TIMING);
    #define PRINT_GREEDY_MSG(TIMING)                                        \
        PRFX_PRINTF(TAB "greedy phase done in "TIMING_FORMAT" seconds\n",   \
                    TIMING);
#else
    #define PRINT_NEXTPRIME_MSG(TIMING) /* intentionally left empty */
    #define PRINT_SQRTM_MSG(TIMING)     /* intentionally left empty */
    #define PRINT_GREEDY_MSG(TIMING)    /* intentionally left empty */
    #define PRINT_FERMAT_FACT_DONE_MSG  /* intentionally left empty */
#endif

//
// Messages pertaining to SQUFOF only
//
#define PRINT_UPDATE_RACE_MSG PRINT_MSG(UPDATE_RACE_STRING)
#define PRINT_FWD_CYCL_MSG    PRINT_MSG(FWD_CYCL_STRING)
#define PRINT_INV_SQRT_MSG    PRINT_MSG(INV_SQRT_STRING)
#define PRINT_REV_CYCL_MSG    PRINT_MSG(REV_CYCL_STRING)
//
// Messages pertaining to trial division only
//
#if __TIMING__ || __VERBOSE__
    #define PRINT_DIVIDING_MSG(X)                              \
        do {                                                   \
            char __tmp_str__[64];                              \
            (void)sprintf(__tmp_str__, DIVIDING_STRING_FORMAT, \
                          (unsigned long int)X);               \
            PRINT_MSG(__tmp_str__);                            \
        } while (0)
#else
    #define PRINT_DIVIDING_MSG(X) /* intentionally left empty */
#endif
//
// Messages pertaining to timings
//
#if __TIMING__
    #define PRINT_NAMED_TIMING(NAME)                      \
            PRINTF("done in " TIMING_FORMAT " seconds\n", \
                   GET_NAMED_TIMING(NAME));
    #define PRINT_NAMED_TOTAL(NAME)                       \
            PRINTF("total   " TIMING_FORMAT " seconds\n", \
                   GET_NAMED_TIMING(NAME));
#else
    #define PRINT_NAMED_TIMING(NAME) /* intentionally left empty */
    #define PRINT_NAMED_TOTAL(NAME)  /* intentionally left empty */
#endif
#define PRINT_TIMING PRINT_NAMED_TIMING()
#define PRINT_TOTAL  PRINT_NAMED_TOTAL()

#ifdef __cplusplus
}
#endif

#endif
