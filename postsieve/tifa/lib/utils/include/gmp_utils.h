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
 * \file    gmp_utils.h
 * \author  Jerome Milan
 * \date    Wed Mar 1 2006
 * \version 1.0
 *
 * \brief Various GMP small utilities.
 *
 * GMP small utilities' definitions should go here.
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_GMP_UTILS_H_)
   /**
    * \def _TIFA_GMP_UTILS_H_
    * Standard include guard.
    */
#define _TIFA_GMP_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>
#include "hashtable.h"

   /**
    * \struct struct_mpz_pair_t array.h lib/utils/include/gmp_utils.h
    * \brief  A pair of \c mpz_t integers.
    *
    * This very simple structure defines a pair of \c mpz_t integers.
    */
struct struct_mpz_pair_t {
        /**
         * The first \c mpz_t integer of the pair.
         */
    mpz_t x;
        /**
         * The second \c mpz_t integer of the pair.
         */
    mpz_t y;
};

   /**
    * \typedef mpz_pair_t
    * \brief Equivalent to <tt>struct struct_mpz_pair_t</tt>.
    */
typedef struct struct_mpz_pair_t mpz_pair_t;

   /**
    * \brief inits a <tt>mpz_pair_t</tt>.
    *
    * Inits a <tt>mpz_pair_t</tt> by initializing each of its \c mpz_t
    * element.
    *
    * \param[in] pair A pointer to the <tt>mpz_pair_t</tt> to init.
    */
inline void init_mpz_pair(mpz_pair_t* pair);

   /**
    * \brief Clears a <tt>mpz_pair_t</tt>.
    *
    * Clears a <tt>mpz_pair_t</tt>.
    *
    * \param[in] pair A pointer to the <tt>mpz_pair_t</tt> to clear.
    */
inline void clear_mpz_pair(mpz_pair_t* pair);

   /**
    * \brief Empties a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s.
    *
    * Empties a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s and clears
    * the memory associated to the keys and their associated data.
    *
    * \param[in] pair A pointer to the <tt>hashtable_t</tt> to empty.
    */
void empty_mpzpair_htable(hashtable_t* const htable);

   /**
    * \brief Clears a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s.
    *
    * Clears a <tt>hashtable_t</tt> holding <tt>mpz_pair_t</tt>'s. It clears
    * the memory associated to the keys, their associated data and the
    * hashtable itself.
    *
    * \param[in] pair A pointer to the <tt>hashtable_t</tt> to clear.
    */
void clear_mpzpair_htable(hashtable_t* const htable);


#ifdef __cplusplus
}
#endif

#endif
