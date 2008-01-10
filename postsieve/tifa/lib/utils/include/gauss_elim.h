//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
// in Computer Science and Control)
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
 * \file    gauss_elim.h
 * \author  Jerome Milan
 * \date    Wed Mar 1 2006
 * \version 1.0
 *
 * \brief Gaussian elimination over GF(2) (from a paper by D. Parkinson and
 * M. Wunderlich).
 *
 * Gaussian elimination over GF(2) as presented in the paper "A compact
 * algorithm for Gaussian elimination over GF(2) implemented on highly
 * parallel computers" written by Dennis Parkinson and Marvin Wunderlich
 * (Parallel Computing 1, 1984).
 *
 * \see "A compact algorithm for Gaussian elimination over GF(2) implemented
 * on highly parallel computers", D. Parkinson and M. Wunderlich,
 * <i>Parallel Computing 1</i>, 1984, pages 65-73.
 */

#if !defined(_TIFA_GAUSS_ELIM_H_)
   /**
    * \def _TIFA_GAUSS_ELIM_H_
    * Standard include guard.
    */
#define _GAUSS_ELIM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include "matrix.h"
#include "x_array_list.h"

   /**
    * \brief Gaussian elimination on a <tt>binary_matrix_t</tt>.
    *
    * Performs a gaussian elimination on a <tt>binary_matrix_t</tt> as
    * described in the paper "A compact algorithm for Gaussian elimination
    * over GF(2) implemented on highly parallel computers", by D. Parkinson
    * and M. Wunderlich (Parallel Computing 1 (1984) 65-73).
    *
    * Solutions (if any) of this linear system are stored in <tt>relations</tt>
    * where each entry is a \c uint32_array_t containing the indexes of the rows
    * (from the \e original matrix) composing a solution. In other words, for
    * each entry, the sum of the indexed rows (from the original matrix) is a
    * nul binary vector.
    *
    * \param[in]  matrix A pointer to the <tt>binary_matrix_t</tt> giving
    *                    the linear system to solve.
    * \param[out] relations A pointer to a \c uint32_array_list_t holding the
    *                       solutions of the system, if any.
    *
    */
void gaussian_elim(uint32_array_list_t* relations,
                   binary_matrix_t* const matrix);

#ifdef __cplusplus
}
#endif

#endif
