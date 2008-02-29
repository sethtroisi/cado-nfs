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
 * \file    matrix.h
 * \author  Jerome Milan
 * \date    Tue Jun 19 2007
 * \version 1.1
 *
 * \brief Matrices over GF(2) and associated functions.
 *
 * Defines binary matrices (i.e. matrices over GF(2)) and their
 * associated functions.
 */

 /*
  *  History:
  *    1.1: Tue Jun 19 2007 by JM:
  *          - Added clone_binary_matrix function.
  *    1.0: Wed Mar 1 2006 by JM:
  *          - Initial version.
  */

#if !defined(_TIFA_MATRIX_H_)
   /**
    * \def _TIFA_MATRIX_H_
    * Standard include guard.
    */
#define _TIFA_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "tifa_config.h"

#include <inttypes.h>

   /**
    * \def NO_SUCH_ROW
    * Value returned by the <tt>first_row_with_one_on_col(col, matrix)</tt>
    * function if no row of the matrix has a bit 1 in its \c col-th column.
    */
#define NO_SUCH_ROW UINT32_MAX

/*
 *-----------------------------------------------------------------------------
 *              binary_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_binary_matrix_t matrix.h lib/utils/include/matrix.h
    * \brief  Defines a matrix of bits.
    *
    * This structure defines a matrix of bits which knows its current
    * dimensions and its allocated memory space.
    *
    * \note Internally, a matrix of bits is represented as a matrix of
    * \c TIFA_BITSTRING_T elements.
    */
struct struct_binary_matrix_t {
       /**
        * Maximum number of rows of the matrix.
        */
    uint32_t nrows_alloced;
       /**
        * Maximum number of columns of the matrix.
        * Since bits are packed in \c TIFA_BITSTRING_T elements, the maximum
        * number of bits per column is 8 * <tt>sizeof(TIFA_BITSTRING_T)</tt> *
        * \c nrows_alloced.
        *
        * Hence the total allocated memory for the \c data field is
        * <tt>nrows_alloced * ncols_alloced</tt> *
        * <tt>sizeof(TIFA_BITSTRING_T)</tt> bytes.
        */
    uint32_t ncols_alloced;
       /**
        * Current number of rows of the matrix.
        */
    uint32_t nrows;
       /**
        * Current number of columns of the matrix.
        */
    uint32_t ncols;
       /**
        * 2D array of \c TIFA_BITSTRING_T elements whose dimensions are given
        * by the \c nrows_alloced and \c ncols_alloced fields.
        */
    TIFA_BITSTRING_T** data;
};
   /**
    * \typedef binary_matrix_t
    * \brief Equivalent to <tt>struct struct_binary_matrix_t</tt>.
    */
typedef struct struct_binary_matrix_t binary_matrix_t;

   /**
    * \brief Allocates and returns a new <tt>binary_matrix_t</tt>.
    *
    * Allocates and returns a new <tt>binary_matrix_t</tt> such that:
    * \li its \c nrows_alloced field is set to nrows.
    * \li its \c ncols_alloced field is set to the minimum number
    *  of \c TIFA_BITSTRING_T integers needed to store ncols bits.
    * \li its \c nrows field is set to nrows.
    * \li its \c ncols field is set to ncols.
    * \li its \c data array of arrays is completely filled with zeroes.
    *
    * \note The behaviour of this alloc function differs from the ones in
    * \link array.h \endlink . This discrepancy will be corrected in later
    * versions.
    *
    * \param[in] nrows The maximum number of rows of the \c binary_matrix_t
    *                    to allocate.
    * \param[in] ncols The maximum number of bits per row of the \c
    *                    binary_matrix_t to allocate.
    * \return A pointer to the newly allocated \c binary_matrix_t structure.
    * Note that this matrix may hold more that \c ncols bits per column if
    * \c ncols is not a multiple of 8 * <tt>sizeof(TIFA_BITSTRING_T)</tt>.
    */
binary_matrix_t* alloc_binary_matrix(uint32_t nrows, uint32_t ncols);

   /**
    * \brief Allocates and returns a cloned <tt>binary_matrix_t</tt>.
    *
    * Allocates and returns a clone of the <tt>binary_matrix_t</tt> 
    * pointed by <tt>matrix</tt>.
    *
    * \param[in] matrix A pointer to the binary matrix to clone.
    * \return A pointer to the newly allocated \c binary_matrix_t clone.
    */
binary_matrix_t* clone_binary_matrix(const binary_matrix_t * const matrix);

   /**
    * \brief Sets a <tt>binary_matrix_t</tt> to zero.
    *
    * Sets the <tt>binary_matrix_t matrix</tt> to the zero matrix.
    *
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt> to reset.
    */
void reset_binary_matrix(binary_matrix_t* const matrix);

   /**
    * \brief Clears a <tt>binary_matrix_t</tt>.
    *
    * Clears a <tt>binary_matrix_t</tt>, or, more precisely, clears the memory
    * space used by the arrays pointed by the \c data field of a
    * <tt>binary_matrix_t</tt>. Also set its \c nrows_alloced ,
    * \c ncols_alloced , \c nrows and \c ncols fields to zero.
    *
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt> to clear.
    */
void clear_binary_matrix(binary_matrix_t* const matrix);

   /**
    * \brief Prints a <tt>binary_matrix_t</tt>.
    *
    * Prints a <tt>binary_matrix_t</tt>'s on the standard output.
    *
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt> to print.
    */
void print_binary_matrix(const binary_matrix_t* const matrix);

   /**
    * \brief Returns the value of a given bit in a <tt>binary_matrix_t</tt>.
    *
    * Returns the value of the bit at the <tt>row</tt>-th row and the
    * <tt>col</tt>-th column of the <tt>binary_array_t</tt> pointed by
    * \c array, as either 0 or 1.
    *
    * \param[in] row The row of the bit to read.
    * \param[in] col The column of the bit to read.
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt>.
    * \return The value of the bit at the (<tt>row</tt>,<tt>col</tt>)
    *         position: either 0 or 1.
    */
uint8_t get_matrix_bit(uint32_t row, uint32_t col,
                       const binary_matrix_t* const matrix);

   /**
    * \brief Sets to one the value of a given bit in a <tt>binary_matrix_t</tt>.
    *
    * Sets to one the value of the bit at the <tt>row</tt>-th row and the
    * <tt>col</tt>-th column of the <tt>binary_array_t</tt> pointed by
    * \c array.
    *
    * \param[in] row The row of the bit to set.
    * \param[in] col The column of the bit to set.
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt>.
    */
void set_matrix_bit_to_one(uint32_t row, uint32_t col,
                           binary_matrix_t* const matrix);

   /**
    * \brief Sets to zero the value of a given bit in a
    * <tt>binary_matrix_t</tt>.
    *
    * Sets to zero the value of the bit at the <tt>row</tt>-th row and the
    * <tt>col</tt>-th column of the <tt>binary_array_t</tt> pointed by
    * \c array.
    *
    * \param[in] row The row of the bit to set.
    * \param[in] col The column of the bit to set.
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt>.
    */
void set_matrix_bit_to_zero(uint32_t row, uint32_t col,
                            binary_matrix_t* const matrix);

   /**
    * \brief Flips the value of a given bit in a <tt>binary_matrix_t</tt>.
    *
    * Flips the value of the bit at the <tt>row</tt>-th row and the
    * <tt>col</tt>-th column of the <tt>binary_array_t</tt> pointed by
    * \c array.
    *
    * \param[in] row The row of the bit to flip.
    * \param[in] col The column of the bit to flip.
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt>.
    */
void flip_matrix_bit(uint32_t row, uint32_t col, binary_matrix_t* const matrix);

   /**
    * \brief Returns the index of the first row of a <tt>binary_matrix_t</tt>
    *        with a one in a given column.
    *
    * Returns the index of the first row of a <tt>binary_matrix_t</tt> which has
    * a one on its <tt>col</tt>-th column. It returns NO_SUCH_ROW if no such row
    * is found.
    *
    * \param[in] col The column of the matrix.
    * \param[in] matrix A pointer to the <tt>binary_matrix_t</tt>.
    * \return
    *  \li An unsigned integer \c row between 0 and <tt>matrix->nrows</tt>-1
    *      such that <tt>(matrix->data[row][col] == 1)</tt>
    *  \li NO_SUCH_ROW if, for all valid \c i, <tt>(matrix->data[i][col] !=
    *      1)</tt>.
    *
    * \note This function is needed in the gaussian elimination algorithm
    * described in the paper "A compact algorithm for Gaussian elimination
    * over GF(2) implemented on highly parallel computers", by D. Parkinson
    * and M. Wunderlich (Parallel Computing 1 (1984) 65-73).
    *
    * It could be argued that such a function should then be declared and
    * implemented in the files relevant to the aforementionned algorithm.
    * However, this would lead to a very inefficient implementation of this
    * function since proprer programming pratices would lead to consider the
    * matrix as some kind of opaque structure. Granted, nothing could have
    * prevented us to implement \c first_row_with_one_on_col exactly as in
    * \c matrix.c, but the future maintainer would have the burden to check
    * and update code scattered around several files should the inner structure
    * of \c binary_matrix_t be modified.
    *
    * This can be seen as a moot point: after all, the TIFA code does not
    * strictly enforce type encapsulation. Indeed, some parts of the code
    * do assume (a minimal!) knowledge of the internal structures of some
    * types (have a look at siqs.c for instance). That does not make it the
    * right thing to do though. Unless when facing a real bottleneck, let's try
    * to keep the "programmer's omniscience" to a manageable level...
    */
uint32_t first_row_with_one_on_col(uint32_t col,
                                   const binary_matrix_t* const matrix);

/*
 *-----------------------------------------------------------------------------
 *              byte_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_byte_matrix_t matrix.h lib/utils/include/matrix.h
    * \brief  Defines a matrix of bytes.
    *
    * This structure defines a matrix of bytes which knows its current
    * dimensions and its allocated memory space.
    *
    * \note Internally, a matrix of bytes is represented as a matrix of
    * \c unsigned char elements.
    */
struct struct_byte_matrix_t {
       /**
        * Maximum number of rows of the matrix.
        */
    uint32_t nrows_alloced;
       /**
        * Maximum number of columns of the matrix.
        */
    uint32_t ncols_alloced;
       /**
        * Current number of rows of the matrix.
        */
    uint32_t nrows;
       /**
        * Current number of columns of the matrix.
        */
    uint32_t ncols;
       /**
        * 2D array of \c unsigned char elements whose dimensions are given
        * by the \c nrows_alloced and \c ncols_alloced fields.
        */
    unsigned char** data;
};

   /**
    * \typedef byte_matrix_t
    * \brief Equivalent to <tt>struct struct_byte_matrix_t</tt>.
    */
typedef struct struct_byte_matrix_t byte_matrix_t;

   /**
    * \brief Allocates and returns a new <tt>byte_matrix_t</tt>.
    *
    * Allocates and returns a new <tt>byte_matrix_t</tt> such that:
    * \li its \c nrows_alloced field is set to nrows.
    * \li its \c ncols_alloced field is set to ncols.
    * \li its \c nrows field is set to nrows.
    * \li its \c ncols field is set to ncols.
    * \li its \c data array of arrays is completely filled with zeroes.
    *
    * \note The behaviour of this alloc function differs from the ones in
    * array.h. This discrepancy will be corrected in later versions.
    *
    * \param[in] nrows The maximum number of rows of the \c byte_matrix_t
    *                  to allocate.
    * \param[in] ncols The maximum number of columns \c byte_matrix_t
    *                  to allocate.
    * \return A pointer to the newly allocated \c byte_matrix_t structure.
    */
byte_matrix_t* alloc_byte_matrix(uint32_t nrows, uint32_t ncols);

   /**
    * \brief Allocates and returns a cloned <tt>byte_matrix_t</tt>.
    *
    * Allocates and returns a clone of the <tt>byte_matrix_t</tt> 
    * pointed by <tt>matrix</tt>.
    *
    * \param[in] matrix A pointer to the byte matrix to clone.
    * \return A pointer to the newly allocated \c byte_matrix_t clone.
    */
byte_matrix_t* clone_byte_matrix(const byte_matrix_t * const matrix);

   /**
    * \brief Sets a <tt>byte_matrix_t</tt> to zero.
    *
    * Sets the <tt>byte_matrix_t matrix</tt> to the zero matrix.
    *
    * \param[in] matrix A pointer to the <tt>byte_matrix_t</tt> to reset.
    */
void reset_byte_matrix(byte_matrix_t* const matrix);

   /**
    * \brief Clears a <tt>byte_matrix_t</tt>.
    *
    * Clears a <tt>byte_matrix_t</tt>, or, more precisely, clears the memory
    * space used by the arrays pointed by the \c data field of a
    * <tt>byte_matrix_t</tt>. Also set its \c nrows_alloced,
    * \c ncols_alloced, \c nrows and \c ncols fields to zero.
    *
    * \param[in] matrix A pointer to the <tt>byte_matrix_t</tt> to clear.
    */
void clear_byte_matrix(byte_matrix_t* const matrix);

   /**
    * \brief Prints a <tt>byte_matrix_t</tt>.
    *
    * Prints a <tt>byte_matrix_t</tt>'s on the standard output.
    *
    * \param[in] matrix A pointer to the <tt>byte_matrix_t</tt> to print.
    */
void print_byte_matrix(const byte_matrix_t* const matrix);

#ifdef __cplusplus
}
#endif

#endif
