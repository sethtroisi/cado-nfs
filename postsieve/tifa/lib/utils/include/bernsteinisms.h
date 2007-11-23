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
 * \file    bernsteinisms.h
 * \author  Jerome Milan
 * \date    Fri Oct 12 2007
 * \version 1.1
 *
 * \brief Algorithms from two D. J. Bernstein's papers on the factorization
 * of small integers.
 *
 * Algorithms from two D. J. Bernstein's papers on the factorization
 * of small integers:
 * \li "How to find small factors of integers", http://cr.yp.to/papers/sf.pdf
 * \li "How to find smooth parts of integers",
 *   http://cr.yp.to/factorization/smoothparts-20040510.pdf
 */

 /*
  *  History:
  *  
  *  1.1:   Fri Oct 12 2007 by JM:
  *         - Added multi-step early abort strategy (see smooth_filter.h)
  *
  *  1.0.2: circa March 2007 by JM:
  *         - Added function bern_21_rt_pairs_siqs.
  *         - Added function bern_21_rt_pairs_lp_siqs.
  *
  *  1.0.1: Mon May 22 2006 by JM:
  *         - Added function bern_21_pairs.
  *         - Added function bern_21_pairs_lp.
  *
  *  1.0.0: Fri Mar 3 2006 by JM:
  *         - Initial version.
  */

#if !defined(_TIFA_BERNSTEINISMS_H_)
   /**
    * \def _TIFA_BERNSTEINISMS_H_
    * Standard include guard.
    */
#define _TIFA_BERNSTEINISMS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include "array.h"
#include "x_array_list.h"
#include "hashtable.h"
#include "smooth_filter.h"

//-----------------------------------------------------------------------------
// Ref. "How to find small factors of integers", Daniel J. Bernstein
//       http://cr.yp.to/papers/sf.pdf
//-----------------------------------------------------------------------------

   /**
    * \brief Daniel J. Bernstein's algorithm 5.1.
    *
    * Given a positive integer \c b and an odd positive integer \c u, returns
    * a non negative integer \c v < 2^\c b such that
    * 1 + <tt>u*v</tt> = 0 (mod 2^\c b).
    *
    * This is the algorithm 5.1 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    * http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] b A positive integer.
    * \param[in] u An odd positive \c mpz_t integer.
    * \return A non negative \c mpz_t integer \c v < 2^\c b such that
    *         1 + <tt>u*v</tt> = 0 (mod 2^\c b).
    */
mpz_t* bern_51(uint32_t b, const mpz_t u);

   /**
    * \brief Daniel J. Bernstein's algorithm 5.3.
    *
    * Given an odd positive integer \c u < 2^\c c and a non negative integer
    * \c x < 2^(\c b + \c c), returns a non negative integer
    * \c r < 2^(\c c + 1) such that \c r*2^\c b = \c x (mod \c u).
    *
    * This is the algorithm 5.3 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    *  http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] b A positive integer.
    * \param[in] u An odd positive \c mpz_t integer.
    * \param[in] x An non negative \c mpz_t integer.
    * \return A non negative \c mpz_t integer \c r such that
    *         \c r*2^\c b = \c x (mod \c u).
    */
mpz_t* bern_53(uint32_t b, const mpz_t u, const mpz_t x);

   /**
    * \brief Daniel J. Bernstein's algorithm 6.3.
    *
    * Given a non negative integer \c x and given the product tree \c tree
    * of a sequence of odd positive integers p_i, returns the integers
    * p_i such that: \c x mod p_i = 0.
    *
    * This is the algorithm 6.3 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    *  http://cr.yp.to/papers/sf.pdf
    *
    * \param[in] x A non negative positive integer.
    * \param[in] tree The product tree of a sequence of odd positive
    *                 integers p_i.
    * \return A pointer to an \c uint32_array_t holding the integers
    *         p_i such that: \c x mod p_i = 0.
    */
uint32_array_t* bern_63(const mpz_t x, const mpz_array_t* const tree);

   /**
    * \brief Daniel J. Bernstein's algorithm 7.1.
    *
    * Given a sequence of odd primes p_j given by \c odd_primes and a
    * set of integers n_i given by \c to_be_factored, determines, for each
    * n_i, the list of odd primes p_j such that (n_i mod p_j = 0) and stores
    * them in <tt>decomp_list</tt>. Each entry in \c decomp_list->data[i]
    * is a pointer to a \c mpz_array_t listing the p_j for the integer
    * <tt>to_be_factored->data[i]</tt>.
    *
    * This is the algorithm 7.1 described in Daniel J. Bernstein's paper: "How
    * to find small factors of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find small factors of integers",
    * http://cr.yp.to/papers/sf.pdf
    *
    * \param[out] decomp_list A pointer to the list of matching p_j for
    *                         each n_i.
    * \param[in] to_be_factored A pointer to the set of integers n_i.
    * \param[in] odd_primes A pointer to the set of integers p_j.
    */
void bern_71(
    uint32_array_list_t* const decomp_list,
    const mpz_array_t* const to_be_factored,
    const uint32_array_t* const odd_primes
);

//-----------------------------------------------------------------------------
// Ref. "How to find smooth parts of integers", Daniel J. Bernstein
//       http://cr.yp.to/factorization/smoothparts-20040510.pdf
//-----------------------------------------------------------------------------

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 (with computation of a
    * remainder tree).
    *
    * Given the prime numbers p_j listed by \c pj and the positive integers
    * x_i listed by <tt>xi</tt>, determines the {p_j}-smooth part of each
    * x_i and stores them in <tt>smooth</tt>, so that \c smooth->data[i] is
    * the {p_j}-smooth part of \c xi->data[i].
    *
    * The function stops when each integer from \c xi has been checked for
    * smoothness or when the \c smooth array is completely filled. It then
    * returns the index of the last integer in \c xi that has been checked
    * for smoothness.
    *
    * This is the algorithm 2.1 described in Daniel J. Bernstein's paper: "How
    * to find smooth parts of integers".
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] smooth A pointer to the {p_j}-smooth parts of each
    *                    x_i integer.
    * \param[in] xi A pointer to the list of the x_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c xi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt(
    mpz_array_t* const smooth,
    const mpz_array_t* const xi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 (without computation of a
    * remainder tree).
    *
    * Given the prime numbers p_j listed by \c pj and the positive integers
    * x_i listed by <tt>xi</tt>, determines the {p_j}-smooth part of each
    * x_i and stores them in <tt>smooth</tt>, so that \c smooth->data[i] is
    * the {p_j}-smooth part of \c xi->data[i].
    *
    * The function stops when each integer from \c xi has been checked for
    * smoothness or when the \c smooth array is completely filled. It then
    * returns the index of the last integer in \c xi that has been checked
    * for smoothness.
    *
    * This is the algorithm 2.1 described in Daniel J. Bernstein's paper: "How
    * to find smooth parts of integers".
    *
    * \note This function differs from \c bern_21_rt only because no
    * remainder tree is computed. This can sometimes be faster than the full
    * fledged version.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] smooth A pointer to the {p_j}-smooth parts of each
    *                    x_i integer.
    * \param[in] xi A pointer to the list of the x_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c xi that has been checked
    *          for smoothness.
    */
uint32_t bern_21(
    mpz_array_t* const smooth,
    const mpz_array_t* const xi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified.
    *
    * Given the prime numbers p_j listed by \c pj and the positive integers
    * y_i listed by <tt>cand_yi</tt>, determines the y_i that are {p_j}-smooth
    * and stores them in <tt>smooth_yi</tt>, so that \c smooth_yi->data[i]
    * is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] xi A pointer to the {x_i} associated to the {p_j}-smooth
    *                y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth y_i integer.
    * \param[in] cand_xi A pointer to the list of the x_i integers.
    * \param[in] cand_yi A pointer to the list of the y_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs(
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified (without computation
    * of a remainder tree).
    *
    * Given the prime numbers p_j listed by \c pj and the positive integers
    * y_i listed by <tt>cand_yi</tt>, determines the y_i that are {p_j}-smooth
    * and stores them in <tt>smooth_yi</tt>, so that \c smooth_yi->data[i]
    * is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem.
    *
    * \note This function is very similar to \c bern_21_rt_pairs. The
    *       only difference is that in \c bern_21_pairs no remainder tree is
    *       computed. This can sometimes be faster than the full fledged
    *       version.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] xi A pointer to the {x_i} associated to the {p_j}-smooth
    *                y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth y_i integer.
    * \param[in] cand_xi A pointer to the list of the x_i integers.
    * \param[in] cand_yi A pointer to the list of the y_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_pairs(
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified, with large
    * primes variation.
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * Moreover, this function implements the so-called large primes variation.
    * If a given y_i is not {p_j}-smooth but is the product of a prime \c p by
    * a {p_j}-smooth number, it is stored in the hashtable <tt>htable</tt>.
    * Subsequently, if another y_j is the product of a {p_j}-smooth number
    * by the same prime number \c p, then <tt>y_i*y_j/(p^2)</tt> is stored in
    * \c smooth_yi and <tt>x_i*x_j*pinv</tt> is stored in \c xi where
    * \c pinv is the inverse of \c p in Z/<tt>c</tt>Z.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[in] n The integer to factor.
    * \param[in] htable A pointer to the hashtable used for the large prime
    *                   variation.
    * \param[out] xi A pointer to the {x_i} associated to the {p_j}-smooth
    *                y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth y_i integer.
    * \param[in] cand_xi A pointer to the list of the x_i integers.
    * \param[in] cand_yi A pointer to the list of the y_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs_lp(
    const mpz_t n,
    hashtable_t* const htable,
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified, with large
    * primes variation (without computation of a remainder tree).
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * Moreover, this function implements the so-called large primes variation.
    * If a given y_i is not {p_j}-smooth but is the product of a prime \c p by
    * a {p_j}-smooth number, it is stored in the hashtable <tt>htable</tt>.
    * Subsequently, if another y_j is the product of a {p_j}-smooth number
    * by the same prime number \c p, then <tt>y_i*y_j/(p^2)</tt> is stored in
    * \c smooth_yi and <tt>x_i*x_j*pinv</tt> is stored in \c xi where
    * \c pinv is the inverse of \c p in Z/<tt>c</tt>Z.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem.
    *
    * \note This function is very similar to \c bern_21_rt_pairs_lp. The only
    *       difference is that here no remainder tree is computed. This can
    *       sometimes be faster than the full fledged version.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[in] n The integer to factor.
    * \param[in] htable A pointer to the hashtable used for the large prime
    *                   variation.
    * \param[out] xi A pointer to the {x_i} associated to the {p_j}-smooth
    *                y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth y_i integer.
    * \param[in] cand_xi A pointer to the list of the x_i integers.
    * \param[in] cand_yi A pointer to the list of the y_i integers.
    * \param[in] z The product of the the p_j prime numbers in the facotr base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_pairs_lp(
    const mpz_t n,
    hashtable_t* const htable,
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified for SIQS (with
    * computation of a remainder tree).
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem, particularly
    * to our SIQS implementation where we need to keep track of additionnal
    * \c a_i integers associated to each \c y_i integers. These extra integers
    * are stored in \c cand_a whereas the a_i associated to the smooth y_i will
    * be stored in the \c a_for_smooth_gx array.
    *
    * \note The \c a_i integers are actually the values of the first parameter
    *       of the polynomials used in the SIQS algorithm.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[out] xi A pointer to the \c x_i associated to the {p_j}-smooth
    *                \c y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth \c y_i integer.
    * \param[out] a_for_smooth_yi A pointer to the array of the
    * \param[in] cand_xi A pointer to the list of the \c x_i integers.
    * \param[in] cand_yi A pointer to the list of the \c y_i integers.
    * \param[in] cand_a A pointer to the list of the \c a_i integers associated
    *                   to the {p_j}-smooth y_i integers.
    * \param[in] z The product of the the \c p_j prime numbers in the factor
    *              base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs_siqs(
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    mpz_array_t* const a_for_smooth_gx,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_array_t* const cand_a,
    const mpz_t z
);

   /**
    * \brief Daniel J. Bernstein's algorithm 2.1 modified for SIQS, with large
    * primes variation (with computation of a remainder tree).
    *
    * Given <tt>z</tt>, the product of prime numbers p_j and the positive
    * integers y_i listed by <tt>cand_yi</tt>, determines the y_i that are
    * {p_j}-smooth and stores them in <tt>smooth_yi</tt>, so that
    * \c smooth_yi->data[i] is indeed {p_j}-smooth.
    *
    * In a typical factorization problem, we other found ourselves in situations
    * where each y_i is associated to another integer x_i. The x_i associated
    * to the {p_j}-smooth y_i are hence stored in <tt>xi</tt>.
    *
    * Moreover, this function implements the so-called large primes variation.
    * If a given y_i is not {p_j}-smooth but is the product of a prime \c p by
    * a {p_j}-smooth number, it is stored in the hashtable <tt>htable</tt>.
    * Subsequently, if another y_j is the product of a {p_j}-smooth number
    * by the same prime number \c p, then <tt>y_i*y_j/(p^2)</tt> is stored in
    * \c smooth_yi and <tt>x_i*x_j*pinv</tt> is stored in \c xi where
    * \c pinv is the inverse of \c p in Z/<tt>c</tt>Z.
    *
    * The function stops when each integer from \c cand_yi has been checked for
    * smoothness or when the \c smooth_yi array is completely filled. It then
    * returns the index of the last integer in \c cand_yi that has been checked
    * for smoothness.
    *
    * This function uses the algorithm 2.1 described in Daniel J. Bernstein's
    * paper: "How to find smooth parts of integers" except that this function
    * has been tailored to better suit the factorization problem, particularly
    * to our SIQS implementation where we need to keep track of additionnal
    * \c a_i integers associated to each \c y_i integers. These extra integers
    * are stored in \c cand_a whereas the a_i associated to the smooth y_i will
    * be stored in the \c a_for_smooth_gx array.
    *
    * \note The \c a_i integers are actually the values of the first parameter
    *       of the polynomials used in the SIQS algorithm.
    *
    * \see Daniel J. Bernstein's paper: "How to find smooth parts of integers",
    * http://cr.yp.to/factorization/smoothparts-20040510.pdf
    *
    * \param[in] n The integer to factor.
    * \param[in] htable A pointer to the hashtable used for the large prime
    *                   variation.
    * \param[out] xi A pointer to the \c x_i associated to the {p_j}-smooth
    *                \c y_i integer.
    * \param[out] smooth_yi A pointer to the {p_j}-smooth \c y_i integer.
    * \param[out] a_for_smooth_yi A pointer to the array of the
    * \param[in] cand_xi A pointer to the list of the \c x_i integers.
    * \param[in] cand_yi A pointer to the list of the \c y_i integers.
    * \param[in] cand_a A pointer to the list of the \c a_i integers associated
    *                   to the {p_j}-smooth y_i integers.
    * \param[in] z The product of the the \c p_j prime numbers in the factor
    *              base.
    * \return  The index of the last integer in \c cand_yi that has been checked
    *          for smoothness.
    */
uint32_t bern_21_rt_pairs_lp_siqs(
    const mpz_t n,
    hashtable_t* const htable,
    mpz_array_t* const xi,
    mpz_array_t* const smooth_yi,
    mpz_array_t* const a_for_smooth_yi,
    const mpz_array_t* const cand_xi,
    const mpz_array_t* const cand_yi,
    const mpz_array_t* const cand_a,
    const mpz_t z
);

  /**
    * \brief Daniel J. Bernstein's algorithm 2.1 adapted to be used with a
    * \c smooth_filter_t.
    *
    * <h3>If \c filter->nsteps == 0 </h3>
    *
    * In such a case, no early abort strategy is performed. The effect of
    * the function is the same as \c bern_21_*_pairs_* called with:
    * \verbatim
    filter->n
    filter->htable
    filter->accepted_xi
    filter->accepted_yi
    filter->accepted_ai
    filter->candidate_xi
    filter->candidate_yi
    filter->candidate_ai
    filter->prod_pj[0]
    \endverbatim
    * 
    * <h3>If \c filter->nsteps != 0 </h3>
    *
    * An early abort strategy is performed.
    *
    * \li If 1 <= \c step < \c filter->nsteps:<br />
    * Relations at step \c step-1 from \c filter,  
    * (\c filter->filtered_*[\c step-1]) are used as "candidate" arrays to
    * populate either \c filter->accepted_* or \c filter->filtered_*[step].
    * 
    * \li If \c step == 0:<br />
    * The candidate relations are taken from \c filter->candidate_*.
    *
    * \li If \c step == \c filter->nsteps:<br />
    * The candidate relations are taken from
    * \c filter->filtered_*[\c filter->nsteps - 1] and 'good' relations
    * will be stored in \c filter->accepted_*.
    *
    * \see The bern_21_* functions.
    *
    * \warning Using \c filter->nsteps != 0 is not recommended. First, it
    * certainly does not make any sense to try to early-abort the batch.
    * Second, even if it is useful (for some weird reasons that I'm not
    * aware of), the cases \c filter->nsteps != 0 have not been tuned / fully
    * debugged.
    *
    * \param filter pointer to the \c smooth_filter_t to use
    * \param step the step number in the early abort strategy
    * \return The number of relations used from the "candidate" arrays.
    */
inline uint32_t djb_batch_rt(
    smooth_filter_t* const filter,
    unsigned long int step
);

#ifdef __cplusplus
}
#endif

#endif
