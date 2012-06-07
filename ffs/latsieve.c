#include <stdlib.h>

#include "latsieve.h"
#include "types.h"
#include "ijvec.h"
#include "sublat.h"
#include "gray.h"


/*
 * Given j is a multiple of p. 
 * Compute the next multiple of p, in lex order. 
 * The output is of degree less than J. If failure because this bound is
 * reached, return 0, otherwise return 1.
 *
 * For that, make us of the given basis of multiples of p, as precomputed
 * in the normalized_echelon_multiples() function.
 */
static MAYBE_UNUSED
int next_projective_j(ij_t rj, ij_t j, ij_t *basis, int degp, int J)
{
    // First use monic_set_next() on the high part of j.
    ij_t jhi, njhi;
    ij_div_ti(jhi, j, degp);
    int rc = ij_monic_set_next(njhi, jhi, J-degp);
    if (!rc)
        return 0;
    // The degree of the difference between in and out of set_next gives
    // the basis-element to add. (magic!)
    ij_diff(njhi, njhi, jhi);
    int d = ij_deg(njhi);
    ij_add(rj, j, basis[d]);
    // There is an adjustment to do in the case where the degree of rj is
    // larger than the degree of j, due to the fact that we deal with
    // monic polynomials.
    if (d >= ij_deg(jhi)) {
        ASSERT(d > ij_deg(jhi)); // monic case
        // Have to kill the bit d-1, which is currently 2, without
        // touching lower bits.
        for (int k = 2; k < FP_CHAR; ++k) {
            if (d > 0)
                ij_add(rj, rj, basis[d-1]);
            if (d > 1)
                ij_add(rj, rj, basis[d-2]);
        }
    }
    return 1;
}


static inline
void sieve_hit(uint8_t *S, uint8_t degp, ijpos_t pos,
        MAYBE_UNUSED fbprime_srcptr p,
        MAYBE_UNUSED fbprime_srcptr r,
        MAYBE_UNUSED ijpos_t pos0)
{
#ifdef TRACE_POS
  if (pos0+pos == TRACE_POS) {
    fprintf(stderr, "TRACE_POS(%lu): ", pos0+pos);
    fbprime_out(stderr, p); fprintf(stderr, " ");
    fbprime_out(stderr, r); fprintf(stderr, "\n");
    fprintf(stderr, "TRACE_POS(%lu): degnorm is now %d\n",
            pos0+pos, S[pos]-degp);
  }
#endif
#ifndef NDEBUG
  if (S[pos] < degp)
    fprintf(stderr, "faulty pos is %lu\n", pos0+pos);
  ASSERT(S[pos] >= degp);
#endif
  S[pos] -= degp;
}


// TODO: This code is completely under-efficient. Fix this!
void sieveFB(uint8_t *S, factor_base_srcptr FB, unsigned I, unsigned J,
             ij_t j0, ijpos_t pos0, ijpos_t size, sublat_ptr sublat)
{
  ijbasis_t basis;
  ijbasis_init(basis, I, J);

  ij_t *tmp = (ij_t *)malloc(J * sizeof(ij_t));
  ASSERT_ALWAYS(tmp != NULL);

  fp_t one, two;
  fp_set_one(one);
  fp_add(two, one, one);

    for (unsigned int ii = 0; ii < FB->n; ++ii) {
        fbideal_ptr gothp = FB->elts[ii];
        int L = gothp->degp;
        if (UNLIKELY(gothp->power))
            L = fbprime_deg(gothp->p);
        // Larger primes are left to the bucket sieve.
        if ((unsigned)L >= I) break;
        // List of cases that are not handled yet:
        if (use_sublat(sublat) && L == 1) continue;

        // projective roots are handled differently
        if (gothp->proj) {
            // FIXME: For the time being, projective-root sieving is disabled.
            // See above for an explanation why.
            // sieve_projective_root(S, gothp, I, J, pos0, size, sublat);
            continue;
        }


        ijbasis_compute(NULL, basis, gothp);

        // Only the first time round.
        if (UNLIKELY(!pos0)) {
          if (use_sublat(sublat)) {
            // In the case of sublattices, compute the starting point for the
            // sieve by gothp for the current sublattice.
            // TODO: Way too expensive!
            // xi and yi have degree 2
            ij_t i0, xi, xip, yip;
            ij_set_16(xi, sublat->lat[sublat->n][0]);
            fbprime_t tmp0;
            ij_t tmp1;
            fbprime_set_16(tmp0, sublat->lat[sublat->n][1]);
            fbprime_mulmod(tmp0, gothp->lambda, tmp0, gothp->p);
            ij_set_fbprime(yip, tmp0);
            ij_add(xip, xi, yip);
            ij_set_16(tmp1, sublat->modulus);
            ij_mulmod(xip, xip, gothp->tildep, tmp1);
            ij_set_fbprime(tmp1, gothp->p);
            ij_mul(i0, xip, tmp1);
            ij_add(i0, i0, yip);

            ij_t ijmod;
            ij_set_16(ijmod, sublat->modulus);
            ij_sub(i0, i0, xi);
            ij_div(gothp->i0, i0, ijmod);
          } else {
            ij_set_zero(gothp->i0);
          }
        }

        for (unsigned i = 0; i < J; ++i)
          ij_smul(tmp[i], basis->v[I-L+i]->i, two);
        for (unsigned i = I-L+1; i < basis->dim; ++i)
          ijvec_add(basis->v[i], basis->v[i], basis->v[i-1]);

        ij_t i, j, jj;
        int rcj = 1, degj, degjj;
        for (ij_set(j, j0), degj = MAX(ij_deg(j), 0); rcj; ) {
          ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
          if (start >= size)
            break;
          ij_set(i, gothp->i0);
          ijpos_t pos = start + ijvec_get_offset(i, I);
          if (LIKELY(pos0 || pos))
            sieve_hit(S, gothp->degp, pos, gothp->p, gothp->r, pos0);

          // Size of Gray codes to use for the inner loop.
#         if   defined(USE_F2)
#           define ENUM_LATTICE_UNROLL 8
#         elif defined(USE_F3)
#           define ENUM_LATTICE_UNROLL 5
#         endif

          // Unrolled p-ary Gray code of size ENUM_LATTICE_UNROLL.
          static const uint8_t gray[] = { GRAY(ENUM_LATTICE_UNROLL) };
          unsigned             ngray  = GRAY_LENGTH(ENUM_LATTICE_UNROLL);

          // Just in case the dimension of the vector space is lower than
          // ENUM_LATTICE_UNROLL.
          unsigned gray_dim = MIN(I-L, ENUM_LATTICE_UNROLL);
          unsigned k0       = ngray - GRAY_LENGTH(gray_dim);

          ij_t s, t;
          ij_set_zero(t);
          int rc = I-L > ENUM_LATTICE_UNROLL;
          do {
            // Inner-level Gray code enumeration: just go through the Gray code
            // array, each time adding the indicated basis vector.
            for (unsigned k = k0; k < ngray; ++k) {
              ij_add(i, i, basis->v[gray[k]]->i);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, gothp->degp, pos, gothp->p, gothp->r, pos0);
            }
            k0 = 0;

            // Outer-level Gray code enumeration: using ij_set_next, the degree
            // of the difference with the previous one indicates which basis
            // vector should be added to the current lattice point.
            // rc is set to 0 when all vectors have been enumerated.
            ij_set(s, t);
            rc = rc && ij_set_next(t, t, I-L-ENUM_LATTICE_UNROLL);
            if (rc) {
              ij_diff(s, s, t);
              ij_add(i, i, basis->v[ij_deg(s)+ENUM_LATTICE_UNROLL]->i);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, gothp->degp, pos, gothp->p, gothp->r, pos0);
            }
          } while (rc);

          ij_set(jj, j);
          rcj = ij_monic_set_next(j, j, J);
          if (rcj) {
            ij_diff(jj, jj, j);
            degjj = ij_deg(jj);
            ij_add(gothp->i0, gothp->i0, basis->v[I-L+degjj]->i);
            if (degjj > degj) {
              ij_sub(gothp->i0, gothp->i0, tmp[degjj-1]);
              degj = degjj;
            }
          }
        }
    }

    free(tmp);
    ijbasis_clear(basis);
}


// TODO: This code is completely under-efficient. Fix this!
void sieveFB2(uint8_t *S, small_factor_base_ptr FB, unsigned I, unsigned J,
             ij_t j0, ijpos_t pos0, ijpos_t size, sublat_ptr sublat)
{
    for (unsigned int ii = 0; ii < FB->n; ++ii) {
        small_fbideal_ptr gothp = FB->elts[ii];
        int L = gothp->degq;

        // Larger primes are left to the bucket sieve.
        ASSERT((unsigned)L < I);

        // List of cases that are not handled yet:
        if (use_sublat(sublat) && L == 1) continue;

        // projective roots are handled differently
        if (gothp->proj) {
          // First time round?
          if (UNLIKELY(!pos0)) {
            // Find the first line to fill
            ij_set_zero(gothp->current);
            if (use_sublat(sublat)) {
              ASSERT_ALWAYS(0);  // XXX Sublat is broken
            }
          }
          ij_t i, j;
          int rcj = 1;
          ij_set(j, gothp->current);
          while (rcj) {
            ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
            if (start >= size)
              break;

            // Sieve the whole line
            int rci = 1;
            for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
              ijpos_t pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, gothp->degp, pos, gothp->q, gothp->r, pos0);
            }

            rcj = next_projective_j(j, j, gothp->projective_basis, L, J);
          }
          ij_set(gothp->current, j); // remember the next line to sieve.
          continue;
        }

        // Only the first time round.
        if (UNLIKELY(!pos0)) {
          if (use_sublat(sublat)) {
#if 0
            // FIXME: sublat are Broken XXX
            
            // In the case of sublattices, compute the starting point for the
            // sieve by gothp for the current sublattice.
            // TODO: Way too expensive!
            // xi and yi have degree 2
            ij_t i0, xi, xip, yip;
            ij_set_16(xi, sublat->lat[sublat->n][0]);
            fbprime_t tmp0;
            ij_t tmp1;
            fbprime_set_16(tmp0, sublat->lat[sublat->n][1]);
            fbprime_mulmod(tmp0, gothp->lambda, tmp0, gothp->p);
            ij_set_fbprime(yip, tmp0);
            ij_add(xip, xi, yip);
            ij_set_16(tmp1, sublat->modulus);
            ij_mulmod(xip, xip, gothp->tildep, tmp1);
            ij_set_fbprime(tmp1, gothp->p);
            ij_mul(i0, xip, tmp1);
            ij_add(i0, i0, yip);

            ij_t ijmod;
            ij_set_16(ijmod, sublat->modulus);
            ij_sub(i0, i0, xi);
            ij_div(gothp->current, i0, ijmod);
#endif
            ij_set_zero(gothp->current);
          } else {
            ij_set_zero(gothp->current);
          }
        }

        ij_t i, j, jj;
        int rcj = 1, degj, degjj;
        for (ij_set(j, j0), degj = MAX(ij_deg(j), 0); rcj; ) {
          ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
          if (start >= size)
            break;
          ij_set(i, gothp->current);
          ijpos_t pos = start + ijvec_get_offset(i, I);
          if (LIKELY(pos0 || pos))
            sieve_hit(S, gothp->degp, pos, gothp->q, gothp->r, pos0);

          // Size of Gray codes to use for the inner loop.
#         if   defined(USE_F2)
#           define ENUM_LATTICE_UNROLL 8
#         elif defined(USE_F3)
#           define ENUM_LATTICE_UNROLL 5
#         endif

          // Unrolled p-ary Gray code of size ENUM_LATTICE_UNROLL.
          static const uint8_t gray[] = { GRAY(ENUM_LATTICE_UNROLL) };
          unsigned             ngray  = GRAY_LENGTH(ENUM_LATTICE_UNROLL);

          // Just in case the dimension of the vector space is lower than
          // ENUM_LATTICE_UNROLL.
          unsigned gray_dim = MIN(I-L, ENUM_LATTICE_UNROLL);
          unsigned k0       = ngray - GRAY_LENGTH(gray_dim);

          ij_t s, t;
          ij_set_zero(t);
          int rc = I-L > ENUM_LATTICE_UNROLL;
          do {
            // Inner-level Gray code enumeration: just go through the Gray code
            // array, each time adding the indicated basis vector.
            for (unsigned k = k0; k < ngray; ++k) {
              ij_add(i, i, gothp->basis->v[gray[k]]->i);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, gothp->degp, pos, gothp->q, gothp->r, pos0);
            }
            k0 = 0;

            // Outer-level Gray code enumeration: using ij_set_next, the degree
            // of the difference with the previous one indicates which basis
            // vector should be added to the current lattice point.
            // rc is set to 0 when all vectors have been enumerated.
            ij_set(s, t);
            rc = rc && ij_set_next(t, t, I-L-ENUM_LATTICE_UNROLL);
            if (rc) {
              ij_diff(s, s, t);
              ij_add(i, i, gothp->basis->v[ij_deg(s)+ENUM_LATTICE_UNROLL]->i);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, gothp->degp, pos, gothp->q, gothp->r, pos0);
            }
          } while (rc);

          ij_set(jj, j);
          rcj = ij_monic_set_next(j, j, J);
          if (rcj) {
            ij_diff(jj, jj, j);
            degjj = ij_deg(jj);
            ij_add(gothp->current, gothp->current,
                    gothp->basis->v[I-L+degjj]->i);
            if (degjj > degj) {
              ij_sub(gothp->current, gothp->current,
                      gothp->adjustment_basis->v[degjj-1]->i);
              degj = degjj;
            }
          }
        }
    }
}
