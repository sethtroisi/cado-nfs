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
 *
 * NB: in the case of sublat, the given j is not really a multiple of p.
 * The computation is still valid.
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


void sieveSFB(uint8_t *S, unsigned int *thr,
    small_factor_base_ptr FB, unsigned I, unsigned J,
    ij_t j0, ijpos_t pos0, ijpos_t size, sublat_ptr sublat)
{
    *thr = 0;
    for (unsigned int ii = 0; ii < FB->n; ++ii) {
        small_fbideal_ptr gothp = FB->elts[ii];
        int L = gothp->degq;
        int degp = gothp->degp;

        // Larger primes are left to the bucket sieve.
        ASSERT((unsigned)L < I);

        // In case of sublat, the primes of degree 1 gives a uniform
        // contribution and it is better to handle them globally using
        // thresholds.
        // TODO: we recompute this for each bucket region, whereas it
        // could be computed once and for all.
        if (use_sublat(sublat) && L == 1) {
          if (!gothp->proj) {
            fppol16_t qq, rr;
            fppol16_set_fbprime(qq, gothp->q);
            fppol16_set_fbprime(rr, gothp->lambda);
            fppol16_mul(rr, rr, sublat->lat[sublat->n][1]);
            fppol16_add(rr, rr, sublat->lat[sublat->n][0]);
            fppol16_rem(rr, rr, qq);
            if (fppol16_is_zero(rr))
              *thr += gothp->degp;
          } else {
            fppol16_t qq;
            fppol16_set_fbprime(qq, gothp->q);
            fppol16_rem(qq, sublat->lat[sublat->n][1], qq);
            if (fppol16_is_zero(qq))
              *thr += gothp->degp;
          }
          continue;
        }

        // projective roots are handled differently
        if (gothp->proj) {
          // First time round?
          if (UNLIKELY(!pos0)) {
            // Find the first line to fill. If no sublat, this is zero.
            // Otherwise, there is a bit of computation.
            ij_set_zero(gothp->current);
            if (use_sublat(sublat)) {
              ij_t tmp0, tmp1, ijmod;
              ij_set_16(ijmod, sublat->modulus);
              ij_set_16(tmp0, sublat->lat[sublat->n][1]);
              ij_mulmod(tmp0, gothp->tildep, tmp0, ijmod);
              ij_set_fbprime(tmp1, gothp->q);
              ij_mul(tmp1, tmp1, tmp0);
              ij_div(gothp->current, tmp1, ijmod);
            }
          }
          ij_t j;
          int rcj = 1;
          ij_set(j, gothp->current);
          while (rcj) {
            ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
            if (start >= size)
              break;

            // Sieve the whole line
#ifndef USE_F2
            // TODO: this should be a big memsub().
            ij_t i;
            int rci = 1;
            for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
              ijpos_t pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, degp, pos, gothp->q, gothp->r, pos0);
            }
#else
            // For GF(2), this becomes so simple (and Gcc does it well).
            for(int i=0; i < 1<<I; ++i)
              S[start+i] -= degp;
#endif

            rcj = next_projective_j(j, j, gothp->projective_basis, L, J);
          }
          ij_set(gothp->current, j); // remember the next line to sieve.
          continue;
        }

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
            fbprime_mulmod(tmp0, gothp->lambda, tmp0, gothp->q);
            ij_set_fbprime(yip, tmp0);
            ij_add(xip, xi, yip);
            ij_set_16(tmp1, sublat->modulus);
            ij_mulmod(xip, xip, gothp->tildep, tmp1);
            ij_set_fbprime(tmp1, gothp->q);
            ij_mul(i0, xip, tmp1);
            ij_add(i0, i0, yip);

            ij_t ijmod;
            ij_set_16(ijmod, sublat->modulus);
            ij_sub(i0, i0, xi);
            ij_div(gothp->current, i0, ijmod);
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
            sieve_hit(S, degp, pos, gothp->q, gothp->r, pos0);

          // Size of Gray codes to use for the inner loop.
#         if   defined(USE_F2)
#           define ENUM_LATTICE_UNROLL 5
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
#if defined(USE_F2) 
            if (k0 == 0 && sizeof(ij_t) == 4) {
#             define DOGRAY(n)                      \
                "xorl     " #n "(%[B]), %k[i]\n\t"  \
                "subb     %[degp], (%[S],%[i])\n\t"
#             define DOALLGRAY2 DOGRAY(0)  DOGRAY(4)  DOGRAY(0)
#             define DOALLGRAY3 DOALLGRAY2 DOGRAY(8)  DOALLGRAY2
#             define DOALLGRAY4 DOALLGRAY3 DOGRAY(12) DOALLGRAY3
#             define DOALLGRAY5 DOALLGRAY4 DOGRAY(16) DOALLGRAY4
#             define DOALLGRAY6 DOALLGRAY5 DOGRAY(20) DOALLGRAY5
#             define DOALLGRAY7 DOALLGRAY6 DOGRAY(24) DOALLGRAY6
#             define DOALLGRAY8 DOALLGRAY7 DOGRAY(28) DOALLGRAY7
#             define DOALLGRAY  CAT(DOALLGRAY, ENUM_LATTICE_UNROLL)
              uint64_t ii = i[0];
              uint8_t  dd = degp;
              __asm volatile( DOALLGRAY
                            : [i]   "+r" (ii)
                            : [S]    "r" (S+start),
                              [B]    "r" (gothp->basis),
                              [degp] "r" (dd)
                            : "memory");
              ASSERT((ii >> 32) == 0);
              i[0] = ii;
            } else 
#endif
            {
                for (unsigned k = k0; k < ngray; ++k) {
                    ij_add(i, i, gothp->basis[gray[k]]);
                    pos = start + ijvec_get_offset(i, I);
                    sieve_hit(S, degp, pos, gothp->q, gothp->r, pos0);
                }
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
              ij_add(i, i, gothp->basis[ij_deg(s)+ENUM_LATTICE_UNROLL]);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, degp, pos, gothp->q, gothp->r, pos0);
            }
          } while (rc);

          ij_set(jj, j);
          rcj = ij_monic_set_next(j, j, J);
          if (rcj) {
            ij_diff(jj, jj, j);
            degjj = ij_deg(jj);
            ij_add(gothp->current, gothp->current,
                   gothp->basis[I-L+degjj]);
            if (degjj > degj) {
              ij_sub(gothp->current, gothp->current,
                     gothp->adjustment_basis[degjj-1]);
              degj = degjj;
            }
          }
        }
    }
}
