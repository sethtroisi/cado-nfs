#include <stdlib.h>

#include "latsieve.h"
#include "types.h"
#include "ijvec.h"
#include "sublat.h"
#include "gray.h"



// FIXME: Multiplication is not compatible with lexicographical order, so there
// is no guarantee that we visit lines by increasing j's.
// However, since the positive-j part of the basis of the GF(p)-vector-space is
// (0, p*t^k) for 0 <= k < J-deg(p), and since this basis is independent from
// the special-Q or the sublattice, we can precompute a reduced row-echelon
// form of this basis for all small primes while reading the factor bases.
// Then, a simple enumeration via a p-ary monic Gray code should nicely do the
// trick (TODO!).
#if 0
void sieve_projective_root(uint8_t *S, fbideal_ptr gothp, unsigned I,
                           unsigned J, ijpos_t pos0, ijpos_t size,
                           sublat_ptr sublat)
{
    ASSERT(gothp->degp < I);
    ij_t i, j, jj, p;
    ij_set_fbprime(p, gothp->p);

    // First time round?
    if (UNLIKELY(!pos0)) {
      // Find the first line to fill
      if (use_sublat(sublat)) {
        // take the first multiple of p that is congruent to y.
        ij_t tmp, y;
        ij_set_16(tmp, sublat->modulus);
        ij_set_16(y, sublat->lat[sublat->n][1]);
        ij_mulmod(gothp->j0, y, gothp->tildep, tmp);
        ij_mul   (gothp->j0, gothp->j0, p);
        // map it into the sublattice.
        ij_sub(gothp->j0, gothp->j0, y);
        ij_div(gothp->j0, gothp->j0, tmp);
      } else {
        ij_set_zero(gothp->j0);
      }
      ij_set_zero(gothp->j);
    }

    int rci, rcj = 1;
    for(ij_set(jj, gothp->j); rcj;
        rcj = ij_monic_set_next(jj, jj, J-gothp->degp)) {
        ij_mul(j, jj, p);
        if (use_sublat(sublat))
            ij_add(j, j, gothp->j0);
        if (start >= size)
          break;
        rci = 1;
        for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
            ijpos_t pos = start + ijvec_get_offset(i, I);
#ifdef TRACE_POS
            if (pos == TRACE_POS) {
                fprintf(stderr, "TRACE_POS(%lu): ", pos);
                fbprime_out(stderr, gothp->p); fprintf(stderr, " ");
                fbprime_out(stderr, gothp->r); fprintf(stderr, "\n");
                fprintf(stderr, "TRACE_POS(%lu): degnorm is now %d\n", pos,
                        S[pos]-gothp->degp);
            }
#endif
#ifndef NDEBUG
            if (pos != 0 && (S[pos] < gothp->degp)) {
                fprintf(stderr, "faulty pos is %lu\n", pos);
            }
            ASSERT(pos == 0 || (S[pos] >= gothp->degp)); 
#endif
            S[pos] -= gothp->degp;

        }
    }
    ij_set(gothp->j, jj);
}
#endif


static inline
void sieve_hit(uint8_t *S, fbideal_srcptr gothp,
               MAYBE_UNUSED ijpos_t pos0, ijpos_t pos)
{
#ifdef TRACE_POS
  if (pos0+pos == TRACE_POS) {
    fprintf(stderr, "TRACE_POS(%lu): ", pos0+pos);
    fbprime_out(stderr, gothp->p); fprintf(stderr, " ");
    fbprime_out(stderr, gothp->r); fprintf(stderr, "\n");
    fprintf(stderr, "TRACE_POS(%lu): degnorm is now %d\n",
            pos0+pos, S[pos]-gothp->degp);
  }
#endif
#ifndef NDEBUG
  if (S[pos] < gothp->degp)
    fprintf(stderr, "faulty pos is %lu\n", pos0+pos);
  ASSERT(S[pos] >= gothp->degp);
#endif
  S[pos] -= gothp->degp;
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
            sieve_hit(S, gothp, pos0, pos);

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
              sieve_hit(S, gothp, pos0, pos);
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
              sieve_hit(S, gothp, pos0, pos);
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
