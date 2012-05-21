#include <stdio.h>
#include <stdlib.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "norm.h"
#include "qlat.h"
#include "ijvec.h"
#include "sublat.h"


#ifdef WANT_NORM_STATS
static uint64_t norm_stats_n[2]   = {0, 0};
static int norm_stats_min[2] = {0, 0};
static uint64_t norm_stats_sum[2] = {0, 0};
static int norm_stats_max[2] = {0, 0};

void norm_stats_print() {
  printf("# Statistics on the sizes of the norms:\n");
  for (int i = 0; i < 2; ++i) {
    printf("#   Side %d: min = %u ; avg = %1.1f ; max = %u\n",
            i,
            norm_stats_min[i],
            (double)norm_stats_sum[i] / (double)norm_stats_n[i],
            norm_stats_max[i]);
  }
}
#endif


/*
  Take the N less significant bits of p and the N less significant bits of q,
  multiply them, and put the N most significant bits of the result in the low
  part of r.
  NB: in fact, it is a good practice to give p and q of degree at most N-1,
  so that all the bits of the input are indeed used.
  Then this function puts zeroes in the high part of r.
  For the moment, we put asserts to check that the input are well formed.
*/
static MAYBE_UNUSED void fppol64_mul_high(fppol64_ptr r,
        fppol64_srcptr p, fppol64_srcptr q,
        unsigned int N)
{
  ASSERT(N <= 32);
  ASSERT(N > 0);
  ASSERT(fppol64_deg(p) < (int)N);
  ASSERT(fppol64_deg(q) < (int)N);
  if (N <= 8) {         // N in [0,8]
    fppol8_t pp, qq;
    fppol16_t rr;
    fppol8_set_64(pp, p);
    fppol8_set_64(qq, q);
    fppol16_mul_8x8(rr, pp, qq);
    fppol16_div_ti(rr, rr, N-1);
    fppol64_set_16(r, rr);
  } else if (N <= 16) { // N in ]8,16]
    fppol16_t pp, qq;
    fppol32_t rr;
    fppol16_set_64(pp, p);
    fppol16_set_64(qq, q);
    fppol32_mul_16x16(rr, pp, qq);
    fppol32_div_ti(rr, rr, N-1);
    fppol64_set_32(r, rr);
  } else {              // N in ]16,32]
    fppol32_t pp, qq;
    fppol64_t rr;
    fppol32_set_64(pp, p);
    fppol32_set_64(qq, q);
    fppol64_mul_32x32(rr, pp, qq);
    fppol64_div_ti(r, rr, N-1);
  }
  ASSERT(fppol64_deg(r) < (int)N);
}


#if 0   // This is no longer used.
/* Function fppol_pow computes the power-th power of a polynomial
   should power be an unsigned int?  
   Should this function survive? 
   Will anyone use it ?  
   Does it have its place here? */
static void fppol_pow(fppol_t res, fppol_t in, int power)
{
  fppol_t tmp;
  fppol_init(tmp);
  fppol_set(tmp, in);
  for (int i = 0; i < power; i++)
    fppol_mul(tmp, tmp, in);
  fppol_set(res, tmp);
  fppol_clear(tmp);
}  
#endif

/* Function computing the norm of ffspol at (a,b)
   norm = b^d * ffspol(a/b), d = deg(ffspol) */

void ffspol_norm(fppol_t norm, ffspol_srcptr ffspol, fppol_t a, fppol_t b)
{
  fppol_t *pow_b;
  fppol_t pow_a;
  fppol_t pol_norm_i;
  fppol_t tmp_norm;
  
  fppol_init(pow_a);
  fppol_init(pol_norm_i);
  fppol_init(tmp_norm);

  fppol_set_zero(pol_norm_i);
  fppol_set_zero(tmp_norm);
  fppol_set_one(pow_a);

  /* pow_b contains b^d, b^{d-1}, ... , b^2, b, 1 */
  pow_b = (fppol_t *)malloc((ffspol->alloc) * sizeof(fppol_t));
  fppol_init(pow_b[ffspol->deg]);
  fppol_set_one(pow_b[ffspol->deg]);

  for (int i = ffspol->deg - 1; i > -1; i--) {
    fppol_init(pow_b[i]);
    fppol_mul(pow_b[i], pow_b[i+1], b);
  }
  for (int i = 0; i < ffspol->deg + 1; i++) {
    fppol_mul(pol_norm_i, ffspol->coeffs[i], pow_b[i]);
    fppol_mul(pol_norm_i, pol_norm_i, pow_a);
    fppol_add(tmp_norm, tmp_norm, pol_norm_i);
    fppol_mul(pow_a, pow_a, a);
  }
  fppol_set(norm, tmp_norm);
  fppol_clear(pow_a);
  fppol_clear(pol_norm_i);
  fppol_clear(tmp_norm);
  for (int i = 0; i <= ffspol->deg; i++)
    fppol_clear(pow_b[i]);
  free(pow_b);
}

/* Function computing ffspol_ij, a polynomial such that
   norm(ffspol,a,b) = norm(ffspol_ij, i, j), i.e.
   it is possible to apply the function norm directly on (i,j)
   with the transformed polynomial ffspol_ij, without having 
   to compute a and b with ij2ab().
   Nevertheless, due to type considerations, it should be necessary
   to make a new function norm taking as imput i and j and using the
   appropriate multiplication function.
*/

/* ffspol_ab = [f_0, ..., f_{d-1}, f_d] 
   ffspol_ij = [h_0, ..., h_{d-1}, h_d] 
   powb_ij = [g_0, ... , g_{d-1}, g_d]
 
   We consider the expression
   f(a, b) = f_d a^d + f_{d-1} a^{d-1} b + ... + f_0 b^d

   We can write it as
   f(a, b) = ff(a,b) a + f_0 b^d.
   with a = a0 i + a1 j and b = b0 i + b1 j, we have:
   f(a0 i + a1 j, b0 i + b1 j) = ff(a0 i + a1 j, b0 i + b1 j) (a0 i + a1 j) + f_0 (b0 i + b1 j)^d

   We use the formula recursively, i.e. iteration number k uses:
   h(i, j) = hh(i, j) (a0 i + a1 j) + f_{d-k} (b0 i + b1 j)^k
   We use a table powb_ij to store the coefficients of (b0 i + b1 j)^k = [g_k, g_{k-1}, ..., g_0].
 
   We have the following formula to compute from step (k-1) to step k:
   For hh(i,j) * (a0 i + a1 j):
   h_{d-k} = a1 * h_{d-k+1}
   ...
   h_{d-l} = a0 * h_{d-l} + a1 * h_{d-l+1} , 0 < l < k
   ...
   h_{d-1} = a0 * h_{d-1} + a1 * h_d
   h_d = a0 * h_d

   For (b0 i + b1 j)^k:
   g_k = b0 * g_{k-1}
   ...
   g_{k-l} = b0 * g_{k-l-1} + b1 * g_{k-l}, 0 < l < k
   ...
   g_0 = b1 * g_0
   We then multiply (b0 i + b1 j)^k by f_{d-k} and add it to hh(i,j) (a0 i + a1 j) we have computed.
   We do it for k = 0 to k = d.
   Warning: l and k for the iterations in the function definition are taken in reverse order 
   compared to this comment.
*/
static void ffspol_2ij(ffspol_ptr ffspol_ij, ffspol_srcptr ffspol_ab,
        qlat_t qlat)
{
  int d = ffspol_ab->deg;
  fppol_t *powb_ij;
  fppol_t tmp1, tmp2;

  fppol_init(tmp1);
  fppol_init(tmp2);
  
  /* Init of powb_ij which contains the coefficients of (b0 i + b1 j)^k */
  powb_ij = (fppol_t *)malloc((d+1) * sizeof(fppol_t));
  for (int k = 0; k <= d; ++k)
    fppol_init(powb_ij[k]);
  
  /* Step 0 */
  fppol_set(ffspol_ij->coeffs[d], ffspol_ab->coeffs[d]);
  ffspol_ij->deg = d;
  fppol_set_one(powb_ij[0]);
  
  for (int k = d - 1; k >= 0; --k) {
    /* For hh(i,j) * (a0 i + a1 j) */
    fppol_mul_ai(ffspol_ij->coeffs[k], ffspol_ij->coeffs[k + 1], qlat->a1);
    for (int l = k + 1; l < d; ++l) {
      fppol_mul_ai(tmp1, ffspol_ij->coeffs[l], qlat->a0);
      fppol_mul_ai(tmp2, ffspol_ij->coeffs[l + 1], qlat->a1);
      fppol_add(ffspol_ij->coeffs[l], tmp1, tmp2);
    }
    fppol_mul_ai(ffspol_ij->coeffs[d], ffspol_ij->coeffs[d], qlat->a0);
  
    /* For (b0 i + b1 j)^{d-k} */
    fppol_mul_ai(powb_ij[d - k], powb_ij[d - k - 1], qlat->b0);
    for (int l = d - k - 1; l > 0; --l) {
      fppol_mul_ai(tmp1, powb_ij[l - 1], qlat->b0);
      fppol_mul_ai(tmp2, powb_ij[l], qlat->b1);
      fppol_add(powb_ij[l], tmp1, tmp2);
    }
    fppol_mul_ai(powb_ij[0], powb_ij[0], qlat->b1);

    /* Multiply (b0 i + b1 j)^{d-k} by f_k and add it to hh(i,j) (a0 i + a1 j) we have computed */
    for (int l = k; l <= d; ++l) {
      fppol_mul(tmp1, powb_ij[l - k], ffspol_ab->coeffs[k]);
      fppol_add(ffspol_ij->coeffs[l], ffspol_ij->coeffs[l], tmp1);
    }
  }
  for (int k = 0; k <= d; ++k)
    fppol_clear(powb_ij[k]);
  free(powb_ij);

  fppol_clear(tmp1);
  fppol_clear(tmp2);
}

/* max_special(prev_max, j, &repeated) returns the maximum of prev_max
   and j
   The value repeated should be initialized to 1
   If there is only one maximum, then repeated is one
   If prev_max == j, the variable (*repeated) is *incremented* 
   This function is completely ad hoc to be used in a for loop
   like max = max_special(max, tab[i], &repeated) so the result
   will be different from max = max_special(tab[i], max, &repeated) */

static int max_special(int prev_max, int j, int *repeated)
{
  if (prev_max == j) {
    (*repeated)++;
    return prev_max;
  } 
  else if (prev_max > j) {
      return prev_max;
  }
  else {
    *repeated = 1;
    return j;
  }
}

static int deg_norm_ij(ffspol_ptr ffspol_ij, ij_t i, ij_t j)
{
  int deg, max_deg = -1;
  int repeated = 1;
  int degi = ij_deg(i);
  int degj = ij_deg(j);

  for (int k = 0; k < ffspol_ij->deg + 1; k++) {
    deg = fppol_deg(ffspol_ij->coeffs[k]) + k * degi + (ffspol_ij->deg - k) * degj;
    max_deg = max_special(max_deg, deg, &repeated);
  }

#if FP_SIZE == 2
  if (repeated & 1u) {
#else
  if (repeated == 1) {
#endif
      return max_deg;
  }
  else {
    /* We should think about a cheaper way to compute this degree
       otherwise */
    fppol_t norm, ii, jj;
    fppol_init(norm);
    fppol_init(ii);
    fppol_init(jj);
    fppol_set_ij(ii, i);
    fppol_set_ij(jj, j);
    ffspol_norm(norm, ffspol_ij, ii, jj);
    deg = fppol_deg(norm);
    fppol_clear(norm);
    fppol_clear(ii);
    fppol_clear(jj);
    return deg;
  }
}

/* Function init_norms 
   For each (i,j), it computes the corresponding (a,b) using ij2ab from
   qlat.h. Then it computes deg_norm(ffspol, a, b).
   In a first approximation, we will assume that it fits in an
   uint8_t.
   
   The sqside parameter is a boolean that tells whether we are on the
   side of the special q. If so, then the degree of q must be subtracted
   from the norm.
   */
void init_norms(uint8_t *S, ffspol_srcptr ffspol, unsigned I, unsigned J,
                ij_t j0, ijpos_t pos0, ijpos_t size, qlat_t qlat,
                int sqside, sublat_ptr sublat, MAYBE_UNUSED int side)
{
  ffspol_t ffspol_ij;
  ffspol_init2(ffspol_ij, ffspol->alloc);

  // TODO: this could be precomputed once for all and stored in qlat
  ffspol_2ij(ffspol_ij, ffspol, qlat);  

  int degq = 0;
  if (sqside)
      degq = sq_deg(qlat->q);

  ij_t i, j;
  ij_t hati, hatj;
  int rci, rcj = 1;
  for (ij_set(j, j0); rcj; rcj = ij_monic_set_next(j, j, J)) {
    ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
    if (start >= size)
      break;
    rci = 1;
    for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
      ijpos_t pos = start + ijvec_get_offset(i, I);
 
      if (S[pos] != 255) {
        // If we have sublattices, have to convert (i,j) to (hat i, hat j)
        ij_convert_sublat(hati, hatj, i, j, sublat);
#ifdef TRACE_POS
        if (pos == TRACE_POS) {
          fprintf(stderr, "TRACE_POS(%d): (hat i, hat j) = (", pos);
          ij_out(stderr, hati); fprintf(stderr, " ");
          ij_out(stderr, hatj); fprintf(stderr, ")\n");
          fprintf(stderr, "TRACE_POS(%d): norm = ", pos);
          fppol_t norm;
          fppol_init(norm);
          ffspol_norm_ij(norm, ffspol_ij, hati, hatj);
          fppol_out(stderr, norm);
          fppol_clear(norm);
          fprintf(stderr, "\n");
          fprintf(stderr, "TRACE_POS(%d): degnorm - deg(sq) = %d\n",
                  pos, fppol_deg(norm)-degq);
        }
#endif
        int deg = deg_norm_ij(ffspol_ij, hati, hatj);
        if (deg > 0) {
          ASSERT_ALWAYS(deg < 255);
          S[pos] = deg - degq;
#ifdef WANT_NORM_STATS
          norm_stats_n[side]++;
          norm_stats_sum[side] += deg;
          if ((deg < norm_stats_min[side]) || (norm_stats_min[side] == 0))
              norm_stats_min[side] = deg;
          if (deg > norm_stats_max[side])
              norm_stats_max[side] = deg;
#endif
        }
        else
          S[pos] = 255;
      }
    }
  }
  ffspol_clear(ffspol_ij); 
}
