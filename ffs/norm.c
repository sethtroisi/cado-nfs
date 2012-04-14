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



/* Function fppol_pow computes the power-th power of a polynomial
   should power be an unsigned int?  
   Should this function survive? 
   Will anyone use it ?  
   Does it have its place here? */
void fppol_pow(fppol_t res, fppol_t in, int power)
{
  fppol_t tmp;
  fppol_init(tmp);
  fppol_set(tmp, in);
  for (int i = 0; i < power; i++)
    fppol_mul(tmp, tmp, in);
  fppol_set(res, tmp);
  fppol_clear(tmp);
}  

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

void ffspol_2ij(ffspol_ptr ffspol_ij, ffspol_srcptr ffspol_ab, qlat_t qlat)
{
  fppol_t *pow_a0;
  fppol_t *pow_a1;
  fppol_t *pow_b0;
  fppol_t *pow_b1;
  fppol_t tmp, tmp1;
  
  /* Table containing binomial coefficients binom[n][k] for k and n
     between 0 and degree of ffspol_ab.
     It uses Lucas' theorem : binom[n][k] = prod binom[n_i][k_i] mod(p)
     where n_i and k_i are the coefficients of the p-adic representation
     of n and k.
     We need the coefficients to be constant polynomials so that
     it would be possible to multiply them with polynomials */

  /* For the moment it only works for characteristic 2 */
  fppol_t **binom;
  binom = (fppol_t **)malloc((ffspol_ab->deg + 1) * sizeof(fppol_t*));
  for (int n = 0; n < ffspol_ab->deg + 1; ++n) {
    binom[n] = (fppol_t *)malloc((n + 1) * sizeof(fppol_t));
    for(int k = 0; k < n + 1; ++k) {
      fppol_init(binom[n][k]);
      if ((n & k) == n)
	fppol_set_one(binom[n][k]);
      else
	fppol_set_zero(binom[n][k]);
    }
  }
	    
  /* 4 tables containing the powers between 0 and the degree of
     ffspol_ab of the basis vector of the q-lattice */
  pow_a0 = (fppol_t *)malloc((ffspol_ab->alloc) *
     sizeof(fppol_t)); fppol_init(pow_a0[0]);
     fppol_set_one(pow_a0[0]);

  pow_a1 = (fppol_t *)malloc((ffspol_ab->alloc) * sizeof(fppol_t));
  fppol_init(pow_a1[0]);
  fppol_set_one(pow_a1[0]);

  pow_b0 = (fppol_t *)malloc((ffspol_ab->alloc) * sizeof(fppol_t));
  fppol_init(pow_b0[0]);
  fppol_set_one(pow_b0[0]);
  
  pow_b1 = (fppol_t *)malloc((ffspol_ab->alloc) * sizeof(fppol_t));
  fppol_init(pow_b1[0]);
  fppol_set_one(pow_b1[0]);
  
  for (int n = 0; n < ffspol_ab->deg + 1; ++n) {
    fppol_init(pow_a0[n]);
    fppol_init(pow_a1[n]);
    fppol_init(pow_b0[n]);
    fppol_init(pow_b1[n]);
    fppol_mul_ai(pow_a0[n+1], pow_a0[n], qlat->a0);
    fppol_mul_ai(pow_a1[n+1], pow_a1[n], qlat->a1);
    fppol_mul_ai(pow_b0[n+1], pow_b0[n], qlat->b0);
    fppol_mul_ai(pow_b1[n+1], pow_b1[n], qlat->b1);
  }

  /* Computation of the transformed polynomial ffspol_ij */
  
  fppol_init(tmp);
  fppol_init(tmp1);
  
  for (int w = 0; w < ffspol_ab->deg + 1; ++w) {
    for (int k = 0; k < ffspol_ab->deg + 1; ++k) {
      fppol_set_zero(tmp);
      for (int u = 0; u < k + 1; ++u) {
	if ((u < k + 1) && (w - u < ffspol_ab->deg - k + 1)) {
	  fppol_set_one(tmp1);
	  fppol_mul(tmp1, tmp1, binom[k][u]);
	  fppol_mul(tmp1, tmp1, binom[ffspol_ab->deg - k][w-u]);
	  fppol_mul(tmp1, tmp1, pow_a0[u]);
	  fppol_mul(tmp1, tmp1, pow_a1[k-u]);
	  fppol_mul(tmp1, tmp1, pow_b0[w-u]);
	  fppol_mul(tmp1, tmp1, pow_b1[ffspol_ab->deg - k-w-u]);
	  fppol_add(tmp, tmp, tmp1);
	}
      }
      fppol_mul(tmp, tmp, ffspol_ab->coeffs[k]);
    }
    fppol_add(ffspol_ij->coeffs[w], ffspol_ij->coeffs[w], tmp);
  } 
 
  /* Freeing everyone */
  for (int n = 0; n < ffspol_ab->deg + 1; ++n) {
    for(int k = 0; k < n + 1; ++k) 
      fppol_clear(binom[n][k]);
    free(binom[n]);
    fppol_clear(pow_a0[n]);
    fppol_clear(pow_a1[n]);
    fppol_clear(pow_b0[n]);
    fppol_clear(pow_b1[n]);
  }
  fppol_clear(tmp1);
  fppol_clear(tmp);
  free(binom);
  free(pow_a0);
  free(pow_a1);
  free(pow_b0);
  free(pow_b1);
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

     
/* deg_norm returns the degree of the norm. 
   If during computation of pol_norm_i, only one pol_norm_i has
   maximal degree then deg_norm is equal to this degree otherwise, we
   have to compute the norm and call the fppol_deg function on it */
int deg_norm(ffspol_srcptr ffspol, fppol_t a, fppol_t b)
{
  int deg, max_deg = -1;
  int repeated = 1;
  int dega = fppol_deg(a);
  int degb = fppol_deg(b);
  static int c_deg = 0;
  static int c_tot = 0;

#if 0
  if ((c_tot & 0xFFF) == 1)
    fprintf(stderr, "deg_norm stat: %d / %d\n", c_deg, c_tot);
#endif
  c_tot++;
  for (int i = 0; i < ffspol->deg + 1; i++) {
    deg = fppol_deg(ffspol->coeffs[i]) + i * dega + (ffspol->deg - i) * degb;
    max_deg = max_special(max_deg, deg, &repeated);
  }

#if FP_SIZE == 2
  if (repeated & 1u) {
#else
  if (repeated == 1) {
#endif
      c_deg++;
      return max_deg;
  }
  else {
    /* We should think about a cheaper way to compute this degree
       otherwise */
    fppol_t norm;
    fppol_init(norm);
    ffspol_norm(norm, ffspol, a, b);
    deg = fppol_deg(norm);
    fppol_clear(norm);
    return deg;
  }
}


/* Function init_norms 
   For each (i,j), it compute the corresponding (a,b) using ij2ab from
   qlat.h. Then it computes deg_norm(ffspol, a, b).
   In a first approximation, we will assume that it fits in an
   uint8_t.
   For the moment i and j are assumed to be unsigned int for their
   limb part 
   Enumerating i, j the following way only works for charac. 2 
   
   The last parameter is a boolean that tells whether we are on the
   side of the special q. If so, then the degree of q must be subtracted
   from the norm.
   */

void init_norms(uint8_t *S, ffspol_srcptr ffspol, unsigned I, unsigned J,
                ij_t j0, ijpos_t pos0, ijpos_t size, qlat_t qlat,
                int sqside, sublat_ptr sublat)
{
  fppol_t a, b;
  fppol_init(a);
  fppol_init(b);
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
          ij2ab(a, b, hati, hatj, qlat);
          ffspol_norm(norm, ffspol, a, b);
          fppol_out(stderr, norm);
          fppol_clear(norm);
          fprintf(stderr, "\n");
          fprintf(stderr, "TRACE_POS(%d): degnorm - deg(sq) = %d\n",
                  pos, fppol_deg(norm)-degq);
        }
#endif
        ij2ab(a, b, hati, hatj, qlat);
        int deg = deg_norm(ffspol, a, b);
        if (deg > 0) {
          ASSERT_ALWAYS(deg < 255);
          S[pos] = deg - degq;
        }
        else
          S[pos] = 255;
      }
    }
  }

  fppol_clear(a);
  fppol_clear(b);
}
