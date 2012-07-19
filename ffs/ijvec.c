#include "ijvec.h"



/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/

// Fill the basis with the vectors (i*t^k, j*t^k) as long as their degrees
// stay below the given bounds I and J, respectively.
static inline
unsigned fill_gap(ijvec_t *v, fbprime_t i, fbprime_t j,
                  int max_degi, int max_degj, unsigned I)
{
  int degi = fbprime_deg(i),  degj = fbprime_deg(j);
  int di   = max_degi - degi, dj   = max_degj - degj;
  int n    = degi < 0 ? dj : MIN(di, dj);
  if (n <= 0) return 0;
  ij_t ii, jj;
  ij_set_fbprime(ii, i);
  ij_set_fbprime(jj, j);
  ijvec_set_i_j(v[0], ii, jj, I);
  for (int k = 1; k < n; ++k)
    ijvec_mul_ti(v[k], v[k-1], 1);
  return n;
}

static inline
unsigned fill_euclid(ijvec_t *v, fbprime_t i, fbprime_t j, 
                     int max_degi, int max_degj, unsigned I)
{
  if ((fbprime_deg(i) < max_degi) && (fbprime_deg(j) < max_degj)) {
    ij_t ii, jj;
    ij_set_fbprime(ii, i);
    ij_set_fbprime(jj, j);
    ijvec_set_i_j(v[0], ii, jj, I);
    return 1;
  }
  return 0;
}

// Compute the (i,j)-basis of a given p-lattice.
// Small case.
void ijbasis_compute_small(ij_t *basis, ij_t *adjustment_basis,
        small_fbideal_srcptr gothp, fbprime_srcptr lambda,
        unsigned I, unsigned J)
{
  unsigned L = gothp->degq;

  // First the canonical basis:
  // Basis is { (     q*t^k,       0  ) : k in [0..I-L-1] } join
  //          { (lambda*t^k mod q, t^k) : k in [0..J-1]   }.
  // The J vectors are not stored, however.
  unsigned k = 0;
  ij_set_fbprime(basis[k], gothp->q);
  while (++k < I-L)
    ij_mul_ti(basis[k], basis[k-1], 1);
  ij_set_fbprime(basis[k], lambda);
  while (++k < I+J-L)
    ij_multmod(basis[k], basis[k-1], basis[0]);

  // Transform the second part into triangular instead of diagonal (on
  // the j part) and construct the adjustment part.
  fp_t one, two;
  fp_set_one(one);
  fp_add(two, one, one);
  for (unsigned k = 0; k < J; ++k)
    ij_smul(adjustment_basis[k], basis[I-L+k], two);
  for (unsigned k = I-L+1; k < I+J-L; ++k)
    ij_add(basis[k], basis[k], basis[k-1]);
}


#ifdef USE_F2
static
void specific_euclid_char2(ijvec_t *basis,  unsigned *basis_dim, 
                           unsigned I,      unsigned  J,
                           ijvec_t *euclid, unsigned *euclid_dim,
                           unsigned hatI,   unsigned  hatJ,
                           fbprime_srcptr alpha0, fbprime_srcptr beta0,
                           fbprime_srcptr alpha1, fbprime_srcptr beta1)
{
    ASSERT(__fbprime_SIZE <= 32);
    fppol64_t v0, v1;
    v0[0] = alpha0[0] | ((uint64_t)beta0[0] << 32);
    v1[0] = alpha1[0] | ((uint64_t)beta1[0] << 32);

    int da0, da1;
    da0 = fbprime_deg(alpha0);
    da1 = fbprime_deg(alpha1);

    while ((unsigned)fppol64_deg(v1) < 32+J && da1 >= 0) {
        // alpha0 = alpha0 mod alpha1  (and betas follow)
        uint64_t mask = ((1U<<(da0+1))-1) ^ ((1U<<da1)-1);
        uint64_t shiftv1 = v1[0] << (da0-da1);
        uint64_t mask1 = -(uint64_t)1;
        uint64_t maskbit = 1U<<da0;
        do {
            v0[0] ^= shiftv1 & mask1;
            if (!(v0[0] & mask))
                break;
            maskbit >>= 1;
            mask1 = (v0[0] & maskbit) ? (-(uint64_t)1) : 0;
            shiftv1 >>= 1;
        } while (1);
        fppol64_swap(v0, v1);
        da0 = da1;
        da1 = fppol32_deg((fppol32_ptr)&v1[0]);
        // fill gap
        fbprime_ptr al1, be1;
        al1 = (fbprime_ptr)&v1[0];
        be1 = (fbprime_ptr)(((fppol32_ptr)&v1[0])+1);
        *basis_dim  += fill_gap   (basis +*basis_dim,  al1, be1,
                                   MIN(da0, (signed)I), J, I);
        *euclid_dim += fill_euclid(euclid+*euclid_dim, al1, be1,
                                   hatI, hatJ, hatI);
    }
}
#endif

// Compute the (i,j)-basis of a given p-lattice.
// Large case.
void ijbasis_compute_large(ijvec_t *basis,  unsigned *basis_dim,
                           unsigned I,      unsigned  J,
                           ijvec_t *euclid, unsigned *euclid_dim,
                           unsigned hatI,   unsigned  hatJ,
                           large_fbideal_srcptr gothp, fbprime_srcptr lambda)
{
  // Basis is obtained from an Euclidian algorithm on
  // (p, 0) and (lambda, 1). See tex file.
  fbprime_t alpha0, beta0, alpha1, beta1;
  fbprime_set(alpha0, gothp->p);
  fbprime_set_zero(beta0);
  fbprime_set(alpha1, lambda);
  fbprime_set_one (beta1);
  *basis_dim  = fill_gap   (basis,  alpha1, beta1, I,    J,    I);
  *euclid_dim = fill_euclid(euclid, alpha1, beta1, hatI, hatJ, hatI);

#if 1 && defined(USE_F2)
  specific_euclid_char2(basis,  basis_dim,  I,    J,
                        euclid, euclid_dim, hatI, hatJ,
                        alpha0, beta0, alpha1, beta1);
#else
  fbprime_t q;
  fp_t lc;
  while ((unsigned)fbprime_deg(beta1) < J && !fbprime_is_zero(alpha1)) {
    fbprime_divrem(q, alpha0, alpha0, alpha1);
    fbprime_submul(beta0, beta1, q);
    fbprime_swap(alpha0, alpha1);
    fbprime_swap(beta0,  beta1);
    // Keep beta1 monic.
    fbprime_get_coeff(lc, beta1, fbprime_deg(beta1));
    fbprime_sdiv(alpha1, alpha1, lc);
    fbprime_sdiv(beta1,  beta1,  lc);
    *basis_dim  += fill_gap   (basis +*basis_dim,  alpha1, beta1,
                               MIN(fbprime_deg(alpha0), (signed)I), J, I);
    *euclid_dim += fill_euclid(euclid+*euclid_dim, alpha1, beta1,
                               hatI, hatJ, hatI);
  }
#endif
}
