#include "ijvec.h"



/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/


// Allocate memory for an (i,j)-basis of degrees bounded by I and J.
void ijbasis_init(ijbasis_ptr basis, unsigned I, unsigned J)
{
  basis->I = I;
  basis->J = J;
  basis->v = (ijvec_t *)malloc((I+J) * sizeof(ijvec_t));
  ASSERT_ALWAYS(basis->v != NULL);
}


// Fill the basis with the vectors (i*t^k, j*t^k) as long as their degrees
// stay below the given bounds I and J, respectively.
static inline
unsigned fill_gap(ijvec_t *v, fbprime_t i, fbprime_t j, unsigned I, unsigned J)
{
  int degi = fbprime_deg(i), degj = fbprime_deg(j);
  int di   = (signed)I-degi, dj   = (signed)J-degj;
  int n    = degi < 0 ? dj : MIN(di, dj);
  if (n <= 0) return 0;
  ij_set_fbprime(v[0]->i, i);
  ij_set_fbprime(v[0]->j, j);
  for (int k = 1; k < n; ++k) {
    ij_mul_ti(v[k]->i, v[k-1]->i, 1);
    ij_mul_ti(v[k]->j, v[k-1]->j, 1);
  }
  return n;
}

static inline
unsigned fill_euclid(ijvec_t *v, fbprime_t i, fbprime_t j, 
        int degI, int degJ)
{
    if ((fbprime_deg(i) < degI) && (fbprime_deg(j) < degJ)) {
        ij_set_fbprime(v[0]->i, i);
        ij_set_fbprime(v[0]->j, j);
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
static void specific_euclid_char2(ijbasis_ptr basis, ijbasis_ptr euclid,
        fbprime_t alpha0, fbprime_t beta0, fbprime_t alpha1, fbprime_t beta1)
{
    ASSERT(__fbprime_SIZE <= 32);
    unsigned I = basis->I;
    unsigned J = basis->J;
    unsigned hatI = euclid->I;
    unsigned hatJ = euclid->J;
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
        basis->dim += fill_gap(basis->v+basis->dim, al1, be1,
                MIN((unsigned)da0, I), J);
        euclid->dim += fill_euclid(euclid->v + euclid->dim,
                al1, be1, hatI, hatJ);
    }
}
#endif

// Compute the (i,j)-basis of a given p-lattice.
// Large case.
void ijbasis_compute_large(MAYBE_UNUSED ijbasis_ptr euclid, ijbasis_ptr basis,
        large_fbideal_srcptr gothp, fbprime_srcptr lambda)
{
  unsigned I = basis->I;
  unsigned J = basis->J;

  // Basis is obtained from an Euclidian algorithm on
  // (p, 0) and (lambda, 1). See tex file.
  unsigned hatI = euclid->I;
  unsigned hatJ = euclid->J;
  fbprime_t alpha0, beta0, alpha1, beta1;
  fbprime_set(alpha0, gothp->p);
  fbprime_set_zero(beta0);
  fbprime_set(alpha1, lambda);
  fbprime_set_one (beta1);
  basis->dim = fill_gap(basis->v, alpha1, beta1, I, J);
  euclid->dim = fill_euclid(euclid->v, alpha1, beta1, hatI, hatJ);

#if 1 && defined(USE_F2)
  specific_euclid_char2(basis, euclid, alpha0, beta0, alpha1, beta1);
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
    basis->dim += fill_gap(basis->v+basis->dim, alpha1, beta1,
        MIN((unsigned)fbprime_deg(alpha0), I), J);
    euclid->dim += fill_euclid(euclid->v + euclid->dim, alpha1, beta1, hatI, hatJ);
  }
#endif
}


// Clean up memory.
void ijbasis_clear(ijbasis_ptr basis)
{
  free(basis->v);
}
