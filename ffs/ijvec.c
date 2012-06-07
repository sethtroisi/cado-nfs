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
void ijbasis_compute_small(ijbasis_ptr basis, ijbasis_ptr adjustment_basis,
        small_fbideal_srcptr gothp, fbprime_srcptr lambda)
{
  unsigned I = basis->I;
  unsigned J = basis->J;
  unsigned L = gothp->degq;

  // First the canonical basis:
  // Basis is { (     q*t^k,       0  ) : k in [0..I-L-1] } join
  //          { (lambda*t^k mod q, t^k) : k in [0..J-1]   }.
  basis->dim = I+J-L;
  unsigned k = 0;
  ij_set_fbprime(basis->v[k]->i, gothp->q);
  ij_set_zero(basis->v[k]->j);
  while (++k < I-L) {
    ij_mul_ti(basis->v[k]->i, basis->v[k-1]->i, 1);
    ij_set_zero(basis->v[k]->j);
  }
  ij_set_fbprime(basis->v[k]->i, lambda);
  ij_set_one(basis->v[k]->j);
  while (++k < I+J-L) {
    ij_multmod(basis->v[k]->i, basis->v[k-1]->i, gothp->q);
    ij_set_ti (basis->v[k]->j, k - (I-L));
  }

  // Transform the second part into triangular instead of diagonal (on
  // the j part) and construct the adjustment part.
  fp_t one, two;
  fp_set_one(one);
  fp_add(two, one, one);
  for (unsigned i = 0; i < J; ++i) {
    ij_smul(adjustment_basis->v[i]->i, basis->v[I-L+i]->i, two);
    ij_smul(adjustment_basis->v[i]->j, basis->v[I-L+i]->j, two);
  }
  for (unsigned i = I-L+1; i < basis->dim; ++i)
    ijvec_add(basis->v[i], basis->v[i], basis->v[i-1]);
  adjustment_basis->dim = J;
}


// Compute the (i,j)-basis of a given p-lattice.
// Large case.
void ijbasis_compute_large(ijbasis_ptr euclid, ijbasis_ptr basis,
        large_fbideal_srcptr gothp, fbprime_srcptr lambda)
{
  unsigned I = basis->I;
  unsigned J = basis->J;

  // Basis is obtained from an Euclidian algorithm on
  // (p, 0) and (lambda, 1). See tex file.
  unsigned hatI = euclid->I;
  unsigned hatJ = euclid->J;
  fp_t lc;
  fbprime_t alpha0, beta0, alpha1, beta1, q;
  fbprime_set(alpha0, gothp->p);
  fbprime_set_zero(beta0);
  fbprime_set(alpha1, lambda);
  fbprime_set_one (beta1);
  basis->dim = fill_gap(basis->v, alpha1, beta1, I, J);
  euclid->dim = fill_euclid(euclid->v, alpha1, beta1, hatI, hatJ);

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
    euclid->dim += fill_euclid(euclid->v + euclid->dim, alpha1, beta1,
        hatI, hatJ);
  }
}


// Clean up memory.
void ijbasis_clear(ijbasis_ptr basis)
{
  free(basis->v);
}
