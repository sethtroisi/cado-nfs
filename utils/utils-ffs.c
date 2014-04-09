/* Contains functions used when doing filter step for FFS instead of NFS */

#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "fppol.h"
#include "fq.h"

#include "types.h"
#include "polyfactor.h"


/* Set r to a/b mod p. Assume b and p are coprime. */
static inline void
ffs_mulinvmod (fppol64_t r, fppol64_t a, fppol64_t b, fppol64_t p)
{
  fppol64_t t;
  int ret = fppol64_invmod (t, b, p);
  ASSERT (ret);
  fppol64_mulmod (r, a, t, p);
}

/* return a/b mod p, and p when gcd(b,p) <> 1: this corresponds to a
   projective root */
p_r_values_t
ffs_relation_compute_r (int64_t a, uint64_t b, p_r_values_t p)
{
  fppol64_t pol_a, pol_b, pol_p;
  fppol64_t pol_r;

  fppol64_set_ui_sparse (pol_a, (uint64_t) a);
  fppol64_set_ui_sparse (pol_b, b);
  fppol64_set_ui_sparse (pol_p, (uint64_t) p);

  fppol64_rem (pol_a, pol_a, pol_p);
  fppol64_rem (pol_b, pol_b, pol_p);
  if (fppol64_is_zero(pol_b))
    return (p_r_values_t) p;
  else
  {
    ffs_mulinvmod (pol_r, pol_a, pol_b, pol_p);
    return (p_r_values_t) fppol64_get_ui_sparse(pol_r);
  }
}

p_r_values_t
ffs_compute_r (int64_t a, uint64_t b, p_r_values_t p, p_r_values_t pk)
{
  fppol64_t pol_a, pol_b, pol_p, pol_pk, pol_rk, pol_tmp;

  fppol64_set_ui_sparse (pol_a, (uint64_t) a);
  fppol64_set_ui_sparse (pol_b, b);
  fppol64_set_ui_sparse (pol_p, (uint64_t) p);
  fppol64_set_ui_sparse (pol_pk, (uint64_t) pk);

  fppol64_rem (pol_tmp, pol_b, pol_p);

  fppol64_rem (pol_a, pol_a, pol_pk);
  fppol64_rem (pol_b, pol_b, pol_pk);

  if (fppol64_is_zero(pol_tmp))
  {
    ffs_mulinvmod (pol_rk, pol_b, pol_a, pol_pk);
    fppol64_add (pol_rk, pol_pk, pol_rk);
    return (p_r_values_t) fppol64_get_ui_sparse (pol_rk);
  }
  else
  {
    ffs_mulinvmod (pol_rk, pol_a, pol_b, pol_pk);
    return (p_r_values_t) fppol64_get_ui_sparse (pol_rk);
  }
}

void
ffs_compute_pk_r (p_r_values_t *pk, p_r_values_t *r, p_r_values_t p,
                  p_r_values_t rk, unsigned int k)
{
  fppol64_t pol_r, pol_rk, pol_p, pol_pk;

  fppol64_set_ui_sparse (pol_p, (uint64_t) p);
  fppol64_set_ui_sparse (pol_pk, (uint64_t) p);
  fppol64_set_ui_sparse (pol_rk, (uint64_t) rk);
  for (unsigned int i = 1; i < k; ++i)
    fppol64_mul (pol_pk, pol_pk, pol_p);
  *pk = (p_r_values_t) fppol64_get_ui_sparse (pol_pk);

  if (fppol64_deg(pol_rk) < fppol64_deg(pol_pk))
  {
    fppol64_rem (pol_r, pol_rk, pol_p);
    *r = (p_r_values_t) fppol64_get_ui_sparse (pol_r);
  }
  else {
    fppol64_t x;
    fppol64_sub (x, pol_rk, pol_pk);
    fppol64_rem (x, x, pol_p);
    fppol64_add (pol_r, pol_p, x);
    *r = (p_r_values_t) fppol64_get_ui_sparse (pol_r);
  }
}

int ffs_poly_set_plist(cado_poly poly, param_list pl)
{
  unsigned long lim;
  int lpb, deg;
  const char *s;

  param_list_parse_ulong(pl, "fbb0", &lim);
  //poly->pols[0]->lim = lim;
  param_list_parse_ulong(pl, "fbb1", &lim);
  //poly->pols[1]->lim = lim;
  
  param_list_parse_int(pl, "lpb0", &lpb);
  //poly->pols[0]->lpb = __FP_BITS + __FP_BITS * lpb;
  param_list_parse_int(pl, "lpb1", &lpb);
  //poly->pols[1]->lpb = __FP_BITS + __FP_BITS * lpb;
  
  s = param_list_lookup_string(pl, "pol0");
  for (deg = 0; *s != '\0'; s++)
    if (*s == ',')
      deg++;
  poly->pols[0]->deg = deg;
  s = param_list_lookup_string(pl, "pol1");
  for (deg = 0; *s != '\0'; s++)
    if (*s == ',')
      deg++;
  poly->pols[1]->deg = deg;

  if (poly->pols[1]->deg == 1)
  {
    poly->rat  = poly->pols[1];
    poly->alg  = poly->pols[0];
  }
  else if (poly->pols[0]->deg == 1)
  {
    poly->rat  = poly->pols[0];
    poly->alg  = poly->pols[1];
  }

  return 1;
}

// returns 0 on failure, 1 on success.
int ffs_poly_read(cado_poly poly, const char *filename)
{
    FILE *file;
    int r;
    param_list pl;

    file = fopen(filename, "r");
    if (file == NULL) 
    {
	    fprintf(stderr, "read_polynomial: could not open %s\n", filename);
      return 0;
    }
    
    param_list_init(pl);
    param_list_read_stream(pl, file);
    r = ffs_poly_set_plist(poly, pl);

    param_list_clear(pl);
    fclose(file);
    return r;
}
