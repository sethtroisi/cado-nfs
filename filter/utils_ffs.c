/* Contains functions used when doing filter step for FFS instead of NFS */

#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "fppol.h"

#include "types.h"
#include "polyfactor.h"



/* return a/b mod p, and p when gcd(b,p) <> 1: this corresponds to a
   projective root */
p_r_values_t
findroot_ffs (int64_t a, uint64_t b, p_r_values_t p)
{
  fppol64_t pol_a, pol_b, pol_p;
  fppol64_t pol_r;

  fppol64_set_ui_sparse (pol_a, (uint64_t) a);
  fppol64_set_ui_sparse (pol_b, b);
  fppol64_set_ui_sparse (pol_p, (uint64_t) p);

#if DEBUG
  if (UNLIKELY(fppol64_deg(pol_p) == -1))
    {
      fprintf(stderr, "p: ul=%lu lx=%lx str=%s\n", p, p, s);
    }
#endif

  fppol64_rem (pol_a, pol_a, pol_p);
  fppol64_rem (pol_b, pol_b, pol_p);
  if (!fppol64_invmod (pol_b, pol_b, pol_p))
      return (index_t) p;

  fppol64_mulmod (pol_r, pol_a, pol_b, pol_p);

  return (index_t) fppol64_get_ui_sparse(pol_r);
}

int sq_is_irreducible(sq_srcptr p)
{
  fppol_t P;
  fppol_init(P);
  fppol_set_sq(P, p);
  int ret = fppol_is_irreducible(P);
  fppol_clear(P);
  return ret;
}

int ffs_poly_set_plist(cado_poly poly, param_list pl)
{
  unsigned long lim;
  int lpb, deg;
  const char *s;

  param_list_parse_ulong(pl, "fbb0", &lim);
  poly->pols[0]->lim = lim;
  param_list_parse_ulong(pl, "fbb1", &lim);
  poly->pols[1]->lim = lim;
  
  param_list_parse_int(pl, "lpb0", &lpb);
  poly->pols[0]->lpb = __FP_BITS + __FP_BITS * lpb;
  param_list_parse_int(pl, "lpb1", &lpb);
  poly->pols[1]->lpb = __FP_BITS + __FP_BITS * lpb;
  
  s = param_list_lookup_string(pl, "pol0");
  for (deg = 0; *s != '\0'; s++)
    if (*s == ',')
      deg++;
  poly->pols[0]->degree = deg;
  s = param_list_lookup_string(pl, "pol1");
  for (deg = 0; *s != '\0'; s++)
    if (*s == ',')
      deg++;
  poly->pols[1]->degree = deg;

  if (poly->pols[1]->degree == 1)
  {
    poly->rat  = poly->pols[1];
    poly->alg  = poly->pols[0];
  }
  else if (poly->pols[0]->degree == 1)
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
	      exit(1);
      }
    
    param_list_init(pl);
    param_list_read_stream(pl, file);
    r = ffs_poly_set_plist(poly, pl);

    param_list_clear(pl);
    fclose(file);
    return r;
}
