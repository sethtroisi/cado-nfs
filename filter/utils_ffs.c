/* Contains functions used when doing filter step for FFS instead of NFS */

#include "fppol.h"
#include "utils.h"
#include "utils_ffs.h"

unsigned long 
findroot_ffs (long a, unsigned long b, unsigned long p)
{
  char s[17];
  fppol64_t pol_a, pol_b, pol_p;
  fppol64_t pol_r;

  sprintf (s, "%lx", (unsigned long) a);
  if (fppol64_set_str (pol_a, s) != 1)
      fprintf(stderr, "Error in findroot_ffs with a=%s\n", s);

  sprintf (s, "%lx", b);
  if (fppol64_set_str (pol_b, s) != 1)
      fprintf(stderr, "Error in findroot_ffs with b=%s\n", s);

  sprintf (s, "%lx", p);
  if (fppol64_set_str (pol_p, s) != 1)
      fprintf(stderr, "Error in findroot_ffs with c=%s\n", s);
      
  fppol64_rem (pol_a, pol_a, pol_p);
  fppol64_rem (pol_b, pol_b, pol_p);
  if (!fppol64_invmod (pol_b, pol_b, pol_p))
      return (unsigned long) -1L;

  fppol64_mulmod (pol_r, pol_a, pol_b, pol_p);

  fppol64_get_str (s, pol_r);
  return strtol (s, NULL, 16);
}

void
computeroots_ffs (relation_t *rel)
{
  unsigned long r;
  int i;

  for (i = 0; i < rel->nb_ap; ++i)
    {
      r = findroot_ffs (rel->a, rel->b, rel->ap[i].p);
      rel->ap[i].r = r;
  }
}

int ffs_poly_set_plist(cado_poly poly, param_list pl)
{
  param_list_parse_ulong(pl, "fbb0", &(poly->rat->lim));
  param_list_parse_int(pl, "lpb0", &(poly->rat->lpb));
  param_list_parse_ulong(pl, "fbb1", &(poly->alg->lim));
  param_list_parse_int(pl, "lpb1", &(poly->alg->lpb));

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

