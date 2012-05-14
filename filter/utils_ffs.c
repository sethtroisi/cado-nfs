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

