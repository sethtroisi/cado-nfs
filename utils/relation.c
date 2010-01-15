#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> /* for isxdigit */
#include <string.h>

#include "cado.h"
#include "mod_ul.h"
#include "relation.h"

/*
  Convention for I/O of rels:
    a and b are printed in decimal
    primes are printed in hexadecimal.
*/

void
copy_rel (relation_t *Rel, relation_t rel)
{
  int i;

  Rel->a = rel.a;
  Rel->b = rel.b;
  Rel->nb_rp = rel.nb_rp;
  Rel->nb_ap = rel.nb_ap;
  Rel->rp = (rat_prime_t*) realloc (Rel->rp, Rel->nb_rp * sizeof(rat_prime_t));
  Rel->ap = (alg_prime_t*) realloc (Rel->ap, Rel->nb_ap * sizeof(alg_prime_t));
  for (i = 0; i < Rel->nb_rp; ++i)
    {
      Rel->rp[i].p = rel.rp[i].p;
      Rel->rp[i].e = rel.rp[i].e;
    }
  for (i = 0; i < Rel->nb_ap; ++i)
    {
      Rel->ap[i].p = rel.ap[i].p;
      Rel->ap[i].r = rel.ap[i].r;
      Rel->ap[i].e = rel.ap[i].e;
    }
}

// Sometimes, we want our space back...!
void
clear_relation (relation_t *rel)
{
    free(rel->rp);
    free(rel->ap);
}

// return 0 on failure.
// The input should be a single line of the form
//   a,b:p1,p2,...:q1,q2,...
// Stores the primes into the relation, by collecting identical primes
// and setting the corresponding exponents. Thus a given prime will
// appear only once in rel (but it can appear with an even exponent).
int
read_relation (relation_t *rel, const char *str)
{
  int i, j, k, ret;
  unsigned long p;
  
  ret = sscanf (str, "%ld,%lu:", &(rel->a), &(rel->b));
  if (UNLIKELY(ret != 2))
    {
      fprintf (stderr, "warning: failed reading a,b in relation '%s'\n", str);
      return 0;
    }
  while (str[0] != ':')
    str++;
  str++;

  /* count number of rational primes by counting commas (have to add one) */
  {
    int cpt = 1;
    const char * pstr = str;

    while (pstr[0] != ':') {
      if (pstr[0] == ',')
	cpt++, pstr += 2; /* there cannot be two consecutive ',' or ':' */
      else
	pstr++;
    }
    rel->nb_rp = cpt;
  }

  /* read rational primes */
  rel->rp = (rat_prime_t*) malloc (rel->nb_rp * sizeof (rat_prime_t));
  /* j is the number of (p,e) pairs already stored in rel */
  for (i = j = 0; i < rel->nb_rp; ++i)
    {
      ret = sscanf (str, "%lx", &p);
      if (UNLIKELY(ret != 1))
        {
          fprintf (stderr, "warning: failed reading rat prime %d\n", i);
          return 0;
        }
      /* check if the prime already appears */
      for (k = j - 1; k >= 0 && rel->rp[k].p != p; k--);
      if (k >= 0) /* p = rel->rp[k].p */
	rel->rp[k].e ++;
      else /* new prime */
        {
          rel->rp[j].p = p;
          rel->rp[j].e = 1;
          j ++;
        }
      while (isxdigit(str[0]))
        str++;
      str++; /* skip ',' or ':' */
    }
  rel->nb_rp = j;
  rel->rp = (rat_prime_t*) realloc (rel->rp, j * sizeof (rat_prime_t));

  /* count number of algebraic primes by counting commas (have to add one) */
  {
    int cpt = 1;
    const char * pstr = str;

    while (pstr[0] != '\n') {
      if (pstr[0] == ',')
	cpt++, pstr += 2; /* there cannot be two consecutive ',' or ':' */
      else
        pstr++;
    }
    rel->nb_ap = cpt;
  }

  /* read algebraic primes */
  rel->ap = (alg_prime_t*) malloc (rel->nb_ap * sizeof (alg_prime_t));
  /* j is the number of (p,e) pairs already stored in rel */
  for (i = j = 0; i < rel->nb_ap; ++i)
    {
      ret = sscanf(str, "%lx", &p);
      /* corresponding root must be computed by the caller */
      if (UNLIKELY(ret != 1))
        {
          fprintf (stderr, "warning: failed reading alg prime %d\n", i);
          return 0;
        }
      /* check if the prime ideal already appears */
      for (k = j - 1; k >= 0 && rel->ap[k].p != p; k--);
      if (k >= 0) /* rel->ap[k].p = k */
	rel->ap[k].e ++;
      else /* new prime */
        {
          rel->ap[j].p = p;
          rel->ap[j].e = 1;
          j ++;
        }
      while (isxdigit(str[0]))
        str++;
      str++; /* skip ',' or ':' */
    }
  rel->nb_ap = j;
  rel->ap = (alg_prime_t*) realloc (rel->ap, j * sizeof (alg_prime_t));

  return 1;
}

/* Read a new relation line from file 'file', and put in in 'str'.
   Return 0 on failure
          1 on success
         -1 on EOF
*/
int
fread_buf (char str[STR_LEN_MAX], FILE *file)
{
  do
    {
      /* fgets returns NULL on error or EOF */
      if (fgets (str, STR_LEN_MAX, file) == NULL)
	return feof (file) ? -1 : 0;
    }
  while (str[0] == '#'); /* skip comments */
  return 1;
}

// return 0 on failure
//        1 on success
//        -1 if EOF
// This reads the next valid (non blank, non commented) line and fill in
// the relation. 
int fread_relation (FILE *file, relation_t *rel)
{
  int ret;
  char str[STR_LEN_MAX];

  ret = fread_buf (str, file);
  if (ret != 1)
    return ret;
  return read_relation (rel, str);
}

void skip_relations_in_file(FILE * f, int n)
{
    char str[STR_LEN_MAX];

    for( ; n-- ; ) {
	fgets(str, STR_LEN_MAX, f);
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Too long string in skip_relations_in_file\n");
	}
    }
}


/* return a/b mod p, and -1 when gcd(b,p) <> 1 */
unsigned long
findroot(long a, unsigned long b, unsigned long p) {
  int sig, inv;
  unsigned long root;
  modulusul_t m;
  residueul_t r, t, pa, br;

  modul_initmod_ul (m, p);
  modul_init (pa, m);

  if (a >= 0) {
    modul_set_ul (pa, (unsigned long)a, m); /* Does reduction mod p */
    sig = 1;
  } else {
    modul_set_ul (pa, (unsigned long)(-a), m); /* Does reduction mod p */
    sig = -1;
  }

  modul_init (r, m);
  modul_init (t, m);
  modul_init (br, m);
  modul_set_ul (br, b, m); /* Does reduction mod p */
  inv = modul_inv(t, br, m);
  if (inv) {
    modul_mul(r, t, pa, m);
    if (sig == -1)
      modul_neg(r, r, m);
    root = modul_get_ul (r, m);
  } else {
    root = (unsigned long) -1L;
  }
  
  modul_clear (pa, m); /* No-ops. Here for the sake of pedantry */
  modul_clear (r, m);
  modul_clear (t, m);
  modul_clear (br, m);
  modul_clearmod (m);
  return root;
}

// root = -1 if we don't know the result (p divides leading coeff)
void
computeroots (relation_t *rel)
{
  unsigned long r;
  int i;

  for (i = 0; i < rel->nb_ap; ++i)
    {
      r = findroot (rel->a, rel->b, rel->ap[i].p);
      rel->ap[i].r = r;
  }
}

void
fprint_relation(FILE *file, relation_t rel) {
  int i;
  fprintf(file, "%ld,%lu:", rel.a, rel.b);
  for (i = 0; i < rel.nb_rp-1; ++i)
    fprintf (file, "%lx,", rel.rp[i].p);
  fprintf (file, "%lx:", rel.rp[rel.nb_rp-1].p);
  for (i = 0; i < rel.nb_ap-1; ++i)
    fprintf (file, "%lx,", rel.ap[i].p);
  fprintf (file, "%lx\n", rel.ap[rel.nb_ap-1].p);
}

/* same as fprint_relation, but exponents > 1 are allowed */
void
fprint_relation_raw (FILE *file, relation_t rel)
{
  int i, j;

  fprintf (file, "%ld,%lu:", rel.a, rel.b);
  for (i = 0; i < rel.nb_rp; ++i)
    {
      ASSERT (rel.rp[i].e >= 1);
      for (j = 0; j < rel.rp[i].e; j++)
        {
          fprintf (file, "%lx", rel.rp[i].p);
          if (i + 1 != rel.nb_rp || j + 1 != rel.rp[i].e)
            fprintf (file, ",");
          else
            fprintf (file, ":");
        }
    }
  for (i = 0; i < rel.nb_ap; ++i)
    {
      ASSERT (rel.ap[i].e >= 1);
      for (j = 0; j < rel.ap[i].e; j++)
        {
          fprintf (file, "%lx", rel.ap[i].p);
          if (i + 1 != rel.nb_ap || j + 1 != rel.ap[i].e)
            fprintf (file, ",");
          else
            fprintf (file, "\n");
        }
    }
}

/* reduces exponents mod 2, and discards primes with even exponent */
void
reduce_exponents_mod2 (relation_t *rel)
{
  int i, j;

  if (rel->nb_rp == 0)
    fprintf (stderr, "WARNING: nb_rp = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_rp; i++)
    {
      rel->rp[i].e &= 1; /* reduce exponent mod 2 */
      if (rel->rp[i].e != 0)
        {
          rel->rp[j].p = rel->rp[i].p;
          rel->rp[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
      fprintf(stderr, "WARNING: j_rp=0 in reduce_exponents_mod2\n");
  else
      rel->rp = (rat_prime_t*) realloc (rel->rp, j * sizeof (rat_prime_t));
  rel->nb_rp = j;

  if(rel->nb_ap == 0)
      fprintf(stderr, "WARNING: nb_ap = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_ap; i++)
    {
      rel->ap[i].e &= 1; /* reduce exponent mod 2 */
      if (rel->ap[i].e != 0)
        {
          rel->ap[j].p = rel->ap[i].p;
          rel->ap[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
      fprintf(stderr, "WARNING: j_ap = 0 in reduce_exponents_mod2\n");
  else
      rel->ap = (alg_prime_t*) realloc (rel->ap, j * sizeof (alg_prime_t));
  rel->nb_ap = j;
}


