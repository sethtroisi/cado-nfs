#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "../mod_ul.c"

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
  if (ret != 2)
    {
      fprintf (stderr, "warning: failed reading a,b\n");
      return 0;
    }
  while (str[0] != ':')
    str++;
  str++;

  // count nb of rp by counting commas (have to add one)
  {
    int cpt = 1;
    const char * pstr = str;

    while (pstr[0] != ':') {
      if (pstr[0] == ',') 
	cpt++;
      pstr++;
    }
    rel->nb_rp = cpt;
  }

  // read rp
  rel->rp = (rat_prime_t*) malloc (rel->nb_rp * sizeof (rat_prime_t));
  for (i = j = 0; i < rel->nb_rp; ++i)
    {
      ret = sscanf (str, "%lx", &p);
      if (ret!=1)
        {
          fprintf (stderr, "warning: failed reading rat prime %d\n", i);
          return 0;
        }
      /* j is the number of (p,e) pairs already stored in rel */
      for (k = 0; k < j && rel->rp[k].p != p; k++);
      if (k < j) /* prime already appeared */
        rel->rp[k].e ++;
      else /* new prime */
        {
          rel->rp[j].p = p;
          rel->rp[j].e = 1;
          j ++;
        }
      while (isxdigit(str[0]))
        str++;
      str++;
    }
  rel->nb_rp = j;
  rel->rp = (rat_prime_t*) realloc (rel->rp, j * sizeof (rat_prime_t));

  // count nb of ap by counting commas (have to add one)
  {
    int cpt = 1;
    const char * pstr = str;

    while (pstr[0] != '\n') {
      if (pstr[0] == ',') 
	cpt++;
      pstr++;
    }
    rel->nb_ap = cpt;
  }

  // read ap
  rel->ap = (alg_prime_t*) malloc (rel->nb_ap * sizeof (alg_prime_t));
  for (i = j = 0; i < rel->nb_ap; ++i)
    {
      ret = sscanf(str, "%lx", &p);
      /* corresponding root must be computed by the caller */
      if (ret != 1)
        {
          fprintf (stderr, "warning: failed reading alg prime %d\n", i);
          return 0;
        }
      /* j is the number of (p,e) pairs already stored in rel */
      for (k = 0; k < j && rel->ap[k].p != p; k++);
      if (k < j) /* prime already appeared */
        rel->ap[k].e ++;
      else /* new prime */
        {
          rel->ap[j].p = p;
          rel->ap[j].e = 1;
          j ++;
        }
      while (isxdigit(str[0]))
        str++;
      str++;
    }
  rel->nb_ap = j;
  rel->ap = (alg_prime_t*) realloc (rel->ap, j * sizeof (alg_prime_t));

  return 1;
}


// return 0 on failure
//        1 on success
//        -1 if EOF
// This reads the next valid (non blank, non commented) line and fill in
// the relation. 
int fread_relation (FILE *file, relation_t *rel)
{
  int c, i;
  char str[1024];

  // skip spaces and commented lines
  do {
    // skip spaces
    do {
      c = fgetc(file);
      if (c == EOF)
	return -1;
    } while (isspace(c));
    // skip commented lines
    if (c == '#') {
      do {
	c = fgetc(file);
	if (c == EOF)
	  return -1;
      } while (c != '\n');
    } else {
      ungetc(c, file);
      break;
    }
  } while (1);
  
  // copy line into str
  i = 0;
  do {
    c = fgetc(file);
    if (c == EOF)
      return -1;
    str[i++] = c;
    if (i == 1024) {
      fprintf(stderr, "warning: line too long\n");
      return 0;
    }
  } while (c != '\n');

  str[i] = '\0';

  return read_relation (rel, str);
}

unsigned long
findroot(long a, unsigned long b, unsigned long p) {
  unsigned long r, t;
  unsigned long pa;
  int sig, inv;

  if (a >= 0) {
    pa = (unsigned long)a;
    sig = 1;
  } else {
    pa = (unsigned long)(-a);
    sig = -1;
  }

  inv = modul_inv(&t, &b, &p);
  if (inv) {
    modul_mul(&r, &t, &pa, &p);
    if (sig == -1)
      modul_neg(&r, &r, &p);
    return r;
  } else {
    return -1;
  }
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
    fprintf(file, "%lx,", rel.rp[i]);
  fprintf(file, "%lx:", rel.rp[rel.nb_rp-1]);
  for (i = 0; i < rel.nb_ap-1; ++i)
    fprintf(file, "%lx,", rel.ap[i]);
  fprintf(file, "%lx\n", rel.ap[rel.nb_ap-1]);
}


