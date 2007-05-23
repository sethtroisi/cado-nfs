#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

/*
  Convention for I/O of rels:
    a and b are printed in decimal
    primes are printed in hexadecimal.
*/



static void
copy_rel(relation_t *Rel, relation_t rel) {
  int i;
  Rel->a = rel.a;
  Rel->b = rel.b;
  Rel->nb_rp = rel.nb_rp;
  Rel->nb_ap = rel.nb_ap;
  Rel->rp = (unsigned long *)realloc(Rel->rp, Rel->nb_rp*sizeof(unsigned long));
  Rel->ap = (unsigned long *)realloc(Rel->ap, Rel->nb_ap*sizeof(unsigned long));
  for (i = 0; i < Rel->nb_rp; ++i)
    Rel->rp[i] = rel.rp[i];
  for (i = 0; i < Rel->nb_ap; ++i)
    Rel->ap[i] = rel.ap[i];
  if (rel.ar != NULL) {
    Rel->ar = (unsigned long *)realloc(Rel->ar, Rel->nb_ap*sizeof(unsigned long));
    for (i = 0; i < Rel->nb_ap; ++i)
      Rel->ar[i] = rel.ar[i];
  }
}



// return 0 on failure.
// The input should be a single line of the form
//   a,b:p1,p2,...:q1,q2,...
static int read_relation(relation_t *rel, const char *str) {
  int i, ret;
  
  ret = sscanf(str, "%ld,%lu:", &(rel->a), &(rel->b));
  if (ret!=2) {
    fprintf(stderr, "warning: failed reading a,b\n");
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
  rel->rp = (unsigned long *) malloc (rel->nb_rp*sizeof(unsigned long));
  for (i = 0; i < rel->nb_rp; ++i) {
    ret = sscanf(str, "%lx", &(rel->rp[i]));
    if (ret!=1) {
      fprintf(stderr, "warning: failed reading rat prime %d\n", i);
      return 0;
    }
    while (isxdigit(str[0]))
      str++;
    str++;
  }

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
  rel->ap = (unsigned long *) malloc (rel->nb_ap*sizeof(unsigned long));
  for (i = 0; i < rel->nb_ap; ++i) {
    ret = sscanf(str, "%lx", &(rel->ap[i]));
    if (ret!=1) {
      fprintf(stderr, "warning: failed reading alg prime %d\n", i);
      return 0;
    }
    while (isxdigit(str[0]))
      str++;
    str++;
  }
  rel->ar = NULL;

  return 1;
}


// return 0 on failure
//        1 on success
//        -1 if EOF
// This reads the next valid (non blank, non commented) line and fill in
// the relation. 
static int fread_relation(FILE *file, relation_t *rel) {
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

  return read_relation(rel, str);
}

static unsigned long
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

// root = -1 if we don;t know the result (p divides leading coeff)
static void
computeroots(relation_t * rel) {
  unsigned long r;
  int i;

  rel->ar = (unsigned long *) realloc (rel->ar, rel->nb_ap*sizeof(unsigned long));

  for (i = 0; i < rel->nb_ap; ++i) {
    r = findroot(rel->a, rel->b, rel->ap[i]);
    rel->ar[i] = r;
  }
}

static void
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


