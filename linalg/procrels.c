#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "relation.h"


typedef struct {
  unsigned long prime;
  unsigned long root;
} rootprime_t;

typedef struct {
  rootprime_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_rootprime_t;

typedef struct {
  unsigned long *tab;
  unsigned long allocated;
  unsigned long length;
} tab_prime_t;

typedef struct {
  long a;
  unsigned long b;
} abpair_t;

typedef struct {
  abpair_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_abpair_t;



// TODO: share this with sieve.c
static unsigned long
iscomposite (const unsigned long n)
{
  unsigned long i, i2;

  if (n % 2 == 0)
    return (n == 2) ? 0 : 2;

  /* (i + 2)^2 = i^2 + 4*i + 4 */
  for (i = 3, i2 = 9; i2 <= n; i2 += (i+1) * 4, i += 2)
    if (n % i == 0)
        return i;

  return 0;
}

void
add_alg_rootprime(tab_rootprime_t *table, tab_prime_t *bad_primes,
    unsigned long p, long a, unsigned long b) {
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

    if (table->allocated == table->length) {
      table->allocated += 100;
      table->tab = (rootprime_t *)realloc((void *)table->tab, (table->allocated)*sizeof(rootprime_t));
      assert (table->tab != NULL);
    }
    table->tab[table->length].prime = p;
    table->tab[table->length].root = r;
    table->length++;
  } else {
    if (bad_primes->allocated == bad_primes->length) {
      bad_primes->allocated += 100;
      bad_primes->tab = (unsigned long *)realloc((void *)bad_primes->tab, (bad_primes->allocated)*sizeof(unsigned long));
      assert (bad_primes->tab != NULL);
    }
    bad_primes->tab[bad_primes->length] = p;
    bad_primes->length++;
  }
}

void 
add_rat_prime(tab_prime_t *table, unsigned long p) {
  if (table->allocated == table->length) {
    table->allocated += 100;
    table->tab = (unsigned long *)realloc(table->tab, (table->allocated)*sizeof(unsigned long));
    assert (table->tab != NULL);
  }
  table->tab[table->length] = p;
  table->length++;
}

// When entering this function, the first colon has been read.
void 
update_one_line(FILE *file, tab_rootprime_t *alg_table, tab_prime_t *rat_table,
    tab_prime_t * bad_primes, long a, unsigned long b) {
  char c;
  unsigned long p;

  do {
    c = fgetc(file);
    if (c == ':') 
      break;
    if (c != ',')
      ungetc(c, file);
    fscanf(file, "%lx", &p);
    add_rat_prime(rat_table, p);
  } while (1);

  do {
    c = fgetc(file);
    if (c == '\n') 
      break;
    if (c != ',')
      ungetc(c, file);
    fscanf(file, "%lx", &p);
    add_alg_rootprime(alg_table, bad_primes, p, a, b);
  } while (1);
}

void
update_tables(FILE *file, tab_rootprime_t *alg_table, tab_prime_t *rat_table,
    tab_prime_t *bad_primes) {
  long a;
  unsigned long b;
  int c;
  int sig;
  int cpt = 0;

  do {
    c = fgetc(file);
    if (c == EOF) 
      break;
    if (c == '#') {
      do {
	c = fgetc(file);
      } while (c != '\n');
      continue;
    }
    if (c == '-')
      sig = -1;
    else {
      sig = 1;
      ungetc(c, file);
    }
    fscanf(file, "%lu,%lu:", &a, &b);
    a *= sig;
    cpt++;
    update_one_line(file, alg_table, rat_table, bad_primes, a, b);
    if ((cpt % 1000) == 0)
      fprintf(stderr, "have read %d rels\n", cpt);
  } while (1);
}

static int
cmpprimes(const void *p, const void *q) {
  unsigned long pp, qq;
  pp = *((unsigned long *)p);
  qq = *((unsigned long *)q);
  if (pp < qq) 
    return -1;
  if (pp > qq)
    return 1;
  return 0;
}

static void
uniqprimes(tab_prime_t * table) {
  int i, j, newlength;
  unsigned long oldp;
  newlength = 0;
  j = 0;
  oldp = 0;
  for (i = 0; i < table->length; ++i) {
    if (table->tab[i] != oldp) {
      oldp = table->tab[i];
      table->tab[j] = oldp;
      newlength++;
      j++;
    }
  }
  table->length=newlength;
}


static int
cmprootprimes(const void *p, const void *q) {
  rootprime_t pp, qq;
  pp = *((rootprime_t *)p);
  qq = *((rootprime_t *)q);
  if ((pp.prime < qq.prime)
      || ((pp.prime == qq.prime) && (pp.root < qq.root)))
    return -1;
  if ((pp.prime > qq.prime)
      || ((pp.prime == qq.prime) && (pp.root > qq.root)))
    return 1;
  return 0;
}

static void
uniqrootprimes(tab_rootprime_t * table) {
  int i, j, newlength;
  unsigned long oldp, oldr;
  newlength = 0;
  j = 0;
  oldp = 0;
  for (i = 0; i < table->length; ++i) {
    if ((table->tab[i].prime != oldp) || (table->tab[i].root != oldr))  {
      oldp = table->tab[i].prime;
      oldr = table->tab[i].root;
      table->tab[j].prime = oldp;
      table->tab[j].root = oldr;
      newlength++;
      j++;
    }
  }
  table->length=newlength;
}

static int 
cmpabpairs(const void *p, const void *q) {
  abpair_t pp, qq;
  pp = *((abpair_t *)p);
  qq = *((abpair_t *)q);
  if ((pp.a < qq.a)
      || ((pp.a == qq.a) && (pp.b < qq.b)))
    return -1;
  if ((pp.a > qq.a)
      || ((pp.a == qq.a) && (pp.b > qq.b)))
    return 1;
  return 0;
}



static int
getindex_rat(tab_prime_t rat_table, unsigned long p) {
  int i, j;
  int found = 0;

  i = 0; j = rat_table.length-1;
  while ((i<j-1) && (rat_table.tab[i] != p) && (rat_table.tab[j] != p)) {
    int m = (i+j)/2;
    if (rat_table.tab[m] > p)
      j = m;
    else
      i = m;
  }
  if (rat_table.tab[i] == p)
    return i;
  if (rat_table.tab[j] == p)
    return j;
  return -1;
}

static int
getindex_abpair(tab_abpair_t table, long a, unsigned long b) {
  int i, j;

  i = 0; j = table.length-1;
  while ((i<j-1) &&
      ( (table.tab[i].a != a) || (table.tab[i].b != b)) &&
      ( (table.tab[j].a != a) || (table.tab[j].b != b))) {
    int m = (i+j)/2;
    if ((table.tab[m].a > a) ||
	((table.tab[m].a == a) && (table.tab[m].b > b)))
      j = m;
    else
      i = m;
  }
  if ((table.tab[i].a == a) && (table.tab[i].b == b))
    return i;
  if ((table.tab[j].a == a) && (table.tab[j].b == b))
    return j;
  return -1;
}


static int
getindex_alg(tab_rootprime_t alg_table, unsigned long p, unsigned long r) {
  int i, j;

  i = 0; j = alg_table.length-1;
  while ((i<j-1) &&
      ( (alg_table.tab[i].prime != p) || (alg_table.tab[i].root != r)) &&
      ( (alg_table.tab[j].prime != p) || (alg_table.tab[j].root != r))) {
    int m = (i+j)/2;
    if ((alg_table.tab[m].prime > p) ||
	((alg_table.tab[m].prime == p) && (alg_table.tab[m].root > r)))
      j = m;
    else
      i = m;
  }
  if ((alg_table.tab[i].prime == p) && (alg_table.tab[i].root == r))
    return i;
  if ((alg_table.tab[j].prime == p) && (alg_table.tab[j].root == r))
    return j;
  return -1;
}

void onepass_singleton_removal(relation_t *rel_table, int *nb_rel, 
    tab_rootprime_t * alg_table, tab_prime_t * rat_table) {
  long l_rat, l_alg;
  int i, j;
  tab_abpair_t ab_single;

  l_rat = rat_table->length;
  l_alg = alg_table->length;

  ab_single.allocated = l_rat+l_alg;
  ab_single.length = l_rat+l_alg;
  ab_single.tab = (abpair_t *) malloc (ab_single.allocated*sizeof(abpair_t));
  for (i = 0; i < ab_single.length; ++i) {
    ab_single.tab[i].a = 0;
    ab_single.tab[i].b = 0;
  }

  for (i = 0; i < *nb_rel; ++i) {
    int j;
    for (j = 0; j < rel_table[i].nb_rp; ++j) {
      int index;
      index = getindex_rat(*rat_table, rel_table[i].rp[j]);
      if (index == -1)
	continue;
      if ((ab_single.tab[index].a == 0) && (ab_single.tab[index].b == 0)) {
	ab_single.tab[index].a = rel_table[i].a;
	ab_single.tab[index].b = rel_table[i].b;
      } else if ((ab_single.tab[index].a != rel_table[i].a) || (ab_single.tab[index].b != rel_table[i].b)) {
	ab_single.tab[index].a = 1;
	ab_single.tab[index].b = 1;
      }
    }
    for (j = 0; j < rel_table[i].nb_ap; ++j) {
      int index;
      index = getindex_alg(*alg_table, rel_table[i].ap[j], rel_table[i].ar[j]);
      if (index == -1)
	continue;
      index += l_rat;
      if ((ab_single.tab[index].a == 0) && (ab_single.tab[index].b == 0)) {
	ab_single.tab[index].a = rel_table[i].a;
	ab_single.tab[index].b = rel_table[i].b;
      } else if ((ab_single.tab[index].a != rel_table[i].a) || (ab_single.tab[index].b != rel_table[i].b)) {
	ab_single.tab[index].a = 1;
	ab_single.tab[index].b = 1;
      }
    }
  }

  // remove singleton rational primes
  {
    unsigned long * ptab = rat_table->tab;
    j = 0;
    for (i = 0; i < l_rat; ++i) {
      if ((ab_single.tab[i].a == 0) && (ab_single.tab[i].b == 0)) {
	fprintf(stderr, "warning: this p is untouched ! %ld\n", i);
	continue;
      }
      if ((ab_single.tab[i].a != 1) || (ab_single.tab[i].b != 1)) {
	// printf("singleton pair: %ld %ld, with rat p= %lx\n", ab_single.tab[i].a,
	//    ab_single.tab[i].b, rat_table->tab[i]);
      } else {
	ptab[j] = rat_table->tab[i];
	j++;
      }
    }
    rat_table->length = j;
  }
  fprintf(stderr, "%d rat primes remain\n", rat_table->length);
  // remove singleton algebraic primes
  {
    rootprime_t *ptab = alg_table->tab;
    j = 0;
    for (i = 0; i < l_alg; ++i) {
      if ((ab_single.tab[l_rat+i].a == 0) && (ab_single.tab[l_rat+i].b == 0)) {
        fprintf(stderr, "warning: this p is untouched ! %ld\n", i);
	continue;
      }
      if ((ab_single.tab[l_rat+i].a != 1) || (ab_single.tab[l_rat+i].b != 1)) {
	//printf("singleton pair: %ld %ld, with alg (p,r)=(%lx,%lx)\n", ab_single.tab[i].a,
	//    ab_single.tab[i].b, alg_table->tab[i-l_rat].prime, alg_table->tab[i-l_rat].root);
      } else {
        ptab[j].prime = alg_table->tab[i].prime; 
        ptab[j].root = alg_table->tab[i].root; 
        j++;
      }
    }
    alg_table->length = j;
  }
  fprintf(stderr, "%d alg primes remain\n", alg_table->length);
  // remove singleton relations (sort bad (a,b) pairs, and then browse
  // through all relations).
  qsort(ab_single.tab, l_alg+l_rat, sizeof(abpair_t), cmpabpairs);
  {
    relation_t * ptab = rel_table;
    j = 0;
    for (i = 0; i < *nb_rel; ++i) {
      if (getindex_abpair(ab_single, rel_table[i].a, rel_table[i].b) == -1) {
	copy_rel(&(ptab[j]), rel_table[i]);
	j++;
      }
    }
    *nb_rel = j;
  }
  fprintf(stderr, "%d relations remain\n", *nb_rel);
}






relation_t *
detect_singleton(FILE *file, tab_rootprime_t * alg_table, tab_prime_t * rat_table, int *nb_rel) {
  long l_rat, l_alg;
  relation_t *rel_table;
  relation_t rel;
  int alloced = 100;
  int ret, i;

  rel_table = (relation_t *) malloc(alloced*sizeof(relation_t));
 
  *nb_rel = 0;
  do {
    if (*nb_rel == (alloced - 1)) {
      alloced += 100;
      rel_table = (relation_t *) realloc(rel_table, alloced*sizeof(relation_t));
    }
    ret = fread_relation(file, &(rel_table[*nb_rel]));
    if (ret == 1)
      (*nb_rel)++;
  } while (ret == 1);

  fprintf(stderr, "loaded %d relations\n", *nb_rel);

  for (i = 0; i < *nb_rel; ++i) 
    computeroots(&(rel_table[i]));

  {
    int old_nb;
    do {
      old_nb = *nb_rel;
      fprintf(stderr, "** Do one pass of singleton removal...\n");
      onepass_singleton_removal(rel_table, nb_rel, alg_table, rat_table);
    } while (*nb_rel != old_nb);
  }

  return rel_table;
}



int main(int argc, char **argv) {
  char filename[] = "rels";
  tab_rootprime_t alg_table;
  tab_prime_t rat_table, bad_primes;
  tab_abpair_t ab_single;
  FILE *file;
  relation_t * rel_table;
  int nb_rel;
  int i;

  if (argc == 1) 
    file = stdin;
  else {
    if (argc != 2) {
      fprintf(stderr, "usage: %s [filename]\n", argv[0]);
      fprintf(stderr, "  if no filename is given, takes input on stdin\n");
      exit(1);
    }
    file = fopen(argv[1], "r");
    if (file == NULL) {
      fprintf(stderr, "problem opening file %s for reading\n", argv[1]);
      exit(1);
    }
  }

  alg_table.allocated = 100;
  alg_table.length = 0;
  alg_table.tab = (rootprime_t *)malloc(alg_table.allocated*sizeof(rootprime_t));
  
  rat_table.allocated = 100;
  rat_table.length = 0;
  rat_table.tab = (unsigned long *)malloc(rat_table.allocated*sizeof(unsigned long));
  
  bad_primes.allocated = 100;
  bad_primes.length = 0;
  bad_primes.tab = (unsigned long *)malloc(bad_primes.allocated*sizeof(unsigned long));

  update_tables(file, &alg_table, &rat_table, &bad_primes);
  
  fprintf(stderr, "sorting...\n");
  qsort((void *)rat_table.tab, rat_table.length, sizeof(unsigned long), cmpprimes);
  uniqprimes(&rat_table);
  
  qsort((void *)bad_primes.tab, bad_primes.length, sizeof(unsigned long), cmpprimes);
  uniqprimes(&bad_primes);

  qsort((void *)alg_table.tab, alg_table.length, sizeof(rootprime_t), cmprootprimes);
  uniqrootprimes(&alg_table);

  rewind(file);


  rel_table = detect_singleton(file, &alg_table, &rat_table, &nb_rel);

  fclose(file);

#if 0
  printf("primes = \n");
  for (i = 0; i < rat_table.length; ++i)
    printf("%lx ", rat_table.tab[i]);
  printf("\nrootprimes = \n");
  for (i = 0; i < alg_table.length; ++i)
    printf("(%lx,%lx) ", alg_table.tab[i].prime, alg_table.tab[i].root);
#endif
  fprintf(stderr, "\nbadprimes = \n");
  for (i = 0; i < bad_primes.length; ++i)
    fprintf(stderr, "%lu ", bad_primes.tab[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "nb of rational primes = %d\n", rat_table.length);
  fprintf(stderr, "nb of algebraic primes = %d\n", alg_table.length);
  fprintf(stderr, "Total number of primes = %d\n",
          rat_table.length + alg_table.length);

  for (i = 0; i < nb_rel; ++i) {
    fprint_relation(stdout, rel_table[i]);
  }

  return 0;
}
