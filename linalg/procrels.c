#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "relation.h"


// Data structure for one algebraic prime and for a table of those.
typedef struct {
  unsigned long prime;
  unsigned long root;
} rootprime_t;

typedef struct {
  rootprime_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_rootprime_t;

// Data structure for a table of rational primes.
typedef struct {
  unsigned long *tab;
  unsigned long allocated;
  unsigned long length;
} tab_prime_t;

// Data strucutre for an (a,b) pair and for a table of those
typedef struct {
  long a;
  unsigned long b;
} abpair_t;

typedef struct {
  abpair_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_abpair_t;

// Table of relations.
typedef struct {
  relation_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_relation_t;



// Append an algebraic prime to a table.
void
add_alg_rootprime(tab_rootprime_t *table, tab_prime_t *bad_primes,
    unsigned long p,  unsigned long r) {
  if (r != -1) {
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

// Append a rational prime to a table.
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


// Create tables with primes that occur in relations.
// one structure for rational side, one for algebraic side, and one for
// bad primes (i.e. those that divide leading coefficient of polynomial).
void
update_tables(tab_relation_t *rel_table, tab_rootprime_t *alg_table, tab_prime_t *rat_table,
    tab_prime_t *bad_primes) {
  long a;
  unsigned long b;
  int i, j;
  relation_t *rel;

  for (i = 0; i < rel_table->length; ++i) {
    rel = &rel_table->tab[i];
    for (j = 0; j < rel->nb_rp; ++j)
      add_rat_prime(rat_table, rel->rp[j]);
    for (j = 0; j < rel->nb_ap; ++j)
      add_alg_rootprime(alg_table, bad_primes, rel->ap[j], rel->ar[j]);
  }
}

// Some utilities (for searching, sorting, etc) on algebraic and rational
// primes.
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


// This function does one pass of singleton removal.
// After having detected which primes occur only once (and the
// corresponding (a,b) pairs), those primes are deleted of the tables and
// the corresponding relations are removed from the table of relation.
void onepass_singleton_removal(tab_relation_t *rel_table, 
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

  for (i = 0; i < rel_table->length; ++i) {
    int j;
    relation_t *rel = &(rel_table->tab[i]);
    for (j = 0; j < rel->nb_rp; ++j) {
      int index;
      index = getindex_rat(*rat_table, rel->rp[j]);
      if (index == -1)
	continue;
      if ((ab_single.tab[index].a == 0) && (ab_single.tab[index].b == 0)) {
	ab_single.tab[index].a = rel->a;
	ab_single.tab[index].b = rel->b;
      } else if ((ab_single.tab[index].a != rel->a) || (ab_single.tab[index].b != rel->b)) {
	ab_single.tab[index].a = 1;
	ab_single.tab[index].b = 1;
      }
    }
    for (j = 0; j < rel->nb_ap; ++j) {
      int index;
      index = getindex_alg(*alg_table, rel->ap[j], rel->ar[j]);
      if (index == -1)
	continue;
      index += l_rat;
      if ((ab_single.tab[index].a == 0) && (ab_single.tab[index].b == 0)) {
	ab_single.tab[index].a = rel->a;
	ab_single.tab[index].b = rel->b;
      } else if ((ab_single.tab[index].a != rel->a) || (ab_single.tab[index].b != rel->b)) {
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
//	fprintf(stderr, "warning: this p is untouched ! %ld\n", i);
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
  //      fprintf(stderr, "warning: this p is untouched ! %ld\n", i);
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
    relation_t * ptab = rel_table->tab;
    j = 0;
    for (i = 0; i < rel_table->length; ++i) {
      if (getindex_abpair(ab_single, rel_table->tab[i].a, rel_table->tab[i].b) == -1) {
	copy_rel(&(ptab[j]), rel_table->tab[i]);
	j++;
      }
    }
    rel_table->length = j;
  }
  fprintf(stderr, "%d relations remain\n", rel_table->length);
  free(ab_single.tab);
}

void
remove_singletons(tab_relation_t *rel_table, tab_rootprime_t * alg_table,
    tab_prime_t * rat_table) {
  int old_nb;
  do {
    old_nb = rel_table->length;
    fprintf(stderr, "** Do one pass of singleton removal...\n");
    onepass_singleton_removal(rel_table, alg_table, rat_table);
  } while (rel_table->length != old_nb);
}

static int 
isBadPrime(unsigned long p, tab_prime_t bad_primes) {
  int i;
  for (i = 0; i < bad_primes.length; ++i) {
    if (p == bad_primes.tab[i])
      return 1;
  }
  return 0;
}


// Print relations in a matrix format:
// don't take into account bad primes and even powers of primes.
void 
fprint_rel_row(FILE *file, relation_t rel, tab_prime_t rat_table, tab_rootprime_t alg_table, tab_prime_t bad_primes) {
  int i;
  int *table_ind;
  int nb_coeff;
  int index, old_index;
  int parity;

  fprintf(file, "%ld %lu ", rel.a, rel.b);
  table_ind = (int*) malloc((rel.nb_rp + rel.nb_ap)*sizeof(int));
  
  nb_coeff = 0;
  index = getindex_rat(rat_table, rel.rp[0]);
  old_index = index;
  parity = 1;
  for (i = 1; i < rel.nb_rp; ++i) {
    index = getindex_rat(rat_table, rel.rp[i]);
    if (index == old_index) {
      parity = 1 - parity;
    } else {
      if (parity == 1) {
	table_ind[nb_coeff++] = old_index;
      }
      old_index = index;
      parity = 1;
    }
  }
  if (parity == 1)
    table_ind[nb_coeff++] = index;

  i = 0;
  while (isBadPrime(rel.ap[i], bad_primes))
    i++;
  index = getindex_alg(alg_table, rel.ap[i], rel.ar[i]);
  old_index = index;
  parity = 1;
  for (;i < rel.nb_ap; ++i) {
    if (isBadPrime(rel.ap[i], bad_primes))
      continue;
    index = getindex_alg(alg_table, rel.ap[i], rel.ar[i]);
    if (index == old_index) {
      parity = 1 - parity;
    } else {
      if (parity == 1) {
	table_ind[nb_coeff++] = old_index;
      }
      old_index = index;
      parity = 1;
    }
  }
  if (parity == 1)
    table_ind[nb_coeff++] = index;

  fprintf(file, "%lu ", nb_coeff);
  for (i = 0; i < nb_coeff; ++i) {
    assert (table_ind[i] != -1);
    fprintf(file, "%lu ", table_ind[i]);
  }

  fprintf(file, "\n");
  free(table_ind);
}

// Read all relations from file.
// NB: the NULL assignments to pointers are to help freeing.
int
fread_relations(FILE *file, tab_relation_t *rel_table) {
  int ret, i;
  rel_table->allocated = 100;
  rel_table->tab = (relation_t *) malloc(rel_table->allocated*sizeof(relation_t));
  rel_table->length = 0;
 
  do {
    if (rel_table->length == (rel_table->allocated - 1)) {
      rel_table->allocated += 100;
      rel_table->tab = (relation_t *) realloc(rel_table->tab, rel_table->allocated*sizeof(relation_t));
      for (i = rel_table->allocated-100; i < rel_table->allocated; ++i) {
	rel_table->tab[i].rp = NULL;
	rel_table->tab[i].ap = NULL;
	rel_table->tab[i].ar = NULL;
      }
    }
    ret = fread_relation(file, &(rel_table->tab[rel_table->length]));
    if (ret == 1)
      (rel_table->length)++;
  } while (ret == 1);

  if (ret == 0) {
    fprintf(stderr, "Warning: error when reading relation nb %d\n", rel_table->length);
  }

  fprintf(stderr, "loaded %d relations\n", rel_table->length);

  for (i = 0; i < rel_table->length; ++i) 
    computeroots(&(rel_table->tab[i]));
  
  return (ret == -1);
}


int main(int argc, char **argv) {
  tab_rootprime_t alg_table;
  tab_prime_t rat_table, bad_primes;
  tab_abpair_t ab_single;
  FILE *file;
  tab_relation_t rel_table;
  int ret;
  int i;

  if (argc == 1) {
    file = stdin;
    fprintf(stderr, "usage: %s [filename]\n", argv[0]);
    fprintf(stderr, "  stdin input is not yet available, sorry.\n");
    exit(1);
  } 
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

  alg_table.allocated = 100;
  alg_table.length = 0;
  alg_table.tab = (rootprime_t *)malloc(alg_table.allocated*sizeof(rootprime_t));
  
  rat_table.allocated = 100;
  rat_table.length = 0;
  rat_table.tab = (unsigned long *)malloc(rat_table.allocated*sizeof(unsigned long));
  
  bad_primes.allocated = 100;
  bad_primes.length = 0;
  bad_primes.tab = (unsigned long *)malloc(bad_primes.allocated*sizeof(unsigned long));

  fprintf(stderr, "reading file of relations...\n");
  ret = fread_relations(file, &rel_table);
  assert (ret);

  fprintf(stderr, "creating tables of primes...\n");
  update_tables(&rel_table, &alg_table, &rat_table, &bad_primes);
  
  fprintf(stderr, "sorting...\n");
  qsort((void *)rat_table.tab, rat_table.length, sizeof(unsigned long), cmpprimes);
  uniqprimes(&rat_table);
  
  qsort((void *)bad_primes.tab, bad_primes.length, sizeof(unsigned long), cmpprimes);
  uniqprimes(&bad_primes);

  qsort((void *)alg_table.tab, alg_table.length, sizeof(rootprime_t), cmprootprimes);
  uniqrootprimes(&alg_table);

  fprintf(stderr, "starting singleton removal...\n");
  remove_singletons(&rel_table, &alg_table, &rat_table);

  fprintf(stderr, "\nbadprimes = \n");
  for (i = 0; i < bad_primes.length; ++i)
    fprintf(stderr, "%lu ", bad_primes.tab[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "nb of rational primes = %d\n", rat_table.length);
  fprintf(stderr, "nb of algebraic primes = %d\n", alg_table.length);
  fprintf(stderr, "Total number of primes = %d\n",
          rat_table.length + alg_table.length);

  printf("%d %d\n", rel_table.length, rat_table.length + alg_table.length);
  for (i = 0; i < rel_table.length; ++i) {
    fprint_rel_row(stdout, rel_table.tab[i], rat_table, alg_table, bad_primes);
  }

  for (i = 0; i < rel_table.allocated; ++i) {
    if (rel_table.tab[i].rp != NULL)
      free(rel_table.tab[i].rp);
    if (rel_table.tab[i].ap != NULL)
      free(rel_table.tab[i].ap);
    if (rel_table.tab[i].ar != NULL)
      free(rel_table.tab[i].ar);
  }
  free(rel_table.tab);
  free(bad_primes.tab);
  free(rat_table.tab);
  free(alg_table.tab);

  return 0;
}
