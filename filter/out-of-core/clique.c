/* This program is an experimental implementation of out-of-core clique
   removal for CADO-NFS.

   Input: a set of relation files (might be compressed), each one on a single
          line in the "filelist" argument. Each file name can be in
          subdirectories like "foo/bar/toto.rels.gz".
   Output: each file with singletons removed will be written in the "new"
          directory (relative to where the program was started). This directory
          must already exist. In addition, if file names contain subdirectories
          those subdirectories must already exist. For example the output for
          the input file "foo/bar/toto.rels.gz" will be stored in
          "new/foo/bar/toto.rels.gz".

   Parameters:
   -M nnn : gives the size nnn of the hash table used; nnn should be 2*p where
            p is a prime. In addition for efficiency nnn should be at least
            two times larger than the maximal number of primes ideals, which
            is about 2^lpba/log(2^lpba) + 2^lpbr/log(2^lpbr).
   -minpa nnn : considers only algebraic ideals >= nnn
   -minpr nnn : considers only rational ideals >= nnn
   -t nnn     : uses nnn threads
   -remove nnn : the number of "cliques" to be removed
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <linux/limits.h>
#include <semaphore.h>
#include <libgen.h>
#include <limits.h> /* for CHAR_BIT */

#define NEW_DIR "new/"
#define MAX_FILES_PER_THREAD 1
#define MAX_STOCK_RAT_PER_LINE 16
#define MAX_STOCK_ALG_PER_LINE 32
#define MAXNAME 1024
#define MAX_THREADS 128
#define PR_TYPE uint32_t

/* thread structure */
typedef struct
{
  unsigned long nlines;
  unsigned long nideal;
  unsigned long nsingl;
  unsigned long ndupli;
  char g[ARG_MAX];
  uint8_t busy; /* 0 = free; 1 = work done; 2 = busy */
  uint8_t thread, nthreads;
} __tab_struct;
typedef __tab_struct tab_t[1];

typedef union
{
  double weight;
  int64_t pointer; /* h+1 if points to H2[h] */
} tab2_t;

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "antebuffer 24 %s|gzip -dc", "gzip -c --fast > %s", },
    { ".bz2", "antebuffer 24 %s|bzip2 -dc", "bzip2 -c --fast > %s", },
    { ".lzma", "lzma -dc  %s", "lzma -c -0 > %s", },
    /* These two have to be present */
    { "", "antebuffer 24 %s", NULL },
    { NULL, NULL, NULL },
};

#ifndef HAVE_SYNC_FETCH
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* mutual exclusion lock */
#endif
long wct0;

uint64_t M = 0, minpa = 0, minpr = 0, Hsize, *H, nrels, remains;
PR_TYPE *PR;
uint64_t target = 0;
tab2_t *H2;
double threshold_weight = 999.0;
sem_t sem_pt;

static const unsigned char ugly[256] = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

void
create_directories (char *filelist)
{
  static char odirg[MAXNAME];
  FILE *f;
  char g[MAXNAME], *pdirg, mkdirg[MAXNAME<<1];
  size_t lg;
  int ret;

  f = fopen (filelist, "r");
  while (!feof (f)) {
    if (fgets (g, MAXNAME, f) && (lg = strlen(g))) {
      g[lg-1] = 0;
      pdirg = dirname(g);
      if (strcmp(odirg,pdirg)) {
	strcpy(odirg,pdirg);
	strcpy(mkdirg, "mkdir -p ");
	strcat(mkdirg, NEW_DIR);
	strcat(mkdirg, pdirg);
	ret = system (mkdirg);
        if (ret == -1)
          {
            fprintf (stderr, "An error occurred during system() call\n");
            exit (1);
          }
      }
    }
  }
  fclose(f);
}

double
compute_delay_stat ()
{
  uint64_t h;
  uint32_t i;

  for (h = (Hsize>>24), i = 0; h; h >>=1, i++);
  return (i) ? (double) (i<<6) : 32.;
}

double
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  return rus.ru_utime.tv_sec + rus.ru_utime.tv_usec / 1000000.0;
}

double
realtime ()
{
  struct timeval tv[1];
  /* struct timeval {
     time_t      tv_sec;     seconds
     suseconds_t tv_usec;    microseconds
     }; */

  gettimeofday (tv, NULL);
  return tv->tv_sec + tv->tv_usec / 1000000.0 ;
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

FILE*
fopen_compressed_r (const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_in) { /* this suffix has an associated deflate command */
            char * command;
            int ret = asprintf(&command, r->pfmt_in, name);
            assert(ret >= 0);
            f = popen(command, "r");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            return f;
        } else {
            f = fopen(name, "r");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE*
fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_out) { /* suffix has an associated compression command */
            char * command;
            int ret = asprintf(&command, r->pfmt_out, name);
            assert (ret >= 0);
            f = popen(command, "w");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            return f;
        } else {
            f = fopen(name, "w");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE * gzip_open(const char * name, const char * mode)
{
    if (strcmp(mode, "r") == 0) {
        return fopen_compressed_r (name, NULL, NULL);
    } else if (strcmp(mode, "w") == 0) {
        return fopen_compressed_w (name, NULL, NULL);
    } else {
        abort();
    }
}

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

void gzip_close(FILE * f, const char * name)
{
    const char * e = name + strlen(name);
    const struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        const char * ne = e - MIN(strlen(name), strlen(r->suffix));
        if (strcmp(r->suffix, ne) != 0)
            continue;
        if (r->pfmt_out) {
            pclose(f);
            return;
        } else {
            fclose(f);
            return;
        }
    }
}

static inline uint64_t
index_hash (uint64_t pr)
{
  PR_TYPE *prh;

  prh = PR + (pr % M);
  while (1) {
    if (*prh == (PR_TYPE) pr)
      return (prh - PR);
    if (!(*prh)) {
      *prh = (PR_TYPE) pr;
      return (prh - PR);
    }
    if (++prh == PR + M)
      prh = PR;
  }
}

#define INDEX_RAT(P) (index_hash ((P) + 1)) /* even for even M */
static inline uint64_t
index_rat (uint64_t p)
{
  return INDEX_RAT(p);
}

#define INDEX_ALG(P,R) (index_hash ((P)+ (M - 2) * (R)))  /* always odd for odd p and even M */
static inline uint64_t
index_alg (uint64_t p, uint64_t r)
{
  return INDEX_ALG(p,r);
}

#define WEIGHT(P) ((H[(P)>>4] >> (((P)&15)<<2)) & 15)
static inline uint64_t
weight (uint64_t h)
{
  return WEIGHT(h);
}

static void
insert (uint64_t h)
{
  uint8_t *hh, d, m, a;
  assert (h < M);
  hh = ((uint8_t*) H) + (h >> 1);
  d = (h & 1) << 2;
  m = 15 << d;
  a = 1 << d;
  if ((*hh & m) != m) *hh += a;
}

static void
delete (uint64_t h)
{
  uint8_t *hh, d, m, a, v;
  assert (h < M);
  static int count = 0;

  hh = ((uint8_t*) H) + (h >> 1);
  d = (h & 1) << 2;
  m = 15 << d;
  a = 1 << d;
  v = *hh & m;
  if (!v)
    {
      if (count++ < 10)
        {
          fprintf (stderr, "Warning: deleting weight-0 ideal 0 at h=%lu\n", h);
        }
    }
  else
    if (v != m)
      *hh -= a; /* remove 1 to bit 4s */
}

static double
clique_weight (uint64_t h)
{
  while (H2[h].pointer > 0)
    h = H2[h].pointer - 1;
  return -H2[h].weight;
}

#if 0
/* Warning: recursive procedure calls might explose the stack */
static double
propagate_weight (uint64_t h)
{
  if (H2[h].pointer > 0)
    H2[h].weight = propagate_weight (H2[h].pointer - 1);
  return H2[h].weight;
}
#else
static void
propagate_weight ()
{
  tab2_t *cH;
  tab2_t **pro_wgt = NULL, **ppro_wgt = NULL, **max_pro_wgt = NULL;
  uint64_t nb_pro_wgt, h;

  for (cH = H2; cH < H2 + M; cH++)
    if (cH->pointer > 0) {
      h = cH - H2;
      do {
	if (ppro_wgt >= max_pro_wgt) {
	  if (!ppro_wgt) {
	    nb_pro_wgt = 512;
	    ppro_wgt = pro_wgt = malloc(nb_pro_wgt * sizeof(*pro_wgt));
	  }
	  else {
	    nb_pro_wgt <<= 1;
	    pro_wgt = realloc(pro_wgt, nb_pro_wgt * sizeof(*pro_wgt));
	  }
	  max_pro_wgt = &(pro_wgt[nb_pro_wgt]);
	}
	*ppro_wgt++ = &(H2[h]);
	h = H2[h].pointer - 1;
      } while (H2[h].pointer > 0);
      while (ppro_wgt > pro_wgt)
	(*--ppro_wgt)->weight = H2[h].weight;
    }
  if (pro_wgt)
    free(pro_wgt);
}
#endif

/* adds the weight w to the connected component of H2[h] */
static uint64_t
insert2 (uint64_t h, double w)
{
  while (H2[h].pointer > 0)
    h = H2[h].pointer - 1;
  H2[h].weight -= w;
  return h;
}

/* connects the connected component H2[h] to H2[h0],
   where H2[h0] stores a weight */
static void
insert2a (uint64_t h, uint64_t h0)
{
  double w0, w;

  while (H2[h].pointer > 0)
    h = H2[h].pointer - 1;
  if (h == h0) /* already connected */
    return;
  w0 = H2[h0].weight;
  w = H2[h].weight;
  if (w0 > 0 || w > 0)
    return; /* avoid an assert which might fail */
  H2[h0].weight = w0 + w;
  H2[h].pointer = h0 + 1;
}

#define SAMPLE 1000       /* number of intervals for component weight */
#define RESOLUTION 100.0  /* 1/width of each interval */

static void
stat2 ()
{
  uint64_t h, nc = 0; /* number of connected components */
  double wmax = 0.0, wsum = 0.0, w;
  uint64_t count[SAMPLE+1], n;

  memset(count, 0, sizeof(*count) * (SAMPLE + 1));
  for (h = 0; h < M; h++)
    if (H2[h].weight < 0)
      {
        nc ++;
        w = -H2[h].weight;
        wsum += w;
        if (w > wmax)
          wmax = w;
        n = (uint64_t) (RESOLUTION * w);
        if (n >= SAMPLE)
          n = SAMPLE - 1;
        count[n] ++;
      }
  fprintf (stderr,
	   "Number of connected components: %lu\n"
	   "Weight: %.2f (avg), %.2f (max)\n", nc, wsum / (double) nc, wmax);
  n = SAMPLE;
  while (n > 0)
    {
      count[n-1] += count[n];
      if (count[n] > 0)
        fprintf (stderr, "Weight >= %.2f: %lu components\n",
                 (double) (n - 1) / RESOLUTION, count[n-1]);
      if (count[n-1] > target)
        break;
      n--;
    }
  threshold_weight = (double) n / RESOLUTION;
  fprintf (stderr, "Will remove %lu components of weight >= %.2f\n",
           count[n], threshold_weight);
}

/* return a*b mod p, assuming 0 <= a, b, p < 2^40 */
uint64_t
mulmod_old (uint64_t a, uint64_t b, uint64_t p)
{
  uint64_t a1, a0, b1, b0, r;

  a1 = a >> 20;
  a0 = a & 1048575;
  b1 = b >> 20;
  b0 = b & 1048575;
  r = (a1 * b1) % p; /* r < 2^40 */
  r = ((r << 20) + a1 * b0 + a0 * b1) % p;
  r = ((r << 20) + a0 * b0) % p;
  return r;
}

/* return a*b mod p, assuming 0 <= a, b, p < 2^40 */
uint64_t
mulmod (uint64_t a, uint64_t b, uint64_t p)
{
  uint64_t r;

  assert (a <= 0x100000000000);
  r = ((b >> 20) * a) % p; /* r < 2^40 */
  return ((r << 20) + (b & 0xfffff) * a) % p;
}

/* return a/b mod p */
static unsigned long
root (long a0, long b0, long p)
{
  long u, w, q, r, t, a, b;

  b = b0 % p;
  assert (b != 0);
  a = b;
  u = 1;
  w = 0;
  b = p;
  /* invariant: a = u*b0 mod p, b = w*b0 mod p */
  while (b) {
      q = a / b;
      r = a - q * b; /* r = a % b; */
      a = b;
      b = r;
      t = u;
      u = w;
      w = (t - q * w) % p;
  }
  /* compute a0*u % p */
  if (u < 0) u += p;
  return (a0 >= 0) ? mulmod (a0, u, p) : p - mulmod (-a0, u, p);
}

/* first pass, return the number of lines read */
static unsigned long
load_pass1 (char *g)
{
  char s[1024], *t;
  FILE *f;
  unsigned long line, b, r;
  long a;
  uint64_t p, p2;
  uint8_t m;


  f = gzip_open (g, "r");
  line = 0;
  while (fgets (s, 1024, f)) {
    line++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = p * 10 + m);
    a *= p;
    assert (a != 0);
    assert (*t == ',');
    t ++;
    for (b = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, b = b * 10 + m);
    assert (*t == ':');
    t ++;
    do {
      for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = (p << 4) + m);
      if (p >= minpr) {
	insert (INDEX_RAT (p));
      }
    } while (*t++ == ',');
    if (b) {
      do {
	for (p = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p = (p << 4) + m);
	if (p >= minpa) {
	  r = root (a, b, p);
	  insert (INDEX_ALG (p, r));
	}
      } while (*t++ == ',');
    }
    else
      if (p >= minpa)
	do {
	  for (p2 = 0; (m = ugly[*((uint8_t *) t)]) != 255; t++, p2 = (p2 << 4) + m);
	  insert (INDEX_ALG (p, p2));
	} while (*t++ == ',');
  }
  gzip_close (f, g);
  return line;
}

/* pass 2: insert in H2 ideals with weight 2 */
static void
load_pass2 (char *g)
{
  uint64_t rp[MAX_STOCK_RAT_PER_LINE], ap[MAX_STOCK_ALG_PER_LINE], ar[MAX_STOCK_ALG_PER_LINE];
  char s[1024], *t;
  FILE *f;
  unsigned long line, output, b, r;
  long a;
  uint64_t p, p2, i, w, h;
  double W;
  unsigned int nr, na;
  uint8_t m;

  f = gzip_open (g, "r");
  line = output = 0;
 next_line:
  while (fgets (s, 1024, f)) {
    line ++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = p * 10 + m);
    a *= p;
    t ++;
    for (b = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, b = b * 10 + m);
    t ++;
    nr = na = 0;
    W = 0.0;
    do {
      for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
      if (p >= minpr) {
	i = INDEX_RAT(p);
	w = WEIGHT (i);
	if (w < 2) goto next_line;
	if (w == 2) {
	  assert(nr < MAX_STOCK_RAT_PER_LINE);
	  rp[nr++] = p;
	}
	else
	  W += 1.0 / (double) w;
      }
    } while (*t++ == ',');
    if (b)
      do {
	for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
	if (p >= minpa) {
	  r = root (a, b, p);
	  i = INDEX_ALG (p, r);
	  w = WEIGHT (i);
	  if (w < 2) goto next_line;
	  if (w == 2) {
	    assert (na < MAX_STOCK_ALG_PER_LINE);
	    ap[na] = p;
	    ar[na++] = r;
	  } else
	    W += 1.0 / (double) w;
	}
      } while (*t++ == ',');
    else
      if (p >= minpa)
	do {
	  for (p2 = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p2 = (p2 << 4) + m);
	  i = INDEX_ALG (p, p2);
	  w = WEIGHT (i);
	  if (w < 2) goto next_line;
	  if (w == 2) {
	    assert (na < MAX_STOCK_ALG_PER_LINE);
	    ap[na] = p;
	    ar[na++] = p2;
	  } else
	    W += 1.0 / (double) w;
	} while (*t++ == ',');
    output++;
    /* Warning: it can be that an ideal appears twice (or more).
       To avoid reducing exponents mod 2, which would be expensive,
       we insert this ideal twice, but then we should avoid for infinite
       loops. */
    if (nr-- > 0) {
      h = insert2 (INDEX_RAT (rp[nr]), W);
      while (nr-- > 0) insert2a (INDEX_RAT (rp[nr]), h);
      while (na-- > 0) insert2a (INDEX_ALG (ap[na], ar[na]), h);
    }
    else if (na-- > 0) {
      h = insert2 (INDEX_ALG (ap[na], ar[na]), W);
      while (na-- > 0) insert2a (INDEX_ALG (ap[na], ar[na]), h);
    }
  }
  gzip_close (f, g);
#ifdef HAVE_SYNC_FETCH
  __sync_fetch_and_add(&remains, output);
  __sync_fetch_and_add(&nrels, line);
#else
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
#endif
}

/* pass 3: writes remaining relations */
static void
load_pass3 (char *g)
{
  uint64_t rp[MAX_STOCK_RAT_PER_LINE], ap[MAX_STOCK_ALG_PER_LINE], ar[MAX_STOCK_ALG_PER_LINE];
  char s[1024], newg[MAXNAME<<1], *t;
  FILE *f, *newf;
  unsigned long line = 0, b, r, output = 0;
  long a;
  uint64_t p, p2, w, h;
  unsigned int nr, na;
  uint8_t m;

  strcpy (newg, NEW_DIR);
  strcat (newg, g);
  f = gzip_open (g, "r");
  newf = gzip_open (newg, "w");
  assert (newf != NULL);
  line = output = 0;
  if (0) {
  next_line: /* update 2-bit counters of removed ideals of the suppress line */
    while (nr-- > 0)
      delete (INDEX_RAT (rp[nr]));
    while (na-- > 0)
      delete (INDEX_ALG (ap[na], ar[na]));
  }
  while (fgets (s, 1024, f)) {
    line ++;
    t = s;
    if (*t != '-')
      a = 1;
    else {
      a = -1;
      t++;
    }
    for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = p * 10 + m);
    a *= p;
    t ++;
    for (b = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, b = b * 10 + m);
    t ++;
    nr = na = 0;
    do {
      for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
      if (p >= minpr) {
	assert(nr < MAX_STOCK_RAT_PER_LINE);
	rp[nr++] = p;
	h = INDEX_RAT (p);
	w = WEIGHT (h);
	if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	  goto next_line;
      }
    } while (*t++ == ',');
    if (b)
      do {
	for (p = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p = (p << 4) + m);
	if (p >= minpa) {
	  assert (na < MAX_STOCK_ALG_PER_LINE);
	  r = root (a, b, p);
	  ap[na] = p;
	  ar[na++] = r;
	  h = INDEX_ALG (p, r);
	  w = WEIGHT (h);
	  if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	    goto next_line;
	}
      } while (*t++ == ',');
    else
      if (p >= minpa)
	do {
	  assert (na < MAX_STOCK_ALG_PER_LINE);
	  for (p2 = 0; ((m = ugly[*((uint8_t *) t)]) != 255); t++, p2 = (p2 << 4) + m);
	  ap[na] = p;
	  ar[na++] = p2;
	  h = INDEX_ALG (p, p2);
	  w = WEIGHT (h);
	  if (w < 2 || ((w == 2) && (clique_weight (h) >= threshold_weight)))
	    goto next_line;
	} while (*t++ == ',');
    fputs (s, newf);
    output ++;
  }
  gzip_close (f, g);
  gzip_close (newf, newg);
  fprintf (stderr, "   new/%s done: remains %lu rels out of %lu\n",
           g, output, line);
#ifdef HAVE_SYNC_FETCH
  __sync_fetch_and_add(&remains, output);
  __sync_fetch_and_add(&nrels, line);
#else
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
#endif
}

void*
pass1_one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  tab[0]->nlines = load_pass1 (tab[0]->g);
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  return NULL;
}

/* for pass2 */
void*
pass2_one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  load_pass2 (tab[0]->g);
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  return NULL;
}

/* for pass3 */
void*
pass3_one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  load_pass3 (tab[0]->g);
  tab[0]->busy = 1;
  sem_post(&sem_pt);
  return NULL;
}

void*
one_thread3 (void* args)
{
  tab_t *tab = (tab_t*) args;
  uint64_t *hd, *he, h, wt[16];

  memset (wt, 0, sizeof(*wt) * 16);
  hd = &(H[(tab[0]->thread * Hsize) / tab[0]->nthreads]);
  he = &(H[((tab[0]->thread + 1) * Hsize) / tab[0]->nthreads]);
  while (hd < he)
    {  /* deal with 64 bits, i.e., 16 entries at a time */
      h = *hd++;
      while (h)
        {
          wt[h & 0xf] ++;
          h >>= 4;
        }
    }
  tab[0]->nideal = tab[0]->nsingl = wt[1];
  tab[0]->ndupli = wt[2];
  for (hd = &(wt[2]); hd < &(wt[16]); tab[0]->nideal += *hd++);
  return NULL;
}

/* given a hash table of n entries, of which ek entries are already filled,
   return an estimate of the number k of insertions in the table.
   We assume that e(k+1) = e(k) + 1 - e(k)/n with e(0) = 1,
   which gives e(k) = n - n*(1-1/n)^k,
   thus k = log(1-e(k)/n)/log(1-1/n) */
static double
estimate_ideals (double n, double ek)
{
  return ek;
}

/* estimate the number of prime ideals below a */
static double
est_excess (double a)
{
  return a / log (a);
}

/* multi-thread version */
static void
stat_mt (int nthreads)
{
  int i;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  unsigned long nideal = 0, nsingl = 0, ndupli = 0;
  long st = cputime (), rt = realtime ();
  double est_ideals;

  T = malloc (nthreads * sizeof (tab_t));
  for (i = 0; i < nthreads; i++)
    {
      T[i]->thread = i;
      T[i]->nthreads = nthreads;
      pthread_create (&tid[i], NULL, one_thread3, (void *) (T + i));
    }
  for (i = 0; i < nthreads; i++)
    {
      pthread_join (tid[i], NULL);
      nideal += T[i]->nideal;
      nsingl += T[i]->nsingl;
      ndupli += T[i]->ndupli;
    }
  free (T);
  fprintf (stderr, "   Stat_mt() took %lds (cpu), %lds (real)\n",
           (unsigned long) (cputime () - st), (unsigned long) (realtime () - rt));
  est_ideals = estimate_ideals ((double) M, (double) nideal);
  fprintf (stderr, "   Read %lu rels (%lu/s), %lu ideals (%.2f%%), %lu singletons, %lu of weight 2\n",
           nrels, (unsigned long) (nrels / (rt - wct0)),
           nideal, (100.0 * nideal) / M, nsingl, ndupli);
  fprintf (stderr, "   Estimated excess %1.0f, need %1.0f\n",
           (double) nrels - est_ideals,
           est_excess ((double) minpr) + est_excess ((double) minpa));
}

/* read all files and fills the hash table */
static void
doit (int nthreads, char *filelist)
{
  FILE *f;
  char g[ARG_MAX], *pg;
  int i, j = 0, nbt = 0, notload, nbf, askstat;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  double st = cputime(), rt = realtime(), crt, art,
    delay_stat = compute_delay_stat();
  size_t lpg;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 1/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  if (f == NULL)
    {
      fprintf (stderr, "Error, cannot open %s\n", filelist);
      exit (1);
    }
  T = malloc (nthreads * sizeof (tab_t));
  memset(T, 0, nthreads * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads - 1);
  crt = art = rt;
  askstat = 0;
  while ((notload = (!feof (f))) || nbt) {
    if (notload && !askstat) {
      pg = g;
      nbf = 0;
      do {
	if (fgets (pg, MAXNAME, f) && (lpg = strlen(pg))) {
	  nbf++;
	  pg += lpg;
	  pg[-1] = ' ';
	}
	else break;
      } while (nbf < MAX_FILES_PER_THREAD && pg < g + ARG_MAX - MAXNAME);
      if (nbf) {
	pg[-1] = 0;
	for (i = 0; T[i]->busy && i < nthreads; i++);
	assert(i < nthreads);
	strcpy (T[i]->g, g);
	T[i]->busy = 2;
	T[i]->thread = i;
	nbt++;
	pthread_create (&tid[i], NULL, pass1_one_thread, (void *) (T + i));
      }
    }
    sem_wait(&sem_pt);
    for (i = 0; i < nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	nrels += T[i]->nlines;
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j + 1 == nthreads) {
      fprintf (stderr, "Pass1; rels: load %lu; krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (cputime() - st),
	       (unsigned long) (art - rt));
      j = 0;
    }
    if ((askstat = ((art - crt) > delay_stat)) && !nbt) {
      askstat = 0;
      stat_mt (nthreads);
      crt = art;
      sem_destroy(&sem_pt);
      sem_init(&sem_pt, 0, nthreads - 1);
    }
  }
  fprintf (stderr, "\n*** Pass 1, final stats: *** \n");
  stat_mt (nthreads);
  fprintf (stderr, "Pass 1 took %lds (cpu), %lds (real)\n",
           (unsigned long) (cputime () - st), (unsigned long) (realtime () - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 1/3 DONE *\n"
	   "*****************\n\n");
}

/* second pass: stores ideals of weight 2 */
static void
pass2 (int nthreads, char *filelist)
{
  FILE *f;
  char g[ARG_MAX], *pg;
  int i, j = 0, notload, nbt = 0, nbf;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  double st = cputime(), rt = realtime(), art;
  size_t lpg;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 2/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  T = malloc (nthreads * sizeof (tab_t));
  memset(T, 0, nthreads * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads - 1);
  art = rt;
  while ((notload = (!feof (f))) || nbt) {
    if (notload) {
      pg = g;
      nbf = 0;
      do {
	if (fgets (pg, MAXNAME, f) && (lpg = strlen(pg))) {
	  nbf++;
	  pg += lpg;
	  pg[-1] = ' ';
	}
	else
	  break;
      } while (nbf < MAX_FILES_PER_THREAD && pg < g + ARG_MAX - MAXNAME);
      if (nbf) {
	pg[-1] = 0;
	for (i = 0; T[i]->busy && i < nthreads; i++);
	assert(i < nthreads);
	strcpy (T[i]->g, g);
	T[i]->busy = 2;
	T[i]->thread = i;
	nbt++;
	pthread_create (&tid[i], NULL, pass2_one_thread, (void *) (T + i));
      }
    }
    sem_wait(&sem_pt);
    for (i = 0; i < nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j + 1 == nthreads) {
      fprintf (stderr, "Pass2; rels: load %lu, used %lu (%.2f%%); "
	       "krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels, remains,
	       (100.0 * remains) / nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (cputime() - st),
	       (unsigned long) (art - rt));
      j = 0;
    }
  }
  art = realtime();
  fprintf (stderr, "\n*** Pass 2, final stats: *** \n"
	   "Relations: load %lu, used %lu (%.2f%%); "
	   "krels/s: %lu; time: %lus cpu, %lus real\n",
	   nrels, remains,
	   (100.0 * remains) / nrels,
	   (unsigned long) (nrels / ((art - rt) * 1000)),
	   (unsigned long) (cputime() - st),
	   (unsigned long) (art - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  stat2 ();
  /* propagate weights to speed up clique_weight */
#if 0
  {
    uint64_t h;
    for (h = 0; h < M; h++)
      propagate_weight (h);
  }
#else
  propagate_weight ();
#endif
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 2/3 DONE *\n"
	   "*****************\n\n");
}

/* third pass: outputs remaining relations */
/* Problem: I cannot agglomerate the files because 1 file load = 1 file written.
   So I have to take them one by one and write them synchronously, one by one too.
   It's a pity for disk access...
*/
static void
pass3 (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME];
  int i, j = 0, k = 0, nbt = 0, notload;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  double st = cputime(), rt = realtime(), art;
  size_t lg;

  fprintf (stderr, "\n"
	   "************\n"
	   "* PASS 3/3 *\n"
	   "************\n\n");
  f = fopen (filelist, "r");
  T = malloc (nthreads * sizeof (tab_t));
  memset(T, 0, nthreads * sizeof (tab_t));
  sem_init(&sem_pt, 0, nthreads - 1);
  art = rt;
  while ((notload = (!feof (f))) || nbt) {
    if (notload && fgets (g, MAXNAME, f) && (lg = strlen(g))) {
      g[lg - 1] = 0;
      for (i = 0; T[i]->busy && i < nthreads; i++);
      assert(i < nthreads);
      strcpy (T[i]->g, g);
      T[i]->busy = 2;
      T[i]->thread = i;
      nbt++;
      pthread_create (&tid[i], NULL, pass3_one_thread, (void *) (T + i));
    }
    sem_wait(&sem_pt);
    for (i = 0; i < nthreads; i++)
      if (T[i]->busy == 1) {
	pthread_join (tid[i], NULL);
	T[i]->busy = 0;
	nbt--;
	j = i;
	break;
      }
    art = realtime();
    if (j + 1 == nthreads && ++k == MAX_FILES_PER_THREAD) {
      fprintf (stderr, "Pass3; rels: load %lu, used %lu (%.2f%%); "
	       "krels/s: %lu; time: %lus cpu, %lus real\n",
	       nrels, remains,
	       (100.0 * remains) / nrels,
	       (unsigned long) (nrels / ((art - rt) * 1000)),
	       (unsigned long) (cputime() - st),
	       (unsigned long) (art - rt));
      j = k = 0;
    }
  }
  fprintf (stderr, "Pass 3 took %lds (cpu), %lds (real)\n",
           (unsigned long) (cputime () - st), (unsigned long) (realtime () - rt));
  sem_destroy(&sem_pt);
  fclose (f);
  free (T);
  fprintf (stderr, "\n"
	   "*****************\n"
	   "* PASS 3/3 DONE *\n"
	   "*****************\n\n");
}

int
isprime (uint64_t n)
{
  uint64_t p;

  if (!(n & 1)) return 0;
  for (p = 3; p * p <= n; p += 3)
    if (!(n % p))
      return 0;
  return 1;
}

int
main (int argc, char *argv[])
{
  int nthreads = 1, i;
  char *filelist = NULL;
  size_t size_malloc;

  /* print command-line arguments */
  fprintf (stderr, "%s", argv[0]);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-M") == 0)
        {
          M = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-minpa") == 0)
        {
          minpa = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-minpr") == 0)
        {
          minpr = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-t") == 0)
        {
          nthreads = atoi (argv[2]);
          assert (nthreads <= MAX_THREADS);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-remove") == 0)
        {
          target = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-filelist") == 0)
        {
          filelist = argv[2];
          argc -= 2;
          argv += 2;
        }
      else if (argv[1][0] == '-')
        {
          fprintf (stderr, "Unknown option: %s\n", argv[1]);
          exit (1);
        }
    }

  if (M == 0)
    {
      fprintf (stderr, "Missing option -M\n");
      exit (1);
    }

  if (minpa == 0)
    {
      fprintf (stderr, "Missing option -minpa\n");
      exit (1);
    }

  if (minpr == 0)
    {
      fprintf (stderr, "Missing option -minpr\n");
      exit (1);
    }

  if (target == 0)
    {
      fprintf (stderr, "Missing option -remove\n");
      exit (1);
    }

  if (filelist == NULL)
    {
      fprintf (stderr, "Missing option -filelist\n");
      exit (1);
    }

  if (argc != 1)
    {
      fprintf (stderr, "Unknown command-line argument: %s\n", argv[1]);
      exit (1);
    }

  /* check that M = 2p where p is prime */
  assert (!(M & 1));
  assert (isprime (M>>1));

  Hsize = 1 + ((M - 1) >> 4); /* each entry uses 4 bits */

  size_malloc = Hsize * sizeof (*H);
  fprintf (stderr, "Creating & clearing hashtable: %lu mB...", size_malloc >> 20);
  H = malloc (size_malloc);
  memset (H, 0, size_malloc);
  fprintf (stderr, " Done. \n");
  size_malloc = M * sizeof (*PR);
  fprintf (stderr, "Creating & clearing (p,r) array: %lu mB (%lu bits per entry)...", size_malloc >> 20, sizeof(PR_TYPE) * CHAR_BIT);
  PR = malloc (size_malloc);
  memset (PR, 0, size_malloc);
  fprintf (stderr, " Done.\n");
  fprintf (stderr, "Creating possible " NEW_DIR "* directories...");
  create_directories(filelist);
  fprintf (stderr, " Done.\n");

  fprintf (stderr, "Each thread processing %d file(s)\n",
           MAX_FILES_PER_THREAD);

  wct0 = realtime ();
  nrels = 0;
  doit (nthreads, filelist);

  size_malloc = M * sizeof (*H2);
  fprintf (stderr, "Creating & clearing weight primes array: %lu mB...", size_malloc >> 20);
  H2 = malloc (size_malloc);
  memset (H2, 0, size_malloc);
  fprintf (stderr, " Done.\n");

  nrels = remains = 0;
  pass2 (nthreads, filelist);

  nrels = remains = 0;
  pass3 (nthreads, filelist);

  free (H);
  free (H2);
  free (PR);

  return 0;
}
