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

#define MAXNAME 1024
#define MAX_THREADS 128
#define EXACT_HASH /* store (p,r) exactly in hash table */

/* thread structure */
typedef struct
{
  char g[MAXNAME];
  int thread;
  unsigned long nlines;
  unsigned long nideal;
  unsigned long nsingl;
  unsigned long ndupli;
  unsigned long nthreads;
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
    { ".gz", "gzip -dc %s", "gzip -c --fast > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    /* These two have to be present */
    { "", NULL, NULL },
    { NULL, NULL, NULL },
};

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* mutual exclusion lock */
long wct0;

uint64_t M = 0, minpa = 0, minpr = 0, Hsize, *H, nrels, remains, *PR;
uint64_t target = 0;
tab2_t *H2;
double threshold_weight = 999.0;

long
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  return rus.ru_utime.tv_sec;
}

long
realtime ()
{
  struct timeval tv[1];
  /* struct timeval {
     time_t      tv_sec;     seconds
     suseconds_t tv_usec;    microseconds
     }; */

  gettimeofday (tv, NULL);
  return tv->tv_sec;
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

#ifndef EXACT_HASH
static inline uint64_t
index_hash (uint64_t pr)
{
  return pr % M;
}
#else
static inline uint64_t
index_hash (uint64_t pr)
{
  uint64_t h = pr % M;

  while (PR[h] != 0 && PR[h] != pr)
    if (++h == M)
      h = 0;
  /* now PR[h] = 0 or PR[h] = pr */
  PR[h] = pr;
  return h;
}
#endif

static inline uint64_t
index_rat (uint64_t p)
{
  uint64_t h;

  h = index_hash (p + 1); /* even for even M */
  return h;
}

static inline uint64_t
index_alg (uint64_t p, uint64_t r)
{
  uint64_t h;

  h = index_hash (p + (M - 2) * r); /* always odd for odd p and even M */
  return h;
}

static inline uint64_t
weight (uint64_t h)
{
  uint64_t r, s;

  r = h >> 4;
  s = h & 15;
  return (H[r] >> (s << 2)) & 15;
}

static void
insert (uint64_t h)
{
  uint64_t r, s, v;
  const uint64_t mask[16] = {0xf, 0xf0, 0xf00, 0xf000, 0xf0000, 0xf00000,
                             0xf000000, 0xf0000000, 0xf00000000, 0xf000000000,
                             0xf0000000000, 0xf00000000000, 0xf000000000000,
                             0xf0000000000000, 0xf00000000000000,
                             0xf000000000000000};

  assert (h < M);
  r = h >> 4;         /* h div 16 */
  s = h & 15;         /* h mod 16 */
  v = H[r] & mask[s];
  if (v != mask[s])
    H[r] += (uint64_t) 1 << (s << 2); /* add 1 to bit 4s */
}

static double
clique_weight (uint64_t h)
{
  while (H2[h].pointer > 0)
    h = H2[h].pointer - 1;
  return -H2[h].weight;
}

static double
propagate_weight (uint64_t h)
{
  double w;

  if (H2[h].pointer > 0)
    {
      w = propagate_weight (H2[h].pointer - 1);
      H2[h].weight = w;
    }
  return H2[h].weight;
}

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
  while (H2[h].pointer > 0)
    h = H2[h].pointer - 1;
  if (h == h0) /* already connected */
    return;
  assert (H2[h0].pointer <= 0);
  H2[h0].weight += H2[h].weight;
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

  for (h = 0; h <= SAMPLE; h++)
    count[h] = 0;
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
  fprintf (stderr, "Number of connected components: %lu\n", nc);
  fprintf (stderr, "Weight: %.2f (avg), %.2f (max)\n", wsum / (double) nc, wmax);
  n = SAMPLE;
  while (n > 0)
    {
      count[n-1] += count[n];
      if (count[n-1] > 0)
        fprintf (stderr, "Weight >= %.2f: %lu components\n",
                 (double) (n - 1) / RESOLUTION, count[n-1]);
      if (count[n-1] > target)
        break;
      n--;
    }
  threshold_weight = (double) n / RESOLUTION;
  fprintf (stderr, "Will remove %lu components of weight >= %.2f\n",
           count[n], threshold_weight);

  /* propagate weights to speed up clique_weight */
  for (h = 0; h < M; h++)
    propagate_weight (h);
}

static void
delete (uint64_t h)
{
  uint64_t r, s, v;
  const uint64_t mask[16] = {0xf, 0xf0, 0xf00, 0xf000, 0xf0000, 0xf00000,
                             0xf000000, 0xf0000000, 0xf00000000, 0xf000000000,
                             0xf0000000000, 0xf00000000000, 0xf000000000000,
                             0xf0000000000000, 0xf00000000000000,
                             0xf000000000000000};

  assert (h < M);
  r = h >> 4;         /* h div 16 */
  s = h & 15;         /* h mod 16 */
  v = H[r] & mask[s];
  if (v == 0)
    fprintf (stderr, "Warning: deleting ideal with weight 0 at h=%lu\n", h);
  else if (v != mask[s])
    H[r] -= (uint64_t) 1 << (s << 2); /* remove 1 to bit 4s */
}

/* return a*b mod p, assuming 0 <= a, b, p < 2^40 */
uint64_t
mulmod (uint64_t a, uint64_t b, uint64_t p)
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

/* return a/b mod p */
static unsigned long
root (long a0, long b0, long p)
{
  long u, w, q, r, t, a, b;

  a = a0 % p;
  b = b0 % p;
  assert (b != 0);
  if (a < 0)
    a += p;

  u = 1;
  w = 0;
  a = b;
  b = p;
  /* invariant: a = u*b0 mod p, b = w*b0 mod p */
  while (b != 0)
    {
      q = a / b;
      r = a % b;
      a = b;
      b = r;
      t = u;
      u = w;
      w = (t - q * w) % p;
    }
  assert (a == 1);
  /* compute a0*u % p */
  if (u < 0)
    u += p;
  assert (u > 0);
  if (a0 >= 0)
    u = mulmod (a0, u, p);
  else
    u = p - mulmod (-a0, u, p);
  assert (u > 0);
  r = ((long) mulmod (u, b0, p) - a0) % p;
  assert (r == 0);
  return u;
}

/* first pass, return the number of lines read */
static unsigned long
foo (char *g)
{
  FILE *f;
  char *ret, s[1024], *t, *endptr;
  unsigned long line = 0, b, r;
  long a;
  uint64_t p, pfree;

  f = gzip_open (g, "r");

  while (1)
    {
      ret = fgets (s, 1024, f);
      if (ret == 0)
        break;
      line ++;
      t = s;
      a = 0;
      a = strtol (t, &endptr, 10);
      assert (a != 0);
      assert (endptr > t);
      t = endptr;
      assert (*t == ',');
      t ++;
      b = strtoul (t, &endptr, 10);
      assert (endptr > t);
      t = endptr;
      assert (*t == ':');
      t ++;
      pfree = 0;
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (b == 0)
            pfree = p;
          if (p >= minpr)
            insert (index_rat (p));
          if (*t++ != ',')
            break; /* end of rational primes */
        }
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (pfree >= minpa)
            insert (index_alg (pfree, p));
          else if (p >= minpa)
            {
              r = root (a, b, p);
              insert (index_alg (p, r));
            }
          if (*t++ != ',')
            break; /* end of algebraic primes */
        }
    }

  gzip_close (f, g);
  return line;
}

/* pass 2: insert in H2 ideals with weight 2 */
static void
bar2 (char *g)
{
  FILE *f;
  char *ret, s[1024], *t, *endptr;
  unsigned long line = 0, b, r, output = 0;
  long a;
  uint64_t p, pfree, w, h;
  int nr, na, keep;
  uint64_t rp[16], ap[16], ar[16];
  double W;

  f = gzip_open (g, "r");

  while (1)
    {
      ret = fgets (s, 1024, f);
      if (ret == 0)
        break;
      line ++;
      t = s;
      a = 0;
      a = strtol (t, &endptr, 10);
      // assert (a != 0);       /* already checked in pass 1 */
      // assert (endptr > t);   /* already checked in pass 1 */
      t = endptr;
      // assert (*t == ',');    /* already checked in pass 1 */
      t ++;
      b = strtoul (t, &endptr, 10);
      // assert (endptr > t);   /* already checked in pass 1 */
      t = endptr;
      // assert (*t == ':');    /* already checked in pass 1 */
      t ++;
      pfree = 0;
      nr = na = 0;
      W = 0.0;
      keep = 1;
      while (keep)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (b == 0)
            pfree = p;
          if (p >= minpr)
            {
              w = weight (index_rat (p));
              if (w <= 1)
                keep = 0;
              else if (w == 2)
                rp[nr++] = p;
              else if (w <= 15)
                W += 1.0 / (double) w;
            }
          if (*t++ != ',')
            break; /* end of rational primes */
        }
      while (keep)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (pfree >= minpa)
            {
              r = p;
              p = pfree;
            }
          else if (p >= minpa) /* non-free case */
            r = root (a, b, p);
          if (p >= minpa)
            {
              w = weight (index_alg (p, r));
              if (w <= 1)
                keep = 0;
              else if (w == 2)
                {
                  ap[na] = p;
                  ar[na++] = r;
                }
              else if (w <= 15)
                W += 1.0 / (double) w;
            }
          if (*t++ != ',')
            break; /* end of algebraic primes */
        }
      if (keep && (nr + na > 0)) /* at least one ideal of weight 2 */
        {
          /* Warning: it can be that an ideal appears twice (or more).
             To avoid reducing exponents mod 2, which would be expensive,
             we insert this ideal twice, but then we should avoid for infinite
             loops. */
          if (nr-- > 0)
            h = insert2 (index_rat (rp[nr]), W);
          else /* na > 0 */
            {
              na --;
              h = insert2 (index_alg (ap[na], ar[na]), W);
            }
          while (nr-- > 0)
            insert2a (index_rat (rp[nr]), h);
          while (na-- > 0)
            insert2a (index_alg (ap[na], ar[na]), h);
        }
      output += keep;
    }

  gzip_close (f, g);

  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
}

/* pass 3: writes remaining relations */
static void
bar (char *g)
{
  FILE *f, *newf;
  char *ret, s[1024], *t, *endptr, newg[1024];
  unsigned long line = 0, b, r, output = 0;
  long a;
  uint64_t p, pfree, w, h;
  int keep, nr, na;
  uint64_t rp[16], ap[16], ar[16];

  f = gzip_open (g, "r");
  strcpy (newg, "new/");
  strcpy (newg + 4, g);
  newf = gzip_open (newg, "w");
  assert (newf != NULL);

  while (1)
    {
      ret = fgets (s, 1024, f);
      if (ret == 0)
        break;
      line ++;
      t = s;
      a = 0;
      a = strtol (t, &endptr, 10);
      // assert (a != 0);          /* already checked in pass 1 */
      // assert (endptr > t);      /* already checked in pass 1 */
      t = endptr;
      // assert (*t == ',');       /* already checked in pass 1 */
      t ++;
      b = strtoul (t, &endptr, 10);
      // assert (endptr > t);      /* already checked in pass 1 */
      t = endptr;
      // assert (*t == ':');       /* already checked in pass 1 */
      t ++;
      pfree = 0;
      keep = 1;
      nr = na = 0;
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (b == 0)
            pfree = p;
          if (p >= minpr)
            {
              rp[nr++] = p;
              h = index_rat (p);
              w = weight (h);
              if (w <= 1)
                keep = 0;
              else if ((w == 2) && (clique_weight (h) >= threshold_weight))
                keep = 0;
            }
          if (*t++ != ',')
            break; /* end of rational primes */
        }
      while (1)
        {
          p = 0;
          while (('0' <= *t && *t <= '9') || ('a' <= *t && *t <= 'f'))
            {
              p = 16 * p;
              if ('0' <= *t && *t <= '9')
                p += *t - '0';
              else
                p += *t - 'a' + 10;
              t++;
            }
          if (pfree >= minpa)
            {
              r = p;
              p = pfree;
            }
          else if (p >= minpa) /* non-free case */
            r = root (a, b, p);
          if (p >= minpa)
            {
              ap[na] = p;
              ar[na++] = r;
              h = index_alg (p, r);
              w = weight (h);
              if (w <= 1)
                keep = 0;
              else if ((w == 2) && (clique_weight (h) >= threshold_weight))
                keep = 0;
            }
          if (*t++ != ',')
            break; /* end of algebraic primes */
        }
      if (keep)
        {
          fprintf (newf, "%s", s);
          output ++;
        }
      else /* update 2-bit counters of removed ideals */
        {
          while (nr-- > 0)
            delete (index_rat (rp[nr]));
          while (na-- > 0)
            delete (index_alg (ap[na], ar[na]));
        }
    }

  gzip_close (f, g);
  gzip_close (newf, newg);

  fprintf (stderr, "   new/%s done: remains %lu rels out of %lu\n",
           g, output, line);
  fflush (stderr);
  pthread_mutex_lock (&lock);
  remains += output;
  nrels += line;
  pthread_mutex_unlock (&lock);
}

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  fprintf (stderr, "1: Thread %d deals with %s\n", tab[0]->thread, tab[0]->g);
  fflush (stderr);
  tab[0]->nlines = foo (tab[0]->g);
  fprintf (stderr, "   %s done (%lu relations)\n", tab[0]->g, tab[0]->nlines);
  fflush (stderr);
  return NULL;
}

/* for pass2 */
void*
pass2_one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  fprintf (stderr, "2: Thread %d deals with %s\n", tab[0]->thread, tab[0]->g);
  fflush (stderr);
  bar2 (tab[0]->g);
  return NULL;
}

/* for pass3 */
void*
pass3_one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  fprintf (stderr, "3: Thread %d deals with %s\n", tab[0]->thread, tab[0]->g);
  fflush (stderr);
  bar (tab[0]->g);
  return NULL;
}

void*
one_thread3 (void* args)
{
  tab_t *tab = (tab_t*) args;
  uint64_t h, r, wt[16];

  for (r = 1; r < 16; r++)
    wt[r] = 0;
  for (r = tab[0]->thread; r < Hsize; r += tab[0]->nthreads)
    {  /* deal with 64 bits, i.e., 16 entries at a time */
      h = H[r];
      while (h != 0)
        {
          wt[h & 0xf] ++;
          h >>= 4;
        }
    }

  for (tab[0]->nideal = 0, r = 1; r < 16; r++)
    tab[0]->nideal += wt[r];
  tab[0]->nsingl = wt[1];
  tab[0]->ndupli = wt[2];

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
#ifdef EXACT_HASH
  return ek;
#else
  return -n * log (1.0 - ek / n);
#endif
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
  fprintf (stderr, "stat_mt() took %lds (cpu), %lds (real)\n",
           cputime () - st, realtime () - rt);
  est_ideals = estimate_ideals ((double) M, (double) nideal);
#ifndef EXACT_HASH
  fprintf (stderr, "Read %lu rels (%lu/s), %lu entries (%.2f%%), est. %1.0f ideals, %lu singletons, %lu of weight 2\n",
           nrels, nrels / (rt - wct0 + (rt == wct0)),
           nideal, 100.0 * (double) nideal / (double) M,
           est_ideals, nsingl, ndupli);
#else
  fprintf (stderr, "Read %lu rels (%lu/s), %lu ideals (%.2f%%), %lu singletons, %lu of weight 2\n",
           nrels, nrels / (rt - wct0 + (rt == wct0)),
           nideal, 100.0 * (double) nideal / (double) M,
           nsingl, ndupli);
#endif
  fprintf (stderr, "Estimated excess %1.0f, need %1.0f\n",
           (double) nrels - est_ideals,
           est_excess ((double) minpr) + est_excess ((double) minpa));
  fflush (stderr);
}

/* read all files and fills the hash table */
static void
doit (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], *ret;
  int i = 0, j;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  long st, rt;

  f = fopen (filelist, "r");
  if (f == NULL)
    {
      fprintf (stderr, "Error, cannot open %s\n", filelist);
      exit (1);
    }
  T = malloc (nthreads * sizeof (tab_t));
  while (!feof (f))
    {
      ret = fgets (g, MAXNAME, f);
      if (ret == NULL)
        break;
      g[strlen(g) - 1] = '\0';
      strcpy (T[i]->g, g);
      T[i]->thread = i;
      i ++;
      if (i == nthreads)
        {
          st = cputime ();
          rt = realtime ();
          for (i = 0; i < nthreads; i++)
            pthread_create (&tid[i], NULL, one_thread, (void *) (T + i));
          for (i = 0; i < nthreads; i++)
            {
              pthread_join (tid[i], NULL);
              nrels += T[i]->nlines;
            }
          fprintf (stderr, "Batch took %lds (cpu), %lds (real)\n",
                   cputime () - st, realtime () - rt);
          i = 0;
          stat_mt (nthreads);
        }
    }
  if (i > 0)
    {
      for (j = 0; j < i; j++)
        pthread_create (&tid[j], NULL, one_thread, (void *) (T + j));
      for (j = 0; j < i; j++)
        {
          pthread_join (tid[j], NULL);
          nrels += T[j]->nlines;
        }
      stat_mt (nthreads);
    }
  fclose (f);
  free (T);
}

/* second pass: stores ideals of weight 2 */
static void
pass2 (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], *ret;
  int i = 0, j;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  long st, rt;

  f = fopen (filelist, "r");
  T = malloc (nthreads * sizeof (tab_t));
  while (!feof (f))
    {
      ret = fgets (g, MAXNAME, f);
      if (ret == NULL)
        break;
      g[strlen(g) - 1] = '\0';
      strcpy (T[i]->g, g);
      T[i]->thread = i;
      i ++;
      if (i == nthreads)
        {
          st = cputime ();
          rt = realtime ();
          for (i = 0; i < nthreads; i++)
            pthread_create (&tid[i], NULL, pass2_one_thread, (void *) (T + i));
          for (i = 0; i < nthreads; i++)
            pthread_join (tid[i], NULL);
          fprintf (stderr, "Batch took %lds (cpu), %lds (real)\n",
                   cputime () - st, realtime () - rt);
          i = 0;
          fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
                   remains, nrels, 100.0 * (double) remains / (double) nrels);
          fflush (stderr);
        }
    }
  if (i > 0)
    {
      for (j = 0; j < i; j++)
        pthread_create (&tid[j], NULL, pass2_one_thread, (void *) (T + j));
      for (j = 0; j < i; j++)
        pthread_join (tid[j], NULL);
      fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
               remains, nrels, 100.0 * (double) remains / (double) nrels);
      fflush (stderr);
    }
  stat2 ();
  fclose (f);
  free (T);
}

/* third pass: outputs remaining relations */
static void
pass3 (int nthreads, char *filelist)
{
  FILE *f;
  char g[MAXNAME], *ret;
  int i = 0, j;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  long st, rt;

  f = fopen (filelist, "r");
  T = malloc (nthreads * sizeof (tab_t));
  while (!feof (f))
    {
      ret = fgets (g, MAXNAME, f);
      if (ret == NULL)
        break;
      g[strlen(g) - 1] = '\0';
      strcpy (T[i]->g, g);
      T[i]->thread = i;
      i ++;
      if (i == nthreads)
        {
          st = cputime ();
          rt = realtime ();
          for (i = 0; i < nthreads; i++)
            pthread_create (&tid[i], NULL, pass3_one_thread, (void *) (T + i));
          for (i = 0; i < nthreads; i++)
            pthread_join (tid[i], NULL);
          fprintf (stderr, "Batch took %lds (cpu), %lds (real)\n",
                   cputime () - st, realtime () - rt);
          i = 0;
          fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
                   remains, nrels, 100.0 * (double) remains / (double) nrels);
          fflush (stderr);
        }
    }
  if (i > 0)
    {
      for (j = 0; j < i; j++)
        pthread_create (&tid[j], NULL, pass3_one_thread, (void *) (T + j));
      for (j = 0; j < i; j++)
        pthread_join (tid[j], NULL);
      fprintf (stderr, "Total %lu rels remaining out of %lu (%.2f%%)\n",
               remains, nrels, 100.0 * (double) remains / (double) nrels);
      fflush (stderr);
    }
  fclose (f);
  free (T);
}

int
isprime (uint64_t n)
{
  uint64_t p;

  if ((n % 2) == 0)
    return 0;
  for (p = 3; p * p <= n; p += 3)
    if ((n % p) == 0)
      return 0;
  return 1;
}

int
main (int argc, char *argv[])
{
  int nthreads = 1, i;
  char *filelist = NULL;

  /* print command-line arguments */
  fprintf (stderr, "%s", argv[0]);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");
  fflush (stderr);

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
  assert ((M % 2) == 0);
  assert (isprime (M/2));

  Hsize = 1 + (M - 1) / 16; /* each entry uses 4 bits */
  H = (uint64_t*) malloc (Hsize * sizeof (uint64_t));
  memset (H, 0, Hsize * sizeof (uint64_t));
#ifdef EXACT_HASH
  PR = (uint64_t*) malloc (M * sizeof (uint64_t));
  memset (PR, 0, M * sizeof (uint64_t));
#endif

  wct0 = realtime ();

  nrels = 0;
  doit (nthreads, filelist);

  H2 = (tab2_t*) malloc (M * sizeof (tab2_t));
  memset (H2, 0, M * sizeof (tab2_t));
  nrels = remains = 0;
  pass2 (nthreads, filelist);

  nrels = remains = 0;
  pass3 (nthreads, filelist);

  free (H);
  free (H2);
#ifdef EXACT_HASH
  free (PR);
#endif

  return 0;
}
