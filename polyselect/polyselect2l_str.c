/* Data struct used for polyselect2l */
#include "polyselect2l_str.h"

void match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
            uint64_t ad, unsigned long d, mpz_t N, unsigned long q,
            mpz_t rq);

void gmp_match (uint32_t p1, uint32_t p2, int64_t i, mpz_t m0,
		uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
		mpz_t rq);

/* LEN_SPECIAL_Q in the header */
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 0 };

//#define LESS_P

static inline uint64_t cputicks()
{
        uint64_t r;
        __asm__ __volatile__(
                "rdtsc\n\t"
                "shlq $32, %%rdx\n\t"
                "orq %%rdx, %%rax\n\t"
                : "=a"(r)
                :
                : "rdx");
        return r;
}


/* init prime array */
unsigned long
initPrimes ( unsigned long P,
             uint32_t **primes )
{
  unsigned long p, nprimes = 0;

#ifdef LESS_P // if impatient for root finding
  unsigned long maxprimes = (unsigned long) (1.2 * (double) P) /
    log (1.2 * (double) P) - (double) P / log ((double) P);
#else
  unsigned long maxprimes = (unsigned long) 2.0 * (double) P /
    log (2.0 * (double) P) - (double) P / log ((double) P);
#endif

  *primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  for (p = 2; p < P; p = getprime (p));

#ifdef LESS_P
  while (p <= (P + P/5)) {
#else
  while (p <= 2 * P) {
#endif
    if (nprimes + 1 >= maxprimes) {
      maxprimes += maxprimes / 10;
      *primes = (uint32_t*) realloc (*primes, maxprimes * sizeof (uint32_t));
      if ( (*primes) == NULL) {
        fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
        exit (1);
      }
    }
    (*primes)[nprimes++] = p;
    p = getprime (p);
  }

  getprime (0); /* free the memory used by getprime */

  *primes = (uint32_t*) realloc (*primes, (nprimes) * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  return nprimes;
}


/* clear prime array */
void
printPrimes ( uint32_t *primes,
              unsigned long size )
{
  unsigned long i;
  for (i = 0; i < size; i++) {
    fprintf (stderr, "(%lu, %"PRIu32") ", i, primes[i]);
    if ((i+1) % 5 == 0)
      fprintf (stderr, "\n");
  }
  fprintf (stderr, "\n");
}


/* clear prime array */
void
clearPrimes ( uint32_t **primes )
{
  free (*primes);
}


/* init the header struct */
void
header_init ( header_t header,
              mpz_t N,
              unsigned long d,
              uint64_t ad )
{
  /* compute Ntilde, m0 */
  mpz_init_set (header->N, N);
  mpz_init (header->Ntilde);
  mpz_init (header->m0);
  header->d = d;
  header->ad = ad;

  /* compute Ntilde, ... from N, ... */
  mpz_set_uint64 (header->Ntilde, header->ad);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_pow_ui (header->Ntilde, header->Ntilde, header->d - 1);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_mul (header->Ntilde, header->Ntilde, header->N); /* d^d * ad^(d-1) * N */
  mpz_root (header->m0, header->Ntilde, header->d);
}


/* clear header struct */
void
header_clear ( header_t header )
{
  mpz_clear (header->m0);
  mpz_clear (header->Ntilde);
  mpz_clear (header->N);
}


/* init proots_t */
void
proots_init ( proots_t R,
              unsigned long size )
{
  R->size = size;

  /* length of nr&roots is known now. lengths of roots[i] are TBD. */
  R->nr = (unsigned int *) malloc (size * sizeof (unsigned int));
  R->roots = (uint64_t **) malloc (size * sizeof (uint64_t*));

  if (R->nr == NULL || R->roots == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in proots_init().\n");
    exit (1);
  }
}


/* add a root to proots_t */
void
proots_add ( proots_t R,
             unsigned long nr,
             uint64_t *roots,
             unsigned long index )
{
  unsigned int i;
  R->nr[index] = nr;

  if (nr != 0) {
    R->roots[index] = (uint64_t *) malloc (nr * sizeof (uint64_t));
    if (R->roots[index] == NULL) {
      fprintf (stderr, "Error, cannot allocate memory in proots_add\n");
      exit (1);
    }

    for (i = 0; i < nr; i++)
      R->roots[index][i] = roots[i];
  }
  else
    R->roots[index] = NULL;
}


/* print roots */
void
proots_print ( proots_t R,
               unsigned long size )
{
  unsigned int i, j;
  for (i = 0; i < size; i++) {
    if (R->nr[i] == 0) {
      fprintf (stderr, "NULL\n");
    }
    else {
      for (j = 0; j < R->nr[i]; j ++)
        fprintf (stderr, "%"PRIu64" ", R->roots[i][j]);
      fprintf (stderr, "\n");
    }
  }
}


/* clear roots */
void
proots_clear ( proots_t R,
               unsigned long size )
{
  unsigned int i;

  free (R->nr);
  for (i = 0; i < size; i++)
    free (R->roots[i]);
  free (R->roots);
}


void
qroots_init (qroots_t R)
{
  R->alloc = 0;
  R->size = 0;
  R->q = NULL;
  R->nr = NULL;
  R->roots = NULL;
}

void
qroots_realloc (qroots_t R, unsigned long newalloc)
{
  assert (newalloc >= R->size);
  R->alloc = newalloc;
  R->q = realloc (R->q, newalloc * sizeof (unsigned int));
  if (R->q == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->nr = realloc (R->nr, newalloc * sizeof (unsigned int));
  if (R->nr == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->roots = realloc (R->roots, newalloc * sizeof (uint64_t*));
  if (R->roots == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
}

/* reorder by nr */
void
qroots_rearrange (qroots_t R)
{
  if (R->size > 1) {
    unsigned int i, j, k, max, tmpq, tmpnr;
    uint64_t *tmpr = malloc (MAX_DEGREE * sizeof (uint64_t));

    for (i = 0; i < R->size; i ++) {
      max = i;
      for (j = i+1; j < R->size; j++) {
        if (R->nr[j] > R->nr[max]) {
          max = j;
        }
      }
      
      tmpq = R->q[i];
      tmpnr = R->nr[i];
      for (k = 0; k < MAX_DEGREE; k ++)
        tmpr[k] = R->roots[i][k];

      R->q[i] = R->q[max];
      R->nr[i] = R->nr[max];
      for (k = 0; k < MAX_DEGREE; k ++)
        R->roots[i][k] = R->roots[max][k];

      R->q[max] = tmpq;
      R->nr[max] = tmpnr;
      for (k = 0; k < MAX_DEGREE; k ++)
        R->roots[max][k] = tmpr[k];
    }
    free (tmpr);
  }
}

void
qroots_add (qroots_t R, unsigned int q, unsigned int nr, uint64_t *roots)
{
  unsigned int i;

  if (nr == 0)
    return;
  if (R->size == R->alloc)
    qroots_realloc (R, R->alloc + R->alloc / 2 + 1);
  R->q[R->size] = q;
  R->nr[R->size] = nr;
  R->roots[R->size] = malloc (MAX_DEGREE * sizeof (uint64_t));
  if (R->roots[R->size] == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in roots_add\n");
    exit (1);
  }
  for (i = 0; i < nr; i++)
    R->roots[R->size][i] = roots[i];
  R->size ++;
}

void
qroots_print (qroots_t R)
{
  unsigned int i, j;
  for (i = 0; i < R->size; i++) {
    fprintf (stderr, "q: %u, r: ", R->q[i]);
    for (j = 0; j < R->nr[i]; j ++)
      fprintf (stderr, "%"PRIu64" ", R->roots[i][j]);
    fprintf (stderr, "\n");
  }
}

void
qroots_clear (qroots_t R)
{
  unsigned int i;

  free (R->q);
  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}


/* init hash table */
void
hash_init (hash_t H, unsigned int init_size)
{
  H->alloc = init_size;
  H->slot = (slot_t*) malloc (H->alloc * sizeof (slot_t));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in hash_init\n");
    exit (1);
  }

  for (unsigned int j = 0; j < H->alloc; j ++) {
    H->slot[j].i = 0;
    H->slot[j].p = 0;
  }
  H->size = 0;
#ifdef DEBUG_HASH_TABLE
  H->coll = 0;
  H->coll_all = 0;
#endif
}

/* init_size is an approximation of the number of entries */
void
shash_init (shash_t H, unsigned int init_size)
{
  int j;
  unsigned int init_size0 = init_size;

  /* round up to multiple of SHASH_NBUCKETS */
  init_size = 1 + (init_size - 1) / SHASH_NBUCKETS;
  init_size += init_size / 10 + 128; /* use 10% margin */
  if (init_size > init_size0)
    init_size = init_size0;
  H->alloc = init_size * (SHASH_NBUCKETS + 1) + 5;
  /* + init_size for guard for the last buckets to avoid seg fault */
  /* + 5 for extreme guard when init_size is too small */
  H->mem = (uint64_t*) malloc (H->alloc * sizeof (uint64_t));
  if (!H->mem)
    {
      fprintf (stderr, "Error, cannot allocate memory in shash_init\n");
      exit (1);
    }
  H->balloc = init_size;
  H->base[0] = H->current[0] = H->mem;
  for (j = 1; j <= SHASH_NBUCKETS; j++)
    {
      H->base[j] = H->current[j] = H->base[j-1] + H->balloc;
    }
  /* Trick for prefetch T in shash_find_collision after the end
     of the last bucket. */
  memset (H->base[SHASH_NBUCKETS], 0, sizeof(**(H->base) * 5));
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
hash_add (hash_t H, unsigned long p, int64_t i, mpz_t m0, uint64_t ad,
          unsigned long d, mpz_t N, unsigned long q, mpz_t rq)
{
  uint32_t h;

  ASSERT(m0 != NULL);
  ASSERT(H->size < H->alloc);

  h = (uint32_t) i % H->alloc;

#ifdef DEBUG_HASH_TABLE
  if (H->slot[h].i != 0)
    H->coll ++;
#endif
  while (H->slot[h].i != 0)
  {
    if (H->slot[h].i == i) /* we cannot have H->slot[h].p = p, since for a
                       given prime p, all (p,i) values entered are different */
      match (H->slot[h].p, p, i, m0, ad, d, N, q, rq);
    if (UNLIKELY(++h == H->alloc))
      h = 0;
#ifdef DEBUG_HASH_TABLE
    H->coll_all ++;
#endif
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}

int
shash_find_collision (shash_t H)
{
  static uint32_t size = 0, mask;
  uint64_t *Hj, *Hjm;
  uint32_t *Th0, *Th1, *Th2, *Th3, *Th4, *T, *Tend;
  uint64_t i0, i1, i2, i3, i4;
  uint32_t k;
  unsigned int key;
  
#define SHASH_RESEARCH(TH,I)				\
  do {							\
    key = ((I) >> 32) + (I);				\
    while (UNLIKELY(*TH)) {				\
      if (UNLIKELY(*TH == key)) { free (T); return 1; }	\
      TH++; if (UNLIKELY(TH == Tend)) TH = T;		\
    }							\
    *TH = key;						\
  } while (0)

#define SHASH_TH_I(TH,I,IND)			\
  do {						\
    I = Hj[IND];				\
    TH = T + ((I >> LN2SHASH_NBUCKETS) & mask); \
    __builtin_prefetch(TH, 1, 3);		\
  } while (0)					\
  
  if (!size) {
    size = H->balloc << 1;
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 1;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = size - 1;
  }
  T = (uint32_t*) malloc (size * sizeof(*T));
  Tend = T + size;
  for (k = 0; k < SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    /* Here, a special guard at the end of shash_init allows
       until Hjm[SHASH_BUCKETS-1] + 5.
       So, it's not needed to test if Hj + 4 < Hjm to avoid prefetch problem. */
    SHASH_TH_I(Th0, i0, 0);
    SHASH_TH_I(Th1, i1, 1);
    SHASH_TH_I(Th2, i2, 2);
    SHASH_TH_I(Th3, i3, 3);
    SHASH_TH_I(Th4, i4, 4);
    Hj += 5;
    while (LIKELY(Hj < Hjm)) {
      __builtin_prefetch(Hj, 0, 3);
      __builtin_prefetch(((void *) Hj) + 32, 0, 3);
      SHASH_RESEARCH(Th0, i0); SHASH_TH_I(Th0, i0, 0);
      SHASH_RESEARCH(Th1, i1); SHASH_TH_I(Th1, i1, 1);
      SHASH_RESEARCH(Th2, i2); SHASH_TH_I(Th2, i2, 2);
      SHASH_RESEARCH(Th3, i3); SHASH_TH_I(Th3, i3, 3);
      SHASH_RESEARCH(Th4, i4); SHASH_TH_I(Th4, i4, 4);
      Hj += 5;
    }
    switch (Hj - Hjm) { /* no break: it's NOT an error! */
    case 0: SHASH_RESEARCH(Th4, i4);
    case 1: SHASH_RESEARCH(Th3, i3);
    case 2: SHASH_RESEARCH(Th2, i2);
    case 3: SHASH_RESEARCH(Th1, i1);
    case 4: SHASH_RESEARCH(Th0, i0);
    }
  }
  free (T);
  return 0;
}
#undef SHASH_TH_I
#undef SHASH_RESEARCH

/* return non-zero iff there is a collision */
#define PREFETCH 16
int
shash_find_collision_old (shash_t H)
{
  static uint32_t size = 0, mask;
  struct {
    uint32_t *Th, key;
  } data[PREFETCH], *pdata, *edata, *ldata;
  uint32_t *Th;
  uint32_t key;
  uint64_t *Hj, *Hjm;
  uint32_t *T, *Tend;
  uint64_t i;
  unsigned int j, k, l;
  
  if (!size) {
    size = H->balloc + (H->balloc >> 1);
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 1;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = size - 1;
  }
  T = (uint32_t*) malloc (size * sizeof(*T));
  Tend = T + size;
  edata = data + PREFETCH;
  for (k = 0; k <SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    pdata = data;
    j = Hjm - Hj;
    if (j > PREFETCH) j = PREFETCH;
    for (l = 0; l < j; l++) {
      i = Hj[l];
      data[l].Th = T + ((i >> LN2SHASH_NBUCKETS) & mask);
      __builtin_prefetch(pdata->Th, 1, 0);
      data[l].key = (i >> 32) + i;
    }
    Hj += j;
    if (LIKELY(j == PREFETCH)) {
      while (LIKELY(Hj != Hjm)) {
	i = *Hj++;
	Th = pdata->Th;
	pdata->Th = T + ((i >> LN2SHASH_NBUCKETS) & mask);
	__builtin_prefetch(pdata->Th, 1, 0);
	key = pdata->key;
	pdata->key = (i >> 32) + i;
	pdata++;
	if (UNLIKELY (pdata == edata)) pdata = data;
	if (LIKELY(!*Th))
	  *Th = key;
	else
	  {
	    do {
	      if (UNLIKELY(*Th == key))
		{
		  free (T);
		  return 1;
		}
	      Th++;
	      if (UNLIKELY(Th == Tend))
		Th = T;
	    } while (*Th);
	    *Th = key;
	  }
      }
      ldata = pdata;
    }
    else
      ldata = data + l;
    do {
      Th = pdata->Th;
      key = pdata->key;
      pdata++;
      if (UNLIKELY (pdata == edata)) pdata = data;
      if (LIKELY(!*Th))
	*Th = key;
      else
        {
          do {
            if (UNLIKELY(*Th == key))
              {
                free (T);
                return 1;
              }
            Th++;
            if (UNLIKELY(Th == Tend))
              Th = T;
          } while (UNLIKELY(*Th));
          *Th = key;
        }
    } while (LIKELY(ldata != pdata));
  }
  free (T);
  return 0;
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
gmp_hash_add (hash_t H, uint32_t p, int64_t i, mpz_t m0, uint64_t ad,
              unsigned long d, mpz_t N, uint64_t q, mpz_t rq)
{
  unsigned long h;

  if (H->size >= H->alloc)
    hash_grow (H);
  if (i >= 0)
    h = ((int)i) % H->alloc;
  else
  {
    h = H->alloc - ( ((int)(-i)) % H->alloc);
    if (h == H->alloc)
      h = 0;
  }

  while (H->slot[h].p != 0)
  {
    if (m0 != NULL && H->slot[h].i== i && H->slot[h].p != p) {
      gmp_match (H->slot[h].p, p, i, m0, ad, d, N, q, rq);
    }
    if (++h == H->alloc)
      h = 0;
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}

void
hash_clear (hash_t H)
{
  free (H->slot);
}

void
shash_clear (shash_t H)
{
  free (H->mem);
}

void
hash_grow (hash_t H)
{
  unsigned long j, old_alloc;
  slot_t *old_slot;
  mpz_t tmp;
  mpz_init (tmp);
  mpz_set_ui (tmp, 0);

  old_alloc = H->alloc;
  old_slot = H->slot;
  H->alloc = 2 * old_alloc;
  H->slot = (slot_t*) malloc (H->alloc * sizeof (slot_t));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in hash_init\n");
    exit (1);
  }
  memset (H->slot, 0, (sizeof(int64_t) + sizeof(uint32_t)) * H->alloc);
  H->size = 0;

  for (j = 0; j < old_alloc; j++)
    if (old_slot[j].p != 0)
      hash_add (H, old_slot[j].p, old_slot[j].i, NULL, 0, 0, NULL, 0, tmp);

  free (old_slot);
  mpz_clear (tmp);
}

#if 0
double
hash_mean_value (hash_t H)
{
  double s = 0;
  unsigned long j;

  for (j = 0; j < H->alloc; j++)
    if (H->p[j] != 0)
      s += fabs ((double) H->i[j]);
  return s / (double) H->size;
}
#endif

