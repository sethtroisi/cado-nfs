#include "cado.h"

#include "fb.h"
#include "utils.h"           /* lots of stuff */
#include "bucket.h"
#include "modredc_ul.h"
#include "modredc_2ul2.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "las-debug.h"
#include "las-arith.h"
#include "las-qlattice.h"
#include "las-fill-in-buckets.h"

#ifdef USE_CACHEBUFFER
#include "cachebuf.h"
#endif

/* For memcpy in fill_in_k_buckets & fill_in_m_buckets.
   When you have to move data whose length is a
   static const uint8_t N <= 16 bytes,
   it is faster to move optimal_move[N]>=N bytes with a memcpy
   (if possible, of course) :
   it is done with one or two instructions */
static const uint8_t optimal_move[] = { 0, 1, 2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16 };

#include "las-plattice.h"

/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/**************************************************************************
 * Global DEFINEs for fill_in_buckets, fill_in_k_buckets, fill_in_m_buckets 
 **************************************************************************/

/* Do we test for co-primality of i,j before writing a bucket update?
   With SKIP_GCD_UV, we test that gcd(u,v) = 1, where (i,j)~ = M_p (u,v)~.
   Without SKIP_GCD_UV, but with SKIP_GCD3, we test whether 3 | gcd(i,j).
   With neither defined, we don't test except for the hard-coded test of
   2 | gcd(i,j) */

// #define SKIP_GCD_UV 1

/* Return if a is divisible by 3 */
#ifdef SKIP_GCD3
static inline int 
is_divisible_3_u32 (uint32_t a)
{
  return (a * (- (UINT32_MAX / 3)) <= UINT32_MAX / 3);
}
#endif

#ifdef SKIP_GCD_UV
    #define INCREASE_UV(uv) do {(uv)++;} while(0)
    #define PROBABLY_COPRIME_IJ(u,v) (gcd_ul(u, v) == 1)
#else
    #define INCREASE_UV(uv) do {} while(0)
    #ifdef SKIP_GCD3
    #define PROBABLY_COPRIME_IJ(u,v)			\
         (!is_divisible_3_u32 (i + I) ||			\
          !is_divisible_3_u32 ((uint32_t) (x >> logI)))			
    #else
    #define PROBABLY_COPRIME_IJ(u,v) 1
    #endif
#endif


#ifdef TRACE_K								
#define FILL_BUCKET_TRACE_K(X) do {					\
    if (trace_on_spot_x(X)) {						\
      WHERE_AM_I_UPDATE(w, N, (X) >> 16);				\
      WHERE_AM_I_UPDATE(w, x, (uint16_t) (X));				\
      fprintf (stderr, "# Pushed factor base entry (%u, %u) (x=%u, %s) to BA[%u]\n",	\
	       (unsigned int) slice_index, (unsigned int) slice_offset, (unsigned int) (uint16_t) (X), sidenames[side],	\
	       (unsigned int) ((X) >> 16));				\
      ASSERT(test_divisible(w));					\
    }									\
  } while(0)
#else
#define FILL_BUCKET_TRACE_K(X)
#endif

#ifdef HAVE_SSE2							
#define FILL_BUCKET_PREFETCH(PT) do {				\
    _mm_prefetch((char *)(PT), _MM_HINT_T0);			\
  } while (0)
#else
#define FILL_BUCKET_PREFETCH(PT)
#endif

#define FILL_BUCKET_INC_X() do {					\
    if (i >= bound1) {							\
      x += inc_a;							\
      INCREASE_UV(u);							\
    }									\
    if (i < bound0) {							\
      x += inc_c;							\
      INCREASE_UV(v);							\
    }									\
  } while (0)
/************************************************************************/

#ifdef USE_CACHEBUFFER
DECLARE_CACHE_BUFFER(bucket_update_t, 256)
#endif

#if 0

/* Returns -(-1)^neg num/den (mod p). If gcd(den, p) > 1, returns p. */
static inline unsigned long
compute_1_root_ul(const unsigned long p, const unsigned long num,
                  const unsigned long den, const int k, const int neg)
{
  modulusredcul_t m;
  residueredcul_t rn, rd;

  modredcul_initmod_ul(m, p);
  modredcul_init_noset0(rn, m);
  modredcul_init_noset0(rd, m);

  modredcul_set_ul(rn, num, m);
  modredcul_set_ul(rd, den, m);
  int ok = modredcul_inv(rd, rd, m);
  if (!ok) {
    /* FIXME: p could be a prime power */
    return p;
  }
  modredcul_mul(rn, rn, rd, m);
  if (!neg)
    modredcul_neg(rn, rn, m);

  for (int i = 0; i < k; i++)
    modredcul_div2(rn, rn, m);
  unsigned long r = modredcul_get_ul(rn, m);

  modredcul_clear(rn, m);
  modredcul_clear(rd, m);
  modredcul_clearmod(m);
  return r;
}

static inline unsigned long
compute_1_root_mpz(const unsigned long p, const mpz_t f0, const mpz_t f1)
{
  return compute_1_root_ul(p, mpz_fdiv_ui(f0, p), mpz_fdiv_ui(f1, p), 0, 0);
}

struct contexts_s {
  modredc2ul2_batch_Q_to_Fp_context_t *context_2ul2;
  unsigned long k;
  unsigned long num_ul, den_ul;
  int neg;
};

static inline void
modredc2ul2_set_mpz(modintredc2ul2_t r, const mpz_t a)
{
  unsigned long t[MODREDC2UL2_SIZE];
  size_t n;
  mpz_export (t, &n, -1, sizeof(unsigned long), 0, 0, a);
  ASSERT_ALWAYS(n <= MODREDC2UL2_SIZE);
  modredc2ul2_intset_uls(r, t, n);
}

static void
init_Q_to_Fp_context(struct contexts_s *contexts, mpz_poly_srcptr f)
{
  const int use_Q_to_Fp = 0;
  contexts->context_2ul2 = NULL;
  contexts->num_ul = contexts->den_ul = 0;

  if (!use_Q_to_Fp || f->deg != 1)
    return;

  contexts->k = mpz_scan1(f->coeff[1], 0);
  mpz_t odd_den_mpz;
  mpz_init(odd_den_mpz);
  mpz_tdiv_q_2exp(odd_den_mpz, f->coeff[1], contexts->k);

  /* If exactly one of the two coefficients is negative, we need to flip
     the sign of the root (mod p). The odd_den and num variables contain
     absolute values. */
  contexts->neg = (mpz_sgn(f->coeff[0]) < 0) != (mpz_sgn(f->coeff[1]) < 0);

  modintredc2ul2_t num, odd_den;
  modredc2ul2_intinit(num);
  modredc2ul2_intinit(odd_den);
  /* Set num = |f->coeff[0]|, odd_den = |odd_den_mpz| */
  modredc2ul2_set_mpz(num, f->coeff[0]);
  modredc2ul2_set_mpz(odd_den, odd_den_mpz);

  if (mpz_fits_ulong_p(f->coeff[0]) && mpz_fits_ulong_p(odd_den_mpz))
  {
    contexts->num_ul = mpz_get_ui(f->coeff[0]);
    contexts->den_ul = mpz_get_ui(odd_den_mpz);
  } else if (mpz_sizeinbase(f->coeff[0], 2) >= MODREDC2UL2_MINBITS &&
             mpz_sizeinbase(f->coeff[0], 2) <= MODREDC2UL2_MAXBITS &&
             mpz_sizeinbase(odd_den_mpz, 2) >= MODREDC2UL2_MINBITS &&
             mpz_sizeinbase(odd_den_mpz, 2) <= MODREDC2UL2_MAXBITS )
  {
    contexts->context_2ul2 = modredc2ul2_batch_Q_to_Fp_init (num, odd_den);
  }

  mpz_clear(odd_den_mpz);
  modredc2ul2_intclear(num);
  modredc2ul2_intclear(odd_den);
}

static void
clear_Q_to_Fp_context(struct contexts_s *contexts)
{
  if (contexts->context_2ul2 != NULL)
    modredc2ul2_batch_Q_to_Fp_clear(contexts->context_2ul2);
  contexts->context_2ul2 = NULL;
  contexts->num_ul = contexts->den_ul = 0;
}


/* Compute up to n roots of the transformed polynomial, either by tranforming
   the root stored in the factor base, or, if the polynomial has degree 1 and
   its coefficients fit into unsigned long, by computing the root from the
   transformed polynomial.
   For each root, the modulus (i.e., factor base prime) and
   the root are stored in p[i] and r[i], resp.
   Returns the number of transformed roots produced, which may be less than n
   if we reach the end of the factor base. */
static inline size_t
transform_n_roots(unsigned long *p, unsigned long *r, fb_iterator t,
                  const size_t n, const int side, sieve_info_srcptr si,
                  struct contexts_s *contexts)
{
  size_t i;
  int ok = 0;
  if (contexts->den_ul == 0 && contexts->context_2ul2 == NULL) {
    /* Do the old transform with separate inverse modulo each p */
    for (i = 0; !fb_iterator_over(t) && i < n; i++, fb_iterator_next(t)) {
      fbprime_t current_p = fb_iterator_get_p(t);
      ASSERT_ALWAYS (current_p & 1);
      fbroot_t R = fb_iterator_get_r(t);
      p[i] = current_p;
      r[i] = fb_root_in_qlattice(current_p, R, t->fb->invp, si);
    }
    return i;
  }

  /* Collect primes from fb, store up to n namy of them in p, and the actual
     number of stored primes in i */
  for (i = 0; !fb_iterator_over(t) && i < n; i++, fb_iterator_next(t)) {
    p[i] = fb_iterator_get_p(t);
  }

  if (contexts->den_ul != 0) {
    ok = modredcul_batch_Q_to_Fp(r, contexts->num_ul, contexts->den_ul, contexts->k, p, i);
    if (ok && !contexts->neg) {
      /* We computed r[i] = |g_0| / |g_1| (mod p[i]). If both g0, g1 have the same sign,
         then the root is -g_0/g_1, and we need a sign flip. */
      for (size_t j = 0; j < i; j++)
        if (r[j] != 0)
          r[j] = p[j] - r[j];
    }
  } else if (contexts->context_2ul2 != NULL) {
    ok = modredc2ul2_batch_Q_to_Fp(r, contexts->context_2ul2, contexts->k, !contexts->neg, p, i);
  }

  if(UNLIKELY(!ok)) {
    /* Modular inverse did not exists, i.e., at least one of the factor base
       entries was not coprime to den. Do factor base entries one at a time,
       handling the non-coprime case properly. This should happen at most a
       few times per special-q. */
    mpz_poly_srcptr f = si->sides[side]->fij;
    for (size_t j = 0; j < i; j++)
      r[j] = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
  }

#if 0
    /* Very expensive but thorough test */
    for (size_t j = 0; j < i; j++) {
      mpz_poly_srcptr f = si->sides[side]->fij;
      unsigned long tr = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
      ASSERT_ALWAYS(tr == r[j]);
    }
#endif

  return i;
}
#endif


static inline
void fill_bucket_heart(bucket_array_t &BA, const uint64_t x, const prime_hint_t hint,
                       const int side MAYBE_UNUSED,
                       const slice_index_t slice_index MAYBE_UNUSED, 
                       where_am_I_ptr w MAYBE_UNUSED)
{
#if 0
  bucket_update_t **pbut = BA.bucket_write + (x >> 16);
  bucket_update_t *but = *pbut;
  FILL_BUCKET_TRACE_K(x);
  WHERE_AM_I_UPDATE(w, N, x >> 16);
  WHERE_AM_I_UPDATE(w, x, (uint16_t) x);
  but->slice_offset = slice_offset;
  but->x = (uint16_t) x;
  *pbut = ++but;
  FILL_BUCKET_PREFETCH(but);
#else
  const int i = (x >> 16);
  bucket_update_t update;
  update.x = (uint16_t) x;
  update.hint = hint;
  push_bucket_update(BA, i, update);
                     
#endif
}


/* {{{ */
void
fill_in_buckets(thread_data_ptr th, const int side,
                const fb_transformed_vector *transformed_vector,
                const fb_slice_interface *slice MAYBE_UNUSED,
                where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  const slice_index_t slice_index = transformed_vector->get_index();
  bucket_array_t BA = th->sides[side]->BA;  /* local copy. Gain a register + use stack */
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  /* Write new set of pointers for the new slice */
  bucket_add_slice_index(&BA, slice_index);

  for (fb_transformed_vector::const_iterator pl_it = transformed_vector->cbegin();
       pl_it != transformed_vector->cend(); pl_it++) {
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */

#if 0
    printf("%s(): sieving side=%d, p=%u, logp = %u, a = (%d, %u), b = (%u, %u)\n", 
           __func__, side, pl_it->det(), (unsigned int) slice->get_logp(), pl_it->get_a0(), pl_it->get_a1(), pl_it->get_b0(), pl_it->get_b1());
#endif
    
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
    const uint32_t maskI = I-1;
    const uint64_t even_mask = (1ULL << logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << logI;
    const prime_hint_t hint = pl_it->hint;

    if (pl_it->get_b0() == 0 || pl_it->get_b1() == 0) {
      /* r == 0 or r == 1/0.
         1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
         this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
         gcd(a,b), for i = 0 we only need to sieve j = 1. 
         So, x = j*I + (i + I/2) = I + I/2.
         2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
         this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
         FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
         (if any) do we sieve? */
      uint64_t x = (pl_it->b1 == 0 ? 1 : I) + (I >> 1);
     /*****************************************************************/
      fill_bucket_heart(BA, x, hint, side, slice_index, w);
      continue;
    }

    const uint32_t bound0 = pl_it->get_bound0(logI), bound1 = pl_it->get_bound1(logI);
    const uint64_t inc_a = pl_it->get_inc_a(logI), inc_c = pl_it->get_inc_c(logI);
    uint64_t x = 1ULL << (logI-1);
    uint32_t i = x;
    FILL_BUCKET_INC_X();

    /* Now, do the real work: the filling of the buckets */
    while (x < IJ){
      /***************************************************************/
#define FILL_BUCKET() do {						\
        unsigned int i = x & maskI;					\
        if (LIKELY(((x & even_mask))		\
                   )) fill_bucket_heart(BA, x, hint, slice_index, side, w);	\
        FILL_BUCKET_INC_X();						\
      } while (0)
      /***************************************************************/
      FILL_BUCKET();
    }
  }
  th->sides[side]->BA = BA;
}

#ifdef HAVE_K_BUCKETS
/* Same than fill_in_buckets, but with 2 passes (k_buckets & buckets). */
void
fill_in_k_buckets(thread_data_ptr th, int side, where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  bucket_array_t BA = th->sides[side]->BA;      /* local copy, gain a register + use stack */
  k_bucket_array_t kBA = th->sides[side]->kBA;
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  fb_iterator t;
  fb_iterator_init_set_fb(t, th->sides[side]->fb_bucket);
  for( ; !fb_iterator_over(t) ; fb_iterator_next(t)) {
    fbprime_t p = t->fb->p;
    ASSERT_ALWAYS (p & 1);
    WHERE_AM_I_UPDATE(w, p, p);
    
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */
    if (UNLIKELY(!mpz_cmp_ui(si->doing->p, p))) continue;
    fbprime_t R = fb_iterator_get_r(t), r = fb_root_in_qlattice(p, R, t->fb->invp, si);
    
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
    const uint32_t maskI = I-1;
    const uint64_t even_mask = (1ULL << logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << logI;

    /* Special cases */
    if (UNLIKELY((!r) || (r >= p))) {
      if (r > p) /* should only happen for lattice-sieved prime powers,
		    which is not possible currently since maxbits < I */
	continue;
      /* r == p or r == 0.
	 1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
	 this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
	 gcd(a,b), for i = 0 we only need to sieve j = 1. 
	 So, x = j*I + (i + I/2) = I + I/2.
	 2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
	 this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
	 FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
	 (if any) do we sieve? */
      uint64_t x = (r ? 1 : I) + (I >> 1);
      prime_hint_t hint = bucket_encode_prime(p);
      /* 1. pkbut must be volatile in BIG_ENDIAN: the write order of prime & x
	 is need because the last byte of x (always 0 because x <<= 8) must
	 be overwritten by the first byte of prime.
	 2. memcpy is good because it's its job to know if it's possible to
	 write more than 1 byte at an odd adress (ODD_ADDRESS_IO_INT).
	 gcc does a optimal job with memcpy & a little + constant length.
      */
      /**************************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define FILL_K_BUCKET_HEART() do {					\
	k_bucket_update_t **pkbut = kBA.bucket_write + (x >> 24);	\
	k_bucket_update_t *kbut = *pkbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	memcpy(kbut, &hint, sizeof(prime_hint_t));			\
	uint32_t i = (uint32_t) x;					\
	memcpy((uint8_t *) kbut + sizeof(prime_hint_t), &i, 4);		\
	*pkbut = ++kbut;						\
	FILL_BUCKET_PREFETCH(kbut);					\
      } while (0)
#else
#define FILL_K_BUCKET_HEART() do {					\
	k_bucket_update_t **pkbut = kBA.bucket_write + (x >> 24);	\
	volatile k_bucket_update_t *kbut = *pkbut;			\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	uint32_t i = (uint32_t) x << 8;					\
	memcpy(kbut, &i, 4);						\
	memcpy((uint8_t *) kbut + 3, &hint, sizeof(prime_hint_t));	\
	*pkbut = ++kbut;						\
	FILL_BUCKET_PREFETCH(kbut);					\
      } while(0)
#endif
      /**************************************************************************/
      FILL_K_BUCKET_HEART();
      continue;
    }
    /* If working with congruence classes, once the loop on the parity goes at the level
       above, this initialization should in fact either be done for each congruence class,
       or saved for later use within the factor base structure. */
    plattice_info_t pli;
    if (UNLIKELY(!reduce_plattice(&pli, p, r, si->I))) {
      verbose_output_print (1, 1, "# fill_in_buckets: reduce_plattice() returned 0 for p = %"
	                    FBPRIME_FORMAT ", r = %" FBPRIME_FORMAT "\n", p, r);
      continue; /* Simply don't consider that (p,r) for now.
                   FIXME: can we find the locations to sieve? */
    }
    /* OK, all special cases are done. */

    const uint32_t bound0 = plattice_bound0(&pli, logI), bound1 = plattice_bound1(&pli, logI);
#if !MOD2_CLASSES_BS
    const uint64_t inc_a = plattice_a(&pli, logI), inc_c = plattice_c(&pli, logI);
    uint64_t x = 1ULL << (logI-1);
    uint32_t i = x;
#ifdef SKIP_GCD_UV
    unsigned long u = 0, v = 0;
#endif
    FILL_BUCKET_INC_X();
    if (x >= IJ) continue;
#else
    for(unsigned int parity = 1 ; parity < 4; parity++) {
      // The sieving point (0,0) is I/2 in x-coordinate
      uint64_t x = plattice_starting_vector(&pli, si, parity);
      if (x >= IJ) continue;
      const uint64_t inc_a = plattice_a(&pli, logI), inc_c = plattice_c(&pli, logI);
#endif
      const prime_hint_t hint = bucket_encode_prime (p);
      
      /* Now, do the real work: the filling of the k-buckets */
      do { 
	/**********************************************************************/
#define FILL_K_BUCKET() do {						\
	  unsigned int i = x & maskI;					\
	  if (LIKELY((MOD2_CLASSES_BS || (x & even_mask))		\
	             && PROBABLY_COPRIME_IJ(u,v)			\
		     )) FILL_K_BUCKET_HEART();				\
	  FILL_BUCKET_INC_X();						\
	} while (0)
	/****************************************************************/
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET(); if (x >= IJ) break;
	FILL_K_BUCKET();
      } while (x < IJ);
#if MOD2_CLASSES_BS
    }
#endif
  }
  th->sides[side]->kBA = kBA;
  
  /* sort : 2nd pass; kBA -> BA */
  bucket_update_t **pbw = BA.bucket_write;
  for (uint32_t kb = 0; kb < kBA.n_bucket; ++kb) {
    uint8_t *kbs = (uint8_t *) (kBA.bucket_start[kb]);
#ifdef USE_CACHEBUFFER
    bucket_update_t_256_cachebuffer cachebuf;
#endif
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) BA.bucket_write + BA.size_b_align - (size_t) pbw;
      if (LIKELY(lg > sizeof(bucket_update_t **) << 8)) lg = sizeof(bucket_update_t **) << 8;
      bucket_update_t **pbl = BA.logp_idx + (kb << 8);
      k_bucket_update_t **pkbl = kBA.logp_idx + kb;
      uint8_t *kbl = (uint8_t *) *pkbl;
      /* There are BA.nr_logp duplicates of all kBA.bucket_write in kBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
#ifdef USE_CACHEBUFFER
        init_bucket_update_t_256_cachebuffer(cachebuf, pbw, MIN(BA.n_bucket, 256));
#endif
	/* Twelve kBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     kbs +  sizeof(k_bucket_update_t)*12 <= kbl;
	     kbs += sizeof(k_bucket_update_t)*12) {
	  /*****************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#ifdef USE_CACHEBUFFER
#define KBA_2_BA(A) do {						\
	    const size_t bucket_idx = kbs[(A) + sizeof(bucket_update_t)]; \
	    add_bucket_update_t_256_to_cachebuffer(cachebuf, bucket_idx, \
	                                           *(bucket_update_t *) (kbs + (A))); \
	  } while (0)
#else
#define KBA_2_BA(A) do {						\
	    bucket_update_t **pbut, *but;				\
	    pbut = pbw + kbs[(A) + sizeof(bucket_update_t)];		\
	    but = *pbut;						\
	    memcpy(but, kbs + (A), optimal_move[sizeof(bucket_update_t)]); \
	    *pbut = ++but;						\
	  } while (0)
#endif /* USE_CACHEBUFFER */
#else
#define KBA_2_BA(A) do {						\
	    bucket_update_t **pbut, *but;				\
	    pbut = pbw + kbs[A];					\
	    but = *pbut;						\
	    memcpy(but, kbs+(A)+1, optimal_move[sizeof(bucket_update_t)]); \
	    *pbut = ++but;						\
	  } while (0)
#endif
	  /*****************************************************************/
	  KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
	  KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
	  KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
	  KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
	  KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
	  KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
	}
	for (; kbs < kbl; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
#ifdef USE_CACHEBUFFER
        flush_bucket_update_t_256_cachebuffer(cachebuf);
#endif
	/* OK, let's duplicate the current (at most) 256 pointers from
	   BA.bucket_write in BA.logp_idx */
	aligned_medium_memcpy(pbl, pbw, lg);
	pbl =    (bucket_update_t **) ((size_t)  pbl +  BA.size_b_align);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	kbl = (uint8_t *) *pkbl;
      }
    }
    /* 2nd part: BA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *kbw = (uint8_t *) (kBA.bucket_write[kb]);
#ifdef USE_CACHEBUFFER
    init_bucket_update_t_256_cachebuffer(cachebuf, pbw, MIN(BA.n_bucket, 256));
#endif
    for (;
	 kbs +  sizeof(k_bucket_update_t)*12 <= kbw;
	 kbs += sizeof(k_bucket_update_t)*12) {
      KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
      KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
      KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
      KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
      KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
      KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
    }
    for (; kbs < kbw; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
#ifdef USE_CACHEBUFFER
    flush_bucket_update_t_256_cachebuffer(cachebuf);
#endif
    pbw += 256;
  }
  th->sides[side]->BA = BA;
}

void
fill_in_m_buckets(thread_data_ptr th, int side, where_am_I_ptr w MAYBE_UNUSED)
{
  WHERE_AM_I_UPDATE(w, side, side);
  sieve_info_srcptr si = th->si;
  bucket_array_t BA = th->sides[side]->BA;  /* local copy; gain a register + use stack */
  k_bucket_array_t kBA = th->sides[side]->kBA;
  m_bucket_array_t mBA = th->sides[side]->mBA;
  // Loop over all primes in the factor base.
  //
  // Note that dispatch_fb already arranged so that all the primes
  // which appear here are >= bucket_thresh and <= pmax (the latter
  // being for the moment unconditionally set to FBPRIME_MAX by the
  // caller of dispatch_fb).
  
  fb_iterator t;
  fb_iterator_init_set_fb(t, th->sides[side]->fb_bucket);
  unsigned char last_logp = 0;
  for( ; !fb_iterator_over(t) ; fb_iterator_next(t)) {
    fbprime_t p = t->fb->p;
    ASSERT_ALWAYS (p & 1);
    WHERE_AM_I_UPDATE(w, p, p);
    
    unsigned char logp = find_logp(th, side, p);

    /* Write new set of pointers if the logp value changed */
    if (UNLIKELY(last_logp != logp)) {
      aligned_medium_memcpy((uint8_t *)mBA.logp_idx + mBA.size_b_align * BA.nr_logp, mBA.bucket_write, mBA.size_b_align);
      BA.logp_val[BA.nr_logp++] = last_logp = logp;
    }
    
    /* If we sieve for special-q's smaller than the factor
       base bound, the prime p might equal the special-q prime q. */
    if (UNLIKELY(mpz_cmp_ui(si->doing->p, p) == 0)) continue;
    fbprime_t R = fb_iterator_get_r(t), r = fb_root_in_qlattice(p, R, t->fb->invp, si);
    
#ifdef SKIP_GCD3
    const uint32_t I = si->I;
    const unsigned int logI = si->conf->logI;
#endif
    const uint32_t maskI = si->I-1;
    const uint64_t even_mask = (1ULL << si->conf->logI) | 1ULL;
    const uint64_t IJ = ((uint64_t) si->J) << si->conf->logI;

    /* Special cases */
    if (UNLIKELY((!r) || (r >= p))) {
      if (r > p) /* should only happen for lattice-sieved prime powers,
		    which is not possible currently since maxbits < I */
	continue;
      /* r == p or r == 0.
	 1. If r == 0 (mod p), this prime hits for i == 0 (mod p), but since p > I,
	 this implies i = 0 or i > I. We don't sieve i > I. Since gcd(i,j) |
	 gcd(a,b), for i = 0 we only need to sieve j = 1. 
	 So, x = j*I + (i + I/2) = I + I/2.
	 2. r == p means root at infinity, which hits for j == 0 (mod p). Since q > I > J,
	 this implies j = 0 or j > J. This means we sieve only (i,j) = (1,0) here.
	 FIXME: what about (-1,0)? It's the same (a,b) as (1,0) but which of these two
	 (if any) do we sieve? */
      uint64_t x = (r ? 1 : si->I) + (si->I >> 1);
      prime_hint_t hint = bucket_encode_prime(p);
      /* memcpy is good because it's its job to know if it's possible to
	 write an int to all even addresses.
	 gcc does a optimal job with memcpy & a little + constant length.
      */
      /**************************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define FILL_M_BUCKET_HEART() do {					\
	m_bucket_update_t **pmbut = mBA.bucket_write + (x >> 32);	\
	m_bucket_update_t *mbut = *pmbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	memcpy(mbut, &hint, sizeof(prime_hint_t));			\
	uint32_t i = (uint32_t) x;					\
	memcpy((uint8_t *) mbut + sizeof(prime_hint_t), &i, 4);	\
	*pmbut = ++mbut;						\
	FILL_BUCKET_PREFETCH(mbut);					\
      } while (0)
#else
#define FILL_M_BUCKET_HEART() do {					\
	m_bucket_update_t **pmbut = mBA.bucket_write + (x >> 32);	\
	m_bucket_update_t *mbut = *pmbut;				\
	FILL_BUCKET_TRACE_K(x);						\
	WHERE_AM_I_UPDATE(w, N, x >> 16);				\
	WHERE_AM_I_UPDATE(w, x, (uint16_t) x);				\
	uint32_t i = (uint32_t) x;					\
	memcpy(mbut, &i, 4);						\
	memcpy((uint8_t *) mbut + 4, &hint, sizeof(prime_hint_t));	\
	*pmbut = ++mbut;						\
	FILL_BUCKET_PREFETCH(mbut);					\
      } while (0)
#endif
      /**************************************************************************/
      FILL_M_BUCKET_HEART();
      continue;
    }
    /* If working with congruence classes, once the loop on the parity goes at the level
       above, this initialization should in fact either be done for each congruence class,
       or saved for later use within the factor base structure. */
    plattice_info_t pli;
    if (UNLIKELY(!reduce_plattice(&pli, p, r, si->I))) {
      verbose_output_print (1, 1, "# fill_in_buckets: reduce_plattice() returned 0 for p = %"
	                    FBPRIME_FORMAT ", r = %" FBPRIME_FORMAT "\n", p, r);
      continue; /* Simply don't consider that (p,r) for now.
                   FIXME: can we find the locations to sieve? */
    }
    /* OK, all special cases are done. */

    const uint32_t bound0 = plattice_bound0(&pli, logI), bound1 = plattice_bound1(&pli, logI);
#if !MOD2_CLASSES_BS
    const uint64_t inc_a = plattice_a(&pli, logI), inc_c = plattice_c(&pli, logI);
    uint64_t x = 1ULL << (si->conf->logI-1);
    uint32_t i = x;
#ifdef SKIP_GCD_UV
    unsigned long u = 0, v = 0;
#endif
    FILL_BUCKET_INC_X();
    if (x >= IJ) continue;
#else
    for(unsigned int parity = 1 ; parity < 4; parity++) {
      // The sieving point (0,0) is I/2 in x-coordinate
      uint64_t x = plattice_starting_vector(&pli, si, parity);
      if (x >= IJ) continue;
      const uint64_t inc_a = plattice_a(&pli, logI), inc_c = plattice_c(&pli, logI);
#endif
      const prime_hint_t hint = bucket_encode_prime (p);

      /* Now, do the real work: the filling of the m-buckets */
      do {
	/***************************************************************/
#define FILL_M_BUCKET() do {						\
	  unsigned int i = x & maskI;					\
	  if (LIKELY((MOD2_CLASSES_BS || (x & even_mask))		\
	             && PROBABLY_COPRIME_IJ(u,v)			\
		     )) FILL_M_BUCKET_HEART();				\
	  FILL_BUCKET_INC_X();						\
	} while (0)
	/***************************************************************/
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET(); if (x >= IJ) break;
	FILL_M_BUCKET();
      } while (x < IJ);
#if MOD2_CLASSES_BS
    }
#endif
  }
  th->sides[side]->mBA = mBA;

  /* sort : 2nd pass; mBA -> kBA */
  k_bucket_update_t **pkbw = kBA.bucket_write;
  for (uint32_t mb = 0; mb < mBA.n_bucket; ++mb) {
    uint8_t *mbs = (uint8_t *) (mBA.bucket_start[mb]);
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) kBA.bucket_write + kBA.size_b_align - (size_t) pkbw;
      if (LIKELY(lg > sizeof(k_bucket_update_t **) << 8)) lg = sizeof(k_bucket_update_t **) << 8;
      k_bucket_update_t **pkbl = kBA.logp_idx + (mb << 8);
      m_bucket_update_t **pmbl = mBA.logp_idx + mb;
      uint8_t *mbl = (uint8_t *) *pmbl;
      /* There are BA.nr_logp duplicates of all mBA.bucket_write in mBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
	/* Ten mBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     mbs +  sizeof(m_bucket_update_t)*10 <= mbl;
	     mbs += sizeof(m_bucket_update_t)*10) {
	  /*****************************************************************/
#ifdef CADO_LITTLE_ENDIAN
#define MBA_2_KBA(A) do {						\
	    k_bucket_update_t **pkbut, *kbut;				\
	    pkbut = pkbw + mbs[(A)+sizeof(k_bucket_update_t)];		\
	    kbut = *pkbut;						\
	    memcpy(kbut, mbs+(A), optimal_move[sizeof(k_bucket_update_t)]); \
	    *pkbut = ++kbut;						\
	  } while (0)
#else
#define MBA_2_KBA(A) do {						\
	    k_bucket_update_t **pkbut, *kbut;				\
	    pkbut = pkbw + mbs[A];					\
	    kbut = *pkbut;						\
	    memcpy(kbut, mbs+(A)+1, optimal_move[sizeof(k_bucket_update_t)]); \
	    *pkbut = ++kbut;						\
	  } while (0)
#endif
	  /*****************************************************************/
	  MBA_2_KBA(0);	                          MBA_2_KBA(sizeof(m_bucket_update_t));
	  MBA_2_KBA(sizeof(m_bucket_update_t)*2); MBA_2_KBA(sizeof(m_bucket_update_t)*3);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*4); MBA_2_KBA(sizeof(m_bucket_update_t)*5);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*6); MBA_2_KBA(sizeof(m_bucket_update_t)*7);
	  MBA_2_KBA(sizeof(m_bucket_update_t)*8); MBA_2_KBA(sizeof(m_bucket_update_t)*9);
	}
	for (; mbs < mbl; mbs += sizeof(m_bucket_update_t)) MBA_2_KBA(0);
	/* OK, let's duplicate the current (at most) 256 pointers in
	   kBA.bucket_write in kBA.logp_idx */
	aligned_medium_memcpy(pkbl, pkbw, lg);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	pmbl = (m_bucket_update_t **) ((size_t) pmbl + mBA.size_b_align);
	mbl = (uint8_t *) *pmbl;
      }
    }
    /* 2nd part: kBA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *mbw = (uint8_t *) (mBA.bucket_write[mb]);
    for (;
	 mbs +  sizeof(m_bucket_update_t)*10 <= mbw;
	 mbs += sizeof(m_bucket_update_t)*10) {
      MBA_2_KBA(0);
      MBA_2_KBA(sizeof(m_bucket_update_t));
      MBA_2_KBA(sizeof(m_bucket_update_t)*2);
      MBA_2_KBA(sizeof(m_bucket_update_t)*3);
      MBA_2_KBA(sizeof(m_bucket_update_t)*4);
      MBA_2_KBA(sizeof(m_bucket_update_t)*5);
      MBA_2_KBA(sizeof(m_bucket_update_t)*6);
      MBA_2_KBA(sizeof(m_bucket_update_t)*7);
      MBA_2_KBA(sizeof(m_bucket_update_t)*8);
      MBA_2_KBA(sizeof(m_bucket_update_t)*9);
    }
    for (; mbs < mbw; mbs += sizeof(m_bucket_update_t)) MBA_2_KBA(0);
    pkbw += 256;
  }
  th->sides[side]->kBA = kBA;

  /* sort : 3th pass; kBA -> BA */
  bucket_update_t **pbw = BA.bucket_write;
#ifdef USE_CACHEBUFFER
  bucket_update_t_256_cachebuffer cachebuf;
#endif
  for (uint32_t kb = 0; kb < kBA.n_bucket; ++kb) {
    uint8_t *kbs = (uint8_t *) (kBA.bucket_start[kb]);
    /* First part: we rewrite 1->256 buckets, and in the same time,
       we have to deal with the rewriting of logp_idx.
       I use a block in these part to see the range of variables */
    {
      size_t lg = (size_t) BA.bucket_write + BA.size_b_align - (size_t) pbw;
      if (LIKELY(lg > sizeof(bucket_update_t **) << 8)) lg = sizeof(bucket_update_t **) << 8;
      bucket_update_t **pbl = BA.logp_idx + (kb << 8);
      k_bucket_update_t **pkbl = kBA.logp_idx + kb;
      uint8_t *kbl = (uint8_t *) *pkbl;
      /* There are BA.nr_logp duplicates of all kBA.bucket_write in kBA.logp_idx. */
      for (uint8_t nr_logp = BA.nr_logp; nr_logp; --nr_logp) {
#ifdef USE_CACHEBUFFER
        init_bucket_update_t_256_cachebuffer(cachebuf, pbw, MIN(BA.n_bucket, 256));
#endif
	/* Twelve kBA records in one time : it's the nearest of a cache line (60 bytes) */
	for (;
	     kbs +  sizeof(k_bucket_update_t)*12 <= kbl;
	     kbs += sizeof(k_bucket_update_t)*12) {
	  /* See the define of KBA_2_BA in fill_in_k_buckets */
	  KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
	  KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
	  KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
	  KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
	  KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
	  KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
	}
	for (; kbs < kbl; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
#ifdef USE_CACHEBUFFER
        flush_bucket_update_t_256_cachebuffer(cachebuf);
#endif
	/* OK, let's duplicate the current (at most) 256 pointers from
	   BA.bucket_write in BA.logp_idx */
	aligned_medium_memcpy(pbl, pbw, lg);
	pbl =    (bucket_update_t **) ((size_t)  pbl +  BA.size_b_align);
	pkbl = (k_bucket_update_t **) ((size_t) pkbl + kBA.size_b_align);
	kbl = (uint8_t *) *pkbl;
      }
    }
    /* 2nd part: BA.logp_idx is rewritten. We finish the rewrite of the current bucket */
    const uint8_t *kbw = (uint8_t *) (kBA.bucket_write[kb]);
#ifdef USE_CACHEBUFFER
    init_bucket_update_t_256_cachebuffer(cachebuf, pbw, MIN(BA.n_bucket, 256));
#endif
    for (;
	 kbs +  sizeof(k_bucket_update_t)*12 <= kbw;
	 kbs += sizeof(k_bucket_update_t)*12) {
      KBA_2_BA(0);                            KBA_2_BA(sizeof(k_bucket_update_t));
      KBA_2_BA(sizeof(k_bucket_update_t)*2);  KBA_2_BA(sizeof(k_bucket_update_t)*3);
      KBA_2_BA(sizeof(k_bucket_update_t)*4);  KBA_2_BA(sizeof(k_bucket_update_t)*5);
      KBA_2_BA(sizeof(k_bucket_update_t)*6);  KBA_2_BA(sizeof(k_bucket_update_t)*7);
      KBA_2_BA(sizeof(k_bucket_update_t)*8);  KBA_2_BA(sizeof(k_bucket_update_t)*9);
      KBA_2_BA(sizeof(k_bucket_update_t)*10); KBA_2_BA(sizeof(k_bucket_update_t)*11);
    }
    for (; kbs < kbw; kbs += sizeof(k_bucket_update_t)) KBA_2_BA(0);
#ifdef USE_CACHEBUFFER
    flush_bucket_update_t_256_cachebuffer(cachebuf);
#endif
    pbw += 256;
  }
  th->sides[side]->BA = BA;
}
#endif


/* TODO: Write loop over all slices in this factor base part, for each slice
   transform loops, call fill_in_buckets() with transformed roots. */

static void
fill_in_buckets_one_side(thread_data_ptr th, const int side)
{
    fb_part *fb = th->sides[side]->fb;
    where_am_I w;
    WHERE_AM_I_UPDATE(w, si, th->si);
    if (th->sides[side]->BA.n_bucket < THRESHOLD_K_BUCKETS) {
      /* Process all slices in this factor base part */
      const fb_slice_interface *slice;
      for (slice_index_t slice_index = fb->get_first_slice_index();
           (slice = fb->get_slice(slice_index)) != NULL;
           slice_index++) {
          const fb_transformed_vector *transformed_vector = slice->make_lattice_bases(th->si->qbasis, th->si->conf->logI);
          fill_in_buckets(th, side, transformed_vector, slice, w);
          delete transformed_vector;
      }
    }
#ifdef HAVE_K_BUCKETS
    else if (th->sides[side]->BA.n_bucket < THRESHOLD_M_BUCKETS)
      fill_in_k_buckets(th, side, w);
    else
      fill_in_m_buckets(th, side, w);
#else
    else abort();
#endif
}

void * fill_in_buckets_both(thread_data_ptr th)
{
    fill_in_buckets_one_side(th, ALGEBRAIC_SIDE);
    fill_in_buckets_one_side(th, RATIONAL_SIDE);
    return NULL;
}
/* }}} */
