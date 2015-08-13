#include "cado.h"

#include "fb.h"
#include "utils.h"           /* lots of stuff */
#include "bucket.h"
#include "modredc_ul.h"
#include "modredc_2ul2.h"
#include "threadpool.h"
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

#include "las-plattice.h"

// FIXME
// Those three lines should go to las-plattice.c, but there is none for
// the moment. And since the present file is the only one that uses
// those, then we put the definition of static variables there.
uint32_t plattice_enumerate_t::maskI;
plattice_x_t plattice_enumerate_t::even_mask;
plattice_x_t plattice_enumerate_t::area;



/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/**************************************************************************
 * Global DEFINEs for fill_in_buckets, fill_in_k_buckets, fill_in_m_buckets 
 **************************************************************************/

#ifdef HAVE_SSE2							
#define FILL_BUCKET_PREFETCH(PT) do {				\
    _mm_prefetch((char *)(PT), _MM_HINT_T0);			\
  } while (0)
#else
#define FILL_BUCKET_PREFETCH(PT)
#endif


#ifdef USE_CACHEBUFFER
DECLARE_CACHE_BUFFER(bucket_update_shorthint_t, 256)
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


/* {{{ */
template <int LEVEL>
void
fill_in_buckets(bucket_array_t<LEVEL, shorthint_t> &orig_BA,
                sieve_info_srcptr const si MAYBE_UNUSED,
                plattices_vector_t *plattices_vector,
                where_am_I_ptr w)
{
  const slice_index_t slice_index = plattices_vector->get_index();
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);
  
  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  for (plattices_vector_t::iterator pl_it = plattices_vector->begin();
       pl_it != plattices_vector->end(); pl_it++) {

    const slice_offset_t hint = pl_it->get_hint();
    WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
    const fb_slice_interface * slice =
        si->sides[w->side]->fb->get_slice(slice_index);
    const fbprime_t p = slice->get_prime(hint); 
    WHERE_AM_I_UPDATE(w, p, p);
#else
    const fbprime_t p = 0;
#endif

    // FIXME
#if 0 // SHOULD be done by one of the callers!
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
      BA.push_update(x, p, hint, slice_index, w);
      continue;
    }
#endif

    /* Now, do the real work: the filling of the buckets */
    while (!pl_it->finished()) {
      if (LIKELY(pl_it->probably_coprime()))
        BA.push_update(pl_it->get_x(), p, hint, slice_index, w);
      pl_it->next();
    }
    pl_it->advance_to_next_area();
  }
  orig_BA.move(BA);
}

class fill_in_buckets_parameters: public task_parameters {
public:
  thread_workspaces &ws;
  const int side;
  sieve_info_srcptr const si;
  const fb_slice_interface * const slice;
  fill_in_buckets_parameters(thread_workspaces &_ws, const int _side, sieve_info_srcptr const _si, const fb_slice_interface *_slice)
  : ws(_ws), side(_side), si(_si), slice(_slice) {}
};

template<int LEVEL>
task_result *
fill_in_buckets_one_slice(const task_parameters *const _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    where_am_I w;
    WHERE_AM_I_UPDATE(w, si, param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());

    /* Do the root transform and lattice basis reduction for this factor base slice */
    plattices_vector_t *plattices_vector =
        param->slice->make_lattice_bases(param->si->qbasis, param->si->conf->logI);
    /* Get an unused bucket array that we can write to */
    bucket_array_t<LEVEL, shorthint_t> &BA = param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
    /* Fill the buckets */
    fill_in_buckets<LEVEL>(BA, param->si, plattices_vector, w);
    /* Release bucket array again */
    param->ws.release_BA(param->side, BA);
    delete plattices_vector;
    delete param;
    return new task_result;
}

static void
fill_in_buckets_one_side(thread_pool &pool, thread_workspaces &ws, const fb_part *fb, sieve_info_srcptr const si, const int side)
{
    /* Process all slices in this factor base part */
    const fb_slice_interface *slice;
    slice_index_t slices_pushed = 0;
    for (slice_index_t slice_index = fb->get_first_slice_index();
         (slice = fb->get_slice(slice_index)) != NULL;
         slice_index++) {
        fill_in_buckets_parameters *param = new fill_in_buckets_parameters(ws, side, si, slice);
        // FIXME: here, this should be toplevel, not 1.
        pool.add_task(fill_in_buckets_one_slice<1>, param, 0);
        slices_pushed++;
    }
    for (slice_index_t slices_completed = 0; slices_completed < slices_pushed; slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
    }
}

void fill_in_buckets_both(thread_pool &pool, thread_workspaces &ws, sieve_info_srcptr si)
{
    int part = si->toplevel;
    plattice_enumerate_t::set_area(si->conf->logI, si->J);
    for (int side = 0; side < 2; ++side)
        fill_in_buckets_one_side(pool, ws, si->sides[side]->fb->get_part(part), si, side);
}
/* }}} */


// first_region0_index is a way to remember where we are in the tree.
// The depth-first is a way to precess all the the region of level 0 in
// increasing order of j-value.
// first_region0_index * nb_lines_per_region0 therefore gives the j-line
// where we are. This is what is called N by WHERE_AM_I and friends.
template <int LEVEL>
void
downsort_tree(uint32_t bucket_index,
        uint32_t first_region0_index,
        thread_workspaces &ws,
        sieve_info_srcptr si,
        typename std::vector<plattices_vector_t *> precomp_plattice[2][FB_MAX_PARTS])
{
    ASSERT_ALWAYS(LEVEL > 0);

    double max_full;

    for (int side = 0; side < 2; ++side) {
        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // FIXME: QUESTION: do we need to reserve???
        bucket_array_t<LEVEL,longhint_t> BAout = 
            ws.reserve_BA<LEVEL,longhint_t>(side);
        // The data that comes from fill-in bucket at level above:
        {
            bucket_array_t<LEVEL+1,shorthint_t> BAin =
                ws.reserve_BA<LEVEL+1,shorthint_t>(side);
            downsort<LEVEL+1>(BAout, BAin, bucket_index);
            ws.release_BA<LEVEL+1,shorthint_t>(side);
        }

        const int toplevel = MAX(si->sides[0]->fb->get_toplevel(),
                si->sides[1]->fb->get_toplevel());                              
        if (LEVEL < toplevel - 1) {
            // What comes from already downsorted data above:
            bucket_array_t<LEVEL+1,longhint_t> BAin =
                ws.reserve_BA<LEVEL+1,longhint_t>(side);
            downsort<LEVEL+1>(BAout, BAin, bucket_index);
            ws.release_BA<LEVEL+1,longhint_t>(side);
        }
        ws.release_BA<LEVEL,longhint_t>(side);
        max_full = std::max(ws.buckets_max_full<LEVEL, longhint_t>());
        ASSERT_ALWAYS(max_full <= 1.0);

        /* SECOND: fill in buckets at this level, for this region. */
        bucket_array_t<LEVEL,shorthint_t> BAin = 
            ws.reserve_BA<LEVEL,shorthint_t>(side);
        for (typename std::vector<plattices_vector_t *>::iterator pl_it =
                precomp_plattice[side][LEVEL].begin();
              pl_it != precomp_plattice[side][LEVEL].end();
              pl_it++) {
            where_am_I_ptr w;
            WHERE_AM_I_UPDATE(w, si, si);
            WHERE_AM_I_UPDATE(w, side, side);
            WHERE_AM_I_UPDATE(w, i, (*pl_it)->get_index());
            fill_in_buckets<LEVEL>(BAin, si, pl_it, w);
        }
        ws.release_BA<LEVEL,shorthint_t>(side);
        max_full = std::max(ws.buckets_max_full<LEVEL, shorthint_t>());
        ASSERT_ALWAYS(max_full <= 1.0);
    }

    /* RECURSE */
    if (LEVEL > 1) {
        for (int i = 0; i < bucket_array_t<LEVEL,shorthint_t>::n_bucket; ++i) {
            uint32_t N = first_region0_index; // FIXME + i*APPROPRIATE_REGION_SIZE.
            downsort_tree<LEVEL-1>(i, N, ws, si, precomp_plattice);
        }
    } else {
        /* PROCESS THE REGIONS AT LEVEL 0 */

        // Should be more or less exactly process_bucket_region(),
        // except that we should also apply the updates of
        //   bucket_array_t<1, longhint_t>
        // coming from downsorting.
    }
}
