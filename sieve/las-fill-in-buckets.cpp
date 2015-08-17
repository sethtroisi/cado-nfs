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
#include "las-norms.h"
#include "las-smallsieve.h"
#include "las-plattice.h"

#ifdef USE_CACHEBUFFER
#include "cachebuf.h"
#endif

// FIXME: this function of las.cpp should be somewhere
template <typename HINT>
void apply_one_bucket (unsigned char *S,
    const bucket_array_t<1, HINT> &BA, const int i,
    const fb_part *fb, where_am_I_ptr w);
void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2);
int factor_survivors (thread_data *th, int N, where_am_I_ptr w MAYBE_UNUSED);

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

template <int LEVEL>
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
        pool.add_task(fill_in_buckets_one_slice<LEVEL>, param, 0);
        slices_pushed++;
    }
    for (slice_index_t slices_completed = 0; slices_completed < slices_pushed; slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
    }
}

void fill_in_buckets_both(thread_pool &pool, thread_workspaces &ws, sieve_info_srcptr si)
{
    plattice_enumerate_t::set_area(si->conf->logI, si->J);
    for (int side = 0; side < 2; ++side) {
        switch (si->toplevel) {
            case 1:
                fill_in_buckets_one_side<1>(pool, ws,
                        si->sides[side]->fb->get_part(si->toplevel), si, side);
                break;
            case 2:
                fill_in_buckets_one_side<2>(pool, ws,
                        si->sides[side]->fb->get_part(si->toplevel), si, side);
                break;
            case 3:
                fill_in_buckets_one_side<3>(pool, ws,
                        si->sides[side]->fb->get_part(si->toplevel), si, side);
                break;
            default:
                ASSERT_ALWAYS(0);
        }
    }
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
    sieve_info_ptr si,
    precomp_plattice_t precomp_plattice,
    thread_data *th)
{
  ASSERT_ALWAYS(LEVEL > 0);

  where_am_I w;
  WHERE_AM_I_UPDATE(w, si, si);
  WHERE_AM_I_UPDATE(w, N, first_region0_index);

  double max_full MAYBE_UNUSED = 0.0;

  for (int side = 0; side < 2; ++side) {
    WHERE_AM_I_UPDATE(w, side, side);
    /* FIRST: Downsort what is coming from the level above, for this
     * bucket index */
    // All these BA are global stuff; see reservation_group.
    // We reserve those where we write, and access the ones for
    // reading without reserving. We require that things at level
    // above is finished before entering here.
    bucket_array_t<LEVEL, longhint_t> & BAout =
      ws.reserve_BA<LEVEL, longhint_t>(side);
    BAout.reset_pointers();
    BAout.add_slice_index(0); // FIXME: This is a fake slice_index!
    // The data that comes from fill-in bucket at level above:
    {
      const bucket_array_t<LEVEL+1,shorthint_t> * BAin
        = ws.cbegin_BA<LEVEL+1,shorthint_t>(side);
      while (BAin != ws.cend_BA<LEVEL+1,shorthint_t>(side)) {
        downsort<LEVEL+1>(BAout, *BAin, bucket_index, w);
        BAin++;
      }
    }

    const int toplevel = si->toplevel;
    if (LEVEL < toplevel - 1) {
      // What comes from already downsorted data above:
      const bucket_array_t<LEVEL+1,longhint_t> * BAin
        = ws.cbegin_BA<LEVEL+1,longhint_t>(side);
      while (BAin != ws.cend_BA<LEVEL+1,longhint_t>(side)) { 
        downsort<LEVEL+1>(BAout, *BAin, bucket_index, w);
        BAin++;
      }
    }
    ws.release_BA<LEVEL,longhint_t>(side, BAout);

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL, longhint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);

    /* SECOND: fill in buckets at this level, for this region. */
    bucket_array_t<LEVEL,shorthint_t> & BAfill =
      ws.reserve_BA<LEVEL, shorthint_t>(side);
    BAfill.reset_pointers();
    for (typename std::vector<plattices_vector_t *>::iterator pl_it =
            precomp_plattice[side][LEVEL].begin();
        pl_it != precomp_plattice[side][LEVEL].end();
        pl_it++) {
      WHERE_AM_I_UPDATE(w, i, (*pl_it)->get_index());
      fill_in_buckets<LEVEL>(BAfill, si, *pl_it, w);
    }
    ws.release_BA<LEVEL,shorthint_t>(side, BAfill);

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL,shorthint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);
  }

  /* RECURSE */
  if (LEVEL > 1) {
    for (unsigned int i = 0; i < si->nb_buckets[LEVEL]; ++i) {
      uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
      uint32_t N = first_region0_index + i*(BRS[LEVEL]/BRS[1]);
      downsort_tree<LEVEL-1>(i, N, ws, si, precomp_plattice, th);
    }
  } else {
    /* PROCESS THE REGIONS AT LEVEL 0 */
    // FIXME: This block duplicates a lot of code with
    // process_bucket_region in las.cpp.
    where_am_I w MAYBE_UNUSED;

    unsigned char * S[2];
    for(int side = 0 ; side < 2 ; side++) {
      thread_side_data &ts = th->sides[side];
      S[side] = ts.bucket_region;
    }

    unsigned char *SS = th->SS;
    memset(SS, 0, BUCKET_REGION_1);
    for (uint32_t ii = 0; ii < si->nb_buckets[LEVEL]; ++ii) {
      uint32_t i = first_region0_index + ii;
      WHERE_AM_I_UPDATE(w, N, i);
      // FIXME: for the descent, could early abort here.

      for (int side = 0; side < 2; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        sieve_side_info_ptr s = si->sides[side];
        thread_side_data &ts = th->sides[side];

        // Init norms
#ifdef SMART_NORM
        init_norms_bucket_region(S[side], i, si, side, 1);
#else
        init_norms_bucket_region(S[side], i, si, side, 0);
#endif
        // FIXME Invalidate the first row except (1,0)

        // Apply buckets that have just been filled-in
        const bucket_array_t<1, shorthint_t> *BA =
          th->ws->cbegin_BA<1, shorthint_t>(side);
        const bucket_array_t<1, shorthint_t> * const BA_end =
          th->ws->cend_BA<1, shorthint_t>(side);
        for (; BA != BA_end; BA++)  {
          apply_one_bucket(SS, *BA, ii, ts.fb->get_part(1), w);
        }

        // Apply downsorted buckets.
        // TODO: there is only one, no ? Why do we iterate here?
        const bucket_array_t<1, longhint_t> *BAd =
          th->ws->cbegin_BA<1, longhint_t>(side);
        const bucket_array_t<1, longhint_t> * const BAd_end =
          th->ws->cend_BA<1, longhint_t>(side);
        for (; BAd != BAd_end; BAd++)  {
          // FIXME: the updates could come from part 3 as well, not only
          // part 2.
          apply_one_bucket(SS, *BAd, ii, ts.fb->get_part(2), w);
        }

        SminusS(S[side], S[side] + BUCKET_REGION_1, SS);

        // Sieve small primes
        sieve_small_bucket_region(SS, i, s->ssd, ts.ssdpos, si, side, w);
        SminusS(S[side], S[side] + BUCKET_REGION_1, SS);
      }

      // Factor survivors
      las_report_ptr rep = th->rep;
      rep->reports += factor_survivors (th, i, w);

      // Reset resiving data
      for(int side = 0 ; side < 2 ; side++) {
        sieve_side_info_ptr s = si->sides[side];
        thread_side_data &ts = th->sides[side];
        int * b = s->fb_parts_x->rs;
        memcpy(ts.rsdpos, ts.ssdpos + b[0], (b[1]-b[0]) * sizeof(int64_t));
      }
    }
  }
}

/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(uint32_t bucket_index MAYBE_UNUSED,
  uint32_t first_region0_index MAYBE_UNUSED,
  thread_workspaces &ws MAYBE_UNUSED,
  sieve_info_ptr si MAYBE_UNUSED,
  precomp_plattice_t precomp_plattice MAYBE_UNUSED,
  thread_data *th MAYBE_UNUSED)
{
    ASSERT_ALWAYS(0);
}

// other fake instances to please level-2 instanciation.
template <>
bucket_array_t<3, longhint_t>::bucket_array_t()
{
    ASSERT_ALWAYS(0);
}

template <>
bucket_array_t<3, longhint_t>::~bucket_array_t()
{
    ASSERT_ALWAYS(0);
}

template <>
void downsort<3>(bucket_array_t<2, longhint_t>&,
        bucket_array_t<3, longhint_t> const&, unsigned int, where_am_I_ptr)
{
    ASSERT_ALWAYS(0);
}

template <>
reservation_array<bucket_array_t<3, longhint_t> > const&
reservation_group::cget<3, longhint_t>() const
{
    ASSERT_ALWAYS(0);
}

// Now the two exported instances

template 
void downsort_tree<1>(uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, sieve_info_ptr si,
  precomp_plattice_t precomp_plattice,
  thread_data *th);

template
void downsort_tree<2>(uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, sieve_info_ptr si,
  precomp_plattice_t precomp_plattice,
  thread_data *th);
