#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>

#include "buckets.h"
#include "ijvec.h"
#include "qlat.h"
#include "macros.h"
#include "gray.h"



/* Bucket updates.
 *
 * Each update contains:
 * - pos:  the position (i, j) of the update, represented using an unsigned
 *         integer of UPDATE_POS_BITS bits.
 * - hint: a prime hint, to quickly recover the irreducible polynomnial of
 *         the corresponding prime ideal. It is used for resieving survivors.
 *         If BUCKET_RESIEVE is not defined, then this hint is nothing,
 *         otherwise, it is the index of the ideal in the factor base
 *         (stored as a difference of index with previous).
 *         The hint is represented using an unsigned integer
 *         of UPDATE_HINT_BITS bits.
 *
 * The actual size in memory of a bucket update is UPDATE_POS_BITS +
 * UPDATE_HINT_BITS, rounded up to a multiple of UPDATE_ALIGN bytes.
 *****************************************************************************/

// Size of the two fields and of a bucket update, in bits.
#define UPDATE_POS_BITS  16
#ifdef BUCKET_RESIEVE
#define UPDATE_HINT_BITS  16
#else
#define UPDATE_HINT_BITS  0
#endif
#define UPDATE_BITS     (UPDATE_POS_BITS + UPDATE_HINT_BITS)

// Memory alignment of bucket updates, in bytes.
#define UPDATE_ALIGN      1

// Smallest unsigned integer in which an update will fit.
#if   UPDATE_BITS <= 16
typedef uint16_t update_t;
#elif UPDATE_BITS <= 32
typedef uint32_t update_t;
#else
typedef uint64_t update_t;
#endif

// In memory, an update will be packed into an array of bytes.
struct __update_packed_struct {
  uint8_t d[((UPDATE_BITS-1) / (8*UPDATE_ALIGN) + 1) * UPDATE_ALIGN];
};

// Internal type for converting between normal and packed representations
// of updates.
typedef union {
  update_t        update;
  update_packed_t packed;
} __update_conv_t;

// Pack an update.
static inline
update_packed_t update_pack(update_t update)
{ __update_conv_t conv = { .update = update };
  return conv.packed; }

// Unpack an update.
static inline
update_t update_unpack(update_packed_t packed)
{ __update_conv_t conv = { .packed = packed };
  return conv.update; }


// Types for positions and hints as exported by the functions.
typedef unsigned pos_t;
typedef unsigned hint_t;

// Construct an update from a prime hint only, leaving space for adjoining the
// position later on.
static inline
update_t update_set_hint(hint_t hint)
{ ASSERT(hint < (hint_t)1<<UPDATE_HINT_BITS);
  return (update_t)hint << UPDATE_POS_BITS; }

// Adjoin a position to a hint-only update.
// /!\ The result is undefined if the update already contains a position.
static inline
update_t update_adj_pos(update_t update, pos_t pos)
{ ASSERT(pos < (pos_t)1<<UPDATE_POS_BITS);
  return update | pos; }

// Construct an update from a position and a prime hint.
static inline
update_t update_set(pos_t pos, hint_t hint)
{ return update_adj_pos(update_set_hint(hint), pos); }

// Retrieve the position from an update.
static inline
pos_t update_get_pos(update_t update)
{ return update & (((pos_t)1<<UPDATE_POS_BITS)-1); }

// Retrieve the prime hint from an update.
static inline
hint_t update_get_hint(update_t update)
{ return (update >> UPDATE_POS_BITS) & (((hint_t)1<<UPDATE_HINT_BITS)-1); }



/* Array of buckets: basic management.
 *****************************************************************************/

// Initialize structure and allocate buckets.
void buckets_init(buckets_ptr buckets, unsigned I, unsigned J,
                  unsigned max_size, unsigned min_degp, unsigned max_degp)
{
  // FIXME: Everything breaks in characteristic 3 if bucket regions are not of
  // size 2^8 or 2^16. This would indeed imply that bucket regions would not
  // be aligned with complete lines of the sieving area. Someday, we'll have to
  // fix this, but for the time being, let's leave it as is.
  // A way to fix this would be to allow for unaligned bucket regions:
  //    _____ ..........
  //   |     |__________: <- line j1
  //   |                |
  //   |  bucket region |
  //   |_________       |
  //   :.........|______| <- line j0
  //
  // Norm initialization and line sieving would be applied to complete lines
  // (overshooting the end of the bucket region, as in line j1), and lines
  // overlapping two bucket regions (such as line j1) will have to be copied
  // and bucket-sieved with the following bucket region (such as what is
  // happening with line j0 here).
  ASSERT(__FP_CHAR == 2 || UPDATE_POS_BITS == 16);

  unsigned n        = 1 + ((ijvec_get_max_pos(I, J)-1) >> UPDATE_POS_BITS);
  buckets->n        = n;
  buckets->max_size = max_size;
  buckets->min_degp = min_degp;
  buckets->max_degp = max_degp;
  buckets->start    =
    (update_packed_t **)malloc(n * sizeof(update_packed_t *));
  ASSERT_ALWAYS(buckets->start != NULL);
  for (unsigned k = 0; k < n; ++k) {
    buckets->start[k] =
      (update_packed_t *)malloc(max_size * sizeof(update_packed_t));
    ASSERT_ALWAYS(buckets->start[k] != NULL);
  }
  unsigned ndegp = max_degp - min_degp;
  buckets->degp_end =
    (update_packed_t **)malloc(ndegp * n * sizeof(update_packed_t *));
  ASSERT_ALWAYS(buckets->degp_end);
#ifdef BUCKET_RESIEVE
  buckets->first_hint   = (unsigned *)malloc(n * sizeof(unsigned));
  buckets->current_hint = (unsigned *)malloc(n * sizeof(unsigned));
  ASSERT_ALWAYS(buckets->first_hint);
  ASSERT_ALWAYS(buckets->current_hint);
  for (unsigned i = 0; i < n; ++i) {
    buckets->first_hint[i] = UINT_MAX;
    buckets->current_hint[i] = UINT_MAX;
  }
#endif
}


// Clean up memory.
void buckets_clear(buckets_ptr buckets)
{
  for (unsigned k = 0; k < buckets->n; ++k)
    free(buckets->start[k]);
  free(buckets->start);
  free(buckets->degp_end);
#ifdef BUCKET_RESIEVE
  free(buckets->first_hint);
  free(buckets->current_hint);
#endif
}


// Return the size of a bucket region.
unsigned bucket_region_size()
{
  return 1u << UPDATE_POS_BITS;
}


// Print information about the buckets.
void print_bucket_info(buckets_srcptr buckets)
{
  printf("# buckets info:\n");
  printf("#   nb of buckets   = %u\n", buckets->n);
  printf("#   size of buckets = %u\n", buckets->max_size);
  printf("#   size of bucket-region (ie, 2^UPDATE_POS_BITS) = %d\n",
         bucket_region_size());
  printf("#   number of bits for the hint = %d\n", UPDATE_HINT_BITS);
  printf("#   bit-size of a bucket-update = %d (rounded to %lu)\n",
         UPDATE_BITS, 8*sizeof(update_packed_t));
}



/* Array of buckets: filling buckets.
 *****************************************************************************/

// In the case of sublattices, compute the starting point for the sieve
// by gothp for the current sublattice.
// If there is no starting point return 0.
// This code is specific to GF(2).
static int compute_starting_point(ijvec_ptr V0,
    MAYBE_UNUSED ijbasis_ptr euclid, sublat_srcptr sublat)
{
  if (!use_sublat(sublat)) {
    ijvec_set_zero(V0);
    return 1;
  }
#ifdef USE_F2  
  int hatI = euclid->I;
  int hatJ = euclid->J;

  fppol8_t i0, j0;
  fppol8_set_16(i0,sublat->lat[sublat->n][0]);
  fppol8_set_16(j0,sublat->lat[sublat->n][1]);

  // case of just one vector in euclid.vec
  if (euclid->dim == 1) {
    fppol8_t rem;
    sublat_mod_ij(rem, euclid->v[0]->i);
    if (fppol8_eq(rem, i0)) {
      sublat_mod_ij(rem, euclid->v[0]->j);
      if (fppol8_eq(rem, j0)) {
        // Got a valid point!
        sublat_div_ij(V0->i, euclid->v[0]->i);
        sublat_div_ij(V0->j, euclid->v[0]->j);
        return 1;
      }
    }
    return 0;
  }

  // Try to combine the first two vectors of the euclidian sequence.
  // If this does not work, the next two should do.
  // For that, we multiply the vector (i0, j0) by the inverse of the 
  // 2x2 matrix formed by the coordinates of the base vectors.
  // since we are working modulo (t^2+t), the determinant of the matrix,
  // which is invertible, is necessarily equal to 1, so its inverse is
  // trivial.
  for (unsigned ind = 0; ind < MIN(2, euclid->dim-1); ++ind) {
    ij_t a,b,c,d;
    ij_set(a, euclid->v[ind]->i);
    ij_set(b, euclid->v[ind]->j);
    ij_set(c, euclid->v[ind+1]->i);
    ij_set(d, euclid->v[ind+1]->j);

    fppol8_t aa, bb, cc, dd;
    sublat_mod_ij(aa, a);
    sublat_mod_ij(bb, b);
    sublat_mod_ij(cc, c);
    sublat_mod_ij(dd, d);

    fppol8_t alpha, beta;
    sublat_mul(alpha, dd, i0);
    sublat_addmul(alpha, alpha, cc, j0);
// TODO: this shortcut is invalid, but one could fix it.
#if 0
    // In the first pass, v0 has max degree, so if alpha has degree 1,
    // it won't work. 
    if ((ind == 0) && (alpha[0] & 2u))
      continue;
#endif
    sublat_mul(beta, bb, i0);
    sublat_addmul(beta, beta, aa, j0);

    ij_t tmp1, tmp2;
    ij_mul_sublat(tmp1, a, alpha);
    ij_mul_sublat(tmp2, c, beta);
    ij_add(tmp1, tmp1, tmp2);
    if (ij_deg(tmp1) >= hatI)
      continue;
    sublat_div_ij(V0->i, tmp1);
    ij_mul_sublat(tmp1, b, alpha);
    ij_mul_sublat(tmp2, d, beta);
    ij_add(tmp1, tmp1, tmp2);
    if (ij_deg(tmp1) >= hatJ)
      continue;
    sublat_div_ij(V0->j, tmp1);
    return 1;
  }
  return 0;
#endif

  // should never go there.
  ASSERT_ALWAYS(0);
}

// Push an update into the corresponding bucket.
static inline
void buckets_push_update(MAYBE_UNUSED buckets_ptr buckets, 
                         update_packed_t **ptr, hint_t hint,
                         ijvec_srcptr v, unsigned I, unsigned J)
{
  ijpos_t  pos = ijvec_get_pos(v, I, J);
  unsigned k   = pos >> UPDATE_POS_BITS;
  pos_t    p   = pos & (((pos_t)1<<UPDATE_POS_BITS)-1);
  ASSERT(ptr[k] - buckets->start[k] < buckets->max_size);
#ifdef BUCKET_RESIEVE
  // TODO: this is a very critical loop. Can we afford this branch???
  if (buckets->first_hint[k] == UINT_MAX) {
      buckets->first_hint[k] = hint;
      buckets->current_hint[k] = hint;
      hint = 0;   // putting 0 allows to add it to the current_hint.
  } else {
      unsigned h = hint;
      hint = hint - buckets->current_hint[k];
      ASSERT(hint < (1U<<UPDATE_HINT_BITS));
      buckets->current_hint[k] = h;
  }
#endif
  *ptr[k]++ = update_pack(update_set(p, hint));
}


// Fill the buckets with updates corresponding to divisibility by elements of
// the factor base.
void buckets_fill(buckets_ptr buckets, large_factor_base_srcptr FB,
        sublat_srcptr sublat, unsigned I, unsigned J, qlat_srcptr qlat)
{
  ijbasis_t basis;
  ijbasis_t euclid;
  unsigned  hatI, hatJ;
  hatI = I + sublat->deg;
  hatJ = J + sublat->deg;

  // The bucket sieve requires all considered ideals to be of degree larger
  // than I.
  ASSERT_ALWAYS(buckets->min_degp >= I);

  ijbasis_init(basis,     I,    J);
  ijbasis_init(euclid, hatI, hatJ);

  // Pointers to the current writing position in each bucket.
  update_packed_t **ptr =
    (update_packed_t **)malloc(buckets->n * sizeof(update_packed_t *));
  ASSERT_ALWAYS(ptr != NULL);
  for (unsigned k = 0; k < buckets->n; ++k)
    ptr[k] = buckets->start[k];
#ifdef BUCKET_RESIEVE
  for (unsigned k = 0; k < buckets->n; ++k)
    buckets->first_hint[k] = UINT_MAX;
#endif

  // Go through the factor base by successive deg(gothp).
  // We should have no small prime in this factor base.
  large_fbideal_srcptr gothp = FB->elts[0];
  ASSERT_ALWAYS(fbideal_deg(gothp) >= buckets->min_degp);
  unsigned i = 0;
  for (unsigned degp = buckets->min_degp; degp < buckets->max_degp; ++degp) {
    ASSERT(degp >= I);

    // Go through each prime ideal of degree degp.
    for (; i < FB->n && fbideal_deg(gothp) == degp; ++i, ++gothp) {
      fbprime_t lambda;
      compute_lambda(lambda, gothp->p, gothp->r, qlat);

      if (fbprime_eq(lambda, gothp->p)) {
          // This is a projective root. For the moment, we skip them.
          continue;
      }
      
      hint_t hint = 0;
#ifdef BUCKET_RESIEVE
      hint = i;
#endif

      ijbasis_compute_large(euclid, basis, gothp, lambda);
      ijvec_t v;
      if (use_sublat(sublat)) {
        int st = compute_starting_point(v, euclid, sublat);
        if (!st)
          continue; // next factor base prime.
        // The first position is to be sieved when we use sublat, but
        // otherwise not, since this is (0,0).
        buckets_push_update(buckets, ptr, hint, v, I, J);
      }
      else
        ijvec_set_zero(v);

      // Size of Gray codes to use for the inner loop.
#     if   defined(USE_F2)
#       define ENUM_LATTICE_UNROLL 8
#     elif defined(USE_F3)
#       define ENUM_LATTICE_UNROLL 5
#     endif

      // Unrolled p-ary Gray code of size ENUM_LATTICE_UNROLL.
      static const uint8_t gray[] = { GRAY(ENUM_LATTICE_UNROLL) };
      unsigned             ngray  = GRAY_LENGTH(ENUM_LATTICE_UNROLL);

      // We only need the "monic" Gray code for the first iteration. Jump
      // directly there in the array.
      unsigned gray_dim = MIN(basis->dim, ENUM_LATTICE_UNROLL);
      unsigned i0       = ngray - GRAY_LENGTH(gray_dim) / (__FP_SIZE-1);

      ij_t s, t;
      ij_set_zero(t);
      int rc = basis->dim > ENUM_LATTICE_UNROLL;
      do {
        // Inner-level Gray code enumeration: just go through the Gray code
        // array, each time adding the indicated basis vector.
        for (unsigned ii = i0; ii < ngray; ++ii) {
          ijvec_add(v, v, basis->v[gray[ii]]);
          buckets_push_update(buckets, ptr, hint, v, I, J);
        }
        i0 = 0;

        // Outer-level Gray code enumeration: using ij_monic_set_next, the
        // degree of the difference with the previous one indicates which basis
        // vector should be added to the current lattice point.
        // rc is set to 0 when all vectors have been enumerated.
        ij_set(s, t);
        rc = rc && ij_monic_set_next(t, t, basis->dim-ENUM_LATTICE_UNROLL);
        if (rc) {
          ij_diff(s, s, t);
          ijvec_add(v, v, basis->v[ij_deg(s)+ENUM_LATTICE_UNROLL]);
          buckets_push_update(buckets, ptr, hint, v, I, J);
        }
      } while (rc);
    }

    // Mark the last position for this degree in the degp_end array.
    for (unsigned i = degp - buckets->min_degp, k = 0; k < buckets->n; ++k)
      buckets->degp_end[i*buckets->n + k] = ptr[k];
  }

  //for (unsigned k = 0; k < buckets->n; ++k)
  //  printf("# #updates[%u] = %u\n", k, ptr[k]-buckets->start[k]);

  free(ptr);
  ijbasis_clear(euclid);
  ijbasis_clear(basis);
}

// Apply all the updates from a given bucket to the sieve region S.
void bucket_apply(uint8_t *S, buckets_srcptr buckets, unsigned k)
{
  // Pointer to the current reading position in the bucket.
  update_packed_t *ptr = buckets->start[k];

  // For debugging purposes.
  MAYBE_UNUSED ijpos_t pos0 = k*bucket_region_size();

  // Iterate through each group of updates having same deg(gothp).
  update_packed_t **end = buckets->degp_end+k;
  for (unsigned degp = buckets->min_degp; degp < buckets->max_degp;
       ++degp, end += buckets->n) {
    while (ptr != *end) {
      update_t update = update_unpack(*ptr++);
      pos_t    pos    = update_get_pos(update);
#ifdef TRACE_POS
      if (pos0+pos == TRACE_POS) {
        fprintf(stderr, "TRACE_POS(%lu): [%u]\n", pos0+pos, degp);
        fprintf(stderr, "TRACE_POS(%lu): degnorm is now %d\n",
                pos0+pos, S[pos]-degp);
      }
#endif
#ifndef NDEBUG
      if (S[pos] < degp)
        fprintf(stderr, "faulty pos is %lu\n", pos0+pos);
      ASSERT(S[pos] >= degp);
#endif
      S[pos] -= degp;
    }
  }
}

#ifdef BUCKET_RESIEVE

static int mycmp(const void *p1, const void *p2) {
  __replayable_update_struct const * r1 = (__replayable_update_struct const *) p1;
  __replayable_update_struct const * r2 = (__replayable_update_struct const *) p2;
  if (r1->pos < r2->pos)
      return -1;
  if (r1->pos > r2->pos)
      return 1;
  return 0;
}

void bucket_prepare_replay(replayable_bucket_ptr bb,
        buckets_srcptr buckets, uint8_t *S, unsigned k)
{
  // Pointer to the current reading position in the bucket.
  update_packed_t *ptr = buckets->start[k];
  unsigned current_hint = buckets->first_hint[k];

  bb->n = 0;

  // Iterate through each group of updates having same deg(gothp).
  update_packed_t **end = buckets->degp_end+k;
  for (unsigned degp = buckets->min_degp; degp < buckets->max_degp;
       ++degp, end += buckets->n) {
    while (ptr != *end) {
      update_t update = update_unpack(*ptr++);
      pos_t    pos    = update_get_pos(update);
      hint_t   hint   = update_get_hint(update);
      current_hint += hint;
      if (S[pos] == 255)
        continue;
      bb->b[bb->n].pos = pos;
      bb->b[bb->n].hint = current_hint;
      bb->n++;
    }
  }

  // sort according to position to facilitate search
  qsort((void *)bb->b, bb->n, sizeof(__replayable_update_struct), mycmp);

  // reset buckets for next turn
  buckets->first_hint[k] = 0;
}

void bucket_apply_at_pos(fppol_ptr norm, ijpos_t pp, 
        replayable_bucket_srcptr buckets, large_factor_base_srcptr FB)
{
  // Find first occurence, if it exists:
  __replayable_update_struct key;
  key.pos = pp;
  __replayable_update_struct *found;
  found = bsearch((void *)&key, (void *)buckets->b, buckets->n,
          sizeof(__replayable_update_struct), mycmp);
  if (found == NULL)
    return;
  // There is no guarantee that bsearch returns the first one, so we have
  // to adjust.
  while (found > buckets->b) {
    __replayable_update_struct *nfound = found -1;
    if (nfound->pos == pp)
      found = nfound;
    else
      break;
  }

  while (found < buckets->b+buckets->n) {
    if (found->pos != pp)
      return;
    fppol_t r;
    fppol_init(r);
    fppol_set_fbprime(r, FB->elts[found->hint]->p);
    fppol_divrem(norm, r, norm, r);
    ASSERT_ALWAYS(fppol_is_zero(r));
    fppol_clear(r);
    found++;
  }
}
#endif
