#include <stdint.h>

#include "buckets.h"
#include "ijvec.h"
#include "macros.h"



/* Bucket updates.
 *
 * Each update contains:
 * - pos:  the position (i, j) of the update, represented using an unsigned
 *         integer of UPDATE_POS_BITS bits.
 * - hint: the most significant coefficients of the prime p corresponding to
 *         the update, represented using an unsigned integer of
 *         UPDATE_HINT_BITS bits.
 *
 * The actual size in memory of a bucket update is UPDATE_POS_BITS +
 * UPDATE_HINT_BITS, rounded up to a multiple of 8 bits.
 *****************************************************************************/

// Size of the two fields and of an update, in bits.
#define UPDATE_POS_BITS  14
#define UPDATE_HINT_BITS  6
#define UPDATE_BITS     (UPDATE_POS_BITS + UPDATE_HINT_BITS)

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
  uint8_t d[(UPDATE_BITS+7) / 8];
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
{ return (update_t)hint << UPDATE_POS_BITS; }

// Adjoin a position to a hint-only update.
// /!\ The result is undefined if this update already contains a position.
static inline
update_t update_adj_pos(update_t update, pos_t pos)
{ return update | pos; }

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

// Return the number of updates in bucket k.
static inline
unsigned bucket_size(buckets_srcptr buckets, unsigned k)
{ return buckets->end[k] - buckets->begin[k]; }

// Return 1 if all updates of bucket k have been processed.
static inline
int bucket_is_done(buckets_srcptr buckets, unsigned k)
{ return buckets->iter[k] == buckets->end[k]; }

// Push an update into bucket k.
static inline
void bucket_push(buckets_ptr buckets, unsigned k, update_t update)
{ *buckets->end[k]++ = update_pack(update); }

// Return the next update in bucket k.
static inline
update_t bucket_next(buckets_ptr buckets, unsigned k)
{ return update_unpack(*buckets->iter[k]++); }


// Initialize structure and allocate buckets.
void buckets_init(buckets_ptr buckets, unsigned I, unsigned J, unsigned size)
{
  unsigned n     = 1 + ((ijvec_get_max_pos(I, J)-1) >> UPDATE_POS_BITS);
  buckets->n     = n;
  buckets->size  = size;
  buckets->begin = (update_packed_t **)malloc(n * sizeof(update_packed_t *));
  buckets->iter  = (update_packed_t **)malloc(n * sizeof(update_packed_t *));
  buckets->end   = (update_packed_t **)malloc(n * sizeof(update_packed_t *));
  ASSERT_ALWAYS(buckets->begin != NULL);
  ASSERT_ALWAYS(buckets->iter  != NULL);
  ASSERT_ALWAYS(buckets->end   != NULL);
  for (unsigned k = 0; k < buckets->n; ++k) {
    buckets->begin[k] = buckets->iter[k] = buckets->end[k] =
      (update_packed_t *)malloc(size * sizeof(update_packed_t));
    ASSERT_ALWAYS(buckets->begin[k] != NULL);
  }
}


// Print information about the buckets
void print_bucket_info(buckets_ptr buckets)
{
  printf("# bucket info:\n");
  printf("#   size of buckets = %u\n", buckets->size);
  printf("#   size of bucket-region (aka 1<<UPDATE_POS_BITS) = %d\n",
          1<<UPDATE_POS_BITS);
  printf("#   nb of buckets   = %u\n", buckets->n);
  printf("#   number of bits for the hint = %d\n", UPDATE_HINT_BITS);
  printf("#   bit-size of a bucket-update = %d\n", UPDATE_BITS);
}


// Re-initialize the bucket pointers so as to start with empty buckets.
void buckets_reset(buckets_ptr buckets)
{
  for (unsigned k = 0; k < buckets->n; ++k)
    buckets->iter[k] = buckets->end[k] = buckets->begin[k];
}


// Clean up memory.
void buckets_clear(buckets_ptr buckets)
{
  for (unsigned k = 0; k < buckets->n; ++k)
    free(buckets->begin[k]);
  free(buckets->begin);
  free(buckets->iter);
  free(buckets->end);
}


// Return the size of a bucket region.
unsigned bucket_region_size()
{ return 1u << UPDATE_POS_BITS; }



/* Array of buckets: filling buckets.
 *****************************************************************************/

// We are going recursive, here, so let's try not to have too many
// parameters to be passed, since they might end-up being pushed/popped
// a lot of times.

typedef struct {
  ijvec_ptr V;
  int I;
  int J;
  ijbasis_srcptr basis;
  buckets_ptr    buckets;
  update_t       update_hint;
} vvs_param_t;


static inline void handle_V(vvs_param_t *vvsp) {
  if (!ij_is_monic(vvsp->V->j) && !ij_is_zero(vvsp->V->j))
    return;
  ijpos_t pos = ijvec_get_pos(vvsp->V, vvsp->I, vvsp->J);
  unsigned b  = pos >> UPDATE_POS_BITS;
  pos_t    p  = pos & (((pos_t)1<<UPDATE_POS_BITS)-1);
  update_t u  = update_adj_pos(vvsp->update_hint, p);
  bucket_push(vvsp->buckets, b, u);
}

// visit the affine vector space V + generated by base vectors [0..k-1]
// V is changed during the execution, but is put back to its value before
// return.
static void visit_vector_space(vvs_param_t *vvsp, unsigned int k) {
    // TODO: is it better to store the orig value to save an add ?
    // In charact 2, certainly not, but in charact 3, probably yes.

    // Unrolled small cases
    if (k == 3) {
        for (int i = 0; i < FP_SIZE; ++i) {
            for (int i = 0; i < FP_SIZE; ++i) {
                for (int i = 0; i < FP_SIZE; ++i) {
                    handle_V(vvsp);
                    ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[0]);
                }
                ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[1]);
            }
            ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[2]);
        }
        return;
    }
    if (k == 2) {
        for (int i = 0; i < FP_SIZE; ++i) {
            for (int i = 0; i < FP_SIZE; ++i) {
                handle_V(vvsp);
                ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[0]);
            }
            ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[1]);
        }
        return;
    }
    if (k == 1) {
        for (int i = 0; i < FP_SIZE; ++i) {
            handle_V(vvsp);
            ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[0]);
        }
        return;
    }
    if (k == 0)
        return;

    // Other cases, using recursivity
    visit_vector_space(vvsp, k-1);
    for (int i = 1; i < FP_SIZE; ++i) {
        ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[k-1]);
        visit_vector_space(vvsp, k-1);
    }
    ijvec_add(vvsp->V, vvsp->V, vvsp->basis->v[k-1]);
}

// In the case of sublattices, compute the starting point for the sieve
// by gothp for the current sublattice.
// If there is no starting point return 0.
static int compute_starting_point(ijvec_ptr V0, ijbasis_ptr euclid,
                                  sublat_srcptr sublat)
{
    if (!use_sublat(sublat)) {
        ijvec_set_zero(V0);
        return 1;
    }
    int hatI = euclid->I;
    int hatJ = euclid->J;
    // TODO: FIXME: WARNING: this whole block is to be rewritten
    // completely!

    // There must be some Thm that says that a combination of
    // the first two or of the second two basis vectors is
    // enough to find a valid starting point (assuming that
    // the basis is sorted in increasing order for deg j).
    // Cf a .tex that is to be written.
    // If p is too large, then the code below is broken (and
    // anyway, this is a weird idea to sieve with such
    // parameters).

    // Try with the first 2 fundamental vectors and if this does not
    // work, try with the second 2 vectors.
    // TODO: question: does it ensures that vectors will be
    // visited in increasing order of j ?
    // FIXME: What if there are only 2 fundamental vectors?
    //
    //

    // Naive approach: check all combinations.
    ij_t ijmod, xi, yi;
    ij_set_16(ijmod, sublat->modulus);
    ij_set_16(xi, sublat->lat[sublat->n][0]);
    ij_set_16(yi, sublat->lat[sublat->n][1]);
    int found = 0;

    // case of just one vector in euclid.vec
    if (euclid->dim == 1) {
      ij_t rem;
      ij_rem(rem, euclid->v[0]->i, ijmod);
      if (ij_eq(rem, xi)) {
        ij_rem(rem, euclid->v[0]->j, ijmod);
        if (ij_eq(rem, yi)) {
          // Got a valid point!
          ij_sub(V0->i, euclid->v[0]->i, xi);
          ij_div(V0->i, V0->i, ijmod);
          ij_sub(V0->j, euclid->v[0]->j, yi);
          ij_div(V0->j, V0->j, ijmod);
          found = 1;
        }
      }
    }

    for (unsigned int ind = 0; (!found) && ind < euclid->dim-1; ++ind) {
      for (int i0 = 0; (!found) && i0 < 2; ++i0)
        for (int i1 = 0; (!found) && i1 < 2; ++i1)
          for (int i2 = 0; (!found) && i2 < 2; ++i2)
            for (int i3 = 0; (!found) && i3 < 2; ++i3){
              ijvec_t W, tmp;
              ijvec_set_zero(W);
              int i01 = i0 ^ i1;
              int i23 = i2 ^ i3;
              if (i0) ijvec_add(W,W,euclid->v[ind]);
              if (i2) ijvec_add(W,W,euclid->v[ind+1]);
              if (i01) {
                ijvec_mul_ti(tmp,euclid->v[ind],1);
                ijvec_add(W, W, tmp);
              }
              if (i23) {
                ijvec_mul_ti(tmp,euclid->v[ind+1],1);
                ijvec_add(W, W, tmp);
              }
              if ((ij_deg(W->i) >= hatI)
                  || (ij_deg(W->j) >= hatJ))
                continue;
              ij_t rem;
              ij_rem(rem, W->i, ijmod);
              if (!ij_eq(rem, xi))
                continue;
              ij_rem(rem, W->j, ijmod);
              if (!ij_eq(rem, yi))
                continue;
              // Got a valid point!
              ij_sub(W->i, W->i, xi);
              ij_div(V0->i, W->i, ijmod);
              ij_sub(W->j, W->j, yi);
              ij_div(V0->j, W->j, ijmod);
              found = 1;
            }
    }
    return found;
}


// Fill the buckets with updates corresponding to divisibility by elements of
// the factor base.
void buckets_fill(buckets_ptr buckets, factor_base_srcptr FB,
                  sublat_srcptr sublat, unsigned I, unsigned J)
{
  ijbasis_t basis;
  ijbasis_t euclid;
  unsigned  hatI, hatJ;
  hatI = I + sublat->deg;
  hatJ = J + sublat->deg;

  ijbasis_init(basis,     I,    J);
  ijbasis_init(euclid, hatI, hatJ);

  unsigned i;
  for (i = 0; i < FB->n && FB->elts[i]->degp < I; ++i);

  for (fbideal_srcptr gothp = FB->elts[i]; i < FB->n; ++i, ++gothp) {
    // List of cases that are not handled yet:
    if (gothp->proj) continue;
    if (use_sublat(sublat) && gothp->degp == 1) continue;

    ijbasis_compute(euclid, basis, gothp);
    ijvec_t V0;
    if (use_sublat(sublat)) {
      int st = compute_starting_point(V0, euclid, sublat);
      if (!st)
        continue; // next factor base prime.
    } else {
      ijvec_set_zero(V0);
    }

    vvs_param_t vvs = { .V = V0,
                        .I = I,
                        .J = J,
                        .basis   = basis,
                        .buckets = buckets,
                        .update_hint = update_set_hint(gothp->degp) };
    visit_vector_space(&vvs, basis->dim);
  }

  //for (unsigned k = 0; k < buckets->n; ++k)
  //  printf("# #updates[%u] = %u\n", k, bucket_size(buckets, k));

  ijbasis_clear(euclid);
  ijbasis_clear(basis);
}


// Apply all the updates from a given bucket to the sieve region S.
void bucket_apply(uint8_t *S, buckets_ptr buckets, unsigned k)
{
  while (!bucket_is_done(buckets, k)) {
    update_t u = bucket_next(buckets, k);
    pos_t    p = update_get_pos(u);
    hint_t   h = update_get_hint(u);
    ASSERT((k == 0 && p == 0) || S[p] >= h);
    S[p] -= h;
  }
}
