#ifndef __FQ_H__
#define __FQ_H__

#include "types.h"

typedef sq_t fq_t;
typedef sq_srcptr fq_srcptr;
typedef sq_ptr fq_ptr;

typedef struct {
    sq_t q;
    unsigned int degq;
    uint64_t order;  // p^degq
} __fq_info;

typedef __fq_info fq_info_t[1];
typedef __fq_info * fq_info_ptr;
typedef const __fq_info * fq_info_srcptr;

///////////////////////////////////
// The type sq2_t is supposed to be enough to contain a polynomial
// of degree twice as large as sq. For the moment we stick to sq_t and
// sq2_t to 64 bits, but this means that fq_t will work only q of degrees
// at most 32.

typedef sq_t sq2_t;
typedef sq_ptr sq2_ptr;
typedef sq_srcptr sq2_srcptr;

static inline void sq2_mul(sq2_ptr r, sq2_srcptr p, sq2_srcptr q)
{ sq_mul(r, p, q); }
static inline void sq2_rem(sq2_ptr r, sq2_srcptr p, sq2_srcptr q)
{ sq_rem(r, p, q); }
static inline void sq2_add(sq2_ptr r, sq2_srcptr p, sq2_srcptr q)
{ sq_add(r, p, q); }

///////////////////////////////////

static inline 
void fq_info_init(fq_info_ptr Fq, sq_srcptr q)
{
    if (1+2*(sq_deg(q)-1) > __sq_SIZE) {
        fprintf(stderr, "Error: the sq2_t is not large enough to allow this size of q\n");
        abort();
    }
    sq_set(Fq->q, q);
    Fq->degq = sq_deg(q);
    Fq->order = 1;
    for (unsigned int i = 0; i < Fq->degq; ++i)
        Fq->order *= FP_SIZE;
}

static inline 
void fq_info_clear(MAYBE_UNUSED fq_info_ptr Fq)
{
    return;
}

static inline 
void fq_out(FILE *file, fq_srcptr p, MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_out(file, p);
}

static inline
void fq_add(fq_ptr r, fq_srcptr p, fq_srcptr q,
        MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_add(r, p, q);
}

static inline
void fq_sub(fq_ptr r, fq_srcptr p, fq_srcptr q,
        MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_sub(r, p, q);
}

static inline 
void fq_opp(fq_ptr r, fq_srcptr p, MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_opp(r, p);
}

static inline 
void fq_mul(fq_ptr r, fq_srcptr p, fq_srcptr q, fq_info_srcptr Fq)
{
    sq_mulmod(r, p, q, Fq->q);
}

static inline 
void fq_sqr(fq_ptr r, fq_srcptr p, fq_info_srcptr Fq)
{
    // TODO: have some sqrmod!
    // sq_sqrmod(r, p, Fq->q);
    sq_mulmod(r, p, p, Fq->q);
}

static inline 
int fq_inv(fq_ptr r, fq_srcptr p, fq_info_srcptr Fq)
{
    return sq_invmod(r, p, Fq->q);
}

static inline 
void fq_mulur(sq2_ptr r, fq_srcptr p, fq_srcptr q,
        MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq2_mul(r, p, q);
}

static inline
void fq_addur(sq2_ptr r, sq2_srcptr p, sq2_srcptr q, 
        MAYBE_UNUSED fq_info_srcptr Fq)
{
        sq2_add(r, p, q);
}       

static inline 
void fq_reduce(fq_ptr r, sq2_srcptr p, fq_info_srcptr Fq)
{
    sq2_rem(r, p,  Fq->q);
}

static inline 
void fq_set(fq_ptr r, fq_srcptr p, MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_set(r, p);
}

static inline
void fq_set_mp(fq_ptr r, fppol_srcptr p, fq_info_srcptr Fq)
{
    fppol_t q, c;
    fppol_init(q);
    fppol_init(c);
    fppol_set_sq(q, Fq->q);
    fppol_rem(c, p, q);
    sq_set_mp(r, c);
    fppol_clear(c);
    fppol_clear(q);
}

static inline
void fq_set_zero(fq_ptr r, MAYBE_UNUSED fq_info_srcptr Fq)
{
    sq_set_zero(r);
}

static inline 
void fq_set_ti(fq_ptr r, unsigned int i, fq_info_srcptr Fq)
{
    if (i >= Fq->degq) {
        fppol_t ti, q;
        fppol_init(ti);
        fppol_init(q);
        fppol_set_ti(ti, i);
        fppol_rem(ti, ti, q);
        sq_set_mp(r, ti);
        fppol_clear(ti);
        fppol_clear(q);
    } else {
        sq_set_ti(r, i);
    }
}


static inline
int fq_is_zero(fq_srcptr r, MAYBE_UNUSED fq_info_srcptr Fq)
{
    return sq_is_zero(r);
}

static inline 
void fq_set_random(fq_ptr rr, MAYBE_UNUSED fq_info_srcptr Fq)
{
    // sq_set_random(r);
#if (RAND_MAX == 2147483647)
#define RAND_BITS  31
#else             
#define RAND_BITS  15 
#endif

    fppol64_t r;

#ifdef USE_F2
    r[0] = 0;
    unsigned int b = 0;
    while ((b*RAND_BITS) < Fq->degq) {
        r[0] |= ((uint64_t)rand()) << (b*RAND_BITS);
        b++;
    }
    r[0] &= (((uint64_t)1) << (Fq->degq))-1;
#elif defined(USE_F3)
    for (int i = 0; i < __FP_BITS; ++i) {
        r[i] = 0;
        unsigned int b = 0;
        while ((b*RAND_BITS) < Fq->degq) {
            r[i] |= ((uint64_t)rand()) << (b*RAND_BITS);
            b++;
        }
        r[i] &= (((uint64_t)1) << (Fq->degq))-1;
    }
    // need to change positions where both bits are on.
    // in that case, we just erase the msb (this is not uniform random).
    uint64_t cy = ~(r[0] & r[1]);
    r[1] &= cy;
#else
#error "random not implemented for this field"
#endif
    sq_set_64(rr, r);
}


#endif   /* __FQ_H__ */
