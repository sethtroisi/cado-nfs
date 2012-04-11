#include "fqpol.h"
#include "string.h"

void fqpol_init(fqpol_ptr r)
{
    r->alloc = 0;
    r->deg = -1;
    r->coeffs = NULL;
}

void fqpol_clear(fqpol_ptr r)
{
    free(r->coeffs);
}

// Reallocate space for n coefficients
static void fqpol_realloc_lazy(fqpol_ptr r, unsigned int n)
{
    r->coeffs = (fq_t *)realloc(r->coeffs, n*sizeof(fq_t));
    ASSERT_ALWAYS(!n || r->coeffs != NULL);
    r->alloc = n;
}

// Update the degree: assume that the given degree is an upper bound, and
// check the nullity of the coefficients.
// /!\ This function works by side-effects.
static void fqpol_update_degree(fqpol_ptr r, fq_info_srcptr Fq)
{
    for (; r->deg >= 0 && fq_is_zero(r->coeffs[r->deg], Fq); --r->deg);
}

void fqpol_set_zero(fqpol_ptr r, MAYBE_UNUSED fq_info_srcptr Fq)
{
    r->deg = -1;
}

int fqpol_is_zero(fqpol_srcptr r, MAYBE_UNUSED fq_info_srcptr Fq)
{
    return (r->deg == -1);
}

void fqpol_set(fqpol_ptr r, fqpol_srcptr p, fq_info_srcptr Fq)
{
    fqpol_realloc_lazy(r, p->deg+1);
    for (int i = 0; i <= p->deg; ++i)
        fq_set(r->coeffs[i], p->coeffs[i], Fq);
    r->deg = p->deg;
}


// Random polynomial of degree at most deg.
void fqpol_set_random(fqpol_ptr r, unsigned int deg, fq_info_srcptr Fq)
{
    fqpol_realloc_lazy(r, deg+1);
    for (unsigned int i = 0; i <= deg; ++i)
        fq_set_random(r->coeffs[i], Fq);
    r->deg = deg;
    fqpol_update_degree(r, Fq);
}

void fqpol_set_ffspol(fqpol_ptr r, ffspol_srcptr f, fq_info_srcptr Fq)
{
    fqpol_realloc_lazy(r, f->deg+1);
    for (int i = 0; i <= f->deg; ++i)
        fq_set_mp(r->coeffs[i], f->coeffs[i], Fq);
    r->deg = f->deg;
    fqpol_update_degree(r, Fq);
}

void fqpol_set_ti(fqpol_ptr r, unsigned int n, fq_info_srcptr Fq)
{
    fqpol_realloc_lazy(r, n+1);
    for (unsigned int i = 0; i < n; ++i)
        fq_set_zero(r->coeffs[i], Fq);
    fq_set_ti(r->coeffs[n], 0, Fq);
    r->deg = n;
}

// This function is for debug.
void fqpol_out(FILE *file, fqpol_srcptr p, fq_info_srcptr Fq)
{
    fprintf(file, "[ ");
    if (p->deg == -1) {
        fprintf(file, "]\n");
        return;
    }
    for (int i = 0; i < p->deg; ++i) {
        fprintf(file, "0x");
        fq_out(file, p->coeffs[i], Fq);
        fprintf(file, ", ");
    }
    fprintf(file, "0x");
    sq_out(file, p->coeffs[p->deg]);
    fprintf(file, " ]\n");
}


void fqpol_add(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq)
{
    int k;
    int need_update = (p->deg == q->deg);

    // Use aliases to have deg(qq) >= deg(pp).
    fqpol_srcptr pp, qq;
    if (p->deg > q->deg) { pp = q; qq = p; }
    else                 { pp = p; qq = q; }

    fqpol_realloc_lazy(r, qq->deg+1);
    for (k = 0; k <= pp->deg; ++k)
        fq_add(r->coeffs[k], pp->coeffs[k], qq->coeffs[k], Fq);
    for (; k <= qq->deg; ++k)
        fq_set(r->coeffs[k], qq->coeffs[k], Fq);
    r->deg = qq->deg;
    if (need_update)
        fqpol_update_degree(r, Fq);
}

void fqpol_sub(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq)
{
    int opp, k;
    int need_update = (p->deg == q->deg);

    // Use aliases to have deg(qq) >= deg(pp).
    fqpol_srcptr pp, qq;
    if (p->deg > q->deg) { pp = q; qq = p; opp = 1; }
    else                 { pp = p; qq = q; opp = 0; }

    fqpol_realloc_lazy(r, qq->deg+1);
    for (k = 0; k <= pp->deg; ++k)
        fq_sub(r->coeffs[k], p->coeffs[k], q->coeffs[k], Fq);
    for (; k <= qq->deg; ++k) {
        if (opp) fq_set(r->coeffs[k], qq->coeffs[k], Fq);
        else     fq_opp(r->coeffs[k], qq->coeffs[k], Fq);
    }
    r->deg = qq->deg;
    if (need_update)
        fqpol_update_degree(r, Fq);
}

void fqpol_mul_ti(fqpol_ptr r, fqpol_srcptr p, unsigned i, fq_info_srcptr Fq)
{
    if (UNLIKELY(fqpol_is_zero(p, Fq))) {
        fqpol_set_zero(r, Fq);
        return;
    }

    fqpol_realloc_lazy(r, p->deg + i + 1);
    for (int k = p->deg; k >= 0; --k)
        fq_set(r->coeffs[k+i], p->coeffs[k], Fq);
    memset(r->coeffs, 0, i*sizeof(fq_t));
    r->deg = p->deg+i;
}

void fqpol_mul(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq)
{
    if (UNLIKELY(fqpol_is_zero(p, Fq) || fqpol_is_zero(q, Fq))) {
        fqpol_set_zero(r, Fq);
        return;
    }

    int kp = p->deg;
    int kq = q->deg;

    sq2_t * r_ur;
    r_ur = (sq2_t *) malloc ((kp+kq+1)*sizeof(sq2_t));
    ASSERT_ALWAYS(r_ur != NULL);
    memset(r_ur, 0, (kp+kq+1) * sizeof(sq2_t));

    for (int i = 0; i <= kp; ++i)
        for (int j = 0; j <= kq; ++j) {
            sq2_t c;
            fq_mulur(c, p->coeffs[i], q->coeffs[j], Fq);
            fq_addur(r_ur[i+j], r_ur[i+j], c, Fq);
        }

    fqpol_realloc_lazy(r, (kp+kq+1));
    for (int i = 0; i <= kp+kq; ++i)
        fq_reduce(r->coeffs[i], r_ur[i], Fq);
    free(r_ur);
    r->deg = p->deg + q->deg;
}

void fqpol_sqr(fqpol_ptr r, fqpol_srcptr p, fq_info_srcptr Fq)
{
    // TODO: do better, specially in characteristic 2.
    fqpol_mul(r, p, p, Fq);
}


void fqpol_smul(fqpol_ptr r, fqpol_srcptr p, fq_srcptr c, fq_info_srcptr Fq)
{
    fqpol_realloc_lazy(r, (p->deg+1));
    for (int i = 0; i <= p->deg; ++i)
        fq_mul(r->coeffs[i], p->coeffs[i], c, Fq);
    r->deg = p->deg;
}

void fqpol_sdiv(fqpol_ptr r, fqpol_srcptr p, fq_srcptr c, fq_info_srcptr Fq)
{
    if (fq_is_zero(c, Fq)) {
        fqpol_set_zero(r, Fq);
        return;
    }

    fq_t ic;
    fq_inv(ic, c, Fq);
    fqpol_realloc_lazy(r, (p->deg+1));
    for (int i = 0; i <= p->deg; ++i)
        fq_mul(r->coeffs[i], p->coeffs[i], ic, Fq);
    r->deg = p->deg;
}

void fqpol_get_coeff(fq_ptr r, fqpol_srcptr p, unsigned i, fq_info_srcptr Fq)
{
    fq_set(r, p->coeffs[i], Fq);
}

void fqpol_set_coeff(fqpol_ptr r, fq_srcptr x, unsigned i, fq_info_srcptr Fq)
{
    int ii = i; // to avoid signedness warning
    if (ii > r->deg && fq_is_zero(x, Fq)) return;
    if (ii > r->deg) {
        fqpol_realloc_lazy(r, ii+1);
        memset(&r->coeffs[r->deg+1], 0, (ii - r->deg)*sizeof(fq_t));
    }
    fq_set(r->coeffs[i], x, Fq);
    if (ii > r->deg)
        r->deg = ii;
    else if (fq_is_zero(x, Fq) && ii == r->deg)
        fqpol_update_degree(r, Fq);
}

int fqpol_deg(fqpol_srcptr r, MAYBE_UNUSED fq_info_srcptr Fq)
{
    return r->deg;
}

int fqpol_divrem(fqpol_ptr q, fqpol_ptr r,
        fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq)
{
    fq_t ilcb; /* Inverse of the leading coeff of b. */
    int  dega, degb;

    /* Trivial cases. */
    degb = fqpol_deg(b, Fq);
    if (UNLIKELY(degb == -1))
        return 0;
    if (UNLIKELY(degb == 0)) {
        fqpol_get_coeff(ilcb, b, 0, Fq);
        fqpol_sdiv(q, a, ilcb, Fq);
        fqpol_set_zero(r, Fq);
        return 1;
    }
    dega = fqpol_deg(a, Fq);
    if (dega < degb) {
        fqpol_set_zero(q, Fq);
        fqpol_set(r, a, Fq);
        return 1;
    }

    /* Check if input alias, and if so, work on a tmp space. */
    fqpol_t   tq, tr;
    fqpol_ptr qq, rr;
    if (q == a || q == b) {
        fqpol_init(tq);
        qq = &tq[0];
    } else
        qq = &q[0];
    if (r == a || r == b) {
        fqpol_init(tr);
        rr = &tr[0];
    } else
        rr = &r[0];

    /* General case. */
    fq_t lcr, cq;
    fqpol_t tmp;
    fqpol_init(tmp);

    /* Invariant: a = b*qq + rr. */
    fqpol_set(rr, a, Fq);
    fqpol_set_zero(qq, Fq);
    fqpol_get_coeff(ilcb, b, degb, Fq);
    fq_inv(ilcb, ilcb, Fq);
    for (int degr = dega; degr >= degb; ) {
        fqpol_get_coeff(lcr, rr, degr, Fq);
        fq_mul(cq, lcr, ilcb, Fq);
        fqpol_set_coeff(qq, cq, degr-degb, Fq);
        fqpol_smul(tmp, b, cq, Fq);
        fqpol_mul_ti(tmp, tmp, degr-degb, Fq);
        fqpol_sub(rr, rr, tmp, Fq);
        degr = fqpol_deg(rr, Fq);
    }
    fqpol_clear(tmp);

    /* If input were aliases, do the final copy. */
    if (q == a || q == b) {
        fqpol_set(q, qq, Fq);
        fqpol_clear(tq);
    }
    if (r == a || r == b) {
        fqpol_set(r, rr, Fq);
        fqpol_clear(tr);
    }
    return 1;
}

int fqpol_div(fqpol_ptr q, fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq)
{
    fqpol_t r;
    fqpol_init(r);
    int ret = fqpol_divrem(q, r, a, b, Fq);
    fqpol_clear(r);
    return ret;
}

int fqpol_rem(fqpol_ptr r, fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq)
{
    fqpol_t q;
    fqpol_init(q);
    int ret = fqpol_divrem(q, r, a, b, Fq);
    fqpol_clear(q);
    return ret;
}

void fqpol_gcd(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq)
{
    if (fqpol_is_zero(p, Fq)) {
        fqpol_set(r, q, Fq);
        return;
    }
    if (fqpol_is_zero(q, Fq)) {
        fqpol_set(r, p, Fq);
        return;
    }
    if ((fqpol_deg(p, Fq) == 0) || (fqpol_deg(q, Fq) == 0)) {
        fqpol_set_ti(r, 0, Fq);
        return;
    }

    // work with copies. Ensure deg(a) >= deg(b).
    fqpol_t a, b;
    fqpol_init(a);
    fqpol_init(b);
    if (p->deg >= q->deg) {
        fqpol_set(a, p, Fq);
        fqpol_set(b, q, Fq);
    } else {
        fqpol_set(a, q, Fq);
        fqpol_set(b, p, Fq);
    }
    fqpol_ptr aa, bb, tmp;
    aa = a; bb = b;

    while (!fqpol_is_zero(bb, Fq)) {
        fqpol_rem(aa, aa, bb, Fq);
        tmp = aa; aa = bb; bb = tmp;
    }
    fq_t lc;
    fqpol_get_coeff(lc, aa, aa->deg, Fq);
    fqpol_sdiv(r, aa, lc, Fq);
    fqpol_clear(a);
    fqpol_clear(b);
}

static void fqpol_powmod(fqpol_ptr r, fqpol_srcptr p,
        uint64_t power, fqpol_srcptr f, fq_info_srcptr Fq)
{
    ASSERT_ALWAYS(fqpol_deg(p, Fq) < fqpol_deg(f, Fq));
    if (power == 0) {
        fqpol_set_ti(r, 0, Fq);
        return;
    }
    if (power == 1) {
        fqpol_set(r, p, Fq);
        return;
    }
    // select msb:
    uint64_t mask = ((uint64_t)1)<<63;
    while ((power & mask) == 0)
        mask >>= 1;
    mask >>= 1;
    // horner:
    fqpol_t res;
    fqpol_init(res);
    fqpol_set(res, p, Fq);
    while (mask != 0) {
        fqpol_sqr(res, res, Fq);
        fqpol_rem(res, res, f, Fq);
        if (power & mask) {
            fqpol_mul(res, res, p, Fq);
            fqpol_rem(res, res, f, Fq);
        }
        mask >>= 1;
    }
    fqpol_set(r, res, Fq);
    fqpol_clear(res);
}

static void fqpol_powerXmod(fqpol_ptr r, uint64_t power,
        fqpol_srcptr f, fq_info_srcptr Fq)
{
    fqpol_t x;
    fqpol_init(x);
    fqpol_set_ti(x, 1, Fq);
    fqpol_powmod(r, x, power, f, Fq);
    fqpol_clear(x);
}

static void fqpol_splitlinear(fq_t * roots, fqpol_srcptr f, fq_info_srcptr Fq)
{
    if (f->deg == 1) {
        fqpol_get_coeff(roots[0], f, 0, Fq);
        fq_opp(roots[0], roots[0], Fq);
        return;
    }

    fqpol_t a, tra;
    fqpol_init(a);
    fqpol_init(tra);

#if (FP_SIZE == 2)
    do {
        fqpol_set_random(a, f->deg - 1, Fq);
        fqpol_set_zero(tra, Fq);
        for (unsigned int i = 0; i < Fq->degq-1; ++i) {
            fqpol_add(tra, tra, a, Fq);
            fqpol_sqr(a, a, Fq);
            fqpol_rem(a, a, f, Fq);
        }
        fqpol_add(tra, tra, a, Fq);

        fqpol_gcd(tra, tra, f, Fq);
    } while (fqpol_deg(tra, Fq) == 0 || fqpol_deg(tra, Fq) == fqpol_deg(f, Fq));
#else
#if ((FP_SIZE & 1) != 1)
#error "this code assume GF(2) or odd characteristic"
#endif
    uint64_t power = (Fq->order - 1) >> 1;
    fqpol_t one;
    fqpol_init(one);
    fqpol_set_ti(one, 0, Fq);
    do {
        fqpol_set_random(a, fqpol_deg(f, Fq)-1, Fq);
        fqpol_powmod(tra, a, power, f, Fq);
        fqpol_sub(tra, tra, one, Fq);
        fqpol_gcd(tra, tra, f, Fq);
    } while (fqpol_deg(tra, Fq) == 0 || fqpol_deg(tra, Fq) == fqpol_deg(f, Fq));
    fqpol_clear(one);
#endif

    if (fqpol_deg(tra, Fq) == 1) {
        fqpol_get_coeff(roots[0], tra, 0, Fq);
        fq_opp(roots[0], roots[0], Fq);
    } else {
        fqpol_splitlinear(roots, tra, Fq);
    }
    roots += tra->deg;

    fqpol_div(tra, f, tra, Fq);

    if (fqpol_deg(tra, Fq) == 1) {
        fqpol_get_coeff(roots[0], tra, 0, Fq);
        fq_opp(roots[0], roots[0], Fq);
    } else {
        fqpol_splitlinear(roots, tra, Fq);
    }
    fqpol_clear(a);
    fqpol_clear(tra);
}

// roots allocated by caller.
// return number of roots.
// NB: if there is a multiple root, it is returned only once.
int fqpol_roots(fq_t * roots, fqpol_srcptr f, fq_info_srcptr Fq)
{
    ASSERT_ALWAYS(!fqpol_is_zero(f, Fq));
    if (f->deg == 0)
        return 0;
    
    // work with a monic copy of f.
    fqpol_t F;
    fqpol_init(F);
    {
        fq_t lcf;
        fqpol_get_coeff(lcf, f, f->deg, Fq);
        fqpol_sdiv(F, f, lcf, Fq);
    }

    if (F->deg == 1) {
        fq_t c;
        fqpol_get_coeff(c, F, 0, Fq);
        fq_opp(roots[0], c, Fq);
        fqpol_clear(F);
        return 1;
    }

    // extract factor of f that contains only linear factors.
    fqpol_t xq, x;
    fqpol_init(xq);
    fqpol_init(x);
    int ret;
    fqpol_powerXmod(xq, Fq->order, F, Fq);
    fqpol_set_ti(x, 1, Fq);
    fqpol_sub(xq, xq, x, Fq);
    fqpol_gcd(xq, xq, F, Fq);
    if (xq->deg <= 0) {
        ret = 0;
    } else {
        ret = xq->deg;
        fqpol_splitlinear(roots, xq, Fq);
    }

    fqpol_clear(xq);
    fqpol_clear(x);
    fqpol_clear(F);
    return ret;
}
