#include "mpz_poly.h"

static void mp_poly_evalz(mpz_t r, mpz_t * poly, int deg, mpz_srcptr a)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
	mpz_mul(r, r, a);
	mpz_add(r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
static void
mp_poly_eval_diffz(mpz_t r, mpz_t * poly, int deg, mpz_srcptr a)
{
    int i;

    mpz_mul_ui(r, poly[deg], (unsigned long) deg);
    for (i = deg - 1; i >= 1; i--) {
	mpz_mul(r, r, a);
	mpz_addmul_ui(r, poly[i], (unsigned long) i);
    }
}

static void mp_poly_eval(mpz_t r, mpz_t * poly, int deg, long a)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
	mpz_mul_si(r, r, a);
	mpz_add(r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
static void mp_poly_eval_diff(mpz_t r, mpz_t * poly, int deg, long a)
{
    int i;

    mpz_mul_ui(r, poly[deg], (unsigned long) deg);
    for (i = deg - 1; i >= 1; i--) {
	mpz_mul_si(r, r, a);
	mpz_addmul_ui(r, poly[i], (unsigned long) i);
    }
}

// assuming that pk is odd, that r is a root mod p^l, and that pk is a
// power of p equal to at most p^(2l), compute a lift of r mod pk.
int lift_root(mpz_t * f, int d, const unsigned long pk, unsigned long *r)
{
    unsigned long res;
    mpz_t aux, aux2, mp_p;

    mpz_init(aux);
    mpz_init(aux2);
    mpz_init_set_ui(mp_p, pk);
    mp_poly_eval_diff(aux, f, d, *r);
    mpz_mod_ui(aux, aux, pk);
    if (!mpz_invert(aux, aux, mp_p)) {
	return 0;
    }
    mp_poly_eval(aux2, f, d, *r);
    mpz_mod(aux2, aux2, mp_p);
    mpz_mul(aux2, aux2, aux);
    mpz_neg(aux2, aux2);
    mpz_add_ui(aux2, aux2, *r);
    mpz_mod(aux2, aux2, mp_p);

    res = mpz_get_ui(aux2);
    mpz_clear(aux);
    mpz_clear(aux2);
    mpz_clear(mp_p);

    *r = res;
    return 1;
}

int lift_rootz(mpz_t * f, int d, mpz_t pk, mpz_t r)
{
    mpz_t aux, aux2;

    mpz_init(aux);
    mpz_init(aux2);
    mp_poly_eval_diffz(aux, f, d, r);
    mpz_mod(aux, aux, pk);
    if (!mpz_invert(aux, aux, pk)) {
	return 0;
    }
    mp_poly_evalz(aux2, f, d, r);
    mpz_mod(aux2, aux2, pk);
    mpz_mul(aux2, aux2, aux);
    mpz_sub(r, r, aux2);
    mpz_mod(r, r, pk);
    mpz_clear(aux);
    mpz_clear(aux2);
    return 1;
}
