#include "types.h"
#include "norm.h"
#include "qlat.h"

// The skewness is a (difference of) degree -> an unsigned int.
// The function returns 1 if it succeeded, i.e. the result fits in the
// bounds for the a_i.
int skewGauss(qlat_t qlat, unsigned int skewness)
{
    if (!qlat->want_longq) {
        sq_t a[2], b[2];

        sq_set(a[0], qlat->q);
        sq_set_zero(b[0]);
        sq_set(a[1], qlat->rho);
        sq_set_ti(b[1], skewness);  // this is one if there is no skewness.

        do {
            sq_t qq;

            sq_div(qq, a[0], a[1]);
            if ((sq_deg(qq) + sq_deg(b[1])) > sq_deg(a[0]))
                break;
            sq_submul(a[0], a[1], qq);
            sq_submul(b[0], b[1], qq);

            sq_div(qq, a[1], a[0]);
            if ((sq_deg(qq) + sq_deg(b[0])) > sq_deg(a[1]))
                break;
            sq_submul(a[1], a[0], qq);
            sq_submul(b[1], b[0], qq);
        } while (sq_deg(a[0]) > sq_deg(b[0]));

        // Compensate for the skewness
        sq_div_ti(b[0], b[0], skewness);
        sq_div_ti(b[1], b[1], skewness);

        // cast and check that the result fits!
        return (ai_set_sq(qlat->a0, a[0]) && ai_set_sq(qlat->a1, a[1])
                && ai_set_sq(qlat->b0, b[0]) && ai_set_sq(qlat->b1, b[1]));
    } else {
        fppol_t a[2], b[2];
        fppol_init(a[0]);
        fppol_init(a[1]);
        fppol_init(b[0]);
        fppol_init(b[1]);

        fppol_set(a[0], qlat->longq);
        fppol_set_zero(b[0]);
        fppol_set(a[1], qlat->longrho);
        fppol_set_ti(b[1], skewness);  // this is one if there is no skewness.

        fppol_t qq;
        fppol_init(qq);
        do {

            fppol_div(qq, a[0], a[1]);
            if ((fppol_deg(qq) + fppol_deg(b[1])) > fppol_deg(a[0]))
                break;
            fppol_submul(a[0], a[1], qq);
            fppol_submul(b[0], b[1], qq);

            fppol_div(qq, a[1], a[0]);
            if ((fppol_deg(qq) + fppol_deg(b[0])) > fppol_deg(a[1]))
                break;
            fppol_submul(a[1], a[0], qq);
            fppol_submul(b[1], b[0], qq);
        } while (fppol_deg(a[0]) > fppol_deg(b[0]));
        fppol_clear(qq);

        // Compensate for the skewness
        fppol_div_ti(b[0], b[0], skewness);
        fppol_div_ti(b[1], b[1], skewness);

        fppol_set(qlat->longa0, a[0]);
        fppol_set(qlat->longa1, a[1]);
        fppol_set(qlat->longb0, b[0]);
        fppol_set(qlat->longb1, b[1]);

        fppol_clear(a[0]);
        fppol_clear(a[1]);
        fppol_clear(b[0]);
        fppol_clear(b[1]);
        return 1;
    }
}

void print_qlat_info(qlat_t qlat)
{
    if (!qlat->want_longq) {
        printf("# q-lattice info:\n");
        printf("#   q = ");   sq_out(stdout, qlat->q); 
        printf(" ; rho = "); sq_out(stdout, qlat->rho); printf("\n");
        printf("#   a0 = ");  ai_out(stdout, qlat->a0);
        printf(" ; a1 = ");  ai_out(stdout, qlat->a1);
        printf(" ; b0 = ");  ai_out(stdout, qlat->b0);
        printf(" ; b1 = ");  ai_out(stdout, qlat->b1); printf("\n");
    } else {
        printf("# q-lattice info:\n");
        printf("#   q = ");  fppol_out(stdout, qlat->longq); 
        printf(" ; rho = "); fppol_out(stdout, qlat->longrho); printf("\n");
        printf("#   a0 = "); fppol_out(stdout, qlat->longa0);
        printf(" ; a1 = ");  fppol_out(stdout, qlat->longa1);
        printf(" ; b0 = ");  fppol_out(stdout, qlat->longb0);
        printf(" ; b1 = ");  fppol_out(stdout, qlat->longb1); printf("\n");
    }
}

int is_valid_sq(qlat_t qlat, ffspol_srcptr F)
{
    // F(rho) = Norm_F(rho, 1)
    fppol_t q, rho, one, norm;
    fppol_init(q);
    fppol_init(rho);
    fppol_init(one);
    fppol_init(norm);
    if (!qlat->want_longq) {
        fppol_set_sq(rho, qlat->rho);
        fppol_set_sq(q, qlat->q);
    } else {
        fppol_set(rho, qlat->longrho);
        fppol_set(q,   qlat->longq);
    }
    fppol_set_ti(one, 0);
    ffspol_norm(norm, F, rho, one);
    fppol_rem(norm, norm, q);
    int ret = fppol_is_zero(norm);
    fppol_clear(norm);
    fppol_clear(q);
    fppol_clear(rho);
    fppol_clear(one);
    return ret;
}



// a = i*a0 + j*a1
// b = i*b0 + j*b1
// TODO: for characteristic 3, maybe normalize (make monic) a or b.
void ij2ab(fppol_t a, fppol_t b, ij_t i, ij_t j, qlat_t qlat)
{
    fppol_t tmp;
    fppol_init(tmp);

    if (!qlat->want_longq) {
        fppol_mul_aixij(tmp, qlat->a0, i);
        fppol_mul_aixij(a, qlat->a1, j);
    } else {
        fppol_mul_ij(tmp, qlat->longa0, i);
        fppol_mul_ij(a, qlat->longa1, j);
    }
    fppol_add(a, a, tmp);

    if (!qlat->want_longq) {
        fppol_mul_aixij(tmp, qlat->b0, i);
        fppol_mul_aixij(b, qlat->b1, j);
    } else {
        fppol_mul_ij(tmp, qlat->longb0, i);
        fppol_mul_ij(b, qlat->longb1, j);
    }
    fppol_add(b, b, tmp);

    fppol_clear(tmp);
}

// i = (a*b1 - b*a1) / q
// j = (-a*b0 + b*a0) / q
// TODO: for characteristic 3, maybe normalize (make monic) i or j.
// The input must be an (a,b) pair in the q-lattice, so that the
// divisions by q are exact (q is the determinant of the base-change
// matrix). If this is not the case, abort.
// In principle, this function will always be called with inputs such
// that i and j fits within their types. So we don't return an error code
// but also abort if this is not the case.
void ab2ij(ij_t i, ij_t j, fppol_t a, fppol_t b, qlat_t qlat)
{
    fppol_t tmp, tmp2, tmpq;
    int fit;
    fppol_init(tmp);
    fppol_init(tmp2);
    fppol_init(tmpq);
    if (!qlat->want_longq)
        fppol_set_sq(tmpq, qlat->q);
    else
        fppol_set(tmpq, qlat->longq);

    if (!qlat->want_longq) {
        fppol_mul_ai(tmp, a, qlat->b1);
        fppol_mul_ai(tmp2, b, qlat->a1);
    } else {
        fppol_mul(tmp, a, qlat->longb1);
        fppol_mul(tmp2, b, qlat->longa1);
    }
    fppol_sub(tmp, tmp, tmp2);
    fppol_divrem(tmp, tmp2, tmp, tmpq);
    ASSERT_ALWAYS(fppol_is_zero(tmp2));
    fit = ij_set_mp(i, tmp);
    ASSERT_ALWAYS(fit);

    if (!qlat->want_longq) {
        fppol_mul_ai(tmp, a, qlat->b0);
        fppol_mul_ai(tmp2, b, qlat->a0);
    } else {
        fppol_mul(tmp, a, qlat->longb0);
        fppol_mul(tmp2, b, qlat->longa0);
    }
    fppol_sub(tmp, tmp2, tmp);
    fppol_divrem(tmp, tmp2, tmp, tmpq);
    ASSERT_ALWAYS(fppol_is_zero(tmp2));
    fit = ij_set_mp(j, tmp);
    ASSERT_ALWAYS(fit);

    fppol_clear(tmp);
    fppol_clear(tmp2);
    fppol_clear(tmpq);
}


// Compute lambda for an element of the factor base.
// If the result is projective, the set lambda to p.
void compute_lambda(fbprime_ptr lambda,
    fbprime_srcptr p, fbprime_srcptr r, qlat_srcptr qlat)
{
    fbprime_t t0, t1;
    fbprime_t a0, a1, b0, b1;

    if (!qlat->want_longq) {
        fbprime_set_ai(a0, qlat->a0); 
        fbprime_set_ai(a1, qlat->a1);
        fbprime_set_ai(b0, qlat->b0);
        fbprime_set_ai(b1, qlat->b1);
    } else {
        fppol_t aux, pp;
        fppol_init(aux);
        fppol_init(pp);
        fppol_set_fbprime(pp, p);
        fppol_rem(aux, qlat->longa0, pp); fbprime_set_mp(a0, aux);
        fppol_rem(aux, qlat->longa1, pp); fbprime_set_mp(a1, aux);
        fppol_rem(aux, qlat->longb0, pp); fbprime_set_mp(b0, aux);
        fppol_rem(aux, qlat->longb1, pp); fbprime_set_mp(b1, aux);
        fppol_clear(aux);
        fppol_clear(pp);
    }

    int was_proj = fbprime_deg(r) == fbprime_deg(p);

    if (was_proj) {
        fbprime_sub(t0, r, p);
        fbprime_mulmod(t0, a0, t0, p);
        fbprime_sub   (t0, t0, b0);
        fbprime_rem   (t0, t0, p);
    } else {
        fbprime_mulmod(t0, b0, r, p);
        fbprime_sub   (t0, a0, t0);
        fbprime_rem   (t0, t0, p);
    }
    if (fbprime_is_zero(t0)) {
      fbprime_set(lambda, p);
      return;
    }
    int ret = fbprime_invmod(t0, t0, p);
    if (!ret) { // this can happen for powers
        // FIXME: is this condition correct ?
        fbprime_set(lambda, p);
        return;
    }
    if (was_proj) {
        fbprime_mulmod(t1, a1, r, p);
        fbprime_sub   (t1, b1, t1);
        fbprime_rem   (t1, t1, p);
    } else {
        fbprime_mulmod(t1, b1, r, p);
        fbprime_sub   (t1, t1, a1);
        fbprime_rem   (t1, t1, p);
    }
    fbprime_mulmod(lambda, t0, t1, p);
}


