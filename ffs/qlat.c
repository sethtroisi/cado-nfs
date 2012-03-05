#include "types.h"
#include "qlat.h"

// The skewness is a (difference of) degree -> an unsigned int.
// The function returns 1 if it succeeded, i.e. the result fits in the
// bounds for the a_i.
int skewGauss(qlat_t qlat, unsigned int skewness)
{
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
}

void print_qlat_info(qlat_t qlat)
{
    printf("# q-lattice info:\n");
    printf("#   q = ");   sq_out(stdout, qlat->q); printf("\n");
    printf("#   rho = "); sq_out(stdout, qlat->rho); printf("\n");
    printf("#   a0 = ");  ai_out(stdout, qlat->a0); printf("\n");
    printf("#   a1 = ");  ai_out(stdout, qlat->a1); printf("\n");
    printf("#   b0 = ");  ai_out(stdout, qlat->b0); printf("\n");
    printf("#   b1 = ");  ai_out(stdout, qlat->b1); printf("\n");
    printf("#   side = %d\n", qlat->side);
}


// a = i*a0 + j*a1
// b = i*b0 + j*b1
// TODO: for characteristic 3, maybe normalize (make monic) a or b.
void ij2ab(fppol_t a, fppol_t b, ij_t i, ij_t j, qlat_t qlat)
{
    fppol_t tmp;
    fppol_init(tmp);

    fppol_mul_aixij(tmp, qlat->a0, i);
    fppol_mul_aixij(a, qlat->a1, j);
    fppol_add(a, a, tmp);

    fppol_mul_aixij(tmp, qlat->b0, i);
    fppol_mul_aixij(b, qlat->b1, j);
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
    fppol_set_sq(tmpq, qlat->q);

    fppol_mul_ai(tmp, a, qlat->b1);
    fppol_mul_ai(tmp2, b, qlat->a1);
    fppol_sub(tmp, tmp, tmp2);
    fppol_divrem(tmp, tmp2, tmp, tmpq);
    ASSERT_ALWAYS(fppol_is_zero(tmp2));
    fit = ij_set_mp(i, tmp);
    ASSERT_ALWAYS(fit);

    fppol_mul_ai(tmp, a, qlat->b0);
    fppol_mul_ai(tmp2, b, qlat->a0);
    fppol_sub(tmp, tmp2, tmp);
    fppol_divrem(tmp, tmp2, tmp, tmpq);
    ASSERT_ALWAYS(fppol_is_zero(tmp2));
    fit = ij_set_mp(j, tmp);
    ASSERT_ALWAYS(fit);

    fppol_clear(tmp);
    fppol_clear(tmp2);
    fppol_clear(tmpq);
}
