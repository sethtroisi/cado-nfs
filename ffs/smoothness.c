#include "fppol.h"
#include "fppol_internal.h"
#include "fppol_facttools.h"


/*
 Iteration formula:
    g <- ( 2*g*t^(deg(g)+deg(f)) - f*g^2 ) div t^xxx
*/
void fppol_msb_preinverse(fppol_ptr invf, fppol_srcptr f, int k)
{
    ASSERT_ALWAYS(fppol_is_monic(f));
    int i = 1;
    fppol_t g, tmp;
    fppol_init(g);
    fppol_init(tmp);
    fppol_set_ti(g, 0);
    while (i <= k) {
        i *= 2;
#ifdef USE_F2
        fppol_qpow(tmp, g);
        fppol_mul(g, tmp, f);
        fppol_div_ti(g, g, fppol_deg(g)-i);
#else
        fppol_mul(tmp, g, g);
        fppol_mul(tmp, tmp, f);
        fppol_add(g, g, g);
        fppol_mul_ti(g, g, fppol_deg(g)+fppol_deg(f));
        fppol_sub(g, g, tmp);
        fppol_div_ti(g, g, fppol_deg(g)-i);
#endif
    }
    fppol_div_ti(invf, g, fppol_deg(g)-k);
    fppol_clear(g);
    fppol_clear(tmp);
}

void fppol_rem_precomp(fppol_ptr r, fppol_srcptr pq, 
        fppol_srcptr m, fppol_srcptr invm)
{
    fppol_t pqt;
    int d = fppol_deg(m);
    ASSERT_ALWAYS(fppol_deg(pq) <= fppol_deg(m) + fppol_deg(invm));
    fppol_init(pqt);
    fppol_div_ti(pqt, pq, d-2);
    fppol_mul(pqt, pqt, invm);
    fppol_div_ti(pqt, pqt, fppol_deg(invm)+2);
    fppol_mul(pqt, pqt, m);
    fppol_sub(r, pq, pqt);
    ASSERT_ALWAYS(fppol_deg(r) < fppol_deg(m));
    fppol_clear(pqt);
}


// The return value is 
//   - 0 if P is not smooth
//   - an upper bound on the degree of the largest factor otherwise
// NB: this algorithm is allowed to fail. If there is a factor of even
// multiplicity that has a degree > B, it is not detected.
int fppol_is_smooth(fppol_srcptr PP, int B)
{
    int B2 = (1+B)>>1;  // Ceiling(B/2)

    if (fppol_deg(PP) <= B)
        return fppol_deg(PP);

    fppol_t P;  // monic version of PP
    fppol_init(P);
    {
        fp_t lc;
        fppol_get_coeff(lc, PP, fppol_deg(PP));
        if (fp_is_one(lc))
            fppol_set(P, PP);
        else
            fppol_sdiv(P, PP, lc);
    }
    fppol_t dP;
    fppol_init(dP);
    fppol_derivative(dP, P);
    
    fppol_t tqi, acc, t, tmp;
    fppol_t invP2, invPq;
    fppol_init(tqi);
    fppol_init(t);
    fppol_init(tmp);
    fppol_init(acc);
    fppol_init(invP2);  // to reduce a product of 2
    fppol_init(invPq);  // to reduce the result of q-power.

    // TODO: the first one can be deduced from the second one.
    fppol_msb_preinverse(invP2, P, fppol_deg(P)-2);
#ifndef USE_F2
    fppol_msb_preinverse(invPq, P, 2*fppol_deg(P)-3);
#else
    fppol_set(invPq, invP2);
#endif

    fppol_set_ti(t, 1);

    int q = FP_SIZE;
    int qi = q;
    int i = 1;
    while (qi < fppol_deg(P)) {
        qi *= q;
        i++;
    }
    fppol_set_ti(tqi, qi);
    fppol_rem_precomp(tqi, tqi, P, invPq);

    while (i < B2) {
        fppol_qpow(tqi, tqi);
        fppol_rem_precomp(tqi, tqi, P, invPq);
        i++;
    }

    int smooth = 0;
    fppol_sub(acc, tqi, t);
    fppol_mul(acc, acc, dP);
    fppol_rem_precomp(acc, acc, P, invP2);
    while (i < B) {
        fppol_qpow(tqi, tqi);
        fppol_rem_precomp(tqi, tqi, P, invPq);
        i++;
        fppol_sub(tmp, tqi, t);
        fppol_mul(acc, acc, tmp);
        fppol_rem_precomp(acc, acc, P, invP2);
        if (fppol_is_zero(acc)) {
            smooth = 1;
            break;
        }
    }

    fppol_clear(P);
    fppol_clear(dP);
    fppol_clear(tqi);
    fppol_clear(t);
    fppol_clear(tmp);
    fppol_clear(acc);
    fppol_clear(invP2);
    fppol_clear(invPq);
    if (!smooth)
        return 0;
    else
        return i;
}

