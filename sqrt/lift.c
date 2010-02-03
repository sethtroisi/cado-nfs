#define _BSD_SOURCE
#define xxxWITH_BARRETT

#include <gmp.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "macros.h"
#include "powers_of_p.h"

double program_starttime;
double print_delay = 0.1;
#define WCT     (wct_seconds() - program_starttime)

struct prime_data {
    unsigned long p;
    // unsigned long * r;
    // unsigned long rj;
    void * powers;      // see .cpp file.

};/* }}} */

int degree;
mpz_t * f_hat;
mpz_t * f_hat_diff;

// above this threshold, we report each multiplication we do.
#define MUL_REPORT_THRESHOLD    1000000

#define REPORT_THIS(na, nb)     \
    (((na) > 1) && ((nb) > 1) && ((na) + (nb) > MUL_REPORT_THRESHOLD))

static void WRAP_mpz_mul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_mul(c,a,b);
    double t1 = seconds();
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f\n", na, nb, (double)na/nb, t1-t0);
    }
}

static void WRAP_mpz_invert(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_invert(c,a,b);
    double t1 = seconds();
    if (REPORT_THIS(na, nb)) {
        printf("inv %zu %zu (%.1f) %.1f\n", na, nb, (double)na/nb, t1-t0);
    }
}

static void WRAP_mpz_addmul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_addmul(c,a,b);
    double t1 = seconds();
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f\n", na, nb, (double)na/nb, t1-t0);
    }
}

static void WRAP_mpz_submul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    double t0 = seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_submul(c,a,b);
    double t1 = seconds();
    if (REPORT_THIS(na, nb)) {
        printf("mul %zu %zu (%.1f) %.1f\n", na, nb, (double)na/nb, t1-t0);
    }
}

static void WRAP_barrett_mod(mpz_ptr c, mpz_srcptr a, mpz_srcptr p, mpz_srcptr q)
{
    double t0 = seconds();
    mp_size_t na = mpz_size(a);
    mp_size_t np = mpz_size(p);
    barrett_mod(c,a,p,q);
    double t1 = seconds();
    if (REPORT_THIS(na, np)) {
        printf("mod %zu %zu (%.1f) %.1f\n", na, np, (double)na/np, t1-t0);
    }
}



// {{{ trivial utility
static const char * size_disp(size_t s, char buf[16])
{
    char * prefixes = "bkMGT";
    double ds = s;
    const char * px = prefixes;
    for( ; px[1] && ds > 500.0 ; ) {
        ds /= 1024.0;
        px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}
// }}}

static void mp_poly_eval_mod(mpz_ptr r, mpz_t * poly, int deg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
        WRAP_mpz_mul(r, r, a);
        // falls back on mpz_mod if qx == NULL
        WRAP_barrett_mod(r, r, q, qx);
        // mpz_mod(r, r, q);
        mpz_add(r, r, poly[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
}

static void mp_2poly_eval_mod(mpz_ptr r, mpz_ptr s, mpz_t * f, mpz_t * g, int degf, int degg, mpz_srcptr a, mpz_srcptr q, mpz_srcptr qx)
{
    int i;

    if (r == NULL) {
        mp_poly_eval_mod(s, g, degg, a, q, qx);
        return;
    }

    if (s == NULL) {
        mp_poly_eval_mod(r, f, degf, a, q, qx);
        return;
    }

    mpz_t w;
    mpz_init(w);
    mpz_set(r,f[0]);
    mpz_set(s,g[0]);
    mpz_set(w, a);
    WRAP_mpz_addmul(r, w, f[1]);
    WRAP_mpz_addmul(s, w, g[1]);
    for(i = 2 ; i <= degf && i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    for( ; i <= degf ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(r, w, f[i]);
    }
    for( ; i <= degg ; i++) {
        WRAP_mpz_mul(w, w, a);
        WRAP_barrett_mod(w, w, q, qx);
        WRAP_mpz_addmul(s, w, g[i]);
    }
    WRAP_barrett_mod(r, r, q, qx);
    WRAP_barrett_mod(s, s, q, qx);
    mpz_clear(w);
}

mpz_ptr my_barrett_init(mpz_srcptr px MAYBE_UNUSED)
{
#ifndef WITH_BARRETT
    return NULL;
#else
    mpz_ptr qx = malloc(sizeof(mpz_t));
    mpz_init(qx);
    barrett_init(qx, px);
    return qx;
#endif
}

void my_barrett_clear(mpz_ptr qx)
{
    if (!qx) return;
    mpz_clear(qx);
    free(qx);
}

#if 0
/* This implements the following iteration */
/*
p:=goodprimes[1];
r1:=GF(p)!goodprimes_lroots[1][1];
z1:=(1/Evaluate(Derivative(PolynomialRing(GF(p))!f),r1));
r:=GF(p)!r1;
z:=GF(p)!z1;
k:=1;
while k lt lift_prec do
k*:=2;
R:=Integers(p^k);
r:=R!Z!r;
z:=R!Z!z;
fr:=Evaluate(PolynomialRing(R)!f, r);
r:=r-fr*z;
fdr:=Evaluate(Derivative(PolynomialRing(R)!f), r);
z:=z-z*(fdr*z-1);
end while;
*/
void root_lift_innerlevels(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/*  */
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, rx, irx, lower);

    if (WCT > w0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n",
                WCT,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_ptr qk = my_barrett_init(pk);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(tb, f_hat, degree, rx, pk, qk);
    WRAP_mpz_mul(ta, tb, irx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(ta), mpz_size(tb), mpz_size(irx));
    mpz_sub(rx, rx, ta);
    WRAP_barrett_mod(rx, rx, pk, qk);

    mp_poly_eval_mod(tb, f_hat_diff, degree-1, rx, pk, qk);
    WRAP_mpz_mul(ta, irx, tb);
    mpz_sub_ui(ta, ta, 1);
    WRAP_barrett_mod(ta, ta, pk, qk);
    WRAP_mpz_mul(tb, ta, irx);
    mpz_sub(ta, irx, tb);
    WRAP_barrett_mod(irx, ta, pk, qk);

    mpz_clear(ta);
    mpz_clear(tb);
    /*
       if (precision < 40) {
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", WCT, rx);
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", WCT, p->invdev_rx);
       }
       */

    my_barrett_clear(qk);
}
void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, rx, irx, lower);

    if (WCT > w0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n",
                WCT,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_ptr qk = my_barrett_init(pk);

    mpz_t ta;
    mpz_init(ta);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(ta, f_hat, degree, rx, pk, qk);
    WRAP_mpz_mul(ta, ta, irx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(ta), mpz_size(tb), mpz_size(p->invdev_rx));
    mpz_sub(rx, rx, ta);
    WRAP_barrett_mod(rx, rx, pk, qk);
    // don't do the lower part (thus irx is not relevant).

    mpz_clear(ta);

    my_barrett_clear(qk);
}/*  */
#else
void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/*  */
{
    double w0 = WCT;
    assert(precision > 0);

    if (precision == 1) {
        // mpz_srcptr pk = power_lookup_const(p->powers, precision);
        // we want f'(r) mod p !
        // mp_2poly_eval_mod(NULL, fprime, f_hat, f_hat_diff, degree, degree-1, rx, pk, NULL);
        // gmp_fprintf(stderr, "f'(r_1) mod p^1 == %Zd\n", fprime);

        return;
    }
    int lower = precision - precision / 2;

    // recurse.
    root_lift(p, rx, irx, lower);
    // gmp_fprintf(stderr, "f'(r_%d) mod p^%d == %Zd\n", lower, lower, fprime_here);

    if (WCT > w0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n",
                WCT,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    mpz_srcptr pl = power_lookup_const(p->powers, lower);
    mpz_ptr qk = my_barrett_init(pk);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    mpz_t fprime;
    mpz_init(fprime);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_2poly_eval_mod(tb, fprime, f_hat, f_hat_diff, degree, degree-1, rx, pk, qk);
    // gmp_fprintf(stderr, "f(r_%d) mod p^%d == %Zd\n", lower, precision, tb);
    // if (fprime) gmp_fprintf(stderr, "f'(r_%d) mod p^%d == %Zd\n", precision, precision, fprime);
    // WRAP_mpz_invert(ta, fprime, pl);
    /* use irx. only one iteration of newton.
     */
    WRAP_barrett_mod(fprime, fprime, pl, NULL); /* FIXME */
    WRAP_mpz_mul(ta, irx, fprime);
    WRAP_barrett_mod(ta, ta, pk, qk);
    mpz_sub_ui(ta,ta,1);
    WRAP_mpz_submul(irx, irx, ta);
    WRAP_barrett_mod(irx, irx, pk, qk);

    mpz_clear(fprime);

    WRAP_mpz_mul(ta, irx, tb);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(ta), mpz_size(tb), mpz_size(irx));
    mpz_sub(rx, rx, ta);
    WRAP_barrett_mod(rx, rx, pk, qk);

    // gmp_fprintf(stderr, "r_%d = %Zd\n", precision, rx);

    mpz_clear(ta);
    mpz_clear(tb);
    /*
       if (precision < 40) {
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", WCT, rx);
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", WCT, p->invdev_rx);
       }
       */

    my_barrett_clear(qk);
}

#endif

int main(int argc, char * argv[])
{
    program_starttime = wct_seconds();

    unsigned long p = 9223372036854799379UL;
    unsigned long r = 8197937682795196680UL;
    unsigned int precision = 1908226;
    const char * coeff_strings[] = {
        "33702859150680615562179587684939483",
        "2824683557949990996869965408188",
        "48702062013735864254260849",
        "-969379411138214076748",
        "-3170645766299012",
        "9108387600",
    };
    degree = sizeof(coeff_strings) / sizeof(coeff_strings[0]) - 1;

    if (argc >= 3 && strcmp(argv[1], "-prec") == 0) {
        precision = atoi(argv[2]);
    }

    fprintf(stderr, "GMP header: %d.%d.%d, using library %s\n",
            __GNU_MP_VERSION,
            __GNU_MP_VERSION_MINOR,
            __GNU_MP_VERSION_PATCHLEVEL,
            gmp_version);
    fprintf(stderr, "GMP compiled with %s, flags %s\n",
            __GMP_CC,
            __GMP_CFLAGS);
#ifdef  WITH_BARRETT
    fprintf(stderr, "This run uses barrett reduction\n");
#else
    fprintf(stderr, "This run does not use barrett reduction\n");
#endif

    struct prime_data prime[1];
    prime->p = p;
    prime->powers = power_lookup_table_init(prime->p);

    f_hat = malloc((degree+1) * sizeof(mpz_t));
    f_hat_diff = malloc((degree) * sizeof(mpz_t));

    for(int i = 0 ; i <= degree ; i++) {
        mpz_init_set_str(f_hat[i], coeff_strings[i], 0);
    }
    {
        mpz_t w;
        mpz_init(w);
        mpz_set_ui(w,1);
        for(int i = degree-2 ; i >= 0 ; i--) {
            mpz_mul(w,w,f_hat[degree]);
            mpz_mul(f_hat[i],f_hat[i],w);
        }
        mpz_set_ui(f_hat[degree],1);
        mpz_clear(w);
    }

    for(int i = 0 ; i <= degree-1 ; i++) {
        mpz_init(f_hat_diff[i]);
        mpz_mul_ui(f_hat_diff[i], f_hat[i+1], i+1);
    }

    {
        char sbuf[32];
        fprintf(stderr, "# [%2.2lf] Lifting to precision l=%d (p^l is approx %s)\n", WCT, precision, size_disp(precision * log(prime->p)/M_LN2 / 8, sbuf));
    }

    fprintf(stderr, "# [%2.2lf] Computing powers of p\n", WCT);
    power_lookup(prime->powers, precision);
    fprintf(stderr, "# [%2.2lf] Computing powers of p: done.\n", WCT);

    mpz_t rx, irx;
    mpz_init(rx);
    mpz_init(irx);

    mpz_set_ui(rx, r);

    mpz_srcptr p1 = power_lookup_const(prime->powers, 1);

    mp_poly_eval_mod(irx, f_hat_diff, degree-1, rx, p1, NULL);
    mpz_invert(irx, irx, p1);

    fprintf(stderr, "# [%2.2lf] Lifting (%lu,x-%lu)\n", WCT, p, r);

    {
        double t0, t1;
        double w0, w1;
        double rate;

        t0 = seconds();
        w0 = WCT;

        root_lift(prime, rx, irx, precision);

        t1 = seconds();
        w1 = WCT;
        rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;
        fprintf(stderr, "# [%2.2lf] lift completed. %.2lf, wct %.2lf, %.1f%%\n",
                WCT, t1-t0, w1-w0, rate);
        fprintf(stderr, "# [%2.2lf] limb0 of lifted root is %lu\n",
                WCT, mpz_getlimbn(rx, 0));
    }

    mpz_clear(irx);
    mpz_clear(rx);

    for(int i = 0 ; i <= degree ; i++) {
        mpz_clear(f_hat[i]);
    }
    for(int i = 0 ; i <= degree-1 ; i++) {
        mpz_clear(f_hat_diff[i]);
    }

    free(f_hat);
    free(f_hat_diff);
}

