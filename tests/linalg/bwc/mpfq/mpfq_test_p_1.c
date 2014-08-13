#include "cado.h"
#define _GNU_SOURCE     /* asprintf */
#define _BSD_SOURCE     /* asprintf sometimes (I think) */
/* test_p_1.c is sed- generated from test.c.meta */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <unistd.h>
#include <gmp.h>
#include <time.h>

#include "mpfq_p_1.h"
#include "mpfq_name_K.h"

Kfield K;

#define DO_ONE_TEST(name, CODE)						\
	do {								\
	for(int i = 0 ; i < ntests ; i++) {				\
		do { CODE } while (0);					\
		if (Kcmp(r1, r2) == 0)					\
			continue;					\
		fprintf(stderr, "Test failed [" name "]/%d\n", i);	\
		fprintf(stderr, "Seed is %lu, nb_tests is %d\n", seed, ntests);\
		abort();						\
	}								\
	if (!quiet) fprintf(stderr, "ok - [" name "], %d times\n", ntests);\
        } while (0)

#define DO_ONE_TEST_VEC(name, CODE)					\
	do {								\
	for(int i = 0 ; i < ntests ; i++) {				\
		do { CODE } while (0);					\
		if (Kvec_cmp(v1, v2, length) == 0)				\
			continue;					\
		fprintf(stderr, "Test failed [" name "]/%d\n", i);	\
		fprintf(stderr, "Seed is %lu, nb_tests is %d\n", seed, ntests);\
		abort();						\
	}								\
	if (!quiet) fprintf(stderr, "ok - [" name "], %d times\n", ntests);\
        } while (0)
#define DO_ONE_TEST_POLY(name, CODE)					\
	do {								\
	for(int i = 0 ; i < ntests ; i++) {				\
		do { CODE } while (0);					\
		if (Kpoly_cmp(p1, p2) == 0)				\
			continue;					\
		fprintf(stderr, "Test failed [" name "]/%d\n", i);	\
		fprintf(stderr, "Seed is %lu, nb_tests is %d\n", seed, ntests);\
		abort();						\
	}								\
	if (!quiet) fprintf(stderr, "ok - [" name "], %d times\n", ntests);\
        } while (0)


void usage() {
    fprintf(stderr, "usage: ./test [-q] [-N <nb_loops>] [-n <nb_tests>] [-s <seed>]\n");
    fprintf(stderr, "  -N 0 yields an infinite loop\n");
    fprintf(stderr, "  -q means quiet\n");
    exit(1);
}


#ifndef FIX_PRIME
void get_random_prime(mpz_t z, mp_bitcnt_t n, int quiet, gmp_randstate_t rnd) {
    do {
        mpz_urandomb(z, rnd, n);
        if (mpz_sizeinbase(z, 2) != n) continue;
    } while (!mpz_probab_prime_p(z, 5));
    if (!quiet) gmp_fprintf(stderr, "Using prime p = %Zd\n", z);
}
#endif





int main(int argc, char * argv[])
{
    int ntests = 100;
    int nloops = 1;
    int quiet = 0;
    Kelt a0, a1, a2, a3, a4, a5;
    Kelt  r1, r2;
    Kvec  v1, v2;
    Kvec  w1, w2, w3, w4;
    Kpoly  p1, p2;
    Kpoly  q1, q2, q3, q4;
    unsigned long seed = (unsigned long)time(NULL);
    const char * prime_str = NULL;

    while (argc > 1 && argv[1][0] == '-') {
        if (argc > 2 && strcmp(argv[1], "-s") == 0) {
            seed = atol(argv[2]);
            argc -= 2;
            argv += 2;
        } else if (argc > 2 && strcmp(argv[1], "-N") == 0) {
            nloops = atol(argv[2]);
            argc -= 2;
            argv += 2;
        } else if (argc > 2 && strcmp(argv[1], "-n") == 0) {
            ntests = atol(argv[2]);
            argc -= 2;
            argv += 2;
#ifndef FIX_PRIME
        } else if (argc > 2 && strcmp(argv[1], "-p") == 0) {
            prime_str = argv[2];
            argc -= 2;
            argv += 2;
#endif
        } else if (argc > 1 && strcmp(argv[1], "-q") == 0) {
            quiet = 1;
            argc--;
            argv++;
        } else
            usage();
    }
    if (argc > 1)
        usage();


    if (!quiet) fprintf(stderr, "--- testing for p_1\n");

    int i = 0;
    while ( (nloops == 0) || (i < nloops) ) {
        gmp_randstate_t  rstate;

        Kfield_init();

        if (!quiet)
            fprintf(stderr, "seeding random generator with %lu\n", seed);
        gmp_randinit_mt(rstate);
        gmp_randseed_ui(rstate, seed);

#ifndef CHAR2
#ifndef FIX_PRIME
#ifndef EXTENSION_OF_GFP
        mpz_t p;
        mpz_init(p);
#ifdef  VARIABLE_SIZE_PRIME
        mp_bitcnt_t size_prime= 12*GMP_NUMB_BITS;
#endif
#ifndef VARIABLE_SIZE_PRIME
        mp_bitcnt_t size_prime=Kimpl_max_characteristic_bits();
#endif
        if (prime_str) {
            mpz_set_str(p, prime_str, 0);
        } else {
            get_random_prime(p, size_prime, quiet, rstate);
        }
        Kfield_specify(MPFQ_PRIME_MPZ, p);
        mpz_clear(p);
#endif
#endif
#endif

#ifdef EXTENSION_OF_GFP
        mpz_t p;
        mpz_init(p);
#ifdef  VARIABLE_SIZE_PRIME
        mp_bitcnt_t size_prime= 12*GMP_NUMB_BITS;
#else
        mp_bitcnt_t size_prime=Kimpl_max_characteristic_bits();
#endif
        get_random_prime(p, size_prime, quiet, rstate);
        Kfield_specify(MPFQ_PRIME_MPZ, p);
        mpz_clear(p);

        int extdeg = 5;
        MPFQ_CREATE_FUNCTION_NAME(BFIELD, poly) defpol;
        MPFQ_CREATE_FUNCTION_NAME(BFIELD, poly_init) (K->kbase, defpol, 0);
        MPFQ_CREATE_FUNCTION_NAME(BFIELD, poly_random) (K->kbase, defpol, extdeg, rstate);
        MPFQ_CREATE_FUNCTION_NAME(BFIELD, poly_setcoeff_ui) (K->kbase, defpol, 1, extdeg);
        MPFQ_CREATE_FUNCTION_NAME(BFIELD, poly_setcoeff_ui) (K->kbase, defpol, 1, 0);
        Kfield_specify(MPFQ_POLYNOMIAL, defpol);
#endif


        Kinit(&a0);
        Kinit(&a1);
        Kinit(&a2);
        Kinit(&a3);
        Kinit(&a4);
        Kinit(&a5);
        Kinit (&r1);
        Kinit (&r2);

        /*-----------------------------------------------------------*/
        /*          Common tests                                     */
        /*-----------------------------------------------------------*/

        DO_ONE_TEST("add commutativity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kadd (r1, a0, a1);
                Kadd (r2, a1, a0);
                });

        DO_ONE_TEST("sub = add o neg", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Ksub (r1, a0, a1);
                Kneg(r2, a1);
                Kadd(r2, r2, a0);
                });

        DO_ONE_TEST("add o sub = id", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Ksub (r1, a0, a1);
                Kadd (r1, r1, a1);
                Kset (r2, a0);
                });

        DO_ONE_TEST("get_ui o sub_ui(x) o set_ui(x+1) == 1", {
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                Kset_ui(a1, x+1);
                Ksub_ui(a2, a1, x);
                unsigned long r = Kget_ui(a2);
                Kset_ui(r1, r);
                Kset_ui(r2, 1);
                });

        DO_ONE_TEST("sub_ui(y) o set_ui(x) == set_mpz(x-y)", {
                mpz_t z;
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                unsigned long y = gmp_urandomb_ui(rstate, 32);
                mpz_init_set_ui(z, x);
                mpz_sub_ui(z, z, y);
                Kset_ui(a1, x);
                Ksub_ui(r1, a1, y);
                Kset_mpz(r2, z);
                mpz_clear(z);
                });

        DO_ONE_TEST("add_ui(y) o neg o set_ui(x) == set_mpz(y-x)", {
                mpz_t z;
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                unsigned long y = gmp_urandomb_ui(rstate, 32);
                mpz_init_set_ui(z, y);
                mpz_sub_ui(z, z, x);
                Kset_ui(a1, x);
                Kneg(a1, a1);
                Kadd_ui(r1, a1, y);
                Kset_mpz(r2, z);
                mpz_clear(z);
                });

        DO_ONE_TEST("add_ui o sub_ui = id", {
                Krandom2 (a0, rstate);
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                Ksub_ui (r1, a0, x);
                Kadd_ui (r1, r1, x);
                Kset (r2, a0);
                });

        DO_ONE_TEST("add o neg = sub", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Ksub (r1, a0, a1);
                Kneg (r2, a1);
                Kadd (r2, r2, a0);
                });

        DO_ONE_TEST("mul commutativity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kmul (r1, a0, a1);
                Kmul (r2, a1, a0);
                });

        DO_ONE_TEST("sqr(x) = mul(x,x)", {
                Krandom2 (a0, rstate);
                Kmul (r1, a0, a0);
                Ksqr (r2, a0);
                });

        DO_ONE_TEST("mul distributivity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Krandom2 (a2, rstate);
                Kadd (a3, a1, a2);
                Kmul (r1, a0, a3);
                Kmul (a4, a0, a1);
                Kmul (a5, a0, a2);
                Kadd (r2, a4, a5);
                });

#ifndef HALF_WORD
        DO_ONE_TEST("mul_ui distributivity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Krandom2 (a2, rstate);
                unsigned long xx;
                xx = ((unsigned long *)(void *)(a0))[0];
                Kadd (a3, a1, a2);
                Kmul_ui (r1, a3, xx);
                Kmul_ui (a4, a1, xx);
                Kmul_ui (a5, a2, xx);
                Kadd (r2, a4, a5);
                });
#endif

        /* we can't be sure that our polynomial is irreducible... */
        DO_ONE_TEST("inversion", {
                do {
                    Krandom2 (r1, rstate);
                    if (Kcmp_ui(r1, 0) == 0) continue;
                } while (!Kinv (a1, r1));
                Kmul (a2, a1, r1);
                Kmul (r2, r1, a2);
                });

        DO_ONE_TEST("reduce o mul_ur = mul", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kelt_ur tmp;
                Kelt_ur_init(&tmp);
                Kmul_ur(tmp, a0, a1);
                Kreduce(r1, tmp);
                Kelt_ur_clear(&tmp);
                Kmul(r2, a0, a1);
                });

        DO_ONE_TEST("reduce o sqr_ur = sqr", {
                Krandom2 (a0, rstate);
                Kelt_ur tmp;
                Kelt_ur_init(&tmp);
                Ksqr_ur(tmp, a0);
                Kreduce(r1, tmp);
                Kelt_ur_clear(&tmp);
                Ksqr(r2, a0);
                });
        /* currently mpfq_pz does not have sqrt implemented */
#ifndef VARIABLE_SIZE_PRIME
        DO_ONE_TEST("sqr o sqrt o sqr = sqr", {
                Krandom2 (a0, rstate);
                /* force testing x==0 at least once ! */
                if (i == 0) Kset_ui (a0, 0);
                Ksqr(r1, a0);
                Ksqrt(r2, r1);
                Ksqr(r2, r2);
                });
#endif
        /*-----------------------------------------------------------*/
        /*          Tests specific to prime fields                   */
        /*-----------------------------------------------------------*/
#ifndef CHAR2
        DO_ONE_TEST("sscan o asprint = id", {
                Krandom2 (a0, rstate);
                /* force testing x==0 at least once ! */
                if (i == 0) Kset_ui (a0, 0);
                char *str;
                Kset(r1, a0);
                Kasprint(&str, a0);
                int ret = Ksscan(r2, str);
                free(str);
                if (!ret) abort();
                });

        {
            /* Now do some I/O tests */
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            fprintf(f, "\n");   /* ensure we properly ignore leading ws */
            DO_ONE_TEST("fprint", {
                    Kset_ui(a0, i);
                    Kinv(a0, a0);
                    Kfprint(f, a0);
                    fprintf(f, "\n");
                    Kset_ui(r1, 1);
                    Kset(r2, r1);
                    });
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            DO_ONE_TEST("fscan", {
                    Kset_ui(a0, i);
                    Kfscan(f, a1);
                    Kinv(r1, a1);
                    Kset(r2, a0);
                    });
            fclose(f);
            unlink(filename);
            free(filename);
        }

        DO_ONE_TEST("mul by 3 = add o add", {
                Krandom2 (a0, rstate);
                Kset_ui(a1, 3);
                Kmul (r1, a0, a1);
                Kadd (r2, a0, a0);
                Kadd (r2, r2, a0);
                });
#ifndef EXTENSION_OF_GFP
        DO_ONE_TEST("Fermat by pow", {
                Krandom2 (a0, rstate);
                Kset(r1, a0);
                Kpowz(r2, a0, K->p);
                });
#endif
        DO_ONE_TEST("x^(q-1) == 1/x", {
                do { Krandom2(a0, rstate); } while (Kcmp_ui(a0,0) == 0);
                mpz_t z;
                mpz_init_set_si(z, -1);
                Kpowz(r1, a0, z);
                Kinv(r2, a0);
                mpz_clear(z);
                });
        DO_ONE_TEST("x^0 == 1", {
                Krandom2 (a0, rstate);
                mpz_t z;
                mpz_init_set_ui(z, 0);
                Kpowz(r1, a0, z);
                Kset_ui(r2, 1);
                mpz_clear(z);
                });
        DO_ONE_TEST("x^(#K^2) == x", {
                Krandom2 (a0, rstate);
                mpz_t z;
                mpz_init(z);
                Kfield_characteristic(z);
                mpz_pow_ui(z, z, Kfield_degree() * 2);
                Kpowz(r1, a0, z);
                Kset(r2, a0);
                mpz_clear(z);
                });
#ifndef VARIABLE_SIZE_PRIME
        /* pz has no Tonelli Shanks for now */
        DO_ONE_TEST("is_sqr o (mul(sqr,nsqr)) = false", {
                do {
                    Krandom2 (a0, rstate);
                } while (Kcmp_ui(a0, 0) == 0);
                Ksqr(r1, a0);
                Kmul(r1, r1, (Ksrc_elt)K->ts_info.z);
                Kset_ui(r2, Kis_sqr(r1));
                if (Kcmp_ui(r1, 0) == 0)
                    Kset_ui(r1, 1);
                else
                    Kset_ui(r1, 0);
                });
#endif
        DO_ONE_TEST("is_sqr o sqr = true", {
                do {
                    Krandom2 (a0, rstate);
                } while (Kcmp_ui(a0, 0) == 0);
                Kset_ui(r1, 1);
                Ksqr(a0, a0);
                Kset_ui(r2, Kis_sqr(a0));
                });
        DO_ONE_TEST("ur_add 500 times and reduce", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kelt_ur tmp0;
                Kelt_ur tmp1;
                Kelt_ur_init(&tmp0);
                Kelt_ur_init(&tmp1);
                Kmul_ur(tmp0, a0, a1);
                Kelt_ur_set_ui(tmp1, 0);
                {
                  int j;
                  for (j = 0; j < 500; ++j)
                    Kelt_ur_add(tmp1, tmp1, tmp0);
                }
                Kreduce(r1, tmp1);
                Kelt_ur_clear(&tmp0);
                Kelt_ur_clear(&tmp1);
                Kmul(r2, a0, a1);
                Kmul_ui(r2, r2, 500);
                });

        DO_ONE_TEST("ur_sub 500 times and reduce", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kelt_ur tmp0;
                Kelt_ur tmp1;
                Kelt_ur_init(&tmp0);
                Kelt_ur_init(&tmp1);
                Kmul_ur(tmp0, a0, a1);
                Kelt_ur_set_ui(tmp1, 0);
                {
                  int j;
                  for (j = 0; j < 500; ++j)
                    Kelt_ur_sub(tmp1, tmp1, tmp0);
                }
                Kreduce(r1, tmp1);
                Kelt_ur_clear(&tmp0);
                Kelt_ur_clear(&tmp1);
                Kneg(r2, a0);
                Kmul(r2, r2, a1);
                Kmul_ui(r2, r2, 500);
                });

        DO_ONE_TEST("ur_neg o ur_sub = ur_add", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Krandom2 (a2, rstate);
                Krandom2 (a3, rstate);
                Kelt_ur tmp0;
                Kelt_ur tmp1;
                Kelt_ur tmp2;
                Kelt_ur_init(&tmp0);
                Kelt_ur_init(&tmp1);
                Kelt_ur_init(&tmp2);
                Kmul_ur(tmp0, a0, a1);
                Kmul_ur(tmp1, a2, a3);
                Kelt_ur_sub(tmp2, tmp1, tmp0);
                Kreduce(r1, tmp2);
                Kelt_ur_neg(tmp0, tmp0);
                Kelt_ur_add(tmp2, tmp1, tmp0);
                Kreduce(r2, tmp2);
                Kelt_ur_clear(&tmp0);
                Kelt_ur_clear(&tmp1);
                Kelt_ur_clear(&tmp2);
                });

#endif

        /*-----------------------------------------------------------*/
        /*          Tests specific to some I/O corner cases          */
        /*-----------------------------------------------------------*/

        {
            int ntests = 1;
            DO_ONE_TEST("sscan(garbage) returns 0", {
                    Kset_ui(r1, 0);
                    Kset_ui(r2, Ksscan(a0, "garbage"));
                    });
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-dummy-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            DO_ONE_TEST("fscan(empty file) returns 0", {
                    Kset_ui(r1, 0);
                    Kset_ui(r2, Kfscan(f, a1));
                    });
            fclose(f);
            unlink(filename);
            free(filename);
        }
        /* This is meant to exert the realloc() feature in fscan */
        DO_ONE_TEST("fscan(long file)", {
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-dummy-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            mpz_t z;
            mpz_init(z);
            Kfield_characteristic(z);
            gmp_fprintf(f, "%Zd", z);
            mpz_clear(z);
            for(int j = 0 ; j < i ; j++) gmp_fprintf(f, "0");
            gmp_fprintf(f, "1");
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            Kset_ui(r1, 1);
            Kset_ui(r2, 0);
            Kfscan(f, r2);
            fclose(f);
            unlink(filename);
            free(filename);
        });



        /*-----------------------------------------------------------*/
        /*          Tests specific to Montgomery representation      */
        /*-----------------------------------------------------------*/
#ifdef MGY
        DO_ONE_TEST("mgy_enc o mgy_dec = id", {
                Krandom2 (a0, rstate);
                Kset(r1, a0);
                Kmgy_enc(r2, a0);
                Kmgy_dec(r2, r2);
                });
#endif
        /*-----------------------------------------------------------*/
        /*          Tests specific to characteristic 2               */
        /*-----------------------------------------------------------*/

#ifdef CHAR2
        DO_ONE_TEST("add_uipoly o sub_uipoly = id", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Ksub_uipoly (r1, a0, a1[0]);
                Kadd_uipoly (r1, r1, a1[0]);
                Kset (r2, a0);
                });

        DO_ONE_TEST("mul_uipoly distributivity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Krandom2 (a2, rstate);
                Kadd (a3, a1, a2);
                Kmul_uipoly(r1, a3, a0[0]);
                Kmul_uipoly(a4, a1, a0[0]);
                Kmul_uipoly(a5, a2, a0[0]);
                Kadd (r2, a4, a5);
                });

        DO_ONE_TEST("squaring period", {
                Krandom2 (r1, rstate);
                Kset (r2, r1);
                for(int j = 0 ; j < Kdegree ; j++) {
                        Ksqr (r2, r2);
                }
        });

        DO_ONE_TEST("inv by pow", {
                // inv of 0 is undefined.
                do { Krandom2(a0, rstate); } while (Kcmp_ui(a0,0) == 0);
                mpz_t zz;
                mpz_init_set_ui(zz, 1);
                mpz_mul_2exp(zz, zz, Kfield_degree());
                mpz_sub_ui(zz, zz, 2);
                Kpowz(r1, a0, zz);
                mpz_clear(zz);
                Kinv(r2, a0);
                });

        DO_ONE_TEST("sqrt linearity", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Kadd (a2, a0, a1);
                Ksqrt(a3, a0);
                Ksqrt(a4, a1);
                Kadd (r1, a3, a4);
                Ksqrt(r2, a2);
        });

        DO_ONE_TEST("artin-schreier equation", {
                Krandom2 (r1, rstate);
                Ksqr(a1, r1);
                Kadd(a1, a1, r1);
                Kas_solve(r2, a1);

                /* FIXME ; users of the library may want to be able to do
                 * this kind of thing */
                r1[0] &= ~1UL;
                r2[0] &= ~1UL;
        });

        DO_ONE_TEST("artin-schreier equation (2)", {
                Krandom2 (a0, rstate);
                if (Ktrace(a0) != 0) continue;
                Krandom2 (a2, rstate);
                if (Ktrace(a2) != 0) continue;

                Kas_solve(a1, a0);
                Ksqr(r1, a1);
                Kadd(r1, r1, a1);
                Kadd(r1, r1, a0);

                Kas_solve(a3, a2);
                Ksqr(r2, a3);
                Kadd(r2, r2, a3);
                Kadd(r2, r2, a2);
        });

        DO_ONE_TEST("ur_add and reduce", {
                Krandom2 (a0, rstate);
                Krandom2 (a1, rstate);
                Krandom2 (a2, rstate);
                Krandom2 (a3, rstate);
                Kelt_ur tmp0;
                Kelt_ur tmp1;
                Kelt_ur_init(&tmp0);
                Kelt_ur_init(&tmp1);
                Kmul_ur(tmp0, a0, a1);
                Kmul_ur(tmp1, a2, a3);
                Kelt_ur_add(tmp1, tmp1, tmp0);
                Kreduce(r1, tmp1);
                Kelt_ur_clear(&tmp0);
                Kelt_ur_clear(&tmp1);
                Kmul(r2, a0, a1);
                Kmul(a4, a2, a3);
                Kadd(r2, r2, a4);
                });

        DO_ONE_TEST("sscan o asprint = id, base 10", {
                Krandom2 (a0, rstate);
                unsigned long base=10;
                Kfield_setopt(MPFQ_IO_TYPE, &base);
                char *str;
                Kset(r1, a0);
                Kasprint(&str, a0);
                int ret = Ksscan(r2, str);
                free(str);
                if (!ret) abort();
                });

        DO_ONE_TEST("sscan o asprint = id, base 2", {
                Krandom2 (a0, rstate);
                unsigned long base=2;
                Kfield_setopt(MPFQ_IO_TYPE, &base);
                char *str;
                Kset(r1, a0);
                Kasprint(&str, a0);
                int ret = Ksscan(r2, str);
                free(str);
                if (!ret) abort();
                });

         DO_ONE_TEST("sscan o asprint = id, base 16", {
                Krandom2 (a0, rstate);
                unsigned long base=16;
                Kfield_setopt(MPFQ_IO_TYPE, &base);
                char *str;
                Kset(r1, a0);
                Kasprint(&str, a0);
                int ret = Ksscan(r2, str);
                free(str);
                if (!ret) abort();
                });
#endif

#ifdef  HAVE_mpfq_p_1_hadamard
         /* This does no checking of course. However it is useful to
          * trigger linking of the function, since when it is done in
          * assembly, it is likely to trigger register allocation
          * problems in some instances */
         Khadamard(a0,a1,a2,a3);
#endif
        /*-----------------------------------------------------------*/
        /*          Tests related to vectors                         */
        /*-----------------------------------------------------------*/
        const int length=7;
//        int test_length;
        Kvec_init(&v1, 2*length);
        Kvec_init(&v2, 2*length);
        Kvec_init(&w1, 2*length);
        Kvec_init(&w2, 2*length);
        Kvec_init(&w3, 2*length);
        Kvec_init(&w4, 2*length);
        DO_ONE_TEST_VEC("vec_add commutativity", {
                Kvec_random2(w1, length, rstate);
                Kvec_random2(w2, length, rstate);
                Kvec_add(v1, w1, w2, length);
                Kvec_add(v2, w2, w1, length);
//                test_length = length;
                });
        DO_ONE_TEST_VEC("vec_add associativity", {
                Kvec_random2(w1, length, rstate);
                Kvec_random2(w2, length, rstate);
                Kvec_random2(w3, length, rstate);
                Kvec_add(v1, w1, w2, length);
                Kvec_add(v1, v1, w3, length);
                Kvec_add(v2, w2, w3, length);
                Kvec_add(v2, v2, w1, length);
//                test_length = length;
                });
        DO_ONE_TEST_VEC("vec linearity", {
                Krandom2 (a0, rstate);
                Kvec_random2(w1, length, rstate);
                Kvec_random2(w2, length, rstate);
                Kvec_scal_mul(v1, w1, a0, length);
                Kvec_scal_mul(v2, w2, a0, length);
                Kvec_add(v1, v1, v2, length);
                Kvec_add(w1, w1, w2, length);
                Kvec_scal_mul(v2, w1, a0, length);
//                test_length = length;
                });
        DO_ONE_TEST_VEC("vec_conv linearity", {
                Kvec_random2(w1, length, rstate);
                Kvec_random2(w2, length, rstate);
                Kvec_random2(w3, length, rstate);
                Kvec_add(w4, w2, w3, length);
                Kvec_conv(v1, w1, length, w4, length);
                Kvec_conv(w4, w1, length, w2, length);
                Kvec_conv(v2, w1, length, w3, length);
                Kvec_add(v2, v2, w4, 2*length-1);
//                test_length = 2*length-1;
                });


        DO_ONE_TEST("vec_sscan o vec_asprint = id", {
                /* Do this for vectors of length i */
                unsigned int cap = i < length ? i : length;
                Kvec_random2(w1, cap, rstate);
                char * str;
                Kvec_asprint(&str, w1, cap);
                Kvec tvec;
                unsigned int tlength = 0;
                Kvec_init(&tvec, tlength);
                Kvec_sscan(&tvec, &tlength, str);
                if (tlength != (unsigned int) cap) abort();
                Kvec_set(v1, w1, cap);
                Kvec_set(v2, tvec, cap);
                Kvec_clear(&tvec, tlength);
                free(str);
                });

        {
            /* Now do some I/O tests */
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-vec-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            DO_ONE_TEST("vec_fprint", {
                    for(unsigned int j = 0 ; j < (unsigned int) length ; j++) {
                        Kset_ui(a0, i * length + j);
                        Kinv(a0, a0);
                        Kset(Kvec_coeff_ptr(v1, j), a0);
                    }
                    Kvec_fprint(f, v1, length);
                    fprintf(f, "\n");
                    });
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            DO_ONE_TEST("vec_fscan", {
                    for(unsigned int j = 0 ; j < (unsigned int) length ; j++) {
                        Kset_ui(a0, i * length + j);
                        Kinv(a0, a0);
                        Kset(Kvec_coeff_ptr(v1, j), a0);
                    }
                    Kvec tvec;
                    unsigned int tlength = 0;
                    Kvec_init(&tvec, tlength);
                    Kvec_fscan(f, &tvec, &tlength);
                    if (tlength != (unsigned int) length) abort();
                    Kvec_set(v2, tvec, length);
                    Kvec_clear(&tvec, tlength);
                    });
            fclose(f);
            unlink(filename);
            free(filename);
        }

        /* Some I/O corner cases */
        {
            int ntests = 1;
            DO_ONE_TEST("vec_sscan(garbage) returns 0", {
                    unsigned int tlength = 0;
                    Kvec tvec;
                    Kvec_init(&tvec, tlength);
                    Kvec_random(v1, length, rstate);
                    Kvec_set(v2, v1, length);
                    Kset_ui(Kvec_coeff_ptr(v1, 0), 0);
                    int r = Kvec_sscan(&tvec, &tlength, "garbage");
                    if (tlength) abort();
                    Kset_ui(Kvec_coeff_ptr(v2, 0), r);
                    Kvec_clear(&tvec, tlength);
                    });
            DO_ONE_TEST("vec_sscan(\"[garbage]\") returns 0", {
                    unsigned int tlength = 0;
                    Kvec tvec;
                    Kvec_init(&tvec, tlength);
                    Kvec_random(v1, length, rstate);
                    Kvec_set(v2, v1, length);
                    Kset_ui(Kvec_coeff_ptr(v1, 0), 0);
                    int r = Kvec_sscan(&tvec, &tlength, "[garbage]");
                    if (tlength) abort();
                    Kset_ui(Kvec_coeff_ptr(v2, 0), r);
                    Kvec_clear(&tvec, tlength);
                    });
            DO_ONE_TEST("vec_sscan(\"[1 2]\") returns 0", {
                    unsigned int tlength = 0;
                    Kvec tvec;
                    Kvec_init(&tvec, tlength);
                    Kvec_random(v1, length, rstate);
                    Kvec_set(v2, v1, length);
                    Kset_ui(Kvec_coeff_ptr(v1, 0), 0);
                    int r = Kvec_sscan(&tvec, &tlength, "[1 2]");
                    if (tlength) abort();
                    Kset_ui(Kvec_coeff_ptr(v2, 0), r);
                    Kvec_clear(&tvec, tlength);
                    });
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-vec_dummy-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            DO_ONE_TEST("vec_fscan(empty file) returns 0", {
                    unsigned int tlength = 0;
                    Kvec tvec;
                    Kvec_init(&tvec, tlength);
                    Kvec_random(v1, length, rstate);
                    Kvec_set(v2, v1, length);
                    Kset_ui(Kvec_coeff_ptr(v1, 0), 0);
                    int r = Kvec_fscan(f, &tvec, &tlength);
                    if (tlength) abort();
                    Kset_ui(Kvec_coeff_ptr(v2, 0), r);
                    Kvec_clear(&tvec, tlength);
                    });
            fclose(f);
            unlink(filename);
            free(filename);
        }
        /* This is meant to exert the realloc() feature in vec_fscan */
        DO_ONE_TEST("vec_fscan(long file)", {
            char * filename;
            int rc = asprintf(&filename, "/tmp/mpfq-dummy-test.%lu", gmp_urandomb_ui(rstate, 32));
            if (rc < 0) abort();
            FILE * f = fopen(filename, "w");
            if (f == NULL) abort();
            fprintf(f, "[");
            for(int j = 0 ; j < i ; j++) fprintf(f, " ");
            fprintf(f, "1]");
            fclose(f);
            f = fopen(filename, "r");
            if (f == NULL) abort();
            unsigned int tlength = 0;
            Kvec tvec;
            Kvec_init(&tvec, tlength);
            Kvec_random(v1, length, rstate);
            Kvec_set(v2, v1, length);
            int r = Kvec_fscan(f, &tvec, &tlength);
            if (!r) abort();
            if (tlength != 1) abort();
            Kset_ui(Kvec_coeff_ptr(v1, 0), 1);
            Kvec_set(v2, tvec, 1);
            Kvec_clear(&tvec, tlength);
            fclose(f);
            unlink(filename);
            free(filename);
        });

        Kvec_clear(&v1, 2*length);
        Kvec_clear(&v2, 2*length);
        Kvec_clear(&w1, 2*length);
        Kvec_clear(&w2, 2*length);
        Kvec_clear(&w3, 2*length);
        Kvec_clear(&w4, 2*length);


        /*-----------------------------------------------------------*/
        /*          Tests related to polynomials                     */
        /*-----------------------------------------------------------*/

        const int deg = length - 1;

        Kpoly_init(p1, 2*deg);
        Kpoly_init(p2, 2*deg);
        Kpoly_init(q1, 2*deg);
        Kpoly_init(q2, 2*deg);
        Kpoly_init(q3, 2*deg);
        Kpoly_init(q4, 2*deg);
        DO_ONE_TEST_POLY("poly_add commutativity", {
                if (!(i & 0xf)) {
                    /* test realloc every once in a while */
                    Kpoly_clear(q1);
                    Kpoly_init(q1, 1);
                }
                Kpoly_random(q1, deg, rstate);
                Kpoly_random(q2, deg, rstate);
                Kpoly_add(p1, q1, q2);
                Kpoly_add(p2, q2, q1);
                });
        DO_ONE_TEST_POLY("poly_sub = poly_add o poly_neg", {
                Kpoly_random2 (q1, deg, rstate);
                Kpoly_random2 (q2, deg, rstate);
                if (!(i & 0xf)) {
                    /* test realloc every once in a while */
                    Kpoly_clear(p1);
                    Kpoly_init(p1, 1);
                }
                Kpoly_sub (p1, q1, q2);
                Kpoly_neg(q4, q2);
                Kpoly_add(p2, q4, q1);
                });
        DO_ONE_TEST_POLY("poly_add associativity", {
                Kpoly_random(q1, deg, rstate);
                Kpoly_random(q2, deg, rstate);
                Kpoly_random(q3, deg, rstate);
                Kpoly_add(p1, q1, q2);
                Kpoly_add(p1, p1, q3);
                Kpoly_add(p2, q2, q3);
                Kpoly_add(p2, p2, q1);
                });
        DO_ONE_TEST("poly_add_ui(y) o poly_neg o poly_set_ui(x) == set(y-x)", {
                mpz_t z;
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                unsigned long y = gmp_urandomb_ui(rstate, 32);
                mpz_init_set_ui(z, y);
                mpz_sub_ui(z, z, x);
                Kset_mpz(a0, z);
                mpz_clear(z);

                p1->size = 1;
                Kpoly_setcoeff(p1, a0, 0);

                Kpoly_set_ui(p1, x);
                Kpoly_neg(p1, p1);
                Kpoly_add_ui(p1, p1, y);
                });

        DO_ONE_TEST("add_ui o sub_ui = id", {
                if (!(i & 0xf)) {
                    /* test realloc every once in a while */
                    Kpoly_clear(p1);
                    Kpoly_init(p1, 1);
                }
                Kpoly_random2 (p1, deg, rstate);
                Kpoly_set(p2, p1);
                unsigned long x = gmp_urandomb_ui(rstate, 32);
                Kpoly_sub_ui (p2, p1, x);
                Kpoly_add_ui (p2, p2, x);
                });

        DO_ONE_TEST_POLY("poly linearity", {
                Krandom2 (a0, rstate);
                if (!(i & 0xf)) {
                    /* test realloc every once in a while */
                    Kpoly_clear(p2);
                    Kpoly_init(p2, 1);
                }

                Kpoly_random(q1, deg, rstate);
                Kpoly_random(q2, deg, rstate);
                Kpoly_scal_mul(p1, q1, a0);
                Kpoly_scal_mul(p2, q2, a0);
                Kpoly_add(p1, p1, p2);
                Kpoly_add(q1, q1, q2);
                Kpoly_scal_mul(p2, q1, a0);
                });


        DO_ONE_TEST_POLY("poly_mul linearity", {
                Kpoly_random(q1, deg, rstate);
                Kpoly_random(q2, deg, rstate);
                Kpoly_random(q3, deg, rstate);
                Kpoly_add(q4, q2, q3);
                Kpoly_mul(p1, q1, q4);
                Kpoly_mul(q4, q1, q2);
                Kpoly_mul(p2, q1, q3);
                Kpoly_add(p2, p2, q4);
                });
        DO_ONE_TEST_POLY("poly_divmod", {
                Kpoly_random(p1, deg, rstate);
                do {
                    Kpoly_random(q2, deg, rstate);
                } while (Kpoly_deg(q2) == 0);
                if (i < 10) Kpoly_setmonic(q2, q2);
                Kpoly_divmod(q3, q4, p1, q2);
                Kpoly_mul(q1, q2, q3);
                Kpoly_add(p2, q1, q4);
                });
        DO_ONE_TEST_POLY("poly_gcd", {
                /* make sure leading coefficients are invertible.
                 * Normally we're over a field, so we don't have to.
                 * Unfortunately our testing is so dumb that we don't
                 * promise that we're over a field...
                 */
                do {
                    Kpoly_random(q1, deg, rstate);
                    Kpoly_random(q2, deg, rstate);
                    Kpoly_getcoeff(a1, q1, deg);
                    if (!Kinv(a0, a1)) continue;
                    Kpoly_getcoeff(a2, q2, deg);
                    if (!Kinv(a0, a2)) continue;
                    Kpoly_gcd(p1, q1, q2);
                } while (Kpoly_deg(p1) != 0);
                do {
                    Kpoly_random(q3, deg,rstate);
                    Kpoly_getcoeff(a3, q3, deg);
                } while (!Kinv(a0, a3));
                Kpoly_mul(q1, q1, q3);
                Kpoly_mul(q2, q2, q3);
                Kpoly_setmonic(p2, q3);
                Kpoly_gcd(p1, q1, q2);
                });
        DO_ONE_TEST_POLY("poly_xgcd", {
                do {
                    Kpoly_random(q1, deg, rstate);
                    Kpoly_random(q2, deg, rstate);
                    Kpoly_gcd(p1, q1, q2);
                } while (Kpoly_deg(p1) != 0);
                Kpoly_random(q3, deg, rstate);
                Kpoly_setmonic(p2, q3);
                Kpoly_mul(q1, q1, q3);
                Kpoly_mul(q2, q2, q3);
                Kpoly_xgcd(p1, q3, q4, q1, q2);
                Kpoly_mul(q3, q3, q1);
                Kpoly_mul(q4, q4, q2);
                Kpoly_add(p1, q3, q4);
                });

#if 0   /* grrr, internal... */
        DO_ONE_TEST_POLY("poly_preinv", {
                do { Kpoly_random(q2, deg, rstate); } while (Kpoly_deg(q2) < 0);
                Kpoly_setcoeff_ui(q2, 0, 1);
                Kpoly_preinv(q1, q2, deg);
                Kpoly_mul(q3, q1, q2);
                /* we expect q3=q1*q2 = 1 + X^deg*something */
                Kpoly_set_ui(q4, 0);
                Kpoly_setcoeff_ui(q4, deg, 0);
                Kpoly_divmod(p1, p2, q3, q4);
                /* p1 should be 1 */
                Kpoly_set_ui(p2, 1);
                });
#endif
        DO_ONE_TEST_POLY("poly_mod_pre", {
                do { Kpoly_random(q2, deg, rstate); } while (Kpoly_deg(q2) < 0);
                Kpoly_random(q1, i==0 ? (deg/2) : (2*Kpoly_deg(q2)-2), rstate);

                Kpoly_setmonic(q2, q2);
                if (i&1) {
                    Kpoly_precomp_mod(q3, q2);
                } else {
                    /* test realloc feature every once in a while */
                    Kpoly toosmall;
                    Kpoly_init(toosmall, 1);
                    Kpoly_precomp_mod(toosmall, q2);
                    Kpoly_set(q3, toosmall);
                    Kpoly_clear(toosmall);
                }


                Kpoly_mod_pre(p1, q1, q2, q3);
                Kpoly_divmod(q4, p2, q1, q2);
                });


        /* corner cases */
        {
            int ntests = 1;
            DO_ONE_TEST_POLY("poly_setmonic(0) == 0", {
                    Kpoly_set_ui(p1, 0);
                    Kpoly_set_ui(q1, 0);
                    Kpoly_setmonic(p2, q1);
                    });
            DO_ONE_TEST_POLY("poly_divmod(x, 0) returns 0", {
                    Kpoly_random(q1, deg, rstate);
                    Kpoly_set_ui(p2, 0);
                    int r = Kpoly_divmod(q2, q3, q1, p2);
                    Kpoly_set_ui(p1, r);
                    });
            DO_ONE_TEST_POLY("poly_divmod(0, x) == 0", {
                    do { Kpoly_random(q2, deg, rstate); } while (Kpoly_deg(q2) < 0);
                    Kpoly_set_ui(p1, 0);
                    int r = Kpoly_divmod(p2, q3, p1, q2);
                    if (!r) abort();
                    });
            /* same but with a denormalized polynomial */
            DO_ONE_TEST_POLY("poly_divmod(x, denormalized 0) returns 0", {
                    Kpoly_random(q1, deg, rstate);
                    Kpoly_setcoeff_ui(p2, 0, 1);
                    int r = Kpoly_divmod(q2, q3, q1, p2);
                    Kpoly_set_ui(p2, 0);
                    Kpoly_set_ui(p1, r);
                    });
            DO_ONE_TEST_POLY("poly_divmod(denormalized 0, x) == 0", {
                    do { Kpoly_random(q2, deg, rstate); } while (Kpoly_deg(q2) < 0);
                    Kpoly_setcoeff_ui(p1, 0, 1);
                    int r = Kpoly_divmod(p2, q3, p1, q2);
                    if (!r) abort();
                    Kpoly_set_ui(p1, 0);
                    });
            DO_ONE_TEST_POLY("poly_divmod(a, larger than a) returns 0+a", {
                    Kpoly_random(p1, deg-1, rstate);
                    do { Kpoly_random(q2, deg, rstate); } while (Kpoly_deg(q2) < 0);
                    Kpoly_divmod(q3, p2, p1, q2);
                    });

            DO_ONE_TEST_POLY("poly_gcd(p, 0) == p", {
                    /* make sure leading coefficients are invertible.
                     * Normally we're over a field, so we don't have to.
                     * Unfortunately our testing is so dumb that we don't
                     * promise that we're over a field...
                     */
                    do { Kpoly_random(p1, deg, rstate); } while (Kpoly_deg(p1) < 0);
                    Kpoly_set_ui(q2, 0);
                    Kpoly_gcd(p2, p1, q2);
                    });
            DO_ONE_TEST_POLY("scal_mul(0)(x) == set_ui(0)", {
                    Kpoly_random(q1, deg, rstate);
                    Kset_ui(a0, 0);
                    Kpoly_scal_mul(p1, q1, a0);
                    Kpoly_set_ui(p2, 0);
                    });

            DO_ONE_TEST_POLY("poly_setcoeff_ui(large index)", {
                    Kpoly toosmall;
                    Kpoly_init(toosmall, 1);
                    Kpoly_setcoeff_ui(toosmall, 1, deg);
                    Kpoly_set(p1, toosmall);
                    p2->size = 0;
                    Kpoly_setcoeff_ui(p2, 1, deg);
                    Kpoly_clear(toosmall);
                    });
            DO_ONE_TEST_POLY("poly_setcoeff(large index)", {
                    Kpoly toosmall;
                    Krandom2 (a0, rstate);
                    Kpoly_init(toosmall, 1);
                    Kpoly_setcoeff(toosmall, a0, deg);
                    Kpoly_set(p1, toosmall);
                    Kpoly_setcoeff(p2, a0, deg);
                    Kpoly_clear(toosmall);
                    });
            DO_ONE_TEST_POLY("poly_getcoeff(too large) == 0", {
                    Kpoly_set_ui(p1, 0);
                    q1->size = 0;
                    Kpoly_getcoeff(a0, q1, 1);
                    p2->size = 1;
                    Kpoly_setcoeff(p2, a0, 0);
                    });
            
            DO_ONE_TEST_POLY("poly_xgcd(0, 0)", {
                    Kpoly_set_ui(q1, 0);
                    Kpoly_set_ui(q2, 0);
                    Kpoly_xgcd(p1, q3, q4, q1, q2);
                    Kpoly_mul(q3, q3, q1);
                    Kpoly_mul(q4, q4, q2);
                    Kpoly_add(p2, q3, q4);
                    });

            DO_ONE_TEST_POLY("poly_xgcd(0, 0)", {
                    do { Kpoly_random(q1, deg, rstate); } while (Kpoly_deg(q1) <= 0);
                    Kpoly_mul(q2, q1, q1);
                    int r = Kpoly_cmp(q2, q1); /* should be deg(q1) */
                    Kpoly_set_ui(p1, r);
                    Kpoly_set_ui(p2, Kpoly_deg(q1));
                    });

            ntests = 10;

            DO_ONE_TEST_POLY("poly_xgcd(p, 0)", {
                    do {
                        Kpoly_random(q1, deg, rstate);
                    } while (Kpoly_deg(q1) < 0);
                    Kpoly_set_ui(q2, 0);
                    Kpoly_xgcd(p1, q3, q4, q1, q2);
                    Kpoly_mul(q3, q3, q1);
                    Kpoly_mul(q4, q4, q2);
                    Kpoly_add(p2, q3, q4);
                    });
            DO_ONE_TEST_POLY("poly_xgcd(0, p)", {
                    do {
                        Kpoly_random(q2, deg, rstate);
                    } while (Kpoly_deg(q2) < 0);
                    Kpoly_set_ui(q1, 0);
                    Kpoly_xgcd(p1, q3, q4, q1, q2);
                    Kpoly_mul(q3, q3, q1);
                    Kpoly_mul(q4, q4, q2);
                    Kpoly_add(p2, q3, q4);
                    });
        }

        Kpoly_clear(p1);
        Kpoly_clear(p2);
        Kpoly_clear(q1);
        Kpoly_clear(q2);
        Kpoly_clear(q3);
        Kpoly_clear(q4);

        Kclear (&a0);
        Kclear (&a1);
        Kclear (&a2);
        Kclear (&a3);
        Kclear (&a4);
        Kclear (&a5);
        Kclear (&r1);
        Kclear (&r2);

        Kfield_clear();

        if (quiet) {
            fprintf(stderr, ".");
            fflush(stderr);
        }

        i++;
        seed = gmp_urandomb_ui(rstate, 64);

        gmp_randclear(rstate);
    }
}

/* vim:set ft=c: */