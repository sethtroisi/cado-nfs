#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "cado_poly.h"
#include "portability.h"

const char * sidenames[2] = {
    [RATIONAL_SIDE] = "rational",
    [ALGEBRAIC_SIDE] = "algebraic", };

void cado_poly_init(cado_poly poly)
{
    int i;

    /* ALL fields are zero upon init, EXCEPT the degree field (which is -1) */
    memset(poly, 0, sizeof(poly[0]));
    poly->rat = poly->pols[RATIONAL_SIDE];
    poly->alg = poly->pols[ALGEBRAIC_SIDE];

    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = poly->pols[side];
        ps->f = (mpz_t *) malloc((MAXDEGREE + 1) * sizeof(mpz_t));
        /* mpzs as well are zero */
        for (i = 0; i < (MAXDEGREE + 1); i++) {
            mpz_init_set_ui(ps->f[i], 0);
            ps->degree = -1;
        }
    }
    mpz_init_set_ui(poly->n, 0);
    mpz_init_set_ui(poly->m, 0);
}

void cado_poly_clear(cado_poly poly)
{
    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = poly->pols[side];
        for (int i = 0; i < (MAXDEGREE + 1); i++) {
            mpz_clear(ps->f[i]);
        }
        free(ps->f);
    }
    mpz_clear(poly->n);
    mpz_clear(poly->m);
    memset(poly, 0, sizeof(poly[0]));
}

/* p <- q */
void
cado_poly_set (cado_poly p, cado_poly q)
{
    mpz_set (p->n, q->n);
    p->skew = q->skew;
    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = p->pols[side];
        cado_poly_side_ptr qs = q->pols[side];
        ps->degree = qs->degree;
        for (int i = 0; i <= qs->degree; i++)
            mpz_set (ps->f[i], qs->f[i]);
    }
    mpz_set (p->m, q->m);
}

int cado_poly_set_plist(cado_poly poly, param_list pl)
{
    int have_n = 0;
    int have_f[2][(MAXDEGREE + 1)] = {{ 0, }, {0,}};
    int have_m = 0;
    int i;
    mpz_t tmp;

    have_n = param_list_parse_mpz(pl, "n", poly->n) || param_list_parse_mpz(pl, NULL, poly->n);
    poly->skew = 0.0; /* to ensure that we get an invalid skewness in case
                         it is not given */
    param_list_parse_double(pl, "skew", &(poly->skew));
    for (i = 0; i < (MAXDEGREE + 1); i++) {
        char tag[4];
        snprintf(tag, sizeof(tag), "c%d", i);
        have_f[ALGEBRAIC_SIDE][i] = param_list_parse_mpz(pl, tag, poly->alg->f[i]);
        if (!have_f[ALGEBRAIC_SIDE][i]) {
            snprintf(tag, sizeof(tag), "X%d", i);
            have_f[ALGEBRAIC_SIDE][i] = param_list_parse_mpz(pl, tag, poly->alg->f[i]);
        }
        snprintf(tag, sizeof(tag), "Y%d", i);
        have_f[RATIONAL_SIDE][i] = param_list_parse_mpz(pl, tag, poly->rat->f[i]);
    }
    param_list_parse_ulong(pl, "rlim", &(poly->rat->lim));
    param_list_parse_int(pl, "lpbr", &(poly->rat->lpb));
    param_list_parse_int(pl, "mfbr", &(poly->rat->mfb));
    param_list_parse_double(pl, "rlambda", &(poly->rat->lambda));

    param_list_parse_ulong(pl, "alim", &(poly->alg->lim));
    param_list_parse_int(pl, "lpba", &(poly->alg->lpb));
    param_list_parse_int(pl, "mfba", &(poly->alg->mfb));
    param_list_parse_double(pl, "alambda", &(poly->alg->lambda));

    mpz_set_ui (poly->m, 0);
    param_list_parse_mpz(pl, "m", poly->m);
    have_m = mpz_cmp_ui (poly->m, 0);

    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = poly->pols[side];
        int d;
        for (d = MAXDEGREE; d >= 0 && !have_f[side][d]; d--) {
            if (have_f[side][d] && mpz_cmp_ui(ps->f[d], 0) != 0)
                break;
        }
        ps->degree = d;
    }

    ASSERT_ALWAYS(have_n);

    // compute m, the common root of f and g mod n
    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = poly->pols[side];
        if (ps->degree != 1) continue;
        mpz_init (tmp);
        mpz_invert (tmp, ps->f[1], poly->n);
        mpz_mul (tmp, tmp, ps->f[0]);
        mpz_mod (tmp, tmp, poly->n);
        mpz_sub (tmp, poly->n, tmp);
        mpz_mod (tmp, tmp, poly->n);

        if (have_m && (mpz_cmp (poly->m, tmp) != 0))
        {
            fprintf (stderr, "m is not a root of g mod N\n");
            exit (EXIT_FAILURE);
        }
        have_m = 1;
        mpz_set (poly->m, tmp);
        mpz_clear(tmp);
    }

    if (!have_m) {
        fprintf (stderr, "must provide m for non-linear polynomials\n");
        exit (EXIT_FAILURE);
    }

    cado_poly_check(poly);

    return 1;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_stream(cado_poly poly, FILE * f)
{
    param_list pl;
    param_list_init(pl);
    param_list_read_stream(pl, f);
    int r = cado_poly_set_plist(poly, pl);
    param_list_clear(pl);
    return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read(cado_poly poly, const char *filename)
{
    FILE *file;
    int r;
    file = fopen(filename, "r");
    if (file == NULL) {
	fprintf(stderr, "read_polynomial: could not open %s\n", filename);
	return 0;
    }
    r = cado_poly_read_stream(poly, file);
    fclose(file);
    return r;
}


void fprint_polynomial(FILE * fp, mpz_t * f, const int d)
{
    int i, s, first = 1;
    mpz_t c;

    mpz_init(c);
    for (i = d; i >= 0; i--) {
	s = mpz_cmp_ui(f[i], 0);
	if (s != 0) {
	    if (s > 0) {
		if (first == 0)
		    gmp_fprintf(fp, " + ");
		gmp_fprintf(fp, "%Zd", f[i]);
	    } else if (s < 0) {
		mpz_abs(c, f[i]);
		gmp_fprintf(fp, " - %Zd", c);
	    }
	    first = 0;
	    if (i >= 2)
		gmp_fprintf(fp, "*x^%d", i);
	    else if (i == 1)
		gmp_fprintf(fp, "*x");
	}
    }
    mpz_clear(c);
    fprintf(fp, "\n");
}

/* check that m is a root of both f and g mod n */
void cado_poly_check(cado_poly cpoly)
{
    mpz_t tmp;

    mpz_init(tmp);

    for(int side = 0 ; side < 2 ; side++) {
        cado_poly_side_ptr ps = cpoly->pols[side];

        /* check m is a root of f mod n */
        mpz_set (tmp, ps->f[ps->degree]);
        for (int i = ps->degree - 1; i >= 0; i--)
        {
            mpz_mul (tmp, tmp, cpoly->m);
            mpz_add (tmp, tmp, ps->f[i]);
        }
        mpz_mod (tmp, tmp, cpoly->n);
        if (mpz_cmp_ui (tmp, 0) != 0)
        {
            fprintf (stderr, "m is not a root of the %s polynomial modulo n\n", sidenames[side]);
            exit (EXIT_FAILURE);
        }
    }

    mpz_clear(tmp);
}
