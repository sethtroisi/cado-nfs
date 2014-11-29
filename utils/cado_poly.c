#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "cado_poly.h"
#include "portability.h"

struct cado_poly_s cado_poly_struct;

const char * sidenames[2] = {
    [RATIONAL_SIDE] = "rational",
    [ALGEBRAIC_SIDE] = "algebraic", };

/* Be conservative and allocate two polynomials by default. */
void cado_poly_init(cado_poly poly)
{
    /* ALL fields are zero upon init, EXCEPT the degree field (which is -1) */
    memset(poly, 0, sizeof(poly[0]));
    poly->rat = poly->pols[RATIONAL_SIDE];
    poly->alg = poly->pols[ALGEBRAIC_SIDE];

    poly->nb_polys = 2;
    for(int side = 0 ; side < poly->nb_polys ; side++)
      mpz_poly_init (poly->pols[side], MAXDEGREE);

    mpz_init_set_ui(poly->n, 0);
}

void cado_poly_clear(cado_poly poly)
{
    for(int side = 0 ; side < poly->nb_polys ; side++)
      mpz_poly_clear (poly->pols[side]);

    mpz_clear(poly->n);
    memset(poly, 0, sizeof(poly[0]));
}

/* p <- q */
void
cado_poly_set (cado_poly p, cado_poly q)
{
    mpz_set (p->n, q->n);
    p->skew = q->skew;
    p->nb_polys = q->nb_polys;
    for(int side = 0 ; side < q->nb_polys ; side++)
      mpz_poly_set (p->pols[side], q->pols[side]);
}

// This function is no longer exported
#define BUF_MAX 10000

static
int cado_poly_set_plist(cado_poly poly, param_list pl)
{
    int have_n = 0;
    int have_f[NB_POLYS_MAX][(MAXDEGREE + 1)];
    int i, nb_polys = 2, new_coding = 0;

    for(i = 0; i < NB_POLYS_MAX; i++)
	have_f[i][0] = 0;

    have_n = param_list_parse_mpz(pl, "n", poly->n) 
	     || param_list_parse_mpz(pl, NULL, poly->n);
    poly->skew = 0.0; /* to ensure that we get an invalid skewness in case
                         it is not given */
    param_list_parse_double(pl, "skew", &(poly->skew));

#if 0 // TODO
    for(i = 0; i < NB_POLYS_MAX; i++){
	char tag[5], buf[BUF_MAX];
	snprintf(tag, sizeof(tag), "poly%d", i);
	if(param_list_parse_string(pl, tag, buf, BUF_MAX)){
	    new_coding = 1;
	    fprintf(stderr, "Read %s\n", buf);
	    nb_polys = i+1;
	}
	else
	    break;
    }
    exit(1);
#endif

    poly->nb_polys = nb_polys;
    if(new_coding == 0){
	/* reading polynomials coefficient by coefficient */
	for (i = 0; i < (MAXDEGREE + 1); i++) {
	    char tag[4];
	    snprintf(tag, sizeof(tag), "c%d", i);
	    have_f[ALGEBRAIC_SIDE][i] = param_list_parse_mpz(pl, tag, poly->alg->coeff[i]);
	    if (!have_f[ALGEBRAIC_SIDE][i]) {
		snprintf(tag, sizeof(tag), "X%d", i);
		have_f[ALGEBRAIC_SIDE][i] = param_list_parse_mpz(pl, tag, poly->alg->coeff[i]);
	    }
	    snprintf(tag, sizeof(tag), "Y%d", i);
	    have_f[RATIONAL_SIDE][i] = param_list_parse_mpz(pl, tag, poly->rat->coeff[i]);
	}
    }

    for(int side = 0 ; side < poly->nb_polys ; side++) {
        mpz_poly_ptr ps = poly->pols[side];
        int d;
        for (d = MAXDEGREE; d >= 0 && !have_f[side][d]; d--) {
            if (have_f[side][d] && mpz_cmp_ui(ps->coeff[d], 0) != 0)
                break;
        }
        ps->deg = d;
    }

    ASSERT_ALWAYS(have_n);

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

/* TODO: adapt for more than 2 polynomials:
 * compute for each pair (0,i) the corresponding common root m_i
 * check all m_i are equal
 */
int cado_poly_getm(mpz_ptr m, cado_poly_ptr cpoly, mpz_ptr N)
{
    // have to work with copies, because pseudo_gcd destroys its input
    mpz_poly_t f[2];

    ASSERT_ALWAYS(cpoly->nb_polys == 2);

    for (int i = 0; i < 2; ++i) {
        mpz_poly_init(f[i], cpoly->pols[i]->alloc);
        mpz_poly_set(f[i], cpoly->pols[i]);
    }
    int ret;
    ret = mpz_poly_pseudogcd_mpz(f[0], f[1], N, m); 
    if (ret) {
        ASSERT_ALWAYS(f[0]->deg == 1);
        mpz_t inv;
        mpz_init(inv);
        int ret2 = mpz_invert(inv, f[0]->coeff[1], N);
        // This inversion should always work.
        // If not, it means that N has a small factor (not sure we want
        // to be robust against that...)
        // Or maybe the polynomial selection was really bogus!
        ASSERT_ALWAYS(ret2);
        mpz_mul(inv, inv, f[0]->coeff[0]);
        mpz_neg(inv, inv);
        mpz_mod(m, inv, N);
        mpz_clear(inv);
    }
    for (int i = 0; i < 2; ++i)
        mpz_poly_clear(f[i]);
    return ret;
}

/* Return the rational side or -1 if two algebraic sides */
int
cado_poly_get_ratside (cado_poly_ptr pol)
{
  if (pol->pols[0]->deg != 1 && pol->pols[1]->deg != 1)
    return -1; /* two algrebraic sides */
  else if (pol->pols[0]->deg == 1)
    return 0;
  else
    return 1;
}
