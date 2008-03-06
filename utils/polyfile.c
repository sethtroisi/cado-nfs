#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>		/* for isdigit */
#include "cado.h"

static int separator(char x)
{
    return x == ':' || x == '=' || isspace(x);
}

/* gives a pointer to the config file info behind the key */
static const char *get_rhs(const char *line, const char *tag, int *have)
{
    unsigned int len = strlen(tag);
    if (strncmp(line, tag, len) != 0)
	return NULL;
    if (!separator(line[len]))
	return NULL;

    if (have != NULL && *have != 0) {
	fprintf(stderr, "parse error in get_rhs: %s appears twice\n", tag);
	exit(EXIT_FAILURE);
    }
    if (have) {
	*have = 1;
    }

    for (; !line[len] || separator(line[len]); len++);

    return line + len;
}

/* parse_string wants the size of the destination string */
static int parse_string(char *dst, unsigned int n,
			const char *line, const char *tag, int *have)
{
    const char *ptr = get_rhs(line, tag, have);
    if (ptr == NULL)
	return 0;
    strncpy(dst, ptr, n);
    return 1;
}

static void forbidden_garbage(const char *p, const char *line)
{
    for (; *p && isspace(*p); p++);
    if (*p) {
	fprintf(stderr, "parse error: garbage at end of line:\n%s\n", line);
	exit(EXIT_FAILURE);
    }
}

static int parse_ulong(unsigned long *dst,
		       const char *line, const char *tag, int *have)
{
    const char *ptr = get_rhs(line, tag, have);
    if (ptr == NULL)
	return 0;
    char *eptr;
    unsigned long res;
    res = strtoul(ptr, &eptr, 0);
    forbidden_garbage(eptr, line);
    if (dst)
	*dst = res;
    return 1;
}
static int parse_long(long *dst, const char *line, const char *tag, int *have)
{
    const char *ptr = get_rhs(line, tag, have);
    if (ptr == NULL)
	return 0;
    char *eptr;
    long res;
    res = strtol(ptr, &eptr, 0);
    forbidden_garbage(eptr, line);
    if (dst)
	*dst = res;
    return 1;
}
static int parse_int(int *dst, const char *line, const char *tag, int *have)
{
    long r;
    if (parse_long(&r, line, tag, have) == 0) {
	return 0;
    }
    if (dst)
	*dst = r;
    return 1;
}
static int parse_double(double *dst,
			const char *line, const char *tag, int *have)
{
    const char *ptr = get_rhs(line, tag, have);
    if (ptr == NULL)
	return 0;
    char *eptr;
    double res;
    res = strtod(ptr, &eptr);
    forbidden_garbage(eptr, line);
    if (dst)
	*dst = res;
    return 1;
}
static int parse_mpz(mpz_ptr dst,
		     const char *line, const char *tag, int *have)
{
    int rc;
    unsigned int nread;
    const char *ptr = get_rhs(line, tag, have);

    if (ptr == NULL)
	return 0;
    if (dst) {
	rc = gmp_sscanf(ptr, "%Zd%n", dst, &nread);
    } else {
	/* scan even when the result is not wanted */
	rc = gmp_sscanf(ptr, "%*Zd%n", &nread);
    }
    if (rc != 1) {
	fprintf(stderr, "parse error in parse_mpz while parsing:\n%s", line);
	exit(EXIT_FAILURE);
    }
    forbidden_garbage(ptr + nread, line);
    return 1;
}

void cado_poly_init(cado_poly poly)
{
    int i;

    memset(poly, 0, sizeof(poly));
    poly->f = (mpz_t *) malloc((MAXDEGREE + 1) * sizeof(mpz_t));
    poly->g = (mpz_t *) malloc((MAXDEGREE + 1) * sizeof(mpz_t));
    for (i = 0; i < (MAXDEGREE + 1); i++) {
	mpz_init_set_ui(poly->f[i], 0);
	mpz_init_set_ui(poly->g[i], 0);
    }
    mpz_init(poly->n);
    poly->name[0] = '\0';
    poly->degree = -1;
    poly->type[0] = '\0';
    mpz_init_set_ui(poly->m, 0);
}

void cado_poly_clear(cado_poly poly)
{
    int i;

    for (i = 0; i < (MAXDEGREE + 1); i++) {
	mpz_clear(poly->f[i]);
	mpz_clear(poly->g[i]);
    }
    mpz_clear(poly->n);
    free(poly->f);
    free(poly->g);
    mpz_clear(poly->m);
    memset(poly, 0, sizeof(poly));
}

// returns 0 on failure, 1 on success.
int read_polynomial(cado_poly poly, char *filename)
{
    FILE *file;
    const int linelen = 512;
    char line[linelen];
    int have_name = 0, have_n = 0;
    int i;

    file = fopen(filename, "r");
    if (file == NULL) {
	fprintf(stderr, "read_polynomial: could not open %s\n", filename);
	return 0;
    }

    /* We used to init the cado structure here, although it's only
     * vaguely the place to do so.
     */
    cado_poly_init(poly);

    int have_f[(MAXDEGREE + 1)] = { 0, };
    int have_g[(MAXDEGREE + 1)] = { 0, };
    int degf, degg;

    while (!feof(file)) {
	int ok = 0;
	if (fgets(line, linelen, file) == NULL)
	    break;
	if (line[0] == '#')
	    continue;

	ok += parse_string(poly->name, sizeof(poly->name),
			   line, "name", &have_name);
	ok += parse_mpz(poly->n, line, "n", &have_n);
	ok += parse_double(&(poly->skew), line, "skew", NULL);
	for (i = 0; i < (MAXDEGREE + 1); i++) {
	    char tag[4];
	    snprintf(tag, sizeof(tag), "c%d", i);
	    ok += parse_mpz(poly->f[i], line, tag, &have_f[i]);
	    snprintf(tag, sizeof(tag), "Y%d", i);
	    ok += parse_mpz(poly->g[i], line, tag, &have_g[i]);
	}
	ok +=
	    parse_string(poly->type, sizeof(poly->type), line, "type", NULL);
	ok += parse_ulong(&(poly->rlim), line, "rlim", NULL);
	ok += parse_ulong(&(poly->alim), line, "alim", NULL);
	ok += parse_int(&(poly->lpbr), line, "lpbr", NULL);
	ok += parse_int(&(poly->lpba), line, "lpba", NULL);
	ok += parse_int(&(poly->mfbr), line, "mfbr", NULL);
	ok += parse_int(&(poly->mfba), line, "mfba", NULL);
	ok += parse_double(&(poly->rlambda), line, "rlambda", NULL);
	ok += parse_double(&(poly->alambda), line, "alambda", NULL);
	ok += parse_int(&(poly->qintsize), line, "qintsize", NULL);
	ok += parse_mpz(poly->m, line, "m", NULL);

	ASSERT_ALWAYS(ok < 2);
	if (ok == 0) {
	    fprintf(stderr,
		    "read_polynomial: Cannot parse line %s\nIgnoring.\n",
		    line);
	    continue;
	}
    }

    for (degf = MAXDEGREE; degf >= 0 && !have_f[degf]; degf--) {
	if (have_f[degf] && mpz_cmp_ui(poly->f[degf], 0) != 0)
	    break;
    }
    for (degg = MAXDEGREE; degg >= 0 && !have_g[degg]; degg--) {
	if (have_g[degg] && mpz_cmp_ui(poly->g[degg], 0) != 0)
	    break;
    }
    poly->degree = degf;
    // compute m, the common root of f and g mod n
    if (degg != -1) {
	mpz_t tmp;
	mpz_init(tmp);
	mpz_invert(tmp, poly->g[1], poly->n);
	mpz_mul(tmp, tmp, poly->g[0]);
	mpz_mod(tmp, tmp, poly->n);
	mpz_sub(tmp, poly->n, tmp);
	mpz_mod(tmp, tmp, poly->n);
	if (mpz_cmp_ui(poly->m, 0) != 0) {
	    if (mpz_cmp(poly->m, tmp) != 0) {
		fprintf(stderr, "m is not a root of g mod N\n");
		exit(EXIT_FAILURE);
	    }
	}
	mpz_set(poly->m, tmp);
	mpz_clear(tmp);
    } else {
	ASSERT_ALWAYS(mpz_cmp_ui(poly->m, 0) != 0);
	mpz_set_ui(poly->g[1], 1);
	mpz_neg(poly->g[0], poly->m);
    }

    ASSERT_ALWAYS(have_n);

    fclose(file);

    return 1;
}

void clear_polynomial(cado_poly poly)
{
    cado_poly_clear(poly);
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

/* check that n divides b^d*f(m/d), where g = b*x-m */
void check_polynomials(cado_poly cpoly)
{
    mpz_t r, q;
    int k;
    int d = cpoly->degree;

    mpz_init_set(r, cpoly->f[d]);
    mpz_init_set_ui(q, 1);

    for (k = d - 1; k >= 0; k--) {
	mpz_mul(q, q, cpoly->g[1]);
	/* invariant: q = b^(d-k) */
	mpz_mul(r, r, cpoly->g[0]);
	mpz_neg(r, r);
	mpz_addmul(r, q, cpoly->f[k]);
    }

    if (mpz_divisible_p(r, cpoly->n) == 0) {
	fprintf(stderr, "Error, n does does divide Res(f,g)\n");
	exit(EXIT_FAILURE);
    }

    mpz_clear(r);
    mpz_clear(q);
}

#undef PARSE_MATCH
#undef PARSE_ERROR
#undef PARSE_NOMATCH
