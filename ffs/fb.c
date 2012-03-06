#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "macros.h"
#include "fb.h"



// Return the last character read, or EOF.
static
int skip_spaces(FILE *file)
{
    int c;
    while (isspace(c = getc(file)));
    return ungetc(c, file);
}


// Initialize a factor base, reading the ideals from a file and computing
// the corresponding lambda using the basis of the given q-lattice.
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename)
{
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error reading factor base");
        return 0;
    }

    FB->n    = 0;
    FB->elts = NULL;
    unsigned alloc = 0;
    for (fbideal_t *elt = FB->elts; skip_spaces(file) != EOF; ++FB->n, ++elt) {
        // Need realloc?
        if (alloc <= FB->n) {
            alloc *= 2;
            FB->elts = (fbideal_t *)realloc(FB->elts, alloc*sizeof(fbideal_t));
            ASSERT_ALWAYS(FB->elts != NULL);
            if (elt == NULL) elt = FB->elts;
        }

        // get p
        if (!fbprime_inp(elt->p, file)) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
            fclose(file);
            return 0;
        }
        elt->degp = fbprime_deg(elt->p);

        // remove spaces
        if (skip_spaces(file) == EOF) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
            fclose(file);
            return 0;
        }

        // read ":"
        if (getc(file) != ':') {
            fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
            fclose(file);
            return 0;
        }

        // get r
        if (!fbprime_inp(elt->r, file)) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %u-th ideal.\n", FB->n);
            fclose(file);
            return 0;
        }
    }

    FB->elts = (fbideal_t *)realloc(FB->elts, FB->n*sizeof(fbideal_t));
    ASSERT_ALWAYS(FB->elts != NULL);
    fclose(file);
    return 1;
}


// Precompute lambda for each element of the factor base.
void factor_base_precomp_lambda(factor_base_ptr FB, qlat_srcptr qlat)
{
    for (fbideal_t *elt = FB->elts; elt != FB->elts+FB->n; ++elt) {
        fbprime_t t0, t1;
        fbprime_mulmod(t0, qlat->b0, elt->r, elt->p);
        fbprime_sub   (t0, qlat->a0, t0);
        fbprime_rem   (t0, t0, elt->p);
        if ((elt->proj = fbprime_is_zero(t0))) {
            // This is a projective root! Yurk!
            // TODO: We'll have to handle these, someday.
            continue;
        }
        fbprime_invmod(t0, t0, elt->p);
        fbprime_mulmod(t1, qlat->b1, elt->r, elt->p);
        fbprime_sub   (t1, t1, qlat->a1);
        fbprime_rem   (t1, t1, elt->p);
        fbprime_mulmod(elt->lambda, t0, t1, elt->p);
    }
}


// Clean up memory.
void factor_base_clear(factor_base_ptr FB)
{
    free(FB->elts);
}
