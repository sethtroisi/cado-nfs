#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>


#include "macros.h"
#include "fb.h"

// return 0 on success (!) and EOF otherwise.
int advance_to_non_space(FILE *file)
{
    int c;
    do {
        c = fgetc(file);
        if (c == EOF) {
            ungetc(c, file);
            return EOF;
        }
    } while (isspace(c));
    c = ungetc(c, file);
    if (c == EOF)
        return EOF;
    else 
        return 0;
}




int fbread(factorbase_t * FB, const char *filename, int maxdeg)
{
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("when reading factor base ");
        return 0;
    }
    int n = 0;
    fbideal_t * elts = NULL;
    while (1) {
        fbprime_t p, r;
        int err;
        int c;
        // move to first non-space char
        err = advance_to_non_space(file);
        if (err == EOF)
            break;

        // get p
        err = fbprime_inp(p, file);
        if (!err) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %d-th ideal.\n", n);
            fclose(file);
            return 0;
        }
        
        // remove whites
        err = advance_to_non_space(file);
        if (err == EOF) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %d-th ideal.\n", n);
            fclose(file);
            return 0;
        }
        // read ":"
        c = fgetc(file);
        if (c != ':') {
            fprintf(stderr, "Error parsing ':' in factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %d-th ideal.\n", n);
            fclose(file);
            return 0;
        }

        // get r
        err = fbprime_inp(r, file);
        if (!err) {
            fprintf(stderr, "Error parsing factor base %s.\n", filename);
            fprintf(stderr, "  The error occured at the %d-th ideal.\n", n);
            fclose(file);
            return 0;
        }
        unsigned char degp = fbprime_deg(p);
        if (maxdeg!=0 && degp > maxdeg) 
            break;
        n++;
        elts = (fbideal_t *) realloc(elts, n*sizeof(fbideal_t));
        ASSERT_ALWAYS(elts != NULL);
        elts[n-1].degp = degp;
        fbprime_set(elts[n-1].p, p);
        fbprime_set(elts[n-1].r, r);
    }

    FB->size = n;
    FB->elts = elts;
    fclose(file);
    return 1;
}
