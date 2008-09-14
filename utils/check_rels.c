#define _POSIX_C_SOURCE 2  /* popen/pclose with -std=c99, -pedantic or -ansi */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>
#include "cado.h"
#include "utils.h"

void norm(mpz_t *f, int deg, mpz_t r, long a, unsigned long b)
{
    int i;

    mpz_t tmp;
    mpz_init_set_ui(tmp, 1);
    mpz_set(r, f[deg]);
    for (i = deg - 1; i >= 0; i--) {
        mpz_mul_si(r, r, a);
        mpz_mul_ui(tmp, tmp, b);
        mpz_addmul(r, tmp, f[i]);
    }
    mpz_abs(r, r);
    mpz_clear(tmp);
}


int is_gzip(const char * s)
{
    unsigned int l = strlen(s);
    return l >= 3 && strcmp(s + l - 3, ".gz") == 0;
}

int check_relation(relation_t *rel, __cado_poly_struct* cpoly) {
    mpz_t no, acc;
    int i;

    mpz_init(no);
    mpz_init(acc);

    // algebraic side
    norm(cpoly->f, cpoly->degree, no, rel->a, rel->b);
    mpz_set_ui(acc, 1);
    for(i = 0; i < rel->nb_ap; ++i) {
        int j;
        for (j = 0; j < (rel->ap[i]).e; ++j) 
            mpz_mul_ui(acc, acc, (rel->ap[i]).p);
    }
    if (mpz_cmp(acc, no) != 0) {
        fprintf(stderr,
                "Wrong algebraic side for (%ld, %lu)\n", rel->a, rel->b);
        mpz_clear(no);
        mpz_clear(acc);
        return 0;
    }

    // rational side
    norm(cpoly->g, 1, no, rel->a, rel->b);
    mpz_set_ui(acc, 1);
    for(i = 0; i < rel->nb_rp; ++i) {
        int j;
        for (j = 0; j < (rel->rp[i]).e; ++j) 
            mpz_mul_ui(acc, acc, (rel->rp[i]).p);
    }
    if (mpz_cmp(acc, no) != 0) {
        fprintf(stderr, "Wrong rational side for (%ld, %lu)\n", rel->a, rel->b);
        mpz_clear(no);
        mpz_clear(acc);
        return 0;
    }
    mpz_clear(no);
    mpz_clear(acc);
    return 1;
}

#define MAX_LENGTH 512

int check_stream(const char *name, FILE * stream, __cado_poly_struct* cpoly)
{
    int lnum;
    int nrels = 0;

    for (lnum = 1;; lnum++) {
        char line[MAX_LENGTH];
        if (fgets(line, sizeof(line), stream) == NULL)
            break;
        if (line[strlen(line)-1] != '\n')
          {
            fprintf (stderr, "Line %d of %s is too long, please increase MAX_LENGTH in check_rels.c\n", lnum, name);
            exit (1);
          }
        if (line[0] == '#')
            continue;

        relation_t rel;
        if (!read_relation(&rel, line)) {
            fprintf(stderr, "Failed on line %d in %s: %s\n", lnum, name, line);
            fprintf(stderr,
                    "File check aborted. Please fix this file and re-run.\n");
            fprintf(stderr, "Continuing with next file\n");
            return 1;
        }
        if (!check_relation(&rel, cpoly)) {
            fprintf(stderr, "Failed on line %d in %s: %s\n", lnum, name, line);
            fprintf(stderr,
                    "File check aborted. Please fix this file and re-run.\n");
            fprintf(stderr, "Continuing with next file\n");
            return 1;
        }

        clear_relation(&rel);
        nrels++;
    }
    fprintf(stderr, "%s : ok %d rels\n", name, nrels);
    return 0;
}


void usage_and_die(char *str) {
    fprintf(stderr, "usage: %s -poly <polyfile> <relfile1> <relfile2> ...\n",
            str);
    exit(1);
}

int main(int argc, char * argv[])
{
    cado_poly cpoly;
    int i, had_error = 0;

    if (argc < 4 || strcmp(argv[1], "-poly") != 0) {
        usage_and_die(argv[0]);
    }
    cado_poly_init(cpoly);
    if (!cado_poly_read(cpoly, argv[2])) 
        return 1;

    for(i = 3 ; i < argc ; i++) {
        FILE * f;
        if (is_gzip(argv[i])) {
            char command[1024];
            snprintf(command, sizeof(command), "gzip -dc %s", argv[i]);
            f = popen(command, "r");
            had_error |= check_stream(argv[i], f, cpoly);
            pclose(f);
        } else {
            f = fopen(argv[i], "r");
            had_error |= check_stream(argv[i], f, cpoly);
            fclose(f);
        }
    }

    cado_poly_clear(cpoly);

    return had_error;
}

