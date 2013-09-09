#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>
#include "utils.h"
#include "mpz_poly.h"
#include "portability.h"

int
check_relation (relation_t *rel, cado_poly_ptr cpoly)
{
    int ok = 1;
    for(int side = 0 ; ok && side < 2 ; side++) {
        cado_poly_side_ptr ps = cpoly->pols[side];
        mpz_t no, acc;
        mpz_init (no);
        mpz_init_set_ui(acc, 1);
        mp_poly_homogeneous_eval_siui(no, ps->f, ps->degree, rel->a, rel->b);

        int n = side == RATIONAL_SIDE ? rel->nb_rp : rel->nb_ap;
        for (int i = 0; i < n; ++i)
        {
            unsigned long p = side == RATIONAL_SIDE ? rel->rp[i].p : rel->ap[i].p;
            int e = side == RATIONAL_SIDE ? rel->rp[i].e : rel->ap[i].e;
            for (int j = 0; j < e; ++j)
                mpz_mul_ui (acc, acc, p);
        }
        if (mpz_cmp (acc, no) != 0) {
            ok = 0;
            if (mpz_divisible_p (no, acc)) {
                mpz_divexact (acc, no, acc);
                gmp_fprintf (stderr, "Missing factor %Zd on %s side for (%ld, %lu)\n", acc, sidenames[side], rel->a, rel->b);
            } else {
                mpz_t g;
                mpz_init (g);
                mpz_gcd (g, acc, no);
                mpz_divexact (acc, acc, g);
                fprintf (stderr, "Wrong %s side for (%" PRId64 ", %" PRIu64 ")\n", sidenames[side], rel->a, rel->b);
                gmp_fprintf (stderr, "Given factor %Zd does not divide norm\n", acc);
                mpz_clear (g);
            }
        }

        mpz_clear (no);
        mpz_clear (acc);
    }
    return ok;
}

void usage_and_die(char *str) {
    fprintf(stderr, "usage: %s [-q] -poly <polyfile> [-f] <relfile1> <relfile2> ...\n",
            str);
    exit(3);
}

int
check_relation_files (char ** files, cado_poly_ptr cpoly, int forced_read,
                      int quiet)
{
    relation_stream rs;
    relation_stream_init(rs);
    unsigned long ok = 0;
    unsigned long bad = 0;
    int had_error = 0;
    unsigned long total = 0, nbfiles = 0;

    for( ; *files ; files++) {
        relation_stream_openfile(rs, *files);
        char line[RELATION_MAX_BYTES];
        unsigned long l0 = rs->lnum;
        unsigned long ok0 = ok;
        unsigned long bad0 = bad;
        int nread;
        nbfiles ++;
        for( ; (nread = relation_stream_get (rs, line, sizeof(line), forced_read, 10, 0)) >= 0 ; )
        {
            unsigned long l = rs->lnum - l0;
            if (nread > 0 && check_relation (&rs->rel, cpoly))
              {
                if (quiet == 0)
                  printf ("%s", line);
                ok++;
                continue;
              }
            fprintf (stderr, "Failed at line %lu in %s: %s",
                    l, *files, line);
            bad++;
            had_error = 1;
        }
        if (bad == bad0) {
            fprintf (stderr, "%s : %lu ok\n", *files, ok-ok0);
            total += ok - ok0;
        } else {
            fprintf(stderr, "%s : %lu ok ; FOUND %lu ERRORS\n",
                    *files, ok-ok0, bad-bad0);
        }
        relation_stream_closefile(rs);
    }
    relation_stream_clear(rs);

    if (nbfiles > 1)
      fprintf (stderr, "Total %lu relations ok\n", total);

    return had_error ? -1 : 0;
}

int main(int argc, char * argv[])
{
    int forced_read = 0;
    cado_poly cpoly;
    int had_error = 0;
    int quiet = 0;

    if (argc >= 2 && strcmp (argv[1], "-q") == 0)
      {
        quiet = 1;
        argc --;
        argv ++;
      }

    if (argc < 4 || strcmp(argv[1], "-poly") != 0) {
        usage_and_die(argv[0]);
    }
    cado_poly_init(cpoly);
    if (!cado_poly_read(cpoly, argv[2])) 
        return 2;

    if (strcmp(argv[3], "-f") == 0) {
        forced_read=1;
        argv++,argc--;
    }

    argv++,argc--;
    argv++,argc--;
    argv++,argc--;

    had_error = check_relation_files (argv, cpoly, forced_read, quiet) < 0;

    cado_poly_clear(cpoly);

    return had_error;
}
