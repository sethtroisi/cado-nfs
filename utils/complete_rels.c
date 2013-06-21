/* complete_rels: same as check_rels, but completes small factors if omitted */
#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>
#include "portability.h"
#include "utils.h"
#include "mpz_poly.h"

int
complete_relation (relation_t *rel, cado_poly_ptr cpoly)
{
  for(int side = 0 ; side < 2 ; side++) {
      cado_poly_side_ptr ps = cpoly->pols[side];
      mpz_t no;
      mpz_init (no);
      mp_poly_homogeneous_eval_siui(no, ps->f, ps->degree, rel->a, rel->b);

      int n = side == RATIONAL_SIDE ? rel->nb_rp : rel->nb_ap;
      for (int i = 0; i < n; ++i)
      {
          unsigned long p = side == RATIONAL_SIDE ? rel->rp[i].p : rel->ap[i].p;
          int e = side == RATIONAL_SIDE ? rel->rp[i].e : rel->ap[i].e;
          for (int j = 0; j < e; ++j)
              if (mpz_fdiv_q_ui (no, no, p) != 0)
              {
                  fprintf (stderr,
                          "Wrong %s side for (%" PRId64 ", %" PRIu64 ")\n",
                          sidenames[side], rel->a, rel->b);
                  fprintf (stderr, "Given factor %lu does not divide norm\n", p);
                  mpz_clear (no);
                  return 0;
              }
      }

      for (unsigned long p = 2; mpz_cmp_ui (no, 1) != 0; p = getprime(p)) {
          while (mpz_divisible_ui_p (no, p))
          {
              relation_add_prime (rel, side, p);
              mpz_divexact_ui (no, no, p);
          }
      }
      getprime(0);
      mpz_clear(no);
  }

  // were printing our output, thus we don't care about tidyness of the
  // internal structure.
  // relation_compress_rat_primes(rel);
  // relation_compress_alg_primes(rel);

  return 1;
}

int complete_relation_files(char ** files, cado_poly_ptr cpoly)
{
    relation_stream rs;
    relation_stream_init(rs);
    unsigned long ok = 0;
    unsigned long bad = 0;
    int had_error = 0;
    for( ; *files ; files++) {
        relation_stream_openfile(rs, *files);
        char line[RELATION_MAX_BYTES];
        unsigned long l0 = rs->lnum;
        unsigned long ok0 = ok;
        unsigned long bad0 = bad;
        for( ; relation_stream_get(rs, line, 0, 10) >= 0 ; ) {
            unsigned long l = rs->lnum - l0;
            if (complete_relation(&rs->rel, cpoly)) {
                fprint_relation (stdout, &rs->rel);
                ok++;
                continue;
            }
            fprintf(stderr, "Failed at line %lu in %s: %s",
                    l, *files, line);
            bad++;
            had_error = 1;
        }
        if (bad == bad0) {
            fprintf(stderr, "%s : %lu ok\n", *files, ok-ok0);
        } else {
            fprintf(stderr, "%s : %lu ok ; FOUND %lu ERRORS\n",
                    *files, ok-ok0, bad-bad0);
        }
        relation_stream_closefile(rs);
    }
    relation_stream_clear(rs);

    return had_error ? -1 : 0;
}

void usage_and_die(char *str) {
    fprintf(stderr, "usage: %s [-f] -poly <polyfile> <relfile1> <relfile2> ...\n",
            str);
    exit(3);
}

int main(int argc, char * argv[])
{
    cado_poly cpoly;
    int had_error = 0;
    int forced_read = 0;

    if (argc >= 2 && strcmp (argv[1], "-f") == 0)
      {
        forced_read = 1;
        fprintf (stderr, "Skipping corrupted lines\n");
        argv ++;
        argc --;
      }

    if (argc < 4 || strcmp(argv[1], "-poly") != 0) {
        usage_and_die(argv[-forced_read]);
    }
    cado_poly_init(cpoly);
    if (!cado_poly_read(cpoly, argv[2])) 
        return 2;

    argv++,argc--;
    argv++,argc--;
    argv++,argc--;

    had_error = complete_relation_files(argv, cpoly) < 0;

    cado_poly_clear(cpoly);

    return had_error;
}
