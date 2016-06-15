#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <math.h> /* for sqrt and floor and log and ceil */
#include <pthread.h>
#include "portability.h"
#include "utils.h"
#include <vector>
#include <algorithm>

using namespace std;

// Structure that contains indices (in the sense of renumber.[ch]) for
// one side, with the possibility to get a randomized ideal on the same
// side, and more or less of the same size.
// Usage:
//  - init()
//  - do many append() (increasing order of indices)
//  - finalize()
//  - can call random_index()

struct indexrange {
    vector <index_t>ind;

    void init() { }
    void finalize() { }

    void append(index_t z) {
        ind.push_back(z);
    }

    index_t random_index(index_t z) {
        // Find position of z
        vector<index_t>::iterator it;
        it = lower_bound (ind.begin(), ind.end(), z);
        uint64_t position = it - ind.begin();
        uint64_t range = uint64_t((double(position))*0.2);
        uint64_t low = MAX((int64_t)(position - range), 0);
        uint64_t high = MIN(position + range, ind.size());
        if (high == low)
            high++;
        return ind[low + uint64_t(random()%(high-low))];
    }
};

// Fill in the indexrange data structure from the renumber table,
// gathering indices in two arrays, one for each side.
void prepare_indexrange(indexrange *Ind, renumber_t ren_tab, 
        cado_poly cpoly) {
    Ind[0].init();
    Ind[1].init();
    for (index_t i = 0; i < ren_tab->size; i++) {
        int side = renumber_get_side_from_index(ren_tab, i, cpoly);
        Ind[side].append(i);
    }
    Ind[0].finalize();
    Ind[1].finalize();
}

#define MAXFACTORS 30  // on each side. Should be enough?

struct fake_rel {
    index_t ind[2][MAXFACTORS];
    int nb_ind[2];
};

void read_rel(fake_rel& rel, uint64_t q, int sqside, const char *str,
        renumber_t ren_tab) {
    rel.nb_ind[0] = 0;
    rel.nb_ind[1] = 0;
    int side = 0;
    // read a and b
    int64_t a;
    uint64_t b;
    const char *pstr = str;
    {
        char *endpstr = NULL;
        a = strtoll(pstr, &endpstr, 10);
        assert (endpstr != pstr);
        pstr = endpstr;
        assert (pstr[0]==',');
        pstr++;
        b = strtoull(pstr, &endpstr, 10);
        assert (endpstr != pstr);
        pstr = endpstr;
        // skip ':'
        assert (pstr[0]==':');
        pstr++;
    }
    while (pstr[0] != '\0' && pstr[0] != '\n') {
        if (pstr[0] == ':') {
            // change side
            side++;
            pstr++;
            continue;
        }
        if (pstr[0] == ',') {
            pstr++;
            continue;
        }
        char *endpstr = NULL;
        uint64_t p = strtoull(pstr, &endpstr, 16);
        assert (endpstr != pstr);
        if (side != sqside || q != p) {
            p_r_values_t r = relation_compute_r(a, b, p);
            index_t index;
            int nb;
            if (renumber_is_bad (&nb, &index, ren_tab, p, r, side)) {
                // bad ideal: just pick a random ideal above this prime
                index += (random() % nb);
            } else {
                index = renumber_get_index_from_p_r(ren_tab, p, r, side);
            }
            rel.ind[side][rel.nb_ind[side]] = index;
            rel.nb_ind[side]++;
            assert (rel.nb_ind[side] <= MAXFACTORS);
        }
        pstr = endpstr;
    }
}

void read_sample_file(vector<unsigned int> &nrels, vector<fake_rel> &rels,
        int sqside, const char *filename, renumber_t ren_tab)
{
    FILE * file;
    file = fopen(filename, "r");
    assert (file != NULL);
    uint64_t q = 0;
    unsigned int nr = 0;

    char line[1024];

    do {
        char * ret = fgets(line, 1024, file);
        if (ret == NULL) {
            nrels.push_back(nr);
            break;
        }
        if (line[0] == '#') {
            char *ptr;
            ptr = strstr(line, "# Sieving side-");
            if (ptr != NULL) {
                // starting a new special-q
                assert (sqside == (ptr[15] - '0'));
                ptr += 19; // skip "# Sieving side-0 q="
                q = strtoull(ptr, NULL, 10);
                nrels.push_back(nr);
                nr = 0;
            }
        } else {
            fake_rel rel;
            read_rel(rel, q, sqside, line, ren_tab);
            rels.push_back(rel);
            nr++;
        }
    } while(1);

    fclose(file);
}


void print_random_fake_rel(renumber_ptr ren_info, uint64_t q,
        uint64_t r, int sqside, vector<fake_rel> rels, indexrange *Ind)
{
    index_t indq = renumber_get_index_from_p_r(ren_info, q, r, sqside);
    int i = random() % rels.size();
    printf("%ld,%lu:", (long)random(), (unsigned long)random());
    vector <index_t> frel;
    frel.push_back(indq);
    for (int side = 0; side < 2; ++side) {
        int np = rels[i].nb_ind[side];
        for (int j = 0; j < np; ++j) {
            index_t ind = Ind[side].random_index(rels[i].ind[side][j]);
            frel.push_back(ind);
        }
    }
    sort(frel.begin(), frel.end());
    for (unsigned int i = 0; i < frel.size(); ++i) {
        printf("%" PRid, frel[i]);
        if (i != frel.size()-1) 
            printf(",");
    }
    printf("\n");
}


void advance_prime_in_fb(int *mult, uint64_t *q, uint64_t *roots,
        cado_poly cpoly, int sqside, prime_info pdata)
{
    int nr;
    unsigned long newp;
    do {
        newp = getprime_mt (pdata);
        assert (newp > *q);
        nr = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], newp);
    } while (nr == 0);

    *mult = nr;
    *q = newp;
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "polynomial file");
    param_list_decl_usage(pl, "lpb0", "factor base bound on side 0");
    param_list_decl_usage(pl, "lpb1", "factor base bound on side 1");
    param_list_decl_usage(pl, "q0", "lower bound of the qrange");
    param_list_decl_usage(pl, "q1", "upper bound of the qrange");
    param_list_decl_usage(pl, "sqside", "side of the special-q");
    param_list_decl_usage(pl, "sample", "file where to find a sample of relations");
    param_list_decl_usage(pl, "renumber", "renumber table");
    verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly cpoly;
  int sqside = -1;
  char *argv0 = argv[0];
  int lpb[2] = {0, 0};
  uint64_t q0 = 0;
  uint64_t q1 = 0;

  param_list_init(pl);
  declare_usage(pl);
  cado_poly_init(cpoly);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      FILE *f;
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f, 0);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }
  verbose_interpret_parameters(pl);
  param_list_print_command_line(stdout, pl);

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "lpb0", &lpb[0]);
  if (lpb[0] == 0) {
      fprintf(stderr, "Error: parameter -lpb0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  param_list_parse_int(pl, "lpb1", &lpb[1]);
  if (lpb[1] == 0) {
      fprintf(stderr, "Error: parameter -lpb1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q0", &q0);
  if (q0 == 0) {
      fprintf(stderr, "Error: parameter -q0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q1", &q1);
  if (q1 == 0) {
      fprintf(stderr, "Error: parameter -q1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  if (!cado_poly_read(cpoly, filename))
    {
      fprintf (stderr, "Error reading polynomial file %s\n", filename);
      exit (EXIT_FAILURE);
    }

  param_list_parse_int(pl, "sqside", &sqside);
  if (sqside == -1 || sqside > 2) {
      fprintf(stderr, "Error: sqside must be 0 or 1\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  const char * renumberfile;
  if ((renumberfile = param_list_lookup_string(pl, "renumber")) == NULL) {
      fprintf(stderr, "Error: parameter -renumber is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  renumber_t ren_table;
  renumber_read_table(ren_table, renumberfile);

  // read sample file
  const char * sample;
  if ((sample = param_list_lookup_string(pl, "sample")) == NULL) {
      fprintf(stderr, "Error: parameter -sample is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  vector<fake_rel> rels;
  vector<unsigned int> nrels;
  read_sample_file(nrels, rels, sqside, sample, ren_table);

  param_list_warn_unused(pl);

  // Two index ranges, one for each side
  indexrange Ind[2];
  prepare_indexrange(Ind, ren_table, cpoly);

  // Loop on [q0-q1]
  prime_info pdata;
  prime_info_init(pdata);
  // fast forward until we reach q0
  uint64_t q = 2;
  while (q < q0) {
      q = getprime_mt(pdata);
  }
  uint64_t roots[MAXDEGREE];
  // loop until q1
  do {
      int mult = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], q); 
      for (int i = 0; i < mult; ++i) {
          int n = random()%nrels.size();
          for (unsigned int j = 0; j < nrels[n]; ++j)
              print_random_fake_rel(ren_table, q, roots[i], sqside, rels, Ind);
      }
      advance_prime_in_fb(&mult, &q, roots, cpoly, sqside, pdata);
  } while (q < q1);

  prime_info_clear(pdata);
  renumber_clear(ren_table);
  cado_poly_clear(cpoly);
  param_list_clear(pl);

  return 0;
}
