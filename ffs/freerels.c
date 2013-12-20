#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "timing.h"
#include "params.h"
#include "smoothness.h"
#include "polyfactor.h"
#include "fq.h"
#include "fqpol.h"

#include "renumber.h"

#define MAX_FFS_DEG 20

int sq_is_split(sq_t * roots, sq_srcptr q, ffspol_srcptr F) {
    fq_info_t Fq;
    fq_info_init(Fq, q);

    fqpol_t f;
    fqpol_init(f);
    fqpol_set_ffspol(f, F, Fq);

    int ret = fqpol_is_split(roots, f, Fq);
    fqpol_clear(f);
    fq_info_clear(Fq);
    return ret;
}


int sq_roots(unsigned long *roots, sq_srcptr q, ffspol_srcptr F)
{
  fq_info_t Fq;
  fq_info_init(Fq, q);
  fq_t r[MAX_FFS_DEG];

  fqpol_t f;
  fqpol_init(f);
  fqpol_set_ffspol(f, F, Fq);

  int nb = fqpol_roots(r, f, Fq);
  fqpol_clear(f);
  fq_info_clear(Fq);

  for (int i = 0; i < nb; i++)
    roots[i] = fppol64_get_ui_sparse (r[i]);
  if (f->deg != F->deg) // There is a projective root
  {
    roots[nb] = fppol64_get_ui_sparse (q);
    nb++;
  }
  return nb;
}

int sq_is_irreducible(sq_srcptr p) {
    fppol_t P;
    fppol_init(P);
    fppol_set_sq(P, p);
    int ret = fppol_is_irreducible(P);
    fppol_clear(P);
    return ret;
}

#define NB_POL_PER_THREAD 2048

typedef unsigned long tabroots_t[MAX_FFS_DEG];

struct thread_info
{
  sq_t *p;
  unsigned int nbp;
  int lpb[2];
  ffspol_srcptr ffspol[2];

  tabroots_t *roots0;
  tabroots_t *roots1;
  int *k0;
  int *k1;
};


void * thread_start (void *arg)
{
  struct thread_info *ti = (struct thread_info *) arg;
  unsigned int i;

  for (i = 0; i < ti->nbp; i++)
  {
    if (sq_deg(ti->p[i]) <= ti->lpb[0])
      ti->k0[i] = sq_roots(ti->roots0[i], ti->p[i], ti->ffspol[0]);
    else
      ti->k0[i] = 0;
    if (sq_deg(ti->p[i]) <= ti->lpb[1])
      ti->k1[i] = sq_roots(ti->roots1[i], ti->p[i], ti->ffspol[1]);
    else
      ti->k1[i] = 0;
  }

  return NULL;
}

uint64_t freerels_mt (int nt, ffspol_t pol[2], int lpb[2], int min_deg,
                      int max_deg, renumber_t renumtab)
{
  // Declare tables containing input and output
  int **tabk0;
  int **tabk1;
  tabroots_t **tabr0;
  tabroots_t **tabr1;
  sq_t **tabp;

  //Allocate memory
  tabk0 = (int **) malloc ( nt * sizeof (int *));
  tabk1 = (int **) malloc ( nt * sizeof (int *));
  tabr0 = (tabroots_t **) malloc ( nt * sizeof (tabroots_t *));
  tabr1 = (tabroots_t **) malloc ( nt * sizeof (tabroots_t *));
  tabp = (sq_t **) malloc ( nt * sizeof (sq_t *));

  size_t ksize =  NB_POL_PER_THREAD * sizeof (int);
  size_t rsize =  NB_POL_PER_THREAD * sizeof (tabroots_t);
  size_t sqsize = NB_POL_PER_THREAD * sizeof (sq_t);
  for (int i = 0; i < nt; ++i)
  {
    tabk0[i] = (int *) malloc ( ksize);
    tabk1[i] = (int *) malloc ( ksize);
    tabr0[i] = (tabroots_t *) malloc ( rsize);
    tabr1[i] = (tabroots_t *) malloc ( rsize);
    tabp[i] = (sq_t *) malloc (sqsize);
  }

  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc(nt*sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.

  // Arguments for threads
  struct thread_info *tis;
  tis = (struct thread_info*) malloc(nt*sizeof(struct thread_info));
  for (int i = 0; i < nt; ++i)
  {
    tis[i].lpb[0] = lpb[0];
    tis[i].lpb[1] = lpb[1];
    tis[i].ffspol[0] = pol[0];
    tis[i].ffspol[1] = pol[1];

    tis[i].k0 = tabk0[i];
    tis[i].k1 = tabk1[i];
    tis[i].roots0 = tabr0[i];
    tis[i].roots1 = tabr1[i];
    tis[i].p = tabp[i];
  }

  // Prepare the main loop
  uint64_t nfreerels = 0;
  sq_t p;
  sq_set_ti(p, 1);   // start with the polynomial t, always irreducible

  // Main loop
  while ((sq_deg(p) <= max_deg) || (active_threads > 0))
  {
    // Start / restart threads as many threads as allowed
    if ((active_threads < nt) && (sq_deg(p) <= max_deg))
    {
      unsigned int nbp = 0;
      while (sq_deg(p) <= max_deg && nbp < NB_POL_PER_THREAD)
      {
        sq_set(tis[threads_head].p[nbp], p);
        nbp++;
        do
        {
          sq_monic_set_next(p, p, 64);
        } while (!sq_is_irreducible(p));
      }
      tis[threads_head].nbp = nbp;

      pthread_create(&threads[threads_head], NULL,
                     &thread_start, (void *)(&tis[threads_head]));
      active_threads++;
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
      continue;
    }

    // Wait for the next thread to finish in order to process results
    pthread_join(threads[threads_head], NULL);
    active_threads--;
    for (unsigned long j = 0; j < tis[threads_head].nbp; j++)
    {
      unsigned long p_int = fppol64_get_ui_sparse(tabp[threads_head][j]);
      int k[2] = { tabk0[threads_head][j], tabk1[threads_head][j]};
      unsigned long *roots[2] = {tabr0[threads_head][j], tabr1[threads_head][j]};
      uint64_t l = renumtab->size;

      renumber_write_p (renumtab, p_int, roots, k);

      if (sq_deg(tabp[threads_head][j]) <= min_deg && k[0] == pol[0]->deg
                                                   && k[1] == pol[1]->deg)
      {
        printf ("%lx,0:%lx", p_int, (unsigned long) l);
        for (l = l + 1; l < renumtab->size; l++)
          printf (",%lx", (unsigned long) l);
        printf ("\n");
        nfreerels++;
      }
    }

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (sq_deg(p) > max_deg)
    {
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
    }
  }

  free(tis);
  free(threads);
  for (int i = 0; i < nt; ++i)
  {
    free(tabk0[i]);
    free(tabk1[i]);
    free(tabr0[i]);
    free(tabr1[i]);
    free(tabp[i]);
  }
  free(tabk0);
  free(tabk1);
  free(tabr0);
  free(tabr1);
  free(tabp);

  return nfreerels;
}

uint64_t freerels_mono (ffspol_t pol[2], int lpb[2], int min_deg, int max_deg,
                        renumber_t renumtab)
{
  uint64_t nrel = 0;
  sq_t p;
  sq_set_ti(p, 1);   // start with the polynomial t, always irreducible
  unsigned long *roots[2];
  int k[2], d[2];
  uint64_t old_table_size = renumtab->size;

  d[0] = pol[0]->deg;
  d[1] = pol[1]->deg;
  roots[0] = (unsigned long*) malloc (d[0] * sizeof (unsigned long));
  roots[1] = (unsigned long*) malloc (d[1] * sizeof (unsigned long));
  fprintf(stderr, "lpb[0]=%d lpb[1]=%d\n", lpb[0], lpb[1]);

  while (sq_deg(p) <= max_deg)
  {
    /* first compute the roots */
    for (int i = 0; i < 2; i++)
    {
      if (sq_deg(p) <= lpb[i])
        k[i] = sq_roots(roots[i], p, pol[i]);
      else
        k[i] = 0;
    }

    unsigned long p_int = fppol64_get_ui_sparse (p);
    renumber_write_p (renumtab, p_int, roots, k);

    if (sq_deg(p) <= min_deg && k[0] == d[0] && k[1] == d[1])
    {
      //print the free rels
      uint64_t l;
      printf ("%lx,0:%lx", p_int, (unsigned long) old_table_size);
      for (l = old_table_size + 1; l < renumtab->size; l++)
        printf (",%lx", (unsigned long) l);
      printf ("\n");
      nrel++;
    }
    old_table_size = renumtab->size;

    do
    {
      sq_monic_set_next(p, p, 64);
    } while (!sq_is_irreducible(p));
  }

  free(roots[0]);
  free(roots[1]);
  return nrel;
}

void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0 *           function field polynomial on side 0\n");
    fprintf(stderr, "  pol1 *           function field polynomial on side 1\n");
    fprintf(stderr, "  lpb0 *           large prime bound on side 0\n");
    fprintf(stderr, "  lpb1 *           large prime bound on side 1\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");
    fprintf(stderr, "  -mt n            number of threads (default 1)\n");

    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

int main(int argc, char **argv)
{
    ffspol_t ffspol[2];
    int lpb[2] = {0, 0}, add_full_col = 0;
    char *argv0 = argv[0];
    int gf = 0;
    const char * renumberfilename = NULL;
    const char * badidealsfilename = NULL;
    int nthreads = 1;

    param_list pl;
    param_list_init(pl);
    param_list_configure_switch(pl, "addfullcol", &add_full_col);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }
    // Update parameter list at least once to register argc/argv pointers.
    param_list_update_cmdline(pl, &argc, &argv);
    param_list_print_command_line(stdout, pl);

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
    }

    // read function field polynomials
    {
        const char * polstr;
        ffspol_init(ffspol[0]);
        ffspol_init(ffspol[1]);
        polstr = param_list_lookup_string(pl, "pol0");
        if (polstr == NULL) usage(argv0, "pol0");
        ffspol_set_str(ffspol[0], polstr);
        polstr = param_list_lookup_string(pl, "pol1");
        if (polstr == NULL) usage(argv0, "pol1");
        ffspol_set_str(ffspol[1], polstr);
    }
    // read various bounds
    param_list_parse_int(pl, "lpb0", &lpb[0]);
    param_list_parse_int(pl, "lpb1", &lpb[1]);
    if (lpb[0] == 0) usage(argv0, "lpb0");
    if (lpb[1] == 0) usage(argv0, "lpb1");

    //read filename
    renumberfilename = param_list_lookup_string(pl, "renumber");
    if (renumberfilename == NULL)
      usage (argv0, "renumber");
    badidealsfilename = param_list_lookup_string(pl, "badideals");

    param_list_parse_int(pl, "mt", &nthreads);


    cado_poly dummy_poly;
    renumber_t tab;
    cado_poly_init(dummy_poly);

    dummy_poly->pols[0]->deg = ffspol[0]->deg;
    dummy_poly->pols[1]->deg = ffspol[1]->deg;
    unsigned long dummy_lpb[2] = { __FP_BITS + __FP_BITS * lpb[0],
                                  __FP_BITS + __FP_BITS * lpb[1]};

    if (dummy_poly->pols[1]->deg == 1)
    {
      dummy_poly->rat  = dummy_poly->pols[1];
      dummy_poly->alg  = dummy_poly->pols[0];
    }
    else if (dummy_poly->pols[0]->deg == 1)
    {
      dummy_poly->rat  = dummy_poly->pols[0];
      dummy_poly->alg  = dummy_poly->pols[1];
    }

    renumber_init(tab, dummy_poly, dummy_lpb);
    renumber_init_write (tab, renumberfilename, badidealsfilename, add_full_col);

    int max_deg = MAX(lpb[0], lpb[1]);
    int min_deg = MIN(lpb[0], lpb[1]);
    fprintf (stderr, "Generating freerels up to degree %d\n", min_deg);
    fprintf (stderr, "Generating renumber up to degree %d\n", max_deg);

    uint64_t nrel;
    if (nthreads > 1)
      nrel = freerels_mt (nthreads, ffspol, lpb, min_deg, max_deg, tab);
    else
      nrel = freerels_mono (ffspol, lpb, min_deg, max_deg, tab);

    fprintf(stdin, "# Computed %" PRIu64 " free relations\n", nrel);
    fprintf(stderr, "# Computed %" PRIu64 " free relations\n", nrel);
    renumber_close_write(tab, renumberfilename);
    param_list_clear(pl);
    cado_poly_clear(dummy_poly);
    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);

    return EXIT_SUCCESS;
}
