/* dup2: 2nd pass

   Usage: dup2 -poly name.poly [-basepath <dir>] [-filelist <fl>]
          -K <K> -renumber xxx file1 ... filen

   Put non-duplicate files (among file1 ... filen only) into:
   <dir>/file1 ... <dir>/filen

   Input files can be given on command line, or via a filelist file.

   In case a filelist is given, the -basepath option enables to tell in which
   directory those files are. Note that input files will be replaced in-place.

   By default, the output will be bzipped/gzipped according to the status
   of the input file.

   Allocates a hash-table of 2^k 32-bit entries.

   The relations are renumbered according to the file xxx, as input we have:

       a,b:p1,p2,...,pj:q1,q2,...,qk                                 (*)

   where p1,p2,...,pj are rational ideals (possibly duplicate), and
   q1,q2,...,qk are algebraic ideals (possibly duplicate). Output is:

       a,b:r1,r2,...,rm                                              (**)

   where each index r1,r2,...,rm appears only once, and refers to either
   a rational or to an algebraic ideal.

   The format of each file is recognized by counting the number of ':' in the
   first line: if two we have the raw format (*), if only one we have the
   renumbered format (**). It is assumed that all files in renumbered format
   come first.

   Algorithm: for each (a,b) pair, we compute h(a,b) = (CA*a+CB*b) % 2^64.

   h has 64 bits. We know bits 2..6 are identical for a given slice, thus
   we remove them, and we obtain a 59-bit value h'.

   Let h' = j * 2^k + i.

   We store j % 2^32 at the next empty cell after index i in the hash-table.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"

#ifdef FOR_FFS
#include "fppol.h"
#include "fq.h"
#include "utils_ffs.h"
#endif
#include "filter_badideals.h"

#define DEBUG 0

char *argv0; /* = argv[0] */

/* Renumbering table to convert from (p,r) to an index */
renumber_t renumber_tab;

static uint32_t *H; /* H contains the hash table */
static unsigned long K = 0; /* Size of the hash table */
static double cost = 0.0; /* Cost to insert all rels in the hash table */
/* Number of duplicates and rels on the current file */
static index_t ndup, nrels;
/* Number of duplicates and rels on all read files */
static index_t ndup_tot = 0, nrels_tot = 0;

/* sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
   and check for hash collisions */
unsigned long sanity_size;
int64_t  *sanity_a;
uint64_t *sanity_b;
unsigned long sanity_checked = 0;
unsigned long sanity_collisions = 0;
/* end sanity check */

static double factor = 1.0;

char ** files, ** files_already_renumbered, ** files_new;
unsigned int nb_files, nb_f_new, nb_f_renumbered;

static int is_for_dl; /* Do we reduce mod 2 or not */

static inline void
sanity_check (uint32_t i, int64_t a, uint64_t b)
{
  sanity_checked++;
  if (sanity_a[i] == 0)
  {
    sanity_a[i] = a;
    sanity_b[i] = b;
  }
  else if (sanity_a[i] != a)
  {
    sanity_collisions++;
    fprintf(stderr, "Collision between (%" PRId64 ",%" PRIu64 ") and "
                    "(%" PRId64 ",%" PRIu64 ")\n", sanity_a[i], sanity_b[i],
                    a, b);
  }
}

static inline void
print_warning_size ()
{
  index_t nodup = nrels_tot - ndup_tot;
  double full_table = 100.0 * (double) nodup / (double) K;
  fprintf(stderr, "Warning, hash table is %1.0f%% full\n", full_table);
  if (full_table >= 99.0)
  {
    fprintf(stderr, "Error, hash table is full\n");
    exit(1);
  }
  factor += 1.0;
}

/* Print the relation 'rel' in a line of the form:
    a,b:h_1,h_2,...,h_k
   with a (signed) and b (unsigned) written in hexa and
   and i_1 ... i_k (hexadecimal) are the indices of the ideals
*/
//TODO take care of bad ideals and add a col of 1 if necessary (added column is
//always 0)
static inline void
print_relation (FILE * file, buf_rel_t * rel)
{
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;

  p = d64toa16(buf, rel->a);
  *p++ = ',';
  p = u64toa16(p, rel->b);
  *p++ = ':';

  for (i = 0; i < rel->nb; i++)
  {
    if (rel->primes[i].e > 0)
    {
      op = p;
      p = u64toa16(p, (uint64_t) rel->primes[i].h);
      *p++ = ',';
      t = p - op;
      for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
        memcpy(p, op, t);
    }
  }
  // Add a column of 1 (it always has index 0) if asked
  if (renumber_tab->add_full_col)
  {
    p = u64toa16(p, (uint64_t) 0);
    *p++ = ',';
  }


  *(--p) = '\n';
  p[1] = 0;
  fputs(buf, file);
}

/* if duplicate is_dup = 1, else is_dup = 0
 * return i for sanity check
 */
static inline uint32_t
insert_relation_in_dup_hashtable (buf_rel_t * rel, unsigned int *is_dup)
{
  uint64_t h;
  uint32_t i, j;

  h = CA_DUP2 * (uint64_t) rel->a + CB_DUP2 * rel->b;
  i = h % K;
  j = (uint32_t) (h >> 32);
/* Note: in the case where K > 2^32, i and j share some bits.
 * The high bits of i are in j. These bits correspond therefore to far away
 * positions in the tables, and keeping them in j can only help.
 * FIXME: TODO: that's wrong!!! it would be better do take i from high
   bits instead!
 */
  while (H[i] != 0 && H[i] != j)
  {
    i++;
    if (UNLIKELY(i == K))
      i = 0;
    cost++;
  }

  if (H[i] == j)
    *is_dup = 1;
  else
  {
    H[i] = j;
    *is_dup = 0;
  }

  return i;
}

static inline void
compute_index_rel (buf_rel_t * rel)
{
  unsigned int i;
  p_r_values_t r;
  prime_t *pr = rel->primes;
  int side;
  weight_t len = rel->nb; // rel->nb can be modified by bad ideals

 /* HACK: a relation is on the form a,b:side0:side1
  * we put the side in primes[i].h
  */
  for (i = 0; i < len; i++)
  {
    if (pr[i].e > 0)
    {
      side = (int) pr[i].h;
      if (side != renumber_tab->rat)
      {
#ifndef FOR_FFS
#if DEBUG >= 1
  // Check for this bug : [#15897] [las] output "ideals" that are not prime
        if (!modul_isprime(&(pr[i].p)))
        {
          fprintf (stderr, "Error, relation with a=%" PRId64" b=%" PRIu64 " "
                           "contains %" PRpr " which is not prime.\nRemove "
                           "this relation from the file and re-run dup2.\n",
                           rel->a, rel->b, pr[i].p);
          abort();
        }
#endif
        r = (p_r_values_t) findroot(rel->a, rel->b, pr[i].p);
#else
        r = (p_r_values_t) findroot_ffs(rel->a, rel->b, pr[i].p);
#endif
      }
      else
        r = 0; // on rational side the value of r is meaningless
      
      int nb; //nb of ideal above the bad ideal
      index_t first_index; // first index of the ideals above a bad ideals
      if (renumber_is_bad (&nb, &first_index, renumber_tab, pr[i].p, r, side))
      {
        int exp_above[RENUMBER_MAX_ABOVE_BADIDEALS];
        handle_bad_ideals (exp_above, rel->a, rel->b, pr[i].p, pr[i].e);
        
        if (rel->nb + nb > rel->nb_alloc)
        {
          if (rel->nb_alloc == NB_PRIMES_OPT)
          {
            rel->nb_alloc += rel->nb_alloc >> 1;
            prime_t *p = rel->primes;
            SMALLOC(rel->primes, rel->nb_alloc, "realloc buffer primes");
            memcpy(rel->primes, p, NB_PRIMES_OPT * sizeof(prime_t));
          }
          else
          {
            rel->nb_alloc += rel->nb_alloc >> 1;
            rel->primes = (prime_t *) 
                          realloc(rel->primes, rel->nb_alloc * sizeof(prime_t));
          }
        }

        pr[i].h = first_index;
        pr[i].e = exp_above[0];
        for (int n = 1; n < nb; n++)
        {
          pr[rel->nb].h = first_index + n;
          pr[rel->nb].e = exp_above[n];
          rel->nb++;
        }
      }
      else
        pr[i].h = renumber_get_index_from_p_r(renumber_tab, pr[i].p, r, side);
    }
  }
}

#if 0
static inline void
print_prime (int64_t a, uint64_t b, unsigned long p, int e, int side,
             char *s)
{
  unsigned long j, k;
  p_r_values_t r;
  int nb;
  index_t first_index;

  if (side == renumber_tab->rat)
    r = 0;
  else
#ifndef FOR_FFS
    r = findroot(a, b, p);
#else
    r = findroot_ffs (a, b, p);
#endif

  // Here we have to take care of bad ideals that
  // generate several columns:
  if (renumber_is_bad (&nb, &first_index, renumber_tab, p, r, side))
  {
    int exp_above[RENUMBER_MAX_ABOVE_BADIDEALS];
    handle_bad_ideals (exp_above, a, b, p, e);
    for (int n = 0; n < nb; n++)
    {
      for (k = 0; k < (unsigned long) exp_above[n]; k++)
        sprintf (s, "%s%"PRxid",", s, first_index);
      first_index++;
    }
  }
  else
  {
    j = renumber_get_index_from_p_r (renumber_tab, p, r, side);
    for (k = 0; k < (unsigned long) e; k++)
      sprintf (s, "%s%lx,", s, j);
  }
}
#endif

static void
get_outfilename_from_infilename (char *infilename, const char *outfmt,
                                 const char *outdir, char **oname,
                                 char **oname_tmp)
{
    const char * suffix_in;
    const char * suffix_out;
    get_suffix_from_filename (infilename, &suffix_in);
    suffix_out = outfmt ? outfmt : suffix_in;

    char * newname = strdup(infilename);
    ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
    newname[strlen(newname)-strlen(suffix_in)]='\0';

    if(outdir != NULL)
    {
      int rc;
      const char * basename = path_basename(newname);
      rc = asprintf(oname_tmp, "%s/%s.tmp%s", outdir, basename, suffix_out);
      ASSERT_ALWAYS(rc >= 0);
      rc = asprintf(oname, "%s/%s%s", outdir, basename, suffix_out);
      ASSERT_ALWAYS(rc >= 0);
    }
    else
    {
      int rc;
      rc = asprintf(oname_tmp, "%s.tmp%s", newname, suffix_out);
      ASSERT_ALWAYS(rc >= 0);
      rc = asprintf(oname, "%s%s", newname, suffix_out);
      ASSERT_ALWAYS(rc >= 0);
    }

#if DEBUG >= 1
  fprintf (stderr, "DEBUG: Input file name: %s,\nDEBUG: temporary output file "
                   "name: %s,\nDEBUG: final output file name: %s\n", infilename,
                   *oname_tmp, *oname);
#endif
  free(newname);
}

static void
dup_print_stat (const char *s, index_t nrels, index_t ndup)
{
  index_t nrem = nrels - ndup;
  double pdup = 100.0 * ((double) ndup) / ((double) nrels);
  fprintf (stderr, "%s: nrels=%" PRid " dup=%" PRid " (%.2f%%) rem=%" PRid "\n",
                   s, nrels, ndup, pdup, nrem);
}

static void *
thread_only_hash (buf_arg_t *arg)
{
  unsigned int j;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;
  uint32_t i;
  unsigned int is_dup;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
    {
      if (!is_finish())
        NANOSLEEP();
      else if (cpt_rel_a == cpy_cpt_rel_b)
        pthread_exit(NULL);
    }

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->rels[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP();

    nrels++;
    nrels_tot++;
    i = insert_relation_in_dup_hashtable (my_rel, &is_dup);
#if DEBUG >= 1
    // They should be no duplicate in already renumbered file
    ASSERT_ALWAYS (is_dup == 0);
#endif
    if (i < sanity_size)
      sanity_check(i, my_rel->a, my_rel->b);
    if (cost >= factor * (double) (nrels_tot - ndup_tot))
      print_warning_size ();

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static void *
thread_dup2 (buf_arg_t *arg)
{
  unsigned int j, is_dup;
  uint32_t i;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *myrel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!is_finish())
        NANOSLEEP();
      else if (cpt_rel_a == cpy_cpt_rel_b)
          pthread_exit(NULL);

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    myrel = &(arg->rels[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP();

    nrels++;
    nrels_tot++;
    i = insert_relation_in_dup_hashtable (myrel, &is_dup);
    if (!is_dup)
    {
      if (i < sanity_size)
        sanity_check(i, myrel->a, myrel->b);
      if (cost >= factor * (double) (nrels_tot - ndup_tot))
        print_warning_size ();

      print_relation (arg->fd[0], myrel);
    }
    else
    {
      ndup++;
      ndup_tot++;
    }

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
  return NULL;
}


void *
thread_root(fr_t *mfr)
{
  unsigned int j;

  for (;;)
  {
    switch(mfr->ok)
    {
      case 0:
        NANOSLEEP();
        break;
      case 1 :
        for (j = mfr->num; j <= mfr->end; j++)
        {
          buf_rel_t *myrel = &(mfr->buf_data[j]);

          if (!is_for_dl) /* Do we reduce mod 2 */
            for (unsigned int i = 0; i < myrel->nb; i++)
              myrel->primes[i].e &= 1;

          compute_index_rel (myrel);
        }
        mfr->ok = 0;
        break;
      case 2:
        mfr->ok = 3;
        pthread_exit(NULL);
    }
  }
}

static void
usage(const char *argv0)
{
    fprintf (stderr, "Usage: %s [options] ", argv0);
    fprintf (stderr, "[ -filelist <fl> [-basepath <dir>] | file1 ... filen ]\n");
    fprintf (stderr, "Mandatory command line options:\n");
    fprintf (stderr, "     -poly xxx     - polynomial file\n");
    fprintf (stderr, "     -renumber xxx - file with renumbering table\n");
    fprintf (stderr, "     -K <K>        - size of the hashtable\n");
    fprintf (stderr, "\nOther command line options:\n");
    fprintf (stderr, "    -outdir dir  - by default input files are overwritten\n");
    fprintf (stderr, "    -outfmt .ext - output is written in .ext files\n");
    fprintf (stderr, "    -path_antebuffer <dir> - where is antebuffer\n");
#ifndef FOR_FFS
    fprintf (stderr, "    -dl          - do not reduce exponents modulo 2\n");
#endif
    exit (1);
}

int
main (int argc, char *argv[])
{
    argv0 = argv[0];
    //TODO remove useless polynomials
    cado_poly cpoly;
    const char *renumberfilename = NULL;
    char **p;

    /* print command line */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (int k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

#ifndef FOR_FFS
    is_for_dl = 0; /* By default we do dup2 for factorization */
    param_list_configure_switch(pl, "dl", &is_for_dl);
#else
    is_for_dl = 1; /* With FFS, not for dl is meaningless */
#endif

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }

    const char * polyfilename = param_list_lookup_string(pl, "poly");
    const char * outfmt = param_list_lookup_string(pl, "outfmt");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * outdir = param_list_lookup_string(pl, "outdir");
    renumberfilename = param_list_lookup_string(pl, "renumber");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    param_list_parse_ulong(pl, "K", &K);

    if (param_list_warn_unused(pl) || polyfilename == NULL)
      usage(argv0);

    cado_poly_init (cpoly);
#ifndef FOR_FFS
    if (!cado_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
#else
    if (!ffs_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
#endif

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (K == 0) {
        fprintf (stderr, "The K parameter is required\n");
        usage(argv0);
    }

    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "output compression format unsupported\n");
        usage(argv0);
    }

    if (renumberfilename == NULL)
      {
        fprintf (stderr, "Missing -renumber option\n");
        exit (1);
      }

    set_antebuffer_path (argv0, path_antebuffer);

    renumber_init (renumber_tab, cpoly, NULL);
    renumber_read_table (renumber_tab, renumberfilename);

  /* sanity check: since we allocate two 64-bit words for each, instead of
     one 32-bit word for the hash table, taking K/100 will use 2.5% extra
     memory */
  sanity_size = 1 + (K / 100);
  fprintf (stderr, "[checking true duplicates on sample of %lu cells]\n",
           sanity_size);
  sanity_a = (int64_t*)  malloc (sanity_size * sizeof (int64_t));
  if (sanity_a == NULL)
    {
      fprintf (stderr, "Error, cannot allocate sanity_a\n");
      exit (1);
    }
  memset (sanity_a, 0, sanity_size * sizeof (int64_t));

  sanity_b = (uint64_t*) malloc (sanity_size * sizeof (uint64_t));
  if (sanity_b == NULL)
    {
      fprintf (stderr, "Error, cannot allocate sanity_b\n");
      exit (1);
    }

  H = (uint32_t*) malloc (K * sizeof (uint32_t));
  if (H == NULL)
    {
      fprintf (stderr, "Error, cannot allocate hash table\n");
      exit (1);
    }
  memset (H, 0, K * sizeof (uint32_t));
  fprintf (stderr, "Allocated hash table of %lu entries (%luMb)\n", K,
           (K * sizeof (uint32_t)) >> 20);


  if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Provide either -filelist or freeform file names\n");
      usage(argv0);
  }


  /* Construct the two filelists : new files and already renumbered files */
  fprintf(stderr, "Constructing the two filelists...\n");
  files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;
  for (p = files, nb_files = 0; *p; p++)
    nb_files++;

  SMALLOC(files_already_renumbered, nb_files + 1, "files_already_renumbered");
  SMALLOC(files_new, nb_files + 1, "files_new");

  /* separate already process files
   * check if f_tmp is in raw format a,b:...:... or
   *            in renumbered format a,b:...
   */
  nb_f_new = 0;
  nb_f_renumbered = 0;
  for (p = files; *p; p++)
  {
    unsigned int count = 0;
    char s[1024];
    FILE *f_tmp = fopen_maybe_compressed (*p, "rb");
    ASSERT_ALWAYS(f_tmp != NULL);

    /* Look for first non-comment line */
    while (1) {
      char *ret = fgets (s, 1024, f_tmp);
      if (ret == NULL)
      {
        fprintf (stderr, "Error while reading %s\n", *p);
        exit (1);
      }
      if (strlen (s) >= 1023)
      {
        fprintf (stderr, "Too long line while reading %s\n", *p);
        exit (1);
      }
      size_t i = 0;
      while (s[i] == ' ')
        i++;
      if (s[i] != '#')
        break;
    }
    for (unsigned int i = 0; i < strlen (s); i++)
      count += s[i] == ':';

    if (count == 1)
      files_already_renumbered[nb_f_renumbered++] = *p;
    else if (count == 2)
      files_new[nb_f_new++] = *p;
    else {
      fprintf (stderr, "Error: invalid line (has %u colons): %s", count, s);
      exit(EXIT_FAILURE);
    }

    fclose_maybe_compressed (f_tmp, *p);
  }
  files_new[nb_f_new] = NULL;
  files_already_renumbered[nb_f_renumbered] = NULL;
  ASSERT_ALWAYS (nb_f_new + nb_f_renumbered == nb_files);
  fprintf (stderr, "%u files (%u new and %u already renumbered)\n", nb_files,
                   nb_f_new, nb_f_renumbered);


  fprintf (stderr, "Reading files already renumbered:\n");
  process_rels (files_already_renumbered, &thread_only_hash, NULL, 0, NULL, NULL,
                STEP_DUP2_PASS1);

  fprintf (stderr, "Reading new files:\n");

  for (char **p = files_new; *p != NULL; p++)
  {
    FILE *outputfile = NULL;
    char * oname = NULL;
    char * oname_tmp = NULL;
    char *fic[2] = { *p, NULL};
    info_mat_t info;
    nrels = ndup = 0;

    get_outfilename_from_infilename (*p, outfmt, outdir, &oname, &oname_tmp);
    outputfile = fopen_maybe_compressed(oname_tmp, "w");
    FILE *outfiles[2] = { outputfile, NULL};
    info = process_rels (fic, &thread_dup2, &thread_root, 0, (void **)outfiles, NULL,
                         STEP_DUP2_PASS2);
    ASSERT_ALWAYS (info.nrels == nrels);

    fclose_maybe_compressed(outputfile, oname_tmp);

#ifdef HAVE_MINGW /* For MinGW, rename cannot overwrite an existing file */
    remove (oname);
#endif
    if (rename(oname_tmp, oname))
    {
      fprintf(stderr, "Error while renaming %s into %s\n", oname_tmp, oname);
      abort();
    }

    // stat for the current file
    dup_print_stat (path_basename(*p), nrels, ndup);
    // stat for all the files already read
    dup_print_stat ("Total so far", nrels_tot, ndup_tot);

    free(oname);
    free(oname_tmp);
  }

  fprintf (stderr, "At the end: %" PRid " remaining relations\n",
                   nrels_tot - ndup_tot);

  fprintf (stderr, "At the end: hash table is %1.2f%% full\n"
                   "            hash table cost: %1.2f per relation\n",
                   100.0 * (double) (nrels_tot - ndup_tot) / (double) K,
                   1.0 + cost / (double) nrels_tot);
  fprintf (stderr, "  [found %lu true duplicates on sample of %lu relations]\n",
           sanity_collisions, sanity_checked);

  free (H);
  free (sanity_a);
  free (sanity_b);
  if (filelist)
    filelist_clear(files);
  free(files_already_renumbered);
  free(files_new);

  param_list_clear(pl);
  renumber_free (renumber_tab);
  cado_poly_clear (cpoly);
  return 0;
}
