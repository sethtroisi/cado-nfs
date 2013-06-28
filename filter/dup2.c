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
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>     // for unlink
#include <inttypes.h>
#include <ctype.h>  // for isspace
#include <fcntl.h>   /* for _O_BINARY */
#include "utils.h"
#include "relation.h"

#include "portability.h"
#include "utils.h"
#include "typedefs.h"

#include "filter_utils.h"

#ifdef FOR_FFS
#include "fppol.h"
#include "fq.h"
#include "utils_ffs.h"
#endif
#include "filter_badideals.h"

char *argv0; /* = argv[0] */

renumber_t renumber_table;

index_t dupl = 0;
index_t nodu = 0;
double cost = 0.0;

unsigned long K = 0;
uint32_t *H;

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
    fprintf(stderr, "Collision between (%"PRId64",%"PRIu64") and "
                    "(%"PRId64",%"PRIu64")\n", sanity_a[i], sanity_b[i], a,b);
  }
}

static inline void 
print_warning_size ( )
{
  double full_table = 100.0 * (double) nodu / (double) K;
  fprintf(stderr, "Warning, hash table is %1.0f%% full\n", full_table);
  if (full_table >= 99) 
  {
    fprintf(stderr, "Error, hash table is full\n");
    exit(1);
  }
  factor += 1.0;
}


static inline void
print_prime (int64_t a, uint64_t b, unsigned long p, int e, int side,
             char *s)
{
  unsigned long j, k;
  p_r_values_t r;
  int nb;
  index_t first_index;

  if (side == renumber_table->rat)
    r = 0;
  else
#ifndef FOR_FFS
    r = findroot(a, b, p);
#else
    r = findroot_ffs (a, b, p);
#endif

  // Here we have to take care of bad ideals that
  // generate several columns:
  if (renumber_is_bad (&nb, &first_index, renumber_table, p, r, side))
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
    j = renumber_get_index_from_p_r (renumber_table, p, r, side);
    for (k = 0; k < (unsigned long) e; k++)
      sprintf (s, "%s%lx,", s, j);
  }
}

/* infile is the input file
   if dirname is NULL, no output is done.
   * read new files only */
unsigned long
remove_dup_in_files (char ** files, const char *dirname, const char * outfmt,
                     int is_for_dl, unsigned int ab_base, 
                     renumber_t renumber_table)
{
    FILE * f_in;
    int p_in;
    const char * suffix_in;
    FILE * f_out;
    int p_out;
    const char * suffix_out;

    relation_stream rs;
    relation_stream_init(rs);

    rs->parse_only_ab = 1;

    for( ; *files ; files++) {
        const char * name = *files;
        index_t nodu0 = nodu, dupl0 = dupl;

        f_in = fopen_maybe_compressed2 (name, "r", &p_in, &suffix_in);
        ASSERT_ALWAYS(f_in != NULL);

        suffix_out = outfmt ? outfmt : suffix_in;
        char * newname = strdup(name);
        ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
        newname[strlen(newname)-strlen(suffix_in)]='\0';

        char * oname = NULL;
        char * oname_tmp = NULL;
        if (dirname) {
            int rc;
            rc = asprintf(&oname_tmp, "%s/pre-%s%s", dirname, path_basename(newname), suffix_out);
            ASSERT_ALWAYS(rc >= 0);
            rc = asprintf(&oname, "%s/%s%s", dirname, path_basename(newname), suffix_out);
            ASSERT_ALWAYS(rc >= 0);

            f_out = fopen_maybe_compressed2(oname_tmp, "w", &p_out, NULL);
            if (f_out == NULL) {
              fprintf (stderr, "Could not open file %s for writing: %s\n", 
                       oname_tmp, strerror(errno));
              exit (EXIT_FAILURE);
            }
        } else {
            f_out = NULL;
        }
        free(newname);

        relation_stream_bind(rs, f_in);

        uint32_t i;

        for (;;) {
            char line[RELATION_MAX_BYTES];
            unsigned int is_dup;
            buf_rel_t tmp;

            rs->parse_only_ab = 0;
            if (relation_stream_get (rs, line, 0, ab_base) < 0)
              break;

            tmp.a = rs->rel.a;
            tmp.b = rs->rel.b;
            i = insert_relation_in_dup_hashtable(H, K, &tmp, &cost, &is_dup);

            if(is_dup) {
                dupl++;
                continue;		/* probably duplicate */
            }
            if (i < sanity_size)
              sanity_check(i, rs->rel.a, rs->rel.b);

            nodu++;

            if (cost >= factor * (double) nodu) 
              print_warning_size ();   

            if (f_out) /* output renumbered relation */
              {
                char bufline[RELATION_MAX_BYTES];
                int buggy_line = 0; //to remove lines with non-prime "ideals"

                if (rs->rel.a < 0)
                  sprintf (bufline, "-%" PRIx64 ",%" PRIx64 ":",
                           (uint64_t) (- rs->rel.a), rs->rel.b);
                else
                  sprintf (bufline, "%" PRIx64 ",%" PRIx64 ":",
                            (uint64_t) rs->rel.a, rs->rel.b);

                if (!is_for_dl)
                  reduce_exponents_mod2 (&(rs->rel));

                /* r stands for side 0 (not rational side)
                   a stands for side 1 (not algebraic side)
                */
                for (int i = 0; i < rs->rel.nb_rp; i++)
                {
                  if (rs->rel.rp[i].e > 0)
                  {
#ifndef FOR_FFS
                    if (!modul_isprime(&(rs->rel.rp[i].p)))
                      buggy_line = 1;
                    else
#endif
                      print_prime (rs->rel.a, rs->rel.b, rs->rel.rp[i].p,
                                   rs->rel.rp[i].e, 0, bufline);
                  }
                }
                for (int i = 0; i < rs->rel.nb_ap; i++)
                {
                  if (rs->rel.ap[i].e > 0)
                  {
#ifndef FOR_FFS
                    if (!modul_isprime(&(rs->rel.ap[i].p)))
                      buggy_line = 1;
                    else
#endif
                      print_prime (rs->rel.a, rs->rel.b, rs->rel.ap[i].p,
                                   rs->rel.ap[i].e, 1, bufline);

                  }
                }

                if (!buggy_line)
                {
                  // bufline end with a extra ,
                  if (renumber_table->add_full_col)
                    fprintf (f_out, "%s0\n", bufline); //added col is always 0
                  else
                  {
                    bufline[strlen(bufline)-1] = '\0';
                    fprintf (f_out, "%s\n", bufline);
                  }
                }
                else
                {
                  fprintf (stderr, "WARNING: line with a non-prime \"ideals\"."
                                   "\nWARNING: %s", line);
                  rs->nrels--;
                  nodu--;
                }
              }


            if (relation_stream_disp_progress_now_p(rs)) {
                fprintf(stderr,
                        "Read %"PRid" relations, %"PRid" duplicates (%1.2f%%)"
                        " in %.1f s -- %.1f MB/s, %.1f rels/s\n",
                        rs->nrels, dupl, 100.0 * (double) dupl / (double) rs->nrels,
                        rs->dt, rs->mb_s, rs->rels_s);
            }
        }
        relation_stream_unbind(rs);
        if (f_out != NULL)
          {
            if (p_out) pclose(f_out); else fclose(f_out);
          }
        if (p_in) pclose(f_in); else fclose(f_in);


        if (dirname) {
            int ret = rename(oname_tmp, oname);
            if (ret) {
                perror("Problem renaming result file");
                fprintf(stderr, "Let's hope that it's ok to continue!\n");
            }
            free(oname);
            free(oname_tmp);
        }
        fprintf (stderr, "%s: %"PRid" relations, %"PRid" duplicates, "
                         "remains %"PRid"\n", path_basename(name),
                 nodu + dupl - (nodu0 + dupl0), dupl - dupl0, nodu - nodu0);
    }
    relation_stream_trigger_disp_progress(rs);
    fprintf(stderr,
            "Read %"PRid" relations, %"PRid" duplicates (%1.2f%%)"
            " in %.1f s -- %.1f MB/s, %.1f rels/s\n",
            rs->nrels, dupl, 100.0 * (double) dupl / (double) rs->nrels,
            rs->dt, rs->mb_s, rs->rels_s);
    unsigned long rread = rs->nrels;
    relation_stream_clear(rs);
    return rread;
}

void *
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
        NANOSLEEP;
      else if (cpt_rel_a == cpy_cpt_rel_b)
        pthread_exit(NULL);
    }

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->buf_data[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP;

    i = insert_relation_in_dup_hashtable (H, K, my_rel, &cost, &is_dup);
#if DEBUG >= 1
    // They should be no duplicate in already renumbered file
    ASSERT_ALWAYS (is_dup == 0);
#endif
    if (i < sanity_size)
      sanity_check(i, my_rel->a, my_rel->b);
    nodu++;
    if (cost >= factor * (double) nodu) 
      print_warning_size ();   

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

MAYBE_UNUSED static void *
thread_print(buf_arg_t *arg)
{
  unsigned int j;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!is_finish())
        NANOSLEEP;
      else if (cpt_rel_a == cpy_cpt_rel_b)
          pthread_exit(NULL);

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->buf_data[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP;

    if (my_rel->nb != 0)
      print_relation (arg->f_deleted, my_rel); //FIXME where do we print

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}


void *
thread_root(fr_t *mfr) 
{
  buf_rel_t *myrel;
  unsigned int i, j;
  unsigned int is_dup;

  //ep = &(H.ht[H.hm]);
  for (;;)
  {
    switch(mfr->ok) 
    {
      case 0:
        NANOSLEEP;
        break;
      case 1 :
        for (j = mfr->num; j <= mfr->end; j++)
        {
          myrel = &(mfr->buf_data[j]);
          
          i = insert_relation_in_dup_hashtable (H, K, myrel, &cost, &is_dup);
          if (is_dup)
          {
           // WARNING FIXME H, cost, factor, nodu, dupl is not volatile...
           // FIX les mettre dans fr_t
            if (i < sanity_size)
              sanity_check(i, myrel->a, myrel->b);
            nodu++;
            if (cost >= factor * (double) nodu) 
              print_warning_size ();   

            compute_index_rel (renumber_table, myrel);
          }
          else
          {
            dupl++;
            myrel->nb = 0; /* meaning: do not print this relation */
          }
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
  buf_arg_t buf_arg;
  buf_rel_t *buf_rel;
  p_r_values_t pmin;
  index_t min_index;
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
    int is_for_dl = 0; /* Be default we do dup2 for factorization */
    param_list_configure_switch(pl, "dl", &is_for_dl);
#else
    int is_for_dl = 1; /* With FFS, not for dl is meaningless */
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
    pmin = ulong_nextprime (MIN (cpoly->pols[0]->lim, cpoly->pols[1]->lim));
#else
    if (!ffs_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
    //compute pmin
    {
      sq_t p;
      sq_set_ti(p, MIN (cpoly->pols[0]->lim, cpoly->pols[1]->lim) + 1);
      do
      {
        sq_monic_set_next(p, p, 64);
      } while (!sq_is_irreducible(p));
      pmin = fppol64_get_ui_sparse (p);
    }
#endif

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (K == 0)
        usage(argv0);

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
  
    renumber_init (renumber_table, cpoly);
    renumber_read_table (renumber_table, renumberfilename);
    // Find the index that corresponds to the min value of alim and rlim (for
    // purge)
    // FIXME not correct in the case of two alg side
    min_index = renumber_get_index_from_p_r (renumber_table, pmin, 0,
                                                           renumber_table->rat);
    fprintf (stderr, "Renumbering struct: min_index=%"PRid"\n", min_index);



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
    FILE *f_tmp = fopen_maybe_compressed (*p, "r");
    ASSERT_ALWAYS(f_tmp != NULL);

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
    for (unsigned int i = 0; i < strlen (s); i++)
      count += s[i] == ':';
    
    if (count == 2)
      files_new[nb_f_new++] = *p;
    else
      files_already_renumbered[nb_f_renumbered++] = *p;

    fclose_maybe_compressed (f_tmp, *p);
  }
  files_new[nb_f_new] = NULL;
  files_already_renumbered[nb_f_renumbered] = NULL;
  ASSERT_ALWAYS (nb_f_new + nb_f_renumbered == nb_files);
  fprintf (stderr, "%u files (%u new and %u already renumbered)\n", nb_files, 
                   nb_f_new, nb_f_renumbered);


 //call prempt_scan_rel 2 times with two diff filelist and two diff callback fct

  SMALLOC(buf_rel,SIZE_BUF_REL, "buf_rel");
  MEMSETZERO(&buf_arg, 1);
  buf_arg.buf_data = buf_rel;
  buf_arg.needed = NEEDED_AB;
  
  fprintf (stderr, "Reading files already renumbered:\n");
  prempt_scan_relations (files_already_renumbered, &thread_only_hash, &buf_arg,
                         NULL);

  fprintf (stderr, "Reading new files:\n");
  index_t rread = 0;
  //buf_arg.needed = NEEDED_ABP;
  //prempt_scan_relations (files_new, &thread_print, &buf_arg, &thread_root);
  /* pass 1: we read new files, remove duplicates, and renumber them */
#ifndef FOR_FFS
  rread += remove_dup_in_files (files_new, basepath, outfmt, is_for_dl, 10,
                                renumber_table);
#else
  rread += remove_dup_in_files (files_new, basepath, outfmt, is_for_dl, 16,
                                renumber_table);
#endif


  fprintf (stderr, "Read %"PRid" relations, %"PRid" duplicates (%1.2f%%)\n",
           rread, dupl, 100.0 * (double) dupl / (double) rread);

  free (H);

  fprintf (stderr, "     %"PRid" remaining relations (hash table %1.2f%% full)\n",
           nodu, 100.0 * (double) nodu / (double) K);
  fprintf (stderr, "     Hash-table cost %1.2f per relation\n",
           1.0 + cost / (double) rread);

  free (sanity_a);
  free (sanity_b);
  fprintf (stderr, "[found %lu true duplicates on sample of %lu relations]\n",
           sanity_collisions, sanity_checked);


  if (filelist)
    filelist_clear(files);
  SFREE(files_already_renumbered);
  SFREE(files_new);
  SFREE(buf_rel);

  param_list_clear(pl);
  renumber_free (renumber_table);
  cado_poly_clear (cpoly);
  return 0;
}
