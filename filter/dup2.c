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
#include <unistd.h>     // for unlink
#include <inttypes.h>
#include <ctype.h>  // for isspace
#include <fcntl.h>   /* for _O_BINARY */
#include "utils.h"
#include "relation.h"

#include "portability.h"
#include "macros.h"

#define CA 271828182845904523UL
#define CB 577215664901532889UL

unsigned long dupl = 0;
unsigned long nodu = 0;
double cost = 0.0;

/* sanity check: we store (a,b) pairs for 0 <= i < sanity_size,
   and check for hash collisions */
unsigned long sanity_size;
int64_t  *sanity_a;
uint64_t *sanity_b;
unsigned long sanity_checked = 0;
unsigned long sanity_collisions = 0;
/* end sanity check */

/* infile is the input file
   if dirname is NULL, no output is done */
unsigned long
remove_dup_in_files (char ** files, const char *dirname, const char * outfmt,
                     uint32_t * H, unsigned long K, unsigned int ab_base,
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
        int raw = 0;

        f_in = fopen_maybe_compressed2 (name, "r", &p_in, &suffix_in);
        ASSERT_ALWAYS(f_in != NULL);

        /* check if f_in is in raw format a,b:...:... or in renumbered format
           a,b:... */
        char s[1024];
        char *ret = fgets (s, 1024, f_in);
        if (ret == NULL)
          {
            fprintf (stderr, "Error while reading %s\n", name);
            exit (1);
          }
        if (strlen (s) >= 1023)
          {
            fprintf (stderr, "Too long line while reading %s\n", name);
            exit (1);
          }
        for (unsigned int i = 0; i < strlen (s); i++)
          raw += s[i] == ':';
        raw = (raw == 2) ? 1 : 0;

        /* close and reopen the file */
        if (p_in) pclose(f_in); else fclose(f_in);
        f_in = fopen_maybe_compressed2 (name, "r", &p_in, &suffix_in);
        ASSERT_ALWAYS(f_in != NULL);

        suffix_out = outfmt ? outfmt : suffix_in;
        char * newname = strdup(name);
        ASSERT_ALWAYS(strlen(suffix_in) <= strlen(newname));
        newname[strlen(newname)-strlen(suffix_in)]='\0';

        char * oname = NULL;
        char * oname_tmp = NULL;
        /* if raw = 0, the input file is already in the new format, no need
           to write it again */
        if (dirname && raw) {
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

        uint64_t h;
        uint32_t i, j;
        static double factor = 1.0;
        double full_table;

        for (;;) {
            char line[RELATION_MAX_BYTES];

            /* if raw = 0, i.e., the file contains lines a,b:..., then it was
               already renumbered, thus we only need to read a,b */
            rs->parse_only_ab = (raw == 0);
            if (relation_stream_get (rs, line, 0, ab_base) < 0)
              break;

            h = CA * (uint64_t) rs->rel.a + CB * rs->rel.b;
            /* dup1 uses the high 5 bits of h to identify the slices
               but now we use a different hash function, so we can keep these
               5 bits
               */
            i = h % K;
            j = (uint32_t) (h >> 32);
            /* Note: in the case where K > 2^32, i and j share some bits.
             * The high bits of i are in j. These bits correspond therefore to
             * far away positions in the tables, and keeping them in j can only
             * help.
             * FIXME:
             * TODO: that's wrong!!! it would be better do take i from high
             * bits instead!
             */
            while (H[i] != 0 && H[i] != j) {
                i++;
                if (UNLIKELY(i == K))
                    i = 0;
                cost++;
            }
            if (H[i] == j) {
                dupl++;
                continue;		/* probably duplicate */
            }

            if (i < sanity_size) {
                sanity_checked++;
                if (sanity_a[i] == 0) {
                    sanity_a[i] = rs->rel.a;
                    sanity_b[i] = rs->rel.b;
                } else if (sanity_a[i] != rs->rel.a) {
                    sanity_collisions++;
                    fprintf(stderr,
                            "Collision between (%" PRId64 ",%" PRIu64 ") and (%"
                            PRId64 ",%" PRIu64 ")\n", sanity_a[i], sanity_b[i], rs->rel.a,
                            rs->rel.b);
                }
            }

            nodu++;
            if (cost >= factor * (double) rs->nrels) {
                full_table = 100.0 * (double) nodu / (double) K;
                fprintf(stderr, "Warning, hash table is %1.0f%% full\n",
                        full_table);
                if (full_table >= 99) {
                    fprintf(stderr, "Error, hash table is full\n");
                    exit(1);
                }
                factor += 1.0;
            }
            /* now H[i] = 0 */
            H[i] = j;

            if (f_out) /* output renumbered relation */
              {
                int first = 1;

                fprintf (f_out, "%" PRId64 ",%" PRIu64 ":",
                         rs->rel.a, rs->rel.b);
#if (!defined FOR_FFS) && (!defined FOR_NFS_DL)
                reduce_exponents_mod2 (&(rs->rel));
#endif
                for (int i = 0; i < rs->rel.nb_rp; i++)
                  {
                    unsigned long j;
                    j = renumber_get_index_from_p_r (renumber_table,
                                                     rs->rel.rp[i].p, 0, 0);
                    for (int k = 0; k < rs->rel.rp[i].e; k++)
                      {
                        if (first)
                          first = 0;
                        else
                          fputc (',', f_out);
                        fprintf (f_out, "%lu", j);
                      }
                  }
                for (int i = 0; i < rs->rel.nb_ap; i++)
                  {
                    unsigned long j;
                    j = renumber_get_index_from_p_r (renumber_table,
                                          rs->rel.ap[i].p, rs->rel.ap[i].r, 1);
                    for (int k = 0; k < rs->rel.ap[i].e; k++)
                      {
                        if (first)
                          first = 0;
                        else
                          fputc (',', f_out);
                        fprintf (f_out, "%lu", j);
                      }
                  }
                fprintf (f_out, "\n");
              }


            if (relation_stream_disp_progress_now_p(rs)) {
                fprintf(stderr,
                        "Read %lu relations, %lu duplicates (%1.2f%%)"
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


        if (dirname && raw) {
            fprintf(stderr, "%s/{%s => %s}\n",
                    dirname, path_basename(oname_tmp), path_basename(oname));
            int ret = rename(oname_tmp, oname);
            if (ret) {
                perror("Problem renaming result file");
                fprintf(stderr, "Let's hope that it's ok to continue!\n");
            }
            free(oname);
            free(oname_tmp);
        }
    }
    relation_stream_trigger_disp_progress(rs);
    fprintf(stderr,
            "Read %lu relations, %lu duplicates (%1.2f%%)"
            " in %.1f s -- %.1f MB/s, %.1f rels/s\n",
            rs->nrels, dupl, 100.0 * (double) dupl / (double) rs->nrels,
            rs->dt, rs->mb_s, rs->rels_s);
    unsigned long rread = rs->nrels;
    relation_stream_clear(rs);
    return rread;
}

void usage()
{
    fprintf (stderr, "Usage: dup2 -poly xxx [-out <dir>] [-basepath <dir>] [-filelist <fl>] -renumber xxx -K <K> file1 ... filen\n");
    exit (1);
}

int
main (int argc, char *argv[])
{
    unsigned long K = 0;
    uint32_t *H;
    const char *renumberfilename = NULL;
    cado_poly cpoly;

    /* print command line */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (int k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int bz = 0;
    int ab_hexa = 0;
    param_list_configure_switch(pl, "bz", &bz);
    param_list_configure_switch(pl, "abhexa", &ab_hexa);

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

    param_list_parse_ulong(pl, "K", &K);

    if (param_list_warn_unused(pl) || polyfilename == NULL)
      usage();

    cado_poly_init (cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    if (basepath && !filelist) {
        fprintf(stderr, "-basepath only valid with -filelist\n");
        exit(1);
    }

    if (K == 0)
        usage();

    if (bz) {
        if (outfmt) {
            fprintf(stderr, "-bz and -outfmt are mutually exclusive");
            usage();
        } else {
            outfmt = ".bz2";
        }
    }
    if (outfmt && !is_supported_compression_format(outfmt)) {
        fprintf(stderr, "output compression format unsupported\n");
        usage();
    }

    if (renumberfilename == NULL)
      {
        fprintf (stderr, "Missing -renumber option\n");
        exit (1);
      }

    renumber_t renumber_table;
    renumber_init (renumber_table, cpoly);
    renumber_read_table (renumber_table, renumberfilename);

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
      usage();
  }

  char ** files = filelist ? filelist_from_file (basepath, filelist, 0) : argv;
  unsigned long rread = remove_dup_in_files (files, basepath, outfmt, H, K,
                                             (ab_hexa)?16:10, renumber_table);
  if (filelist) filelist_clear(files);


  fprintf (stderr, "Read %lu relations, %lu duplicates (%1.2f%%)\n",
           rread, dupl, 100.0 * (double) dupl / (double) rread);

  free (H);

  fprintf (stderr, "     %lu remaining relations (hash table %1.2f%% full)\n",
           nodu, 100.0 * (double) nodu / (double) K);
  fprintf (stderr, "     Hash-table cost %1.2f per relation\n",
           1.0 + cost / (double) rread);

  free (sanity_a);
  free (sanity_b);
  fprintf (stderr, "[found %lu true duplicates on sample of %lu relations]\n",
           sanity_collisions, sanity_checked);

  param_list_clear(pl);
  renumber_free (renumber_table);
  cado_poly_clear (cpoly);
  return 0;
}
