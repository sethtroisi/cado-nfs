/* dup2: 2nd pass

   Usage: dup2 [-out <dir>] [-rm] [-bz] [-filelist <fl>] -K <K> file1 ... filen

   Puts non-duplicate files (among file1 ... filen only) into:
   <dir>/file1 ... <dir>/filen

   If -out <dir> is missing, no output is done.

   Input files can be given on command line, or via a filelist file.
   (it is possible to do both)

   If -rm is given, the input files are removed after having been treated.
   This is also the case if -out . is given. 

   By default, the output will be gzipped, but if -bz is passed, then bzip2
   is used.

   Allocates a hash-table of 2^k 32-bit entries.

   Algorithm: for each (a,b) pair, we compute h(a,b) = (CA*a+CB*b) % 2^64.

   h has 64 bits. We know bits 2..6 are identical for a given slice, thus
   we remove them, and we obtain a 59-bit value h'.

   Let h' = j * 2^k + i.

   We store j % 2^32 at the next empty cell after index i in the hash-table.
*/

#define _GNU_SOURCE  /* in ordre to have asprintf */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>     // for unlink
#include <inttypes.h>
#include <ctype.h>  // for isspace

#define LIKELY(x)   __builtin_expect(x,1)
#define UNLIKELY(x) __builtin_expect(x,0)

#define CA 271828182845904523UL
#define CB 577215664901532889UL

#define LINE_SIZE 512
#define COMMAND_SIZE 1024

unsigned long rread = 0;
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

int bz = 0;
int rm = 0;

int is_gzip (const char * s) {
    unsigned int l = strlen (s);
    return l >= 3 && strcmp (s + l - 3, ".gz") == 0;
}

int is_bzip2 (const char * s) {
    unsigned int l = strlen (s);
    return l >= 4 && strcmp (s + l - 4, ".bz2") == 0;
}


char * gz2bz(const char * s)
{
    unsigned int l = strlen (s);
    char *ns = malloc(l+2); // 1 for additional char and 1 for \0
    strcpy(ns, s);
    strcpy(ns+l-3, ".bz2");
    return ns;
}

char * bz2gz(const char * s)
{
    unsigned int l = strlen (s);
    char *ns = strdup(s);
    strcpy(ns+l-4, ".gz");
    return ns;
}


/* infile is the input file
   if dirname is NULL, no output is done */
void
remove_dup (char *infile, char *dirname, uint32_t *H, unsigned long K)
{
  char line[LINE_SIZE];
  int64_t a;
  uint64_t b;
  uint64_t h;
  uint32_t i, j;
  int ret;
  FILE *f, *ofile = NULL;
  char command[COMMAND_SIZE];
  static double factor = 1.0;
  char *outfile = NULL;
  static int filenb=1;
  double full_table;

  fprintf (stderr, "Reading %d-th file: %s\n", filenb++, infile);
  if (is_gzip(infile))
      snprintf (command, COMMAND_SIZE, "gzip -dc %s", infile);
  else if (is_bzip2(infile))
      snprintf (command, COMMAND_SIZE, "bzip2 -dc %s", infile);
  else {
      fprintf(stderr, "Sorry, don't understand format of %s\n", infile);
      fprintf(stderr, "Let's skip it\n");
      return;
  }
  f = popen (command, "r");
  if (dirname != NULL) {
      if (bz) {
          if (is_gzip(infile))
              outfile = gz2bz(infile);
          else
              outfile = strdup(infile);
          snprintf(command, COMMAND_SIZE, "bzip2 -c --fast > %s/%s.part",
                  dirname, basename(outfile));
      } else {
          if (is_gzip(infile))
              outfile = strdup(infile);
          else
              outfile = bz2gz(infile);
          snprintf(command, COMMAND_SIZE, "gzip -c --best > %s/%s.part",
                  dirname, basename(outfile));
      }
      ofile = popen (command, "w");
  }

  for (;;)
    {
      if (rread % 1000000 == 0)
        fprintf (stderr, "Read %lu relations, %lu duplicates (%1.2f%%)\r",
                 rread, dupl, 100.0 * (double) dupl / (double) rread);
      
      if (fgets (line, LINE_SIZE, f) == NULL)
        break; /* end of file */


      ret = sscanf (line, "%ld,%lu", (int64_t *) &a,
                    (uint64_t *) &b);
      if (UNLIKELY(ret != 2))
        {
          fprintf (stderr, "Error, could not parse line\n");
          exit (1);
        }

      rread ++;

      h = CA * (uint64_t) a + CB * b;
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
      while (H[i] != 0 && H[i] != j)
        {
          i ++;
          if (UNLIKELY(i == K))
            i = 0;
          cost ++;
        }
      if (H[i] == j)
        {
          dupl ++;
          continue; /* probably duplicate */
        }

      if (i < sanity_size)
        {
          sanity_checked ++;
          if (sanity_a[i] == 0)
            {
              sanity_a[i] = a;
              sanity_b[i] = b;
            }
          else if (sanity_a[i] != a)
            {
              sanity_collisions ++;
              fprintf (stderr, "Collision between (%" PRId64 ",%" PRIu64 ") and (%" PRId64 ",%" PRIu64 ")\n",
                       sanity_a[i], sanity_b[i], a, b);
            }
        }

      nodu ++;
      if (cost >= factor * (double) rread)
        {
		  full_table = 100.0 * (double) nodu / (double) K;
          fprintf (stderr, "Warning, hash table is %1.0f%% full\n", full_table);
		  if (full_table >= 99) {
			fprintf(stderr, "Error, hash table is full\n");
			exit(1);
		  }
          factor += 1.0;
        }
      /* now H[i] = 0 */
      H[i] = j;
      if (dirname != NULL)
          fprintf (ofile, "%s", line);
    }
  pclose (f);
  fprintf (stderr, "Read %lu relations, %lu duplicates (%1.2f%%)\n",
                 rread, dupl, 100.0 * (double) dupl / (double) rread);
  if (dirname != NULL) {
      pclose (ofile);
      if (rm || (strcmp(dirname, ".") == 0)) {
          fprintf (stderr, "Removing old file %s\n", infile);
          int ret = unlink(infile);
          if (ret) {
              perror("Problem removing file");
              fprintf(stderr, "Let's hope that it's ok to continue!\n");
          }
      }
      char * s1, * s2;
      asprintf(&s1, "%s/%s.part", dirname, basename(outfile));
      asprintf(&s2, "%s/%s", dirname, basename(outfile));
      fprintf (stderr, "Renaming result file %s to %s\n", s1, s2);
      int ret = rename(s1, s2);
      if (ret) {
          perror("Problem renaming result file");
          fprintf(stderr, "Let's hope that it's ok to continue!\n");
      }
      free(s1);
      free(s2);
      free(outfile);
  }
}

int
main (int argc, char *argv[])
{
  char *dirname = NULL;
  unsigned long K = 4294967296; /* 2^32 */
  uint32_t *H;
  char *filelist = NULL;

  while (argc > 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-out") == 0)
        {
          dirname = argv[2];
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-rm") == 0)
        {
          rm = 1;
          argv += 1;
          argc -= 1;
        }
      else if (strcmp (argv[1], "-bz") == 0)
        {
          bz = 1;
          argv += 1;
          argc -= 1;
        }

      else if (strcmp (argv[1], "-filelist") == 0)
        {
          filelist = argv[2];
          argv += 2;
          argc -= 2;
        }

      else if (strcmp (argv[1], "-K") == 0)
        {
          K = (unsigned long) atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  if (K == 0)
    {
      fprintf (stderr, "Usage: dup2 [-rm] [-bz] [-out <dir>] [-filelist <fl>] -K <K> file1 ... filen\n");
      exit (1);
    }

  /* sanity check: since we allocate two 64-bit words for each, instead of
     one 32-bit word for the hash table, taking K/100 will use 2.5% extra
     memory */
  sanity_size = 1 + (K / 100);
  fprintf (stderr, "[checking true duplicates on sample of %lu cells]\n",
           sanity_size);
  sanity_a = (int64_t*)  malloc (sanity_size * sizeof (int64_t));
  memset (sanity_a, 0, sanity_size * sizeof (int64_t));
  sanity_b = (uint64_t*) malloc (sanity_size * sizeof (uint64_t));

  H = (uint32_t*) malloc (K * sizeof (uint32_t));
  if (H == NULL)
    {
      fprintf (stderr, "Error, cannot allocate hash table\n");
      exit (1);
    }
  memset (H, 0, K * sizeof (uint32_t));
  fprintf (stderr, "Allocated hash table of %lu entries (%zuMb)\n", K,
           (K * sizeof (uint32_t)) >> 20);

  if (filelist != NULL) {
      FILE *f;
      f = fopen(filelist, "r");
      if (f == NULL) {
          perror("Problem opening filelist");
          exit(1);
      }
      char relfile[LINE_SIZE];
      while (fgets(relfile, LINE_SIZE, f) != NULL) {
          // skip leading blanks
          char *rfile = relfile;
          while (isspace(rfile[0]))
              rfile++;
          // if empty line or comment line, continue
          if ((rfile[0] == '#') || (rfile[0] == '\0') || (rfile[0] == '\n'))
              continue;
          // get rid of newline char
          char *str = rfile;
          while (str[0] != '\0') {
              str++;
              if (str[0] == '\n')
                  str[0] = '\0';
          }
          remove_dup(rfile, dirname, H, K);
      }
      fclose(f);
  }

  while (argc > 1) {
      remove_dup (argv[1], dirname, H, K);
      argv ++;
      argc --;
  }

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

  return 0;
}
