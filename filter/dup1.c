/* dup1: 1st duplicate pass, split relation files into NSLICES slices
   (adapted from check).

   Usage:
   dup1 [-bz] -out <dir> file1 ... filen

   Files file1 ... filen are split into NSLICES slices in
   <dir>/0/filej ... <dir>/31/filej.

   If option -bz is given, then the output is compressed with bzip2
   instead of gzip.
   Input can be in gzipped or bzipped format.
*/

#define _GNU_SOURCE
#define NSLICES_LOG 2
#define NSLICES (1 << NSLICES_LOG)
#define CA 314159265358979323UL
#define CB 271828182845904523UL

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <limits.h> /* for CHAR_BIT */
#include <unistd.h>
#include <assert.h>

clock_t start;

static int checked=0;

int bz;

    
int
is_gzip (const char * s)
{
  unsigned int l = strlen (s);
  return l >= 3 && strcmp (s + l - 3, ".gz") == 0;
}

int is_bzip2 (const char * s)
{
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


#define LINE_SIZE 512
#define COMMAND_SIZE 1024

/* output relations in dirname/0/name, ..., dirname/31/name */
int
check_stream (const char *name, FILE * stream, const char *dirname, char *Only)
{
    int lnum, i, ok;
    FILE *ofile[NSLICES];
    unsigned long count[NSLICES];
	char line[LINE_SIZE];
    int64_t a;
    uint64_t b;
    uint64_t h;
    char command[COMMAND_SIZE];
    char *newname;

    if (bz && is_gzip(name))
        newname = gz2bz(name);
    else if ((!bz) && is_bzip2(name))
        newname = bz2gz(name);
    else
        newname = strdup(name);

    for (i = 0; i < NSLICES; i++)
      if (Only[i])
        {
          if (!bz)
              snprintf (command, COMMAND_SIZE, "gzip -c --best > %s/%d/%s",
                      dirname, i, basename(newname));
          else 
              snprintf (command, COMMAND_SIZE, "bzip2 -c --fast > %s/%d/%s",
                      dirname, i, basename(newname));
          ofile[i] = popen (command, "w");
          if (ofile[i] == NULL)
            {
              fprintf (stderr, "Error, popen call failed\n");
              exit (1);
            }
          count[i] = 0;
        }

    for (lnum = 0;; lnum++)
      {
	if (fgets (line, LINE_SIZE, stream) == NULL)
          break; /* end of file */

	if (line[0] == '#')
            continue; /* comment line */

	if (sscanf (line, "%ld,%lu", (int64_t *) &a,
                    (uint64_t *) &b) != 2)
          {
	    fprintf (stderr, "Failed on line %d in %s: %s\n", lnum, name,
                     line);
            exit (1); /* errors are fatal here */
        }

	ok = 1;

        h = CA * (uint64_t) a + CB * b;
        /* Using the low bit of h is not a good idea, since then
           odd values of i are twice more likely. The second low bit
           also gives a small bias with RSA768 (but not for random
           coprime a, b). We use here the NSLICES_LOG high bits.
        */
        i = h >> (64 - NSLICES_LOG);

        /* print relation */
        if (Only[i])
          {
            if (fprintf (ofile[i], "%s", line) < 0 ) {
                fprintf (stderr, "Problem writing to file %d\n", i);
                exit(1);
            }
          }
        count[i] ++;

        if (++checked % 1000000 == 0)
          fprintf (stderr, "Split %d relations in %fs\r", checked,
                   (clock()-start)*1.0/CLOCKS_PER_SEC);
      }

    for (i = 0; i < NSLICES; i++)
      if (Only[i])
        {
          pclose (ofile[i]);
          fprintf (stderr, "%d:%lu ", i, count[i]);
        }
    fprintf (stderr, "\n");
    free(newname);
    return 0;
}

int
main (int argc, char * argv[])
{
    int had_error = 0;
    char *dirname = NULL;
    char Only[NSLICES];
    int i, only = -1;
    bz = 0;

    while (argc > 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-out") == 0)
          {
            dirname = argv[2];
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-bz") == 0)
          {
            argv += 1;
            argc -= 1;
            bz = 1;
          }
        else if (strcmp (argv[1], "-only") == 0)
          {
            only = atoi (argv[2]);
            assert (0 <= only && only < NSLICES);
            argv += 2;
            argc -= 2;
          }
        else
          {
            fprintf (stderr, "Unknown option %s\n", argv[1]);
            exit (1);
          }
      }

    if (dirname == NULL)
      {
        fprintf (stderr, "Usage: dup1 [-bz] [-only i] -out <dir> file1 ... filen\n");
        exit (1);
      }

    if (only == -1) /* split all slices */
      {
        for (i = 0; i < NSLICES; i++)
          Only[i] = 1;
      }
    else /* split only slide i */
      {
        for (i = 0; i < NSLICES; i++)
          Only[i] = i == only;
      }

    start = clock();

    if (argc == 1)
      {
        fprintf (stderr, "Error, stdin input is not allowed\n");
        exit (1);
      }
    else
      {
        int i = 0;
        for (i = 1 ; i < argc ; i++) {
            FILE * f;
            if (strcmp(argv[i], "-") == 0) {
              fprintf (stderr, "Error, stdin input is not allowed\n");
              exit (1);
            } else if (is_gzip(argv[i]) || is_bzip2(argv[i])) {
                char command[1024];
                if (is_gzip(argv[i])) 
                    snprintf(command, sizeof(command), "gzip -dc %s", argv[i]);
                else
                    snprintf(command, sizeof(command), "bzip2 -dc %s", argv[i]);
                f = popen(command, "r");
                had_error |= check_stream (argv[i], f, dirname, Only);
                pclose(f);
            } else {
                fprintf (stderr, "Error, compressed input expected\n");
                exit (1);
            }
        }
    }

    fprintf (stderr, "Split %d relations in %fs\n", checked,
             (clock()-start)*1.0/CLOCKS_PER_SEC);

    if (had_error)
      return 1;
    else
      return 0;
}
