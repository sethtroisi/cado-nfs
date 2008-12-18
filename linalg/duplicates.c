/* duplicates --- remove duplicate relations

Copyright 2008 Francois Morain, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>

/* CADO-NFS header files */
#include "utils/utils.h"
#include "hashpair.h"
#include "files.h"
#include "gzip.h"

/* We use a hash table with open adressing and linear probing, each entry
   containing a 64-bit value identifying uniquely the corresponding (a,b)
   pair, at least for |a|, |b| < 2^53. */
typedef struct {
  unsigned long hashmod;
  uint64_t     *hashtab;
} Hashtable_t;

/* Search the first entry in the hash table Hab that is either 0 or H,
   starting from index h. If H was not already in the hash table, add it.
   Return non-zero iff H was not already in the table.
   Assumes h < Hab->hashmod.
*/
static int
is_ab_new (Hashtable_t *Hab, unsigned long h, uint64_t H)
{
  ASSERT(h < Hab->hashmod);

  while (Hab->hashtab[h] != 0 && Hab->hashtab[h] != H)
    h = (h + 1 == Hab->hashmod) ? 0 : h + 1; /* should avoid branching */

  if (Hab->hashtab[h] == 0)
    {
      Hab->hashtab[h] = H;
      return 1; /* new entry */
    }
  else
    return 0; /* duplicate entry */
}

/* Remove duplicates from file 'file' and print them to 'out'.
   irel    - total number of relations read so far
   nrels   - total number of remaining relations (is incremented)
   Hab     - hash table (for current slice)
   slice   - number of slices is 2^slice
   slice0  - index of current slice, 0 <= slice0 < 2^slice
   file    - input file
   verbose - verbosity level
*/
static int
remove_duplicates_from_file (FILE *out, unsigned long *irel,
			     unsigned long *nrels, Hashtable_t *Hab, int slice,
			     int slice0, FILE *file, int verbose)
{
    int ret;
    long a;
    unsigned long b;
    unsigned long file_duplicates = 0;         /* duplicates in this file */
    static unsigned long total_duplicates = 0; /* duplicates in all files */
    char str[STR_LEN_MAX];
    uint64_t HAB;
    unsigned int hab, mask = (((unsigned)1)<<slice)-1;

    while(1)
      {
	ret = fread_buf (str, file);
	if (ret != 1)
	  break; /* wrong line or EOF */

	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "    %lu rel. read (%lu duplicates) at %2.2lf\n",
		    *irel, file_duplicates, seconds ());

	if (sscanf (str, "%ld,%lu", &a, &b) != 2)
	  {
	    fprintf (stderr, "Invalid relation: %s\n", str);
	    continue;
	  }

        HAB = getInitialAddress (a, b, MAGIC_HC0, MAGIC_HC1);
	hab = HAB % Hab->hashmod;

	/* FIXME: with that setting, the low 'slice' bits of hab are
	 *always* equal to 'slice0', thus we only hit 1/2^slice of the
	 hash table at the first time */
	if((slice > 0) && ((hab & mask) != (unsigned) slice0))
	  continue; /* this relation is not in that slice */

	if (is_ab_new (Hab, hab, HAB)){
	    *nrels += 1;
	    fprintf(out, "%s", str);
	}
	else{
	    if (file_duplicates ++ < 10 && verbose)
              fprintf(stderr, "(%ld, %lu) appears more than once\n", a, b);
	}
    }
    total_duplicates += file_duplicates;
    fprintf(stderr, "  Found %lu duplicates in this file (total %lu, %1.0f%%)\n",
	    file_duplicates, total_duplicates,
            100.0 * (double) total_duplicates / (double) *irel);
    return ret;
}

/* Treat all relations from files ficname[0], ..., ficname[nbfic-1]
   for slice 'slice0', and put non-duplicates in 'out'.
   nrels: number of remaining relations (is incremented)
*/
static int
remove_duplicates (FILE *out, char *ficname[], int nbfic, unsigned long *nrels,
		   Hashtable_t *Hab, int slice, int slice0, int verbose)
{
    FILE *relfile;
    int ret = 0;
    unsigned long irel = 0; /* number of relations read in that slice */
    int i;
    
    ASSERT(nbfic > 0);
    for(i = 0; i < nbfic; i++){
        relfile = gzip_open (ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "  Adding file %s\n", ficname[i]);
	ret = remove_duplicates_from_file (out, &irel, nrels, Hab, slice,
					   slice0, relfile, verbose);
	if(ret == 0) {
	    fprintf(stderr, "Warning: error when reading file %s\n", ficname[i]);
	    break;
	}
        gzip_close (relfile, ficname[i]);
    }
    fprintf (stderr, "Scanned %lu relations\n", irel);
    
    return (ret == -1);
}

/* print usage and exit */
static void
usage (void)
{
  fprintf (stderr, "Usage: duplicates -nrels n [options] file1 ... filen\n");
  fprintf (stderr, "       -v            verbose\n");
  fprintf (stderr, "       -nrels n      input files have <= n relations\n");
  fprintf (stderr, "       -slice k      split input in 2^k parts\n");
  fprintf (stderr, "       -out   f      output file is f (if ends in .gz, use gzip)\n");
  exit(1);
}

#define MAX_FILE_NAME 1024

int
main (int argc, char **argv)
{
    FILE *outfile = NULL;
    Hashtable_t Hab;
    char **fic;
    char *outname = NULL;
    unsigned int nfic;
    int ret, k, slice = 0, slice0;
    unsigned int nrelsmax = 0, Hsize;
    unsigned long nrels = 0; /* number of remaining relations */
    int verbose = 0;
    
    /* print command-line arguments */
    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    while (argc > 1 && argv[1][0] == '-'){
	if(argc > 1 && strcmp(argv[1], "-v") == 0){
          verbose ++;
          argc -= 1;
          argv += 1;
	}
	else if(argc > 2 && strcmp(argv[1], "-nrels") == 0){
	    nrelsmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp(argv[1], "-slice") == 0){
	    slice = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp(argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else
	  {
	    fprintf (stderr, "Error, invalid option: %s\n", argv[1]);
	    usage ();
	  }
    }

    if (nrelsmax == 0)
      {
        fprintf (stderr, "Error, missing -nrels option\n");
        usage ();
      }

    if (slice > 0)
	Hsize = nrelsmax >> slice;
    else
	Hsize = nrelsmax;

    fic = argv+1;
    nfic = argc-1;
    fprintf (stderr, "Maximal number of relations is %u\n", nrelsmax);
    Hab.hashmod = getHashMod (Hsize + Hsize / 2, 1);
    Hab.hashtab = (uint64_t*) malloc (Hab.hashmod * sizeof(uint64_t));

    if (outname != NULL)
      outfile = gzip_open(outname, "w");
    else
      outfile = stdout;
    fprintf(stderr, "reading files of relations...\n");
    for (slice0 = 0; slice0 < (1<<slice); slice0++)
      {
	/* (re)set hash table to zero */
	memset (Hab.hashtab, 0, Hab.hashmod * sizeof(uint64_t));
	if (slice)
	  fprintf (stderr, "Performing slice %d\n", slice0);
	ret = remove_duplicates (outfile, fic, nfic, &nrels, &Hab, slice,
				 slice0, verbose);
      }
    if (outname != NULL)
      gzip_close(outfile, outname);
    fprintf(stderr, "Number of relations left: %lu\n", nrels);

    free (Hab.hashtab);

    return 0;
}
