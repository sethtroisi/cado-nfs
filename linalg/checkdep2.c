/* program to check raw dependencies (ker_raw) against input to linalg
   (matrix.small) */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

void
check (char *matname, char *kername, int depnum, int verbose)
{
  FILE *mat, *ker;
  unsigned long nrows, ncols, i, j, n, depr, k, nonzero, totnonzero = 0;
  unsigned long l;
  uint64_t *vec, depw;
  int ret;

  mat = fopen (matname, "r");
  if (mat == NULL)
    {
      fprintf (stderr, "Error, cannot open file %s\n", matname);
      exit (1);
    }

  ker = fopen (kername, "r");
  if (ker == NULL)
    {
      fprintf (stderr, "Error, cannot open file %s\n", kername);
      exit (1);
    }

  ret = fscanf (mat, "%lu %lu\n", &nrows, &ncols);
  if (ret != 2)
    {
      fprintf (stderr, "Error, cannot read nrows/ncols\n");
      exit (1);
    }
  if (verbose)
    printf ("Matrix has %lu rows and %lu columns\n", nrows, ncols);

  /* go to depnum dependency (one per line) */
  if (verbose)
    printf ("Checking dependency %d\n", depnum);
  while (depnum > 0)
    {
      while (getc (ker) != '\n');
      depnum --;
    }

  /* initialize vector with zeroes */
  n = (ncols - 1) / 64 + 1;
  vec = (uint64_t*) malloc (n * sizeof (uint64_t));
  for (i = 0; i < n; i++)
    vec[i] = 0;

  depr = 0;
  for (l = 0; l < nrows; l++)
    {
      /* read next bit from dependency */
      if (depr == 0) /* read next 64-bit word */
        {
          ret = fscanf (ker, "%lx", &depw);
          if (ret == 0)
            {
              fprintf (stderr, "Error, cannot read next dependency word\n");
              exit (1);
            }
          depr = 64;
        }
      
      /* read next row from matrix */
      ret = fscanf (mat, "%lu", &k); /* number of elements */
      if (ret == 0)
        {
          fprintf (stderr, "Error, cannot read length of next matrix row\n");
          exit (1);
        }
      for (i = 0; i < k; i++)
        {
          ret = fscanf (mat, "%lu", &j);
          if (ret == 0)
            {
              fprintf (stderr, "Error, cannot read next matrix element\n");
              exit (1);
            }
          if (depw & 1)
            vec[j >> 6] ^= (uint64_t) 1 << (j & 63);
        }
      fscanf (mat, "\n");
      
      depr --;
      depw >>= 1;
    }

  for (i = 0; i < n; i++)
    if (vec[i] != 0)
      {
        nonzero = 0;
        for (j = 0; j < 64; j++)
          {
            if ((vec[i] >> j) & 1)
              printf ("column %lu is non zero\n", 64 * i + j);
            nonzero += (vec[i] >> j) & 1;
          }
        printf ("word %lu of vector is non zero (%lu bits)\n", i, nonzero);
        totnonzero += nonzero;
      }
  if (totnonzero != 0)
    printf ("Total %lu non-zero bits\n", totnonzero);

  free (vec);
  fclose (mat);
  fclose (ker);
}

int
main (int argc, char *argv[])
{
  char *matname = NULL;
  char *kername = NULL;
  int verbose = 0;
  int depnum = 0;
  
  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-mat") == 0)
        {
          matname = argv[2];
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-ker") == 0)
        {
          kername = argv[2];
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-dep") == 0)
        {
          depnum = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 1 && strcmp (argv[1], "-v") == 0)
        {
          verbose ++;
          argc -= 1;
          argv += 1;
        }
      else
        {
          fprintf (stderr, "Invalid option: %s\n", argv[1]);
          exit (1);
        }
    }

  check (matname, kername, depnum, verbose);

  return 0;
}
