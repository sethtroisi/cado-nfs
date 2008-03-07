/* Simple program to find optimal addition chains and generate an unrolled
   addition chain for stage 1 of P+1 or ECM. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/* We seach for Lucas chains of length at most maxlen for integers less 
   than maxn */
static int maxlen = 10;
static unsigned long maxn = 128;
char *len; /* len[i] is the length of the best chain found for i so far */
unsigned long *prev; /* Stores next-to-last element of chain for i */
static unsigned long chain[255];
static unsigned long pracbest[255];

/* Table of multipliers for PRAC. prac_mul[i], i>0, has continued fraction 
   sequence of all ones but with a 2 in the i-th place, and 
   prac_mul[0] is all ones, i.e. the golden ratio */
#define PRAC_NR 10
static const double prac_mul[PRAC_NR] = 
  {1.61803398874989484820 /* 0 */, 1.38196601125010515179 /* 1 */, 
   1.72360679774997896964 /* 2 */, 1.58017872829546410471 /* 3 */, 
   1.63283980608870628543 /* 4 */, 1.61242994950949500192 /* 5 */,
   1.62018198080741576482 /* 6 */, 1.61721461653440386266 /* 7 */, 
   1.61834711965622805798 /* 8 */, 1.61791440652881789386 /* 9 */};

/* Simple primality test by trial division */
int isprime (unsigned long n)
{
  unsigned long i;

  if (n < 2)
    return 0;
  if (n % 2UL == 0UL)
    return (n == 2UL);
  for (i = 3; i*i <= n; i += 2UL)
    if (n % i == 0UL)
      return 0;

  return 1;
}

/* The cost of the binary addition chain */
int bincost (unsigned long e)
{
  int b, k = 0;
  
  while (e % 2UL == 0UL)
    {
      e >>= 1;
      k++;
    }
  
  if (e == 1UL)
    return k;

  b = (int) floor (log ((double) e) / log (2.));
  if (e < 3UL*(1UL << (b-1)))
    return 2*b + k - 1;
  else
    return 2*b + k;
}


/***********************************************************************
   Generating Lucas chains with Montgomery's PRAC algorithm. Code taken 
   from GMP-ECM, mostly written by Paul Zimmermann, and slightly 
   modified here
************************************************************************/

#define ADD 1
#define DUP 1
static unsigned int
lucas_cost_pp1 (unsigned long n, double v)
{
  unsigned int c;
  unsigned long d, e, r;

  d = n;
  r = (unsigned long) ((double) d / v + 0.5);
  if (r >= n)
    return (ADD * n);
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD; /* initial duplicate and final addition */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
        }
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          c += 3 * ADD; /* 3 additions */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      else if (d <= 4 * e)
        { /* condition 3 */
          d -= e;
          c += ADD; /* one addition */
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          c += ADD + DUP; /* one addition, one duplicate */
        } 
      /* now d+e is odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      /* now d is odd and e even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else /* necessarily e is even */
        { /* condition 9 */
          e /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
    }
  
  return c;
}


/* Returns the length of the shortest Lucas chain for n found by PRAC 
   using the first m multiplers from prac_mul. If mul is not NULL,
   stores the index of the first multiplier that produced such a short
   chain */

static unsigned int prac_best (unsigned long n, int m, int *mul)
{
  int i, bestmul;
  unsigned int bestcost, plaincost, c;

  if (m > PRAC_NR)
    m = PRAC_NR;
  bestcost = plaincost = lucas_cost_pp1 (n, prac_mul[0]);;
  bestmul = 0;
  for (i = 1; i < m; i++)
    {
      c = lucas_cost_pp1 (n, prac_mul[i]);
      if (c < bestcost)
	{
	  bestcost = c;
	  bestmul = i;
	}
      /* If this multiplier is an improvement over the Golden Ratio, 
	 increase its count */
      if (c < plaincost)
	pracbest[i]++;
    }
  
  /* If no multiplier could improve over the Golden Ratio, increase
     the GR's count */
  if (bestcost == plaincost)
    pracbest[0]++;

  if (mul != NULL)
    *mul = bestmul;

  return bestcost;
}


/* Print an addition chain for k. Mostly copied from GMP-ECM */
static void prac_printchain (unsigned long k)
{
  unsigned long d, e, r;
  int m;

  if (k == 2)
    {
      printf ("dup (A, A) /* Start for k = 2 */\n");
      return;
    }

  assert (k % 2 == 1);

  /* Find the best multiplier for this k */
  prac_best (k, PRAC_NR, &m);
  
  d = k;
  r = (unsigned long) ((double) d / prac_mul[m] + 0.5);

  d = k - r;
  e = 2 * r - k;
  printf ("set (B, A) /* Start for k=%lu, multiplier %f: init and double */\n",
	  k, prac_mul[m]);
  printf ("set (C, A)\n");
  printf ("dup (A, A)\n");
  
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
	  printf ("swap (A, B) /* Rule 0: swap */\n");
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
	  printf ("add (t, A, B, C) /* Rule 1 */\n");
	  printf ("add (t2, t, A, B)\n");
	  printf ("add (B, B, t, A)\n");
	  printf ("set (A, t2)\n");
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
	  printf ("add (B, A, B, C) /* Rule 2 */\n");
	  printf ("dup (A, A)\n");
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
	  printf ("adds (C, B, A, C) /* Rule 3 */\n");
	  printf ("swap (B, C)\n");
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
	  printf ("add (B, B, A, C) /* Rule 4 */\n");
	  printf ("dup (A, A)\n");
        }
      else if (d % 2 == 0) /* d+e is now odd */
        { /* condition 5 */
          d /= 2;
	  printf ("add (C, C, A, B) /* Rule 5 */\n");
	  printf ("dup (A, A)\n");
        }
      else if (d % 3 == 0) /* d is odd, e even */
        { /* condition 6 */
          d = d / 3 - e;
	  printf ("dup (t, A) /* Rule 6 */\n");
	  printf ("add (t2, A, B, C)\n");
	  printf ("adds (A,  t, A, A)\n");
	  printf ("adds (C,  t, t2, C)\n");
	  printf ("swap (B, C)\n");
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
	  printf ("add (t, A, B, C) /* Rule 7 */\n");
	  printf ("adds (B, t, A, B)\n");
	  printf ("dup (t, A)\n");
	  printf ("adds (A, A, t, A)\n");
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
	  printf ("add (T, A, B, C) /* Rule 8 */\n");
	  printf ("add (C, C, A, B)\n");
	  printf ("swap (B, t)\n");
	  printf ("dup (t, A)\n");
	  printf ("adds (A, A, T, A)\n");
        }
      else /* necessarily e is even */
        { /* condition 9: never happens? */
          e /= 2;
	  printf ("add (C, C, B, A)\n");
	  printf ("dup (B, B)\n");
        }
    }
  
  printf ("add (A, A, B, C) /* Final add */\n");
  assert (d == 1);
}



/***************************************************************
   Code for finding optimal addition chains by exaustive search 
****************************************************************/


/* At function entry, chain[0 ... curlen] contains a Lucas chain of
   length curlen. If this is a new optimum length for the end value
   of this chain (i.e., chain[curlen]), then this new optimum length 
   is stored in len[].
   If curlen < maxlen, we try to extend this chain by doubling or 
   adding elements of the chain, where for an addition of two elements
   their difference must be in the chain as well */

void chain_extend (int curlen)
{
  int i, j;
  unsigned long s, d;

  /* If the sequence so far is a new optimum for its current end value,
     remember the new optimal length for that value */
  s = chain[curlen];
  if (len[s] == 0 || len[s] > curlen)
    {
      len[s] = curlen;
      if (curlen >= 1)
        prev[s] = chain[curlen - 1];
    }

  if (curlen >= maxlen)
    return;

#if 0
  /* Try doubling the third-to-last sequence element. 
     FIXME: Is this ever worthwhile? */
  if (curlen > 1 && 2 * chain[curlen - 2] > chain[curlen])
    {
      s = 2 * chain[curlen - 2];
      if (s < maxn)
        {
          chain[curlen + 1] = s;
          chain_extend (curlen + 1);
        }
    }
#endif

  /* Try doubling the second-to-last sequence element */
  if (curlen > 0 && 2 * chain[curlen - 1] > chain[curlen])
    {
      s = 2 * chain[curlen - 1];
      if (s < maxn)
        {
          chain[curlen + 1] = s;
          chain_extend (curlen + 1);
        }
    }

  /* Try doubling the last sequence element */
  s = 2 * chain[curlen];
  if (s < maxn)
    {
      chain[curlen + 1] = s;
      chain_extend (curlen + 1);
    }


  /* Look for i so that chain[curlen] - chain[i] is in the chain and if
     such an i is found, extend the chain by chain[curlen] + chain[i].
     FIXME: Do we need to try adding chain[curlen-1] and smaller elements? */
  for (i = 0, j = curlen; i < curlen; i++)
    {
      s = chain[curlen] + chain[i];
      if (s >= maxn)
        break;
      assert (chain[i] < chain[curlen]);
      d = chain[curlen] - chain[i];
      /* See if d is in the chain */
      while (j > 0 && chain[j] > d)
	j--;
      if (chain[j] == d)
	{
	  chain[curlen + 1] = s;
	  chain_extend (curlen + 1);
	}
    }
}


void find_opt_chain ()
{
  unsigned long i, uninit = 0, nr_primes, sum_prac_primes, 
    sum_optimal_primes;
  int j, k, s;
  
  len = malloc (maxn * sizeof (char));
  prev = malloc (maxn * sizeof (unsigned long));
  chain[0] = 1;
  chain[1] = 2;
  len[1] = 0;
  len[2] = 1;
  for (i = 3; i < maxn; i++)
    len[i] = 0;
  chain_extend (1);
  
  nr_primes = 0; /* The number of primes in the search interval */
  sum_prac_primes = 0; /* Sum of lengths of PRAC chains for all primes */
  sum_optimal_primes = 0; /* Sum of lengths of optimal chains for all primes */
  
  /* For each multiplier other than the GR, we remember how many times
     that multiplier beat the GR. We also remember how often GR was optimal.
     Init these counts to 0 */
  for (i = 0; i < 255; i++)
    pracbest[i] = 0;
  
  if (1)
    for (i = 3; i < maxn; i += 2)
      {
	if (len[i] != 0)
	  {
	    if (isprime (i))
	      {
		unsigned int prac_cost = prac_best (i, PRAC_NR, NULL);
		nr_primes++;
		sum_prac_primes += prac_cost;
		sum_optimal_primes += len[i];
		
		assert ((unsigned int) len[i] <= prac_cost);

		/* For each k print k, the length l of the optimal chain 
		   for k, k^(1/l) (i.e. the average multiplier in each 
		   step of the chain), the cost of the binary chain and the
		   cost of the best chain found by PRAC. If the best PRAC 
		   chain is worse than the optimal chain, print an asterisk */
		   
		printf ("%lu p: %d, r=%lu, radix %f, bin: %d, prac: %u%s\n", 
			i, (int) len[i], prev[i], 
			exp(log((double) (i)) / (double) (len[i])),
			bincost (i), prac_cost,
			(prac_cost > (unsigned int) len[i]) ? "*":"");
	      }
	    else if (0)
	      {
		printf ("%lu c: %d, r=%lu, radix %f, binary cost: %d\n", 
			i, (int) len[i], prev[i], 
			exp(log((double) (i)) / (double) (len[i])),
			bincost (i));
	      }
	  }
	else
	  if (uninit == 0)
	    uninit = i;
      }
  
  printf ("%lu primes, sum of optimal lengths %lu, sum of PRAC lengths %lu\n",
	  nr_primes, sum_optimal_primes, sum_prac_primes);
  printf ("Histogram of optimal prac multiplier: ");
  for (i = 0; i < 255; i++)
    if (pracbest[i] != 0)
      printf ("%d:%d ", (int) i, (int) pracbest[i]);
  printf ("\n");

  if (uninit == 0)
    uninit = maxn;
  printf ("Smallest number for which no chain was found: %lu\n", uninit);


  /* Look for cases n=p*q where the concatenated chains of p and q are longer
     than the optimal chain for n */
  s = (int) sqrt ((double) maxn);
  for (k = 3; k <= s; k += 2)
    {
      int p;
      for (j = k, p = k * j; (unsigned long) p < maxn; j += 2, p += 2*k)
	if (len[p] > 0 && (int) len[p] < (int) len[k] + (int) len [j])
	  {
	    if ((int) len[k] + (int) len[j] > (int) len[p] + 1)
	      printf ("%d:%d beats %d*%d (%d+%d), difference %d\n", 
		      p, (int) len[p], k, j, (int) len[k], (int) len[j],
		      (int) len[k] + (int) len[j] - (int) len[p]);
	    else
	      printf ("%d:%d beats %d*%d (%d+%d)\n", 
		      p, (int) len[p], k, j, (int) len[k], (int) len[j]);
	  }
    }
  free (len);
  len = NULL;
  free (prev);
  prev = NULL;
}


/* Print pseudo-code for the Lucas chain for a stage 1 of P+1 or ECM.
   With preprocessor defines for add(), adds(), dup(), set() and swap()
   which implement the arithmetic, this pseudo-code can be compiled into 
   a fully unrolled P+1 or ECM stage 1.

   TBD: this uses only addition chains found by PRAC for the individual 
   primes so far. Use optimal chains where they are shorter than PRAC
   (happens rarely), and take into account chains for products of two 
   primes where the chain for that composite is shorted than the 
   concantenated chain for the two primes (could save a few percent)
*/

void 
generate_stage1 (unsigned long B1, unsigned long oldB1)
{
  unsigned long p, pp;
  
  p = 2;
  while (p <= B1)
    {
      for (pp = p; pp <= B1; pp *= p)
	if (pp > oldB1)
	  prac_printchain (p);

      for (p++; ! isprime(p); p++); /* Find next prime (slow) */
    }
}


void usage()
{
  printf ("-o [maxlen [maxn]]  Find optimal addition chains up to length maxlen for\n"
          "                    integers up to maxn\n");
  printf ("-p [B1 [oldB1]]     Print pseudo-code for addition chain for stage 1 with\n"
	  "B1 and possibly oldB1 (i.e. for extending stage 1)\n");
}


int main (int argc, char ** argv)
{
  if (argc < 2)
    {
      usage();
      exit (EXIT_FAILURE);
    }

  if (strcmp (argv[1], "-o") == 0)
    {
      /* Find optimal chains */
      if (argc > 2)
	{
	  maxlen = atoi (argv[2]);
	  maxn = 1UL << maxlen;
	}
      
      if (argc > 3)
	{
	  maxn = atoi (argv[3]);
	}
      find_opt_chain ();
    }
  else if (strcmp (argv[1], "-p") == 0)
    {
      /* Print code for stage 1 */
      unsigned long B1 = 100, oldB1 = 1;

      if (argc > 2)
	  B1 = atoi (argv[2]);
      if (argc > 3)
	  oldB1 = atoi (argv[3]);

      generate_stage1 (B1, oldB1);
    }
  else
    {
      usage();
      exit (EXIT_FAILURE);
    }

  return 0;
}
