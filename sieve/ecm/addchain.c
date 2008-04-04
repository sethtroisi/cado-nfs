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
static int beat_prac = 0;
static int firstcode = 1;

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

/* One more than the highest code number the byte code generator can produce */
#define MAXCODE 32
typedef char literal_t;
int lastcode = 255;
int compress = 0;
/* Stores in codehist2[prev*MAXCODE+i] the frequency of the "prev, i" byte 
   combination */
unsigned char codehist[MAXCODE];
unsigned char codehist2[MAXCODE*MAXCODE];
#define CODER_HISTLEN 10
literal_t coder_history[CODER_HISTLEN];
int coder_nrstored = 0;
int nr_literals = 0, nr_outputcodes = 0;

#if 1
#define DICT_NRENTRIES 9
const int dict_len[DICT_NRENTRIES] = {2, 2, 2, 3, 3, 3, 4, 4, 6};
const literal_t *dict_entry[DICT_NRENTRIES] = 
  {"\xB\xA", "\x3\x0", "\x3\x3", "\xB\xA\x3", "\x0\x3\x0","\x3\x3\x0",
   "\x3\x0\x3\x0", "\x3\xB\xA\x3", "\x3\x0\x3\x0\x3\x0"};
const int dict_code[DICT_NRENTRIES] = {13, 14, 15, 16, 17, 18, 19, 20, 21};
#else
#define DICT_NRENTRIES 0
const int dict_len[DICT_NRENTRIES] = {};
const literal_t *dict_entry[DICT_NRENTRIES] = {};
const int dict_code[DICT_NRENTRIES] = {};
#endif


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

void print_chain (unsigned long *chain, int len)
{
  int i, j, k, l;
  unsigned long d;

  if (len == 0)
    return;

  printf ("%lu ", chain[0]);

  for (i = 1; i <= len; i++)
    {
      /* First see if we can get chain[i] by a doubling */
      for (j = 0; j < i; j++)
	if (2 * chain[j] == chain[i])
	  {
	    printf ("(%d, %d, -) %lu ", i-1-j, i-1-j, chain[i]);
	    goto cont_i;
	  }
      /* Now we look for j, k, l so that chain[i] == chain[j] + chain[k]
	 and chain[j] - chain[k] == chain[l] */
      for (j = 0; j < i; j++)
	for (k = 0; k < i; k++)
	  {
	    if (chain[i] != chain[j] + chain[k])
	      continue;
	    assert (j != k);
	    d = chain[j] - chain[k];
	    for (l = 0; l < i; l++)
	      if (chain[l] == d)
		{
		  printf ("(%d, %d, %d) %lu ", i-1-j, i-1-k, i-1-l, chain[i]);
		  goto cont_i;
		}
	  }
      abort (); /* Could not generate chain[i] ?!? */
    cont_i: ;
    }
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

static void 
bytecoder_output (int c)
{
  assert (c < MAXCODE);
  if (firstcode)
    printf ("%d", c);
  else
    printf (",%d", c);
  firstcode = 0;
  codehist[c]++;
  if (lastcode != 255)
    codehist2[lastcode*MAXCODE+c]++;
  lastcode = c;
  nr_outputcodes++;
}


/* Returns the length of the matching part of "buf" (with "len" valid bytes)
   and the dictionary entry "entry" */

static int 
dict_matchlen (const literal_t *buf, const int len, const int entry)
{
  int i;

  if (entry == -1)
    return 0;

  assert (0 <= entry && entry < DICT_NRENTRIES);

  for (i = 0; i < len && i < dict_len[entry]; i++)
    if (buf[i] != dict_entry[entry][i])
      break;
  
  return i;
}


/* Returns the index of the dictionary entry that has the longest match with
   "buf", where "buf" contains "len" valid bytes. If "partial" is != 0, 
   then matches that don't match the complete dictionary entry are accepted,
   otherwise only matches with complete dictionary entries are accepted.
   Return -1 if nothing matched at all */

static int
dict_longestmatch (const literal_t *buf, const int len, const int partial)
{
  int i, t, matchlen = 0, matchidx = -1;
  
  for (i = 0; i < DICT_NRENTRIES; i++)
    {
      t = dict_matchlen (buf, len, i);
      if (t > matchlen && (partial || t == dict_len[i]))
	{
	  matchlen = t;
	  matchidx = i;
	}
    }
  
  return matchidx;
}


/* Adds a literal to the coder_history[] FIFO buffer */

static void 
coder_histadd (const int c)
{
  coder_history[coder_nrstored++] = (literal_t) c;
  assert (coder_nrstored < CODER_HISTLEN);
  nr_literals++;
}


/* Removes n literals from the coder_history[] FIFO buffer */

static void
coder_histremove (int n)
{
  int i;
  if (n > coder_nrstored)
    abort();

  for (i = 0; i + n < coder_nrstored; i++)
    coder_history[i] = coder_history[i + n];

  coder_nrstored -= n;
}


/* Output the best dictionary match we have at the moment, or a literal
   if coder_history[] matches no complete dictionary entry */

static void
coder_outputbest ()
{
  int best, bestlen;
  
  assert (coder_nrstored > 0);
  
  /* Find longest complete match */
  best = dict_longestmatch (coder_history, coder_nrstored, 0);
  bestlen = dict_matchlen (coder_history, coder_nrstored, best);
  
  /* If there was a complete dictionary match, output that, 
     otherwise output the first entry in coder_history[] as a literal */
  if (best != -1 && bestlen == dict_len[best])
    {
      bytecoder_output (dict_code[best]);
      coder_histremove (bestlen);
    }
  else
    {
      bytecoder_output ((int) (coder_history[0]));
      coder_histremove (1);
    }
}

static void 
bytecoder (const int c)
{
  /* At this point, the first coder_nrstored literals (posibly 0!) of
     coder_history[] agree with some dictionary entry */
  
  /* Now add the new literal */
  coder_histadd (c);
  
  if (compress)
    {
      /* See if coder_history[] still matches one of the dictionary entries */
      
      int best, bestlen = 0;
      
      /* If all coder_history[] matches some dictionary entry, we do nothing */
      /* Otherwise we repeatedly output the longest complete dictionary match 
	 or literal code until all of coder_history[] (then possibly empty!) 
	 matches some dictionary entry again */
      
      while (1)
	{
	  best = dict_longestmatch (coder_history, coder_nrstored, 1);
	  bestlen = dict_matchlen (coder_history, coder_nrstored, best);
	  if (bestlen == coder_nrstored)
	    break;
	  coder_outputbest ();
	}
    }
  else
    {
      /* No compression: flush the stored code */
      bytecoder_output (coder_history[0]);
      coder_nrstored = 0;
    }
}

static void bytecoder_flush ()
{
  while (coder_nrstored > 0)
    coder_outputbest ();
}

/* Print an addition chain for k. Mostly copied from GMP-ECM */
static void prac_printchain (unsigned long k, int bytecode)
{
  unsigned long d, e, r;
  int m;
  
  if (k == 2)
    {
      if (bytecode)
	bytecoder ((literal_t) 12);
      else
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
  if (bytecode)
    bytecoder ((literal_t) 10);
  else
    {
      printf ("set (B, A) /* Start for k=%lu, multiplier %f: init and "
	      "double */\n", k, prac_mul[m]);
      printf ("set (C, A)\n");
      printf ("dup (A, A)\n");
    }
  
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
	  if (bytecode)
	    bytecoder ((literal_t) 0);
	  else
	    printf ("swap (A, B) /* Rule 0: swap */\n");
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
	  if (bytecode)
	    bytecoder ((literal_t) 1);
	  else
	    {
	      printf ("add (t, A, B, C) /* Rule 1 */\n");
	      printf ("add (t2, t, A, B)\n");
	      printf ("add (B, B, t, A)\n");
	      printf ("set (A, t2)\n");
	    }
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
	  if (bytecode)
	    bytecoder ((literal_t) 2);
	  else
	    {
	      printf ("add (B, A, B, C) /* Rule 2 */\n");
	      printf ("dup (A, A)\n");
	    }
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
	  if (bytecode)
	    bytecoder ((literal_t) 3);
	  else
	    {
	      printf ("adds (C, B, A, C) /* Rule 3 */\n");
	      printf ("swap (B, C)\n");
	    }
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
	  if (bytecode)
	    bytecoder ((literal_t) 4);
	  else
	    {
	      printf ("add (B, B, A, C) /* Rule 4 */\n");
	      printf ("dup (A, A)\n");
	    }
        }
      else if (d % 2 == 0) /* d+e is now odd */
        { /* condition 5 */
          d /= 2;
	  if (bytecode)
	    bytecoder ((literal_t) 5);
	  else
	    {
	      printf ("add (C, C, A, B) /* Rule 5 */\n");
	      printf ("dup (A, A)\n");
	    }
        }
      else if (d % 3 == 0) /* d is odd, e even */
        { /* condition 6 */
          d = d / 3 - e;
	  if (bytecode)
	    bytecoder ((literal_t) 6);
	  else
	    {
	      printf ("dup (t, A) /* Rule 6 */\n");
	      printf ("add (t2, A, B, C)\n");
	      printf ("adds (A,  t, A, A)\n");
	      printf ("adds (C,  t, t2, C)\n");
	      printf ("swap (B, C)\n");
	    }
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
	  if (bytecode)
	    bytecoder ((literal_t) 7);
	  else
	    {
	      printf ("add (t, A, B, C) /* Rule 7 */\n");
	      printf ("adds (B, t, A, B)\n");
	      printf ("dup (t, A)\n");
	      printf ("adds (A, A, t, A)\n");
	    }
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
	  if (bytecode)
	    bytecoder ((literal_t) 8);
	  else
	    {
	      printf ("add (t, A, B, C) /* Rule 8 */\n");
	      printf ("add (C, C, A, B)\n");
	      printf ("swap (B, t)\n");
	      printf ("dup (t, A)\n");
	      printf ("adds (A, A, t, A)\n");
	    }
        }
      else /* necessarily e is even */
        { /* condition 9: never happens? */
          e /= 2;
	  if (bytecode)
	    bytecoder ((literal_t) 9);
	  else
	    {
	      printf ("add (C, C, B, A)\n");
	      printf ("dup (B, B)\n");
	    }
        }
    }
  
  if (bytecode)
    bytecoder ((literal_t) 11);
  else
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
  if (beat_prac)
    {
      if (len[s] > curlen)
	{
	  printf ("Optimal chain for %lu beats PRAC:", s);
	  for (i = 0; i <= curlen; i++)
	    printf (" %lu", chain[i]);
	  printf ("\n");
	  print_chain (chain, curlen);
	  printf ("\n");
	}
    }
  else if (len[s] == 0 || len[s] > curlen)
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
  chain[2] = 3; /* If chain[2] == 4, then all following chain elements are 
		   even, and want chains only for odd multipliers - even
		   ones are trivially reduced */
  len[1] = 0;
  len[2] = 1;
  for (i = 3; i < maxn; i++)
    len[i] = 0;

  /* If we look for cases where the optimal chain beats PRAC, init the
     lengths to the best length PRAC can find */
  if (beat_prac)
    {
      for (i = 3; i < maxn; i += 2)
	if (isprime (i))
	  len[i] = prac_best (i, PRAC_NR, NULL);
      printf ("Histogram of optimal prac multiplier over all primes < %lu: ",
	      maxn);
      for (i = 0; i < 255; i++)
	if (pracbest[i] != 0)
	  printf ("%d:%d ", (int) i, (int) pracbest[i]);
      printf ("\n");
    }

  chain_extend (2);
  
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
generate_stage1 (unsigned long B1, unsigned long oldB1, int bytecode)
{
  unsigned long p, pp;
  
  for (p = 0; p < MAXCODE; p++)
    codehist[p] = 0;
  for (p = 0; p < MAXCODE*MAXCODE; p++)
    codehist2[p] = 0;

  p = 2;
  while (p <= B1)
    {
      for (pp = p; pp <= B1; pp *= p)
	if (pp > oldB1)
	  prac_printchain (p, bytecode);

      for (p++; ! isprime(p); p++); /* Find next prime (slow) */
    }

  if (bytecode)
    {
      bytecoder_flush ();
      printf ("\n");
      for (p = 0; p < MAXCODE; p++)
	printf ("%lu: %d\n", p, (int) codehist[p]);
      for (p = 0; p < MAXCODE*MAXCODE; p++)
	if (codehist2[p] != 0)
	  printf ("%lu,%lu: %d\n", p / MAXCODE, p % MAXCODE, (int)codehist2[p]);
      printf ("%d literals, %d output codes\n", nr_literals, nr_outputcodes);
    }
}


void usage()
{
  printf ("-o [maxlen [maxn]]  Find optimal addition chains up to length maxlen for\n"
          "                    integers up to maxn\n");
  printf ("-p [B1 [oldB1]]     Print pseudo-code for addition chain for stage 1 with\n"
	  "                    B1 and possibly oldB1 (i.e. for extending stage 1)\n");
  printf ("-pb [B1 [oldB1]]    Print bytecode for addition chain\n");
  printf ("-pbc [B1 [oldB1]]   Print bytecode with dictionary compression\n");
  printf ("-b                  Find optimal chains that beat PRAC\n");
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
  else if (strcmp (argv[1], "-b") == 0)
    {
      /* Try to beat PRAC */
      beat_prac = 1;
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

      generate_stage1 (B1, oldB1, 0);
    }
  else if (strcmp (argv[1], "-pb") == 0)
    {
      /* Print bytecode for stage 1 */
      unsigned long B1 = 100, oldB1 = 1;

      if (argc > 2)
	  B1 = atoi (argv[2]);
      if (argc > 3)
	  oldB1 = atoi (argv[3]);

      generate_stage1 (B1, oldB1, 1);
    }
  else if (strcmp (argv[1], "-pbc") == 0)
    {
      /* Print bytecode for stage 1 */
      unsigned long B1 = 100, oldB1 = 1;

      if (argc > 2)
	  B1 = atoi (argv[2]);
      if (argc > 3)
	  oldB1 = atoi (argv[3]);
      compress = 1;
      generate_stage1 (B1, oldB1, 1);
    }
  else
    {
      usage();
      exit (EXIT_FAILURE);
    }

  return 0;
}
