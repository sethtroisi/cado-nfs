#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "prac_bc.h"

#define CODER_HISTLEN 10
static int compress = 0;
static literal_t coder_history[CODER_HISTLEN];
static int coder_nrstored = 0;
static char *buffer;
static size_t bufalloc, buffull; /* Allocated size and current number of 
				    codes in buffer */

#if 0
#define DICT_NRENTRIES 9
const int dict_len[DICT_NRENTRIES] = {2, 2, 2, 3, 3, 3, 4, 4, 6};
const literal_t *dict_entry[DICT_NRENTRIES] = 
  {"\xB\xA", "\x3\x0", "\x3\x3", "\xB\xA\x3", "\x0\x3\x0","\x3\x3\x0",
   "\x3\x0\x3\x0", "\x3\xB\xA\x3", "\x3\x0\x3\x0\x3\x0"};
const int dict_code[DICT_NRENTRIES] = {13, 14, 15, 16, 17, 18, 19, 20, 21};
#else
#define DICT_NRENTRIES 4
const int dict_len[DICT_NRENTRIES] = {2, 2, 3, 4};
const literal_t *dict_entry[DICT_NRENTRIES] = 
  {"\xB\xA", "\x3\x0", "\x3\xB\xA", "\x3\x0\x3\x0"};
const int dict_code[DICT_NRENTRIES] = {13, 14, 15, 16};
#endif

/* Table of multipliers for PRAC. prac_mul[i], i>0, has continued fraction 
   sequence of all ones but with a 2 in the i-th place, and 
   prac_mul[0] is all ones, i.e. the golden ratio */
static const double prac_mul[PRAC_NR_MULTIPLIERS] = 
  {1.61803398874989484820 /* 0 */, 1.38196601125010515179 /* 1 */, 
   1.72360679774997896964 /* 2 */, 1.58017872829546410471 /* 3 */, 
   1.63283980608870628543 /* 4 */, 1.61242994950949500192 /* 5 */,
   1.62018198080741576482 /* 6 */, 1.61721461653440386266 /* 7 */, 
   1.61834711965622805798 /* 8 */, 1.61791440652881789386 /* 9 */};



static void 
bytecoder_output (int c)
{
  assert (c < MAXCODE);
  if (buffull == bufalloc)
    {
      char *newbuf = realloc (buffer, bufalloc * 2);
      if (newbuf == NULL)
	abort ();
      buffer = newbuf;
      bufalloc *= 2;
    }
  buffer[buffull++] = c;
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


/* If want_compress is non-zero, bytecoder will do dictionary compression,
   otherwise it just passes the literals through to the output buffer */

void
bytecoder_init (int want_compress)
{
  compress = want_compress;
  bufalloc = 10 * sizeof (char);
  buffer = malloc (bufalloc);
  buffull = 0;
  return;
}


void bytecoder_flush ()
{
  while (coder_nrstored > 0)
    coder_outputbest ();
}


unsigned int 
bytecoder_size ()
{
  return buffull;
}


void 
bytecoder_read (char *dst)
{
  memmove(dst, buffer, buffull);
  buffull = 0;
}

void
bytecoder_clear ()
{
  bytecoder_flush ();
  free (buffer);
  buffer = NULL;
  bufalloc = 0;
  buffull = 0;
  return;
}


void 
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


/***********************************************************************
   Generating Lucas chains with Montgomery's PRAC algorithm. Code taken 
   from GMP-ECM, mostly written by Paul Zimmermann, and slightly 
   modified here
************************************************************************/


/* Returns the cost of a PRAC chain with initial multiplier v.
   The cost of an addition is addcost, the cost a doubling is doublecost.
*/

static unsigned long
lucas_cost (const unsigned long n, const double v, const unsigned int addcost,
	    const unsigned int doublecost)
{
  unsigned long c, d, e, r;

  d = n;
  r = (unsigned long) ((double) d / v + 0.5);
  if (r >= n)
    return (addcost * n);
  d = n - r;
  e = 2 * r - n;
  c = doublecost + addcost; /* initial duplicate and final addition */
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
          c += 3U * addcost; /* 3 additions */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          c += addcost + doublecost; /* one addition, one duplicate */
        }
      else if (d <= 4 * e)
        { /* condition 3 */
          d -= e;
          c += addcost; /* one addition */
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          c += addcost + doublecost; /* one addition, one duplicate */
        } 
      /* now d+e is odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          c += addcost + doublecost; /* one addition, one duplicate */
        }
      /* now d is odd and e even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          c += 3U * addcost + doublecost; /* three additions, one duplicate */
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          c += 3U * addcost + doublecost; /* three additions, one duplicate */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          c += 3U * addcost + doublecost; /* three additions, one duplicate */
        }
      else /* necessarily e is even */
        { /* condition 9 */
          e /= 2;
          c += addcost + doublecost; /* one addition, one duplicate */
        }
    }
  
  return c;
}


/* Returns the cost of the cheapest Lucas chain for n found by PRAC 
   using the first m multiplers from prac_mul. If mul is not NULL,
   stores the index of the first multiplier that produced such a short
   chain */

unsigned long 
prac_best (double *mul, const unsigned long n, const int m_parm, 
	   const unsigned int addcost, const unsigned int doublecost)
{
  int i, bestmul, m = m_parm;
  unsigned long bestcost, c;

  if (m > PRAC_NR_MULTIPLIERS)
    m = PRAC_NR_MULTIPLIERS;
  bestcost = lucas_cost (n, prac_mul[0], addcost, doublecost);
  bestmul = 0;
  for (i = 1; i < m; i++)
    {
      c = lucas_cost (n, prac_mul[i], addcost, doublecost);
      if (c < bestcost)
	{
	  bestcost = c;
	  bestmul = i;
	}
    }
  
  if (mul != NULL)
    *mul = prac_mul[bestmul];

  return bestcost;
}


/* Write bytecode for an addition chain for k. */
void 
prac_bytecode (const unsigned long k, const unsigned int addcost, 
	       const unsigned int doublecost)
{
  unsigned long d, e, r;
  double m;
  
  if (k == 2)
    {
      bytecoder ((literal_t) 12);
      return;
    }
  
  assert (k % 2 == 1);
  
  /* Find the best multiplier for this k */
  prac_best (&m, k, PRAC_NR_MULTIPLIERS, addcost, doublecost);
  
  d = k;
  r = (unsigned long) ((double) d / m + 0.5);
  
  d = k - r;
  e = 2 * r - k;
  bytecoder ((literal_t) 10);
  
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
	  bytecoder ((literal_t) 0);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
	  bytecoder ((literal_t) 1);
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
	  bytecoder ((literal_t) 2);
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
	  bytecoder ((literal_t) 3);
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
	  bytecoder ((literal_t) 4);
        }
      else if (d % 2 == 0) /* d+e is now odd */
        { /* condition 5 */
          d /= 2;
	  bytecoder ((literal_t) 5);
        }
      else if (d % 3 == 0) /* d is odd, e even */
        { /* condition 6 */
          d = d / 3 - e;
	  bytecoder ((literal_t) 6);
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
	  bytecoder ((literal_t) 7);
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
	  bytecoder ((literal_t) 8);
        }
      else /* necessarily e is even */
        { /* condition 9: never happens? */
          e /= 2;
	  bytecoder ((literal_t) 9);
        }
    }
  
  bytecoder ((literal_t) 11);
  
  assert (d == 1);
}
