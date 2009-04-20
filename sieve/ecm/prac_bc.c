#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "prac_bc.h"

/* Table of multipliers for PRAC. prac_mul[i], i>0, has continued fraction 
   sequence of all ones but with a 2 in the i-th place, and 
   prac_mul[0] is all ones, i.e. the golden ratio */
static const double prac_mul[PRAC_NR_MULTIPLIERS] = 
  {1.61803398874989484820 /* 0 */, 1.38196601125010515179 /* 1 */, 
   1.72360679774997896964 /* 2 */, 1.58017872829546410471 /* 3 */, 
   1.63283980608870628543 /* 4 */, 1.61242994950949500192 /* 5 */,
   1.62018198080741576482 /* 6 */, 1.61721461653440386266 /* 7 */, 
   1.61834711965622805798 /* 8 */, 1.61791440652881789386 /* 9 */};


static code_t
bytecoder_literal_to_code (const literal_t l)
{
  /* Maybe a non-trivial mapping is desired sometime */
  return (code_t) l;
}

/* Output the code c to the output buffer. Reallocate memory for the 
   output buffer if necessary. */
static void 
bytecoder_output (const code_t c, bc_state_t *state)
{
  assert (c < BC_MAXCODE);
  if (state->buffull == state->bufalloc)
    {
      assert (state->bufalloc > 0);
      state->bufalloc *= 2;
      code_t *newbuf = realloc (state->buffer, 
				state->bufalloc * sizeof (code_t));
      if (newbuf == NULL)
	abort ();
      state->buffer = newbuf;
    }
  state->buffer[state->buffull++] = c;
}


/* Returns the length of the matching part of "buf" (with "len" valid bytes)
   and the dictionary entry "entry" */

static size_t 
dict_matchlen (const bc_state_t *state, const int entry)
{
  size_t i;
  const bc_dict_t *dict = state->dict;

  assert (0 <= entry && entry < dict->nr_entries);

  for (i = 0; i < state->nrstored && i < dict->key_len[entry]; i++)
    if (state->history[i] != dict->key[entry][i])
      break;
  
  return i;
}


/* Returns the index of the dictionary key that has the longest match with
   the current history window. If "partial" is != 0, then matches that 
   don't match the complete dictionary key are accepted,
   otherwise only matches with complete dictionary keys are accepted.
   Returns -1 if nothing matched at all.
   If len is non-NULL, length of the match found is stored in *len */

static int
dict_longestmatch (size_t *len, const bc_state_t *state, const int partial)
{
  int i, matchidx = -1;
  size_t t, matchlen = 0;
  const bc_dict_t *dict = state->dict;
  
  for (i = 0; i < dict->nr_entries; i++)
    {
      t = dict_matchlen (state, i);
      if (t > matchlen && (partial || t == dict->key_len[i]))
	{
	  matchlen = t;
	  matchidx = i;
	}
    }
  
  if (len != NULL)
    *len = matchlen;
  
  return matchidx;
}


/* Adds a literal to the history window */

static void 
coder_histadd (const literal_t c, bc_state_t *state)
{
  state->history[state->nrstored++] = c;
  assert (state->nrstored <= state->histlen);
}


/* Removes n literals from the history window */

static void
coder_histremove (const size_t n, bc_state_t *state)
{
  size_t i;
  if (n > state->nrstored)
    abort();

  for (i = 0; i + n < state->nrstored; i++)
    state->history[i] = state->history[i + n];

  state->nrstored -= n;
}


/* Output the best dictionary match we have at the moment, or a literal
   if history[] matches no complete dictionary entry */

static void
coder_outputbest (bc_state_t *state)
{
  int best;
  size_t bestlen;
  const bc_dict_t *dict = state->dict;
  
  assert (state->nrstored > 0);
  
  /* Find longest complete match */
  best = dict_longestmatch (&bestlen, state, 0);
  
  /* If there was a complete dictionary match, output that, 
     otherwise output the first entry in history[] as a literal */
  if (best != -1)
    {
      /* Dictonary code 0 is special and means: no output */
      if (dict->code[best] != 0)
	bytecoder_output (bytecoder_literal_to_code (dict->code[best]), state);
      coder_histremove (bestlen, state);
    }
  else
    {
      bytecoder_output (bytecoder_literal_to_code (state->history[0]), state);
      coder_histremove (1, state);
    }
}


/* If dict is non-NULL, bytecoder will do dictionary compression,
   otherwise it just passes the literals through to the output buffer */

bc_state_t *
bytecoder_init (const bc_dict_t *dict)
{
  int i;
  bc_state_t *state;
  state = malloc (sizeof(bc_state_t));

  if (dict != NULL)
    {
      /* Window size is equal to longest key size + 1 */
      state->histlen = 1;
      for (i = 0; i < dict->nr_entries; i++)
        if (state->histlen < dict->key_len[i] + 1)
          state->histlen = dict->key_len[i] + 1;
      state->history = malloc(state->histlen * sizeof(literal_t));
      if (state->history == NULL)
        abort ();
      state->dict = dict;
    }
  else
    {
      state->histlen = 0;
      state->history = NULL;
      state->dict = NULL;
    }
  state->nrstored = 0;
  state->bufalloc = 16;
  state->buffer = malloc (state->bufalloc * sizeof (code_t));
  if (state->buffer == NULL)
    abort ();
  state->buffull = 0;

  return state;
}


void bytecoder_flush (bc_state_t *state)
{
  while (state->nrstored > 0)
    coder_outputbest (state);
}


/* Return the size in bytes of the data in the output buffer */
size_t 
bytecoder_size (const bc_state_t *state)
{
  return state->buffull * sizeof (code_t);
}


/* Store the data from the output buffer in "dst" and set buffer to empty */
void 
bytecoder_read (code_t *dst, bc_state_t *state)
{
  memmove(dst, state->buffer, state->buffull * sizeof (code_t));
  state->buffull = 0;
}

void
bytecoder_clear (bc_state_t *state)
{
  bytecoder_flush (state);
  free (state->history);
  state->history = NULL;
  state->histlen = 0;
  free (state->buffer);
  state->buffer = NULL;
  state->bufalloc = 0;
  state->buffull = 0;
  state->dict = NULL;
  free (state);
  return;
}


void 
bytecoder (const literal_t c, bc_state_t *state)
{
  /* At this point, the first nrstored literals (possibly 0!) of
     history[] agree with some dictionary entry */
  
  if (state->dict != NULL)
    {
      /* Now add the new literal */
      coder_histadd (c, state);
      
      /* See if history[] still matches one of the dictionary entries */
      
      int best;
      size_t bestlen = 0;
      
      /* If all history[] matches some dictionary entry (possibly partially), 
         we do nothing */
      /* Otherwise, we repeatedly output the longest complete dictionary 
         match or literal code until all of history[] matches some 
         dictionary entry again, or is empty */
      
      while (1)
	{
	  best = dict_longestmatch (&bestlen, state, 1);
	  if (bestlen == state->nrstored)
	    break;
	  coder_outputbest (state);
	}
    }
  else
    {
      /* No compression: simply forward the literal to output */
      bytecoder_output (bytecoder_literal_to_code (c), state);
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
   stores the first multiplier that produced such a short chain */

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
	       const unsigned int doublecost, bc_state_t *state)
{
  unsigned long d, e, r;
  double m;
  
  if (k == 2)
    {
      bytecoder ((literal_t) 12, state);
      return;
    }
  
  assert (k % 2 == 1);
  
  /* Find the best multiplier for this k */
  prac_best (&m, k, PRAC_NR_MULTIPLIERS, addcost, doublecost);
  
  d = k;
  r = (unsigned long) ((double) d / m + 0.5);
  
  d = k - r;
  e = 2 * r - k;
  bytecoder ((literal_t) 10, state);
  
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
	  bytecoder ((literal_t) 0, state);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
	  bytecoder ((literal_t) 1, state);
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
	  bytecoder ((literal_t) 2, state);
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
	  bytecoder ((literal_t) 3, state);
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
	  bytecoder ((literal_t) 4, state);
        }
      else if (d % 2 == 0) /* d+e is now odd */
        { /* condition 5 */
          d /= 2;
	  bytecoder ((literal_t) 5, state);
        }
      else if (d % 3 == 0) /* d is odd, e even */
        { /* condition 6 */
          d = d / 3 - e;
	  bytecoder ((literal_t) 6, state);
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
	  bytecoder ((literal_t) 7, state);
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
	  bytecoder ((literal_t) 8, state);
        }
      else /* necessarily e is even */
        { /* condition 9: never happens? */
          e /= 2;
	  bytecoder ((literal_t) 9, state);
        }
    }
  
  bytecoder ((literal_t) 11, state);
  
  assert (d == 1);
}
