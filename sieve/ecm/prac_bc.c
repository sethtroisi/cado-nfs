#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "prac_bc.h"


/* Table of multipliers for PRAC. prac_mul[i], 1<=<i<=9, has continued 
   fraction sequence of all ones but with a 2 in the (i+1)-st place, 
   prac_mul[i], 10<=<i<=17, has continued fraction sequence of all ones 
   but with a 2 in second and in the (i-7)-th place
   and prac_mul[0] is all ones, i.e. the golden ratio. */
static const double prac_mul[PRAC_NR_MULTIPLIERS] = 
  {1.61803398874989484820 /* 0 */, 1.38196601125010515179 /* 1 */, 
   1.72360679774997896964 /* 2 */, 1.58017872829546410471 /* 3 */, 
   1.63283980608870628543 /* 4 */, 1.61242994950949500192 /* 5 */,
   1.62018198080741576482 /* 6 */, 1.61721461653440386266 /* 7 */, 
   1.61834711965622805798 /* 8 */, 1.61791440652881789386 /* 9 */
#if 0
   ,
   1.41982127170453589529, 1.36716019391129371457,
   1.38757005049050499808, 1.37981801919258423518,
   1.38278538346559613734, 1.38165288034377194202,
   1.38208559347118210614, 1.38192033153010418805,
   1.70431400059211079746};
#else
   };
#endif


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


/* Produce a PRAC chain with initial multiplier v. 
   Returns its arithmetic cost.
   The cost of an addition is addcost, the cost a doubling is doublecost. */

static double
prac_chain (const unsigned long n, const double v, const double addcost,
	    const double doublecost, bc_state_t *state)
{
  unsigned long d, e, r;
  double cost;

  d = n;
  r = (unsigned long) ((double) d / v + 0.5);
  if (r >= n)
    return (addcost * n);
  d = n - r;
  e = 2 * r - n;
  cost = doublecost + addcost; /* initial doubling and final addition */

  bytecoder ((literal_t) 10, state); /* initial doubling (subchain init) */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
	  if (state != NULL) bytecoder ((literal_t) 0, state);
        }
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          cost += 3 * addcost; /* 3 additions */
	  if (state != NULL) bytecoder ((literal_t) 1, state);
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          cost += addcost + doublecost; /* one addition, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 2, state);
        }
      else if (d <= 4 * e)
        { /* condition 3 */
          d -= e;
          cost += addcost; /* one addition */
	  if (state != NULL) bytecoder ((literal_t) 3, state);
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          cost += addcost + doublecost; /* one addition, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 4, state);
        } 
      /* now d+e is odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          cost += addcost + doublecost; /* one addition, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 5, state);
        }
      /* now d is odd and e even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          cost += 3 * addcost + doublecost; /* three additions, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 6, state);
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          cost += 3 * addcost + doublecost; /* three additions, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 7, state);
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          cost += 3 * addcost + doublecost; /* three additions, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 8, state);
        }
      else /* necessarily e is even */
        { /* condition 9 */
          e /= 2;
          cost += addcost + doublecost; /* one addition, one doubling */
	  if (state != NULL) bytecoder ((literal_t) 9, state);
        }
    }
  bytecoder ((literal_t) 11, state); /* final add (subchain end) */
  
  /* Here d = gcd(n, r). If d != 1, then the sequence cannot be
     reduced below d, i.e., the chain does not start at 1. 
     It would be the end of a concatenated chain instead.
     Return 0 in this case. */
  return (d == 1) ? cost : 0.;
}


/* Returns the cost of the cheapest Lucas chain for n found by PRAC 
   using the first m multiplers from prac_mul. If mul is not NULL,
   stores the index of the first multiplier that produced such a cheap
   chain. The cost of an addition, doubling, code byte and code byte change
   are passed via parameters.
   If no valid chain could be found because n was composite and 
   no starting value was coprime to n, return 0. */

double 
prac_best (double *mul, const unsigned long n, const int m_parm, 
	   const double addcost, const double doublecost,
	   const double bytecost, const double changecost,
	   const bc_dict_t *dict)
{
  int i, bestmul = 0, m = m_parm;
  double bestcost = 0., cost;
  bc_state_t *state;

  if (m > PRAC_NR_MULTIPLIERS)
    m = PRAC_NR_MULTIPLIERS;

  state = bytecoder_init (dict);

  for (i = 0; i < m; i++)
    {
      cost = prac_chain (n, prac_mul[i], addcost, doublecost, state);
      if (cost > 0.)
	{
	  size_t j;
	  
	  /* Add the cost for byte codes and byte code changes */
	  bytecoder_flush (state);
	  cost += bytecost * state->buffull;

	  for (j = 1; j < state->buffull; j++)
	    if (state->buffer[j-1] != state->buffer[j])
	      cost += changecost;

	  if (bestcost == 0 || cost < bestcost)
	    {
	      bestcost = cost;
	      bestmul = i;
	    }
	  state->buffull = 0;
	}
    }
  
  if (bestcost > 0. && mul != NULL)
    *mul = prac_mul[bestmul];

  bytecoder_clear (state);

  return bestcost;
}


/* Write bytecode for an addition chain for odd k, and return its cost */
double  
prac_bytecode (const unsigned long k, const double addcost, 
	       const double doublecost, const double bytecost, 
	       const double changecost, bc_state_t *state)
{
  double d, m = 0.;
  
  assert (k % 2 == 1);
  
  /* Find the best multiplier for this k */
  d = prac_best (&m, k, PRAC_NR_MULTIPLIERS, addcost, doublecost, bytecost, 
		 changecost, state->dict);
  if (d == 0.)
    {
      /* This k is composite and prac cannot make a valid chain for it.
	 We could try to factor k and make a composite chain from the
	 prime factors. For now, we bail out - caller shouldn't give
	 a composite k that doesn't have a prac chain we can find. */
      abort ();
    }

  /* Return the cost of modular arithmetic, not including
     length of compressed code or number of code changes */
  d = prac_chain (k, m, addcost, doublecost, state);
  return d;
}
