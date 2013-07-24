#ifndef RENUMBER_H_
#define RENUMBER_H_

#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "hashpair.h"
#include "cado_poly.h"
#include "rootfinder.h"
#include "timing.h"
#include "math.h"
#include "typedefs.h"

#define MAX_LOG_CACHED 20
#define RENUMBER_MAXLINE 1024
#define RENUMBER_MAX_ABOVE_BADIDEALS 8
#define RENUMBER_SPECIAL_VALUE ((p_r_values_t) -1)


struct __bad_ideals_t 
{
  int n; // number of p_r_values that correspond to more than one ideals.
  p_r_values_t * p;  // these p_r_values in two separate tables
  p_r_values_t * r;  // 
  int * side;
  int * nb;               // the number of ideals for each
};

struct __renumber_t
{
  FILE * file;           // file containing the renumbering table
  struct __bad_ideals_t bad_ideals;  // the bad ideals
  p_r_values_t * table;  //renumbering table
  index_t size;          //number of elements in the renumbering table
  uint8_t nb_bits;  // number of bits taken by an index in the file
                    // 32 or 64

  int rat;       // if one poly has degree 1, rat = 0 or 1 depending which one,
                 // else -1 if no rational side
                 // if rat = -1, we add p+1 to roots on side 1
                 // else we add p+1 to roots on rat side
  index_t *cached; // We cached the index for primes < 2^MAX_LOG_CACHED
  index_t first_not_cached;
  int add_full_col; //do we add a col of 1 to all relations
};
typedef struct __renumber_t renumber_t[1];

void renumber_init (renumber_t, cado_poly);
int renumber_is_bad(int *, index_t*,renumber_t, p_r_values_t, p_r_values_t, int);
void renumber_free (renumber_t);
void renumber_init_write (renumber_t, const char *, const char *, int);
void renumber_close_write (renumber_t, const char *);
void renumber_read_table (renumber_t, const char *);
void renumber_write_p (renumber_t, unsigned long, unsigned long * [2], int [2]);
index_t renumber_get_index_from_p_r (renumber_t, p_r_values_t, p_r_values_t,int);
void renumber_get_p_r_from_index (renumber_t, p_r_values_t *, p_r_values_t *,
                                                    int *, index_t, cado_poly);
//for DEBUG, should be remove later
void renumber_debug_print_tab (FILE *, const char *, cado_poly);
#endif
