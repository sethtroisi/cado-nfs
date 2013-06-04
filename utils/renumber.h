#ifndef RENUMBER_H_
#define RENUMBER_H_

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "hashpair.h"
#include "cado_poly.h"
#include "rootfinder.h"
#include "misc.h"
#include "timing.h"
#include "math.h"

struct __renumber_t
{
  FILE * file;           // file containing the renumbering table
  p_r_values_t * table;  //renumbering table
  index_t size;          //number of elements in the renumbering table
  uint8_t nb_bytes;  // number of bytes taken by an index in the file
  
  int rat;       // if one poly has degree 1, rat = 0 or 1 depending which one, 
                 // else -1 if no rational side
                 // if rat = -1, we add p+1 to roots on side 1 
                 // else we add p+1 to roots on rat side
};
typedef struct __renumber_t renumber_t[1];

#define renumber_write_one(r,a) fwrite (&(a), (r)->nb_bytes, 1, (r)->file)
#define renumber_read_one(r,i) fread(&(r->table[i]), r->nb_bytes, 1, r->file)

void renumber_init (renumber_t, cado_poly);
void renumber_print_info (FILE *, renumber_t);
void renumber_free (renumber_t);
void renumber_init_write (renumber_t, const char *);
void renumber_close_write (renumber_t);
void renumber_read_table (renumber_t, const char *);
void renumber_write_p (renumber_t, unsigned long, unsigned long * [2], int [2]);
index_t renumber_get_index_from_p_r (renumber_t, p_r_values_t, p_r_values_t,int);
void renumber_get_p_r_from_index (renumber_t, p_r_values_t *, p_r_values_t *,
                                                            index_t, cado_poly);
//for DEBUG, should be remove later 
void renumber_debug_print_tab (FILE *, const char *, cado_poly);
#endif
