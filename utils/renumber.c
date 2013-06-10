#include "cado.h"
#include "renumber.h"



void 
renumber_init (renumber_t renumber_info, cado_poly pol)
{
  int max_nb_bits = MAX(pol->pols[0]->lpb, pol->pols[1]->lpb);

  if (pol->pols[0]->degree != 1 && pol->pols[1]->degree != 1)
    renumber_info->rat = -1;
  else if (pol->pols[0]->degree == 1)
    renumber_info->rat = 0;
  else
    renumber_info->rat = 1;

  if (renumber_info->rat == -1)
    max_nb_bits++;

  if (max_nb_bits <= 32)
    renumber_info->nb_bytes = 4;
  else
    renumber_info->nb_bytes = (max_nb_bits / 8) + 1;
  ASSERT_ALWAYS (renumber_info->nb_bytes <= sizeof(p_r_values_t));


  renumber_info->size = 0;
  renumber_info->table = NULL;
  renumber_info->file = NULL;
}

void
renumber_print_info (FILE * f, renumber_t renumber_info)
{
  fprintf (f, "Renumbering struct: size=%"PRid", nb_bytes=%"PRIu8"," 
              "sizeof(*table)=%zu, rat=%d\n", renumber_info->size,
              renumber_info->nb_bytes, sizeof(*(renumber_info->table)),
              renumber_info->rat);
}

void 
renumber_free (renumber_t renumber_info)
{
  if (renumber_info->table != NULL)
  {
    free(renumber_info->table);
    renumber_info->table = NULL;
  }
}

/* The renumber_t struct MUST have been initialized before */
void
renumber_init_write (renumber_t renumber_info, const char * filename)
{
  fprintf (stderr, "Opening %s to write the renumbering table\n", filename);
  renumber_info->file = fopen(filename, "wb");
  ASSERT_ALWAYS(renumber_info->file != NULL);

  uint64_t tmp_size = (uint64_t) renumber_info->size;
  fwrite (&(tmp_size), sizeof(uint64_t), 1, renumber_info->file);
  fwrite (&(renumber_info->nb_bytes), sizeof(renumber_info->nb_bytes), 1, 
          renumber_info->file);

  renumber_print_info (stderr, renumber_info);
}

void
renumber_close_write (renumber_t renumber_info)
{
  rewind (renumber_info->file);
  fwrite (&(renumber_info->size), sizeof(renumber_info->size), 1, 
          renumber_info->file);
  fprintf (stderr, "Renumbering struct: size=%"PRid"\n", renumber_info->size);
  
  fclose(renumber_info->file);
}


/* The renumber_t struct MUST have been initialized before */
/* renumber_free MUST be called to free the table afterwards */
void
renumber_read_table (renumber_t renumber_info, const char * filename)
{
  index_t i, ret, report = 0;
  double dt, t = wct_seconds ();
  double mb_s = 0;
  uint8_t old_nb_bytes = renumber_info->nb_bytes;

  //open file for reading
  fprintf (stderr, "Opening %s to read the renumbering table\n", filename);
  renumber_info->file = fopen(filename, "rb");
  ASSERT_ALWAYS(renumber_info->file != NULL);

  // read size of renumbering table
  uint64_t tmp_size = (uint64_t) renumber_info->size;
  ret = fread (&(tmp_size), sizeof(uint64_t), 1, renumber_info->file);
  ASSERT(ret == 1);
  renumber_info->size = (index_t) tmp_size;
  ASSERT_ALWAYS (renumber_info->size != 0);
  ASSERT_ALWAYS ((uint64_t) renumber_info->size == tmp_size);

  // read nb_bytes
  ret = fread (&(renumber_info->nb_bytes), sizeof(renumber_info->nb_bytes), 1, 
             renumber_info->file);
  ASSERT(ret == 1);
  ASSERT_ALWAYS (renumber_info->nb_bytes <= sizeof(p_r_values_t));
  ASSERT_ALWAYS (renumber_info->nb_bytes >= 4);
  if (old_nb_bytes != renumber_info->nb_bytes)
  {
    fprintf (stderr, "Warning, computed value of nb_bytes (%d) is different "
                     "from the read value of nb_bytes (%d).\n", old_nb_bytes,
                     renumber_info->nb_bytes);
  }

  // Allocate the renumbering table
  SMALLOC (renumber_info->table, renumber_info->size, "Renumbering table");

  renumber_print_info (stderr, renumber_info);

  // Begin to read the renumbering table
  ret = 0;
  for (i = 0; i < renumber_info->size; i++)
  {
    ret += renumber_read_one(renumber_info,i);

    if ((ret >> 23) != (report >> 23))
    {
      report = ret;
      dt = wct_seconds () - t;
      if (dt > 0.01)
        mb_s = ((double) (ret * renumber_info->nb_bytes) / (double) dt) * 1.0e-6;
      else
        mb_s = 0.0;
      fprintf(stderr, "Renumbering table: read %"PRid" values from file "
                      "in %.1fs -- %.1f MB/s\n", ret, dt, mb_s);
    }
  }
  
  dt = wct_seconds () - t;
  if (dt > 0.01)
    mb_s = ((double) (ret * renumber_info->nb_bytes) / (double) dt) * 1.0e-6;
  else
    mb_s = 0.0;
  fprintf(stderr, "Renumbering table: end of read. Read %"PRid" values from "
                  "file in %.1fs -- %.1f MB/s\n", ret, dt, mb_s);
  ASSERT_ALWAYS(ret == renumber_info->size);

  fclose(renumber_info->file);
}

/* sort in decreasing order */
static void
sort (unsigned long *r, int n)
{
  int i, j;
  unsigned long v;

  for (i = 1; i < n; i++)
  {
    v = r[i];
    for (j = i; j > 0 && v > r[j - 1]; j--)
      r[j] = r[j - 1];
    r[j] = v;
  }
}


void 
renumber_write_p (renumber_t renumber_info, unsigned long p, unsigned long*r[2], 
                  int k[2])
{
  unsigned long add_value = p+1; 
  int i;

  // We sort roots by decreasing order
  sort(r[0], k[0]);
  sort(r[1], k[1]);    

  if (renumber_info->rat == -1) // two algebraic sides
  {
    // Add p + 1 on side 1
    for (i = 0; i < k[1]; i++)
      r[1][i] += add_value;

    if (k[1] != 0)
      r[1][0] = (p < 1) + 1; // The largest roots become 2p+1
    else
      r[0][0] = (p < 1) + 1; // The largest roots become 2p+1

    for (i = 0; i < k[1]; i++)
      renumber_write_one(renumber_info, r[1][i]);
    for (i = 0; i < k[0]; i++)
      renumber_write_one(renumber_info, r[0][i]);

    renumber_info->size += k[0] + k[1];
  }
  else
  {
    int rat = renumber_info->rat, alg = 1 - rat;

    if(k[rat] == 0) // lpbr < p < lpba
    {
      r[alg][0] = add_value;// We put p+1 instead of the largest roots
    }
    else // p < lpbr and lpba <= lpbr
    {
      ASSERT_ALWAYS (k[rat] == 1);
    
      //If there is a rational side, we put p+1 in the renumbering table.
      renumber_write_one(renumber_info, add_value);
      renumber_info->size++;
    }

    for (i = 0; i < k[alg]; i++)
      renumber_write_one(renumber_info, r[alg][i]);

    renumber_info->size += k[alg];
  }
}

/* side is 0 if (p,r) corresponds to the left part in the relation,
   side is 1 if (p,r) corresponds to the right part.
   For NFS the rational part (if any) corresponds to side 0, and we can
   give any value of r. */
index_t 
renumber_get_index_from_p_r (renumber_t renumber_info, p_r_values_t p, 
                             p_r_values_t r, int side)
{
  float hint = 2.0 * (((float) p) / logf ((float) p));
  index_t i = (index_t) hint;
  // TODO Better int for small prime (up to 2^16?) store the i corresponding to
  // vp in redirection table

  index_t max = renumber_info->size - 1, min = 0;

  p_r_values_t *tab = renumber_info->table;
  p_r_values_t vr, vp; /* values of r and p as they are stored in the table*/


  if (renumber_info->rat == -1)
  {
    vp = (p << 1) + 1;
    vr = (side == 1) ? (p + 1 + r) : r; 
  }
  else
  {
    vp = p + 1;
    vr = (side == renumber_info->rat) ? vp : r; 
  }
  
  /* Looking for vp: the values of vp are ordered in increasing order and are
     always at the beginning of a decreasing sequence */
  while (1)
  {
    index_t old_i = i;

    while (i > 0 && tab[i-1] > tab[i])
      i--;

    if (tab[i] == vp)
      break;
    else if (tab[i] < vp)
    {
      min = old_i;
      i = (old_i + max)/2;
    }
    else
    {
      max = old_i;
      i = (old_i + min)/2;
    }
  }

  if (side != renumber_info->rat)
  {
    while (i < renumber_info->size)
    {
      if (vr < tab[i])
      {
        if (i == renumber_info->size || vr > tab[i+1])
          break;
        else
          i++;
      }
      else
        break;
    }
  }

  return i;
}

// TODO should also return side in the case where there are two alg sides
void
renumber_get_p_r_from_index (renumber_t renumber_info, p_r_values_t *p, 
                             p_r_values_t * r, index_t i, cado_poly pol)
{
  index_t j;
  p_r_values_t *tab = renumber_info->table;

  for (j = i; j>0 && tab[j-1] > tab[j];)
    j--;

  if (renumber_info->rat == -1)
  {
    *p = (tab[j] - 1) > 1;
    if (i == j)
    {
      // TODO Case with two alg sides
      // If nb of roots over p pol[1] is 0 then
      // *r = largest roots of pol[0]
      //else
      // *r = largest roots of pol[1]
    }
    else
    {
      if (tab[i] <= *p)
        *r = tab[i];
      else
        *r = tab[i] - *p - 1; 
    }
  }
  else
  {
    *p = tab[j] - 1;
    if (*p > (1UL << pol->rat->lpb) && i == j)
    {
      // Case where there is only alg side (p >= lpbr) and we are on the largest
      // root on alg side (i == j)
      int k, d = pol->alg->degree;
      unsigned long *roots;

      // if there is a proj root, this the largest (r = p by convention)
      if (mpz_divisible_ui_p (pol->alg->f[d], *p)) 
        *r = *p;
      else
      {
        roots = (unsigned long*) malloc (d * sizeof (unsigned long));
        k = poly_roots_ulong(roots, pol->alg->f, d, *p);
        ASSERT_ALWAYS (k > 0);
        sort(roots, k);
        *r = roots[0];
      }
    }
    else
      *r = tab[i];
  }
}

//for DEBUG, should be remove later 
void renumber_debug_print_tab (FILE *output, const char *filename, 
                                    cado_poly pol)
{
  renumber_t tab;
  index_t i;
  p_r_values_t p, r;

  renumber_init (tab, pol);
  renumber_read_table (tab, filename);
    
  for (i = 0; i < tab->size; i++)
  {
    renumber_get_p_r_from_index (tab, &p, &r, i, pol);
    fprintf (output, "i=%u\ttab[i]=%u\tp=%u\t", i, tab->table[i], p);
    if (r == p + 1)
      fprintf (output, "rat side\n");
    else if (r == p)
      fprintf (output, "alg side proj\n");
    else
      fprintf (output, "alg side r=%u \n", r);
  }
  
  
  renumber_free(tab);
}
