#include "cado.h"
#include "renumber.h"
#include <ctype.h>

#define MAXLINE 1024
#define MAXBADS 20
static void
read_bad_ideals(struct __bad_ideals_t * bad, const char * filename)
{
  FILE *file = fopen(filename, "r");
  ASSERT_ALWAYS(file != NULL);
  char str[MAXLINE], *ret;
  bad->n = 0;

  // Remark: the parsing is done with long, whereas we could have 64-bit
  // exceptional ideals, and that would fail on 32-bit machines.
  // This is probably never occur, so we leave it like that.
  do {
      int c = getc(file);
      if (c == EOF) 
          break;
      if (c == '#') {
          // skip line
          c = ungetc(c, file);
          ASSERT_ALWAYS(c != EOF);
          fgets(str, MAXLINE, file);
          continue;
      }
      c = ungetc(c, file);
      ASSERT_ALWAYS(c != EOF);
      ret = fgets(str, MAXLINE, file);
      ASSERT_ALWAYS(ret != NULL);
      size_t n = strnlen(str, MAXLINE);  
      ASSERT_ALWAYS(n != MAXLINE);
      char * sstr;
      errno = 0;
      long p = strtol(str, &sstr, 10);
      ASSERT_ALWAYS(errno == 0);
      ASSERT_ALWAYS(sstr[0] == ',');
      sstr++;
      long r = strtol(sstr, &sstr, 10);
      ASSERT_ALWAYS(errno == 0);
      ASSERT_ALWAYS(sstr[0] == ':');
      sstr++;
      long s = strtol(sstr, &sstr, 10);
      ASSERT_ALWAYS(errno == 0);
      ASSERT_ALWAYS(sstr[0] == ':');
      sstr++;
      ASSERT_ALWAYS(s == 0 || s == 1);
      while (isspace(sstr[0]))
          sstr++;
      errno = 0;
      long nb = strtol(sstr, &sstr, 10);
      ASSERT_ALWAYS(errno == 0);

      // TODO:
      ASSERT_ALWAYS(s == 1);  // Two algebraic sides not implemented

      ASSERT_ALWAYS(bad->n < MAXBADS);
      bad->p[bad->n] = p;
      bad->r[bad->n] = r;
      bad->nb[bad->n] = nb;
      bad->n++;
  } while (1);
  fclose(file);
}
#undef MAXLINE

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
  renumber_info->bad_ideals.n = 0;
  renumber_info->bad_ideals.p = (p_r_values_t *) malloc(MAXBADS * 
          sizeof(p_r_values_t));
  renumber_info->bad_ideals.r = (p_r_values_t *) malloc(MAXBADS * 
          sizeof(p_r_values_t));
  renumber_info->bad_ideals.nb = (int *) malloc(MAXBADS*sizeof(int));
  ASSERT_ALWAYS(renumber_info->bad_ideals.p != NULL);
  ASSERT_ALWAYS(renumber_info->bad_ideals.r != NULL);
  ASSERT_ALWAYS(renumber_info->bad_ideals.nb != NULL);
}
#undef MAXBADS

void renumber_read_badideals (renumber_t renumber_info, const char * badfile)
{
  if (badfile != NULL) 
    read_bad_ideals(&renumber_info->bad_ideals, badfile);
  else 
    renumber_info->bad_ideals.n = 0;
}

void
renumber_print_info (FILE * f, renumber_t renumber_info)
{
  fprintf (f, "Renumbering struct: size=%"PRid", nb_bytes=%"PRIu8", "
              "sizeof(*table)=%zu, rat=%d\n", renumber_info->size,
              renumber_info->nb_bytes, sizeof(*(renumber_info->table)),
              renumber_info->rat);
}

void
renumber_free (renumber_t renumber_info)
{
  if (renumber_info->table != NULL)
    free(renumber_info->table);
  if (renumber_info->cached != NULL)
    free(renumber_info->cached);

  renumber_info->table = NULL;
  renumber_info->cached = NULL;

  if (renumber_info->bad_ideals.p != NULL)
    free(renumber_info->bad_ideals.p);
  if (renumber_info->bad_ideals.r != NULL)
    free(renumber_info->bad_ideals.r);
  if (renumber_info->bad_ideals.nb != NULL)
    free(renumber_info->bad_ideals.nb);

  renumber_info->bad_ideals.r = NULL;
  renumber_info->bad_ideals.p = NULL;
  renumber_info->bad_ideals.nb = NULL;
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
  // Write first the bad ideals information at the end of the file
  uint64_t n =  (uint64_t) renumber_info->bad_ideals.n;
  fwrite (&(n), sizeof(uint64_t), 1, renumber_info->file);
  for (unsigned int i = 0; i < n; ++i) {
    unsigned long p = renumber_info->bad_ideals.p[i];
    renumber_write_one(renumber_info, p);
    unsigned long r = renumber_info->bad_ideals.r[i];
    renumber_write_one(renumber_info, r);
    unsigned long nb = renumber_info->bad_ideals.nb[i];
    renumber_write_one(renumber_info, nb);
  }

  // Put the right value of size.
  rewind (renumber_info->file);
  fwrite (&(renumber_info->size), sizeof(renumber_info->size), 1,
          renumber_info->file);
  
  index_t offset = 0;
  for (int j = 0; j < renumber_info->bad_ideals.n; ++j)
    offset += renumber_info->bad_ideals.nb[j];
  fprintf (stderr, "Renumbering struct: nprimes=%"PRid"\n", 
          offset + renumber_info->size);

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
  p_r_values_t v, prev_v = 0;

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
  renumber_info->table = (index_t*) malloc(renumber_info->size*sizeof(index_t));
  ASSERT_ALWAYS (renumber_info->table != NULL);
  // Allocate the cached table
  renumber_info->cached = (index_t*) malloc((2<<MAX_LOG_CACHED)*sizeof(index_t));
  ASSERT_ALWAYS (renumber_info->cached != NULL);
  memset (renumber_info->cached, 0, (2 << MAX_LOG_CACHED)*sizeof(index_t));

  // Begin to read the renumbering table
  ret = 0;
  i = 0;

  // we cached value below 2^MAX_LOG_CACHED
  // Remark: the indices in the cache do not take into account 
  // the possible offset due to bad ideals.
  while (i < renumber_info->size)
  {
    ret += renumber_read_one(renumber_info,i);
    v = renumber_info->table[i];
    if ((v >> MAX_LOG_CACHED))
    {
      i++;
      break;
    }
    if (v > prev_v)
      renumber_info->cached[v] = i;

    prev_v = v;
    i++;
  }

  renumber_info->first_not_cached = i;

  for (; i < renumber_info->size; i++)
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

  // Read bad ideal information
  {
    uint64_t n;
    ret = fread (&(n), sizeof(uint64_t), 1, renumber_info->file);
    renumber_info->bad_ideals.n = n;
    for (unsigned int i = 0; i < n; ++i) {
      unsigned long p, r, nb;
      ret = fread (&p, renumber_info->nb_bytes, 1, renumber_info->file);
      ret = fread (&r, renumber_info->nb_bytes, 1, renumber_info->file);
      ret = fread (&nb, renumber_info->nb_bytes, 1, renumber_info->file);
      renumber_info->bad_ideals.p[i] = p;
      renumber_info->bad_ideals.r[i] = r;
      renumber_info->bad_ideals.nb[i] = nb;
    }
  }

  renumber_print_info (stderr, renumber_info);
  fprintf(stderr, "Renumbering struct: first_not_cached=%"PRid"\n",
                  renumber_info->first_not_cached);

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

int renumber_is_bad(renumber_t rn, p_r_values_t p, p_r_values_t r)
{
  int bad = 0;
  for (int i = 0; i < rn->bad_ideals.n; ++i)
    if (p == rn->bad_ideals.p[i] && r == rn->bad_ideals.r[i]) {
      bad = 1;
      break;
    }
  return bad;
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
      r[1][0] = (p << 1) + 1; // The largest roots become 2p+1
    else
      r[0][0] = (p << 1) + 1; // The largest roots become 2p+1

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
  index_t i;
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

  if ((p >> MAX_LOG_CACHED)) // p is not cached
  {
    index_t max = renumber_info->size - 1;
    index_t min = renumber_info->first_not_cached;

    float hint = 2.0 * (((float) p) / logf ((float) p));
    i = (index_t) hint;
    if (i < min)
      i = min;
    else if (i > max)
      i = max;

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
        if (i == old_i)
          i++;
      }
      else
      {
        max = old_i;
        i = (old_i + min)/2;
        if (i == old_i)
          i--;
      }
    }
  }
  else //p is cached
    i = renumber_info->cached[vp];

  /* now i points at the beginning of a decreasing sequence of values of vr */
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

  // In the case where there are bad_ideals, we shift the index by the 
  // number of such ideals
  index_t offset = 0;
  for (int j = 0; j < renumber_info->bad_ideals.n; ++j)
    offset += renumber_info->bad_ideals.nb[j];

  return offset + i;
}

// TODO should also return side in the case where there are two alg sides
void
renumber_get_p_r_from_index (renumber_t renumber_info, p_r_values_t *p,
                             p_r_values_t * r, index_t i, cado_poly pol)
{
  index_t j;
  p_r_values_t *tab = renumber_info->table;
  
  // In the case where there are bad_ideals, we shift the index by the 
  // number of such ideals
  index_t offset = 0;
  for (int j = 0; j < renumber_info->bad_ideals.n; ++j)
    offset += renumber_info->bad_ideals.nb[j];
  i -= offset;

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

  // In the case where there are bad_ideals, we shift the index by the 
  // number of such ideals
  index_t offset = 0;
  for (int j = 0; j < tab->bad_ideals.n; ++j)
    offset += tab->bad_ideals.nb[j];

  for (i = 0; i < tab->size; i++)
  {
    renumber_get_p_r_from_index (tab, &p, &r, i+offset, pol);
    fprintf (output, "i=%u\ttab[i]=%u\tp=%u\t", i, tab->table[i], p);
    if (r == p + 1)
      fprintf (output, "rat side\n");
    else if (r == p)
      fprintf (output, "alg side proj\n");
    else
      fprintf (output, "alg side r=%u \n", r);
  }

  if (tab->bad_ideals.n != 0) {
    fprintf (output, "Bad ideals:\n");
    for (int i = 0; i < tab->bad_ideals.n; ++i) {
      unsigned long p = tab->bad_ideals.p[i];
      unsigned long r = tab->bad_ideals.r[i];
      unsigned long nb = tab->bad_ideals.nb[i];
      fprintf(output, "p=%lu r=%lu nb=%lu\n", p, r, nb);
    }
  }

  renumber_free(tab);
}
