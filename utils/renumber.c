#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <inttypes.h>
#include <stdio.h>
#include "portability.h"
#include "renumber.h"
#include "gzip.h" /* for fopen_maybe_compress */
#include <ctype.h>


/********************** internal functions *****************************/

static const int ugly[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

static int
copy_bad_ideals(const char * filename, FILE *out, index_t *size)
{
  FILE *file = fopen(filename, "r"); /* never compressed, always small */
  ASSERT_ALWAYS(file != NULL);
  char str[RENUMBER_MAXLINE], *ret;
  int n = 0, nb, count, t;

  do
  {
    int c = getc(file);
    if (c == EOF)
      break;
    if (c == '#')
    {
      // skip line
      c = ungetc(c, file);
      ASSERT_ALWAYS(c != EOF);
      ret = fgets(str, RENUMBER_MAXLINE, file);
      continue;
    }
    c = ungetc(c, file);
    ASSERT_ALWAYS(c != EOF);

    ret = fgets(str, RENUMBER_MAXLINE, file);
    ASSERT_ALWAYS(ret != NULL);
    size_t check = strnlen(str, RENUMBER_MAXLINE);
    ASSERT_ALWAYS(check != RENUMBER_MAXLINE);

    const char *ptr = str;
    for (count = 0 ; *ptr != '\n'; ptr++)
      if (*ptr == ':')
      {
        if (count == 0)
          count = 1;
        else
          break;
      }

    ASSERT_ALWAYS(*ptr == ':');
    ptr++;
    while (isspace(ptr[0]))
      ptr++;
    for (nb = 0 ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
      nb = (nb * 10) + t;
    ASSERT_ALWAYS(*ptr == '\n');
    *size += nb;

    int r = fputs (str, out);
    ASSERT_ALWAYS (r >= 0);
    n++;
  } while (1);

  fclose(file);

  return n;
}

static void
parse_one_line_bad_ideals(struct __bad_ideals_t * bad, const char * str, int k)
{
  int t;
  const char *ptr = str;
  p_r_values_t p = 0, r = 0;
  int side = 0, nb = 0;

  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    p = (p << 4) + t;
  ASSERT_ALWAYS(*ptr == ',');
  ptr++;

  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    r = (r << 4) + t;
  ASSERT_ALWAYS(*ptr == ':');
  ptr++;

  side = ugly[(unsigned char) *ptr];
  ptr++;
  ASSERT_ALWAYS(side == 0 || side == 1);
  ASSERT_ALWAYS(*ptr == ':');

  ptr++;
  while (isspace(ptr[0]))
    ptr++;
  for ( ; (t=ugly[(unsigned char) *ptr]) >= 0; ptr++)
    nb = (nb * 10) + t;
  ASSERT_ALWAYS(*ptr == '\n');
  ASSERT_ALWAYS(nb <= RENUMBER_MAX_ABOVE_BADIDEALS);

  bad->p[k] = p;
  bad->r[k] = r;
  bad->nb[k] = nb;
  bad->side[k] = side;
}

static inline p_r_values_t
parse_one_line (char * str)
{
  p_r_values_t v = 0;
  int t;
  char *p;

  for (p = str; (t=ugly[(unsigned char) *p]) >= 0; p++)
    v = (v << 4) + t;

  ASSERT_ALWAYS(*p == '\n');
  return v;
}

static void
print_info (FILE * f, renumber_t renumber_info)
{
  fprintf (f, "Renumbering struct: nb_bits=%"PRIu8", sizeof(*table)=%zu, "
              "rat=%d nb_badideals=%d add_full_col=%d\n",
              renumber_info->nb_bits, sizeof(*(renumber_info->table)),
              renumber_info->rat, renumber_info->bad_ideals.n,
              renumber_info->add_full_col);
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

static void
renumber_alloc (renumber_t tab)
{
  // Allocate the renumbering table
  tab->table = (p_r_values_t *) malloc (tab->size * sizeof(p_r_values_t));
  ASSERT_ALWAYS (tab->table != NULL);
  // Allocate the cached table
  tab->cached = (index_t*) malloc ((2<<MAX_LOG_CACHED) * sizeof(index_t));
  ASSERT_ALWAYS (tab->cached != NULL);
  memset (tab->cached, 0, (2 << MAX_LOG_CACHED) * sizeof(index_t));
  // Allocate memory for bad ideals
  tab->bad_ideals.p = (p_r_values_t *) malloc(tab->bad_ideals.n *
                                                    sizeof(p_r_values_t));
  tab->bad_ideals.r = (p_r_values_t *) malloc(tab->bad_ideals.n *
                                                    sizeof(p_r_values_t));
  tab->bad_ideals.side = (int *) malloc(tab->bad_ideals.n * sizeof(int));
  tab->bad_ideals.nb = (int *) malloc(tab->bad_ideals.n * sizeof(int));
  ASSERT_ALWAYS(tab->bad_ideals.p != NULL);
  ASSERT_ALWAYS(tab->bad_ideals.r != NULL);
  ASSERT_ALWAYS(tab->bad_ideals.nb != NULL);
  ASSERT_ALWAYS(tab->bad_ideals.side != NULL);

}

static inline size_t
get_one_line (FILE *f, char *s)
{
  char *rets;
  size_t n;
  rets = fgets(s, RENUMBER_MAXLINE, f);
  ASSERT_ALWAYS(rets != NULL);
  n = strnlen(s, RENUMBER_MAXLINE);
  ASSERT_ALWAYS(n != RENUMBER_MAXLINE);
  return n;
}

static void
report_read (double t, index_t nread, size_t bytes_read, int end)
{
  double mb_s, dt = wct_seconds () - t;
  if (dt > 0.01)
    mb_s = (double) ((double) bytes_read / (double) dt) * 1.0e-6;
  else
    mb_s = 0.0;

  if (!end)
    fprintf(stderr, "Renumbering table: read %"PRid" values from file "
                    "in %.1fs -- %.1f MB/s\n", nread, dt, mb_s);
  else
    fprintf(stderr, "Renumbering table: end of read. Read %"PRid" values from "
                    "file in %.1fs -- %.1f MB/s\n", nread, dt, mb_s);

}

/* return zero if no roots mod p, else non-zero */
static int
get_largest_root_mod_p (p_r_values_t *r, mpz_t *pol, int deg, p_r_values_t p)
{
  int k;
  unsigned long *roots = NULL;

  // if there is a proj root, this the largest (r = p by convention)
  if (mpz_divisible_ui_p (pol[deg], p))
  {
    *r = p;
    return 1;
  }

  roots = (unsigned long*) malloc (deg * sizeof (unsigned long));
  k = poly_roots_ulong(roots, pol, deg, p);
  if (k > 0)
  {
    sort(roots, k);
    *r = roots[0];
    free(roots);
    return 2;
  }

  free(roots);
  return 0;
}

/*********************** End internal functions  ******************************/


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
    renumber_info->nb_bits = 32;
  else
    renumber_info->nb_bits = 64;
  ASSERT_ALWAYS (renumber_info->nb_bits <= 8 * sizeof(p_r_values_t));

  renumber_info->add_full_col = 0;
  renumber_info->size = 0;
  renumber_info->table = NULL;
  renumber_info->file = NULL;
  renumber_info->cached = NULL;
  renumber_info->bad_ideals.n = 0;
  renumber_info->bad_ideals.p = NULL;
  renumber_info->bad_ideals.r = NULL;
  renumber_info->bad_ideals.nb = NULL;
  renumber_info->bad_ideals.side = NULL;
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
  if (renumber_info->bad_ideals.side != NULL)
    free(renumber_info->bad_ideals.side);

  renumber_info->bad_ideals.r = NULL;
  renumber_info->bad_ideals.p = NULL;
  renumber_info->bad_ideals.nb = NULL;
  renumber_info->bad_ideals.side = NULL;
}

/* The renumber_t struct MUST have been initialized before */
void
renumber_init_write (renumber_t tab, const char *tablefile, const char *badfile,
                     int add_full_col)
{
  // Open a file with no extension. Compression is done later if needed
  char * tablefile_tmp = NULL;
  int rc;
  rc = asprintf(&tablefile_tmp, "%s-tmp", tablefile);
  ASSERT_ALWAYS(rc >= 0);
  fprintf (stderr, "Opening %s to write the temporary renumbering table\n",
                   tablefile_tmp);
  tab->file = fopen(tablefile_tmp, "w");
  ASSERT_ALWAYS(tab->file != NULL);

  tab->add_full_col = (add_full_col) ? 1 : 0;
  tab->size = (add_full_col) ? 1 : 0;

  // Write first the bad ideals information at the beginning of file
  if (badfile != NULL)
    tab->bad_ideals.n = copy_bad_ideals (badfile, tab->file, &tab->size);
  else
    tab->bad_ideals.n = 0;

  print_info (stderr, tab);
  free(tablefile_tmp);
}

void
renumber_close_write (renumber_t tab, const char *tablefile)
{
  char * tablefile_tmp = NULL;
  int rc;
  rc = asprintf(&tablefile_tmp, "%s-tmp", tablefile);
  ASSERT_ALWAYS(rc >= 0);


  // Compression is now if needed
  FILE *final = NULL;
  fclose(tab->file);
  fprintf (stderr, "Opening %s to read the temporary renumbering table\n",
                   tablefile_tmp);
  tab->file = fopen(tablefile_tmp, "r");
  ASSERT_ALWAYS(tab->file != NULL);
  fprintf (stderr, "Opening %s to write the final renumbering table\n",
                   tablefile);
  final = fopen_maybe_compressed (tablefile, "w");
  ASSERT_ALWAYS (final != NULL);

  // First we put some data about the renumbering table.
  fprintf (final, "%u %" PRid " %d %d\n", tab->nb_bits, tab->size,
                                     tab->bad_ideals.n, tab->add_full_col);

  char buffer[128] = "" , *retc;
  int ret;
  do
  {
    ret = fputs (buffer, final);
    ASSERT (ret >= 0);

    retc = fgets (buffer, 128, tab->file);
  } while (retc != NULL);

  ASSERT_ALWAYS (feof(tab->file));

  fclose_maybe_compressed (final, tablefile);
  fclose(tab->file);

  //remove tmp file
  ASSERT_ALWAYS (remove(tablefile_tmp) == 0);

  free (tablefile_tmp);
  fprintf(stderr, "Renumbering struct: nprimes=%"PRid"\n", tab->size);
  renumber_free(tab);
}


/* The renumber_t struct MUST have been initialized before */
/* renumber_free MUST be called to free the table afterwards */
void
renumber_read_table (renumber_t tab, const char * filename)
{
  index_t i, ret, report = 0;
  double t = wct_seconds ();
  uint8_t old_nb_bits = tab->nb_bits;
  char s[RENUMBER_MAXLINE];
  size_t bytes_read = 0;
  p_r_values_t v, prev_v = 0;

  //open file for reading
  fprintf (stderr, "Opening %s to read the renumbering table\n", filename);
  tab->file = fopen_maybe_compressed (filename, "r");
  ASSERT_ALWAYS(tab->file != NULL);

  // read size of renumbering table
  uint64_t tmp_size;
  ret = fscanf (tab->file, "%"SCNu8" %"SCNu64" %d %d\n", &tab->nb_bits,
                        &tmp_size, &tab->bad_ideals.n, &tab->add_full_col);
  ASSERT_ALWAYS (ret == 4);

  ASSERT_ALWAYS (tab->nb_bits <= 8 * sizeof(p_r_values_t));
  ASSERT_ALWAYS (tab->nb_bits == 32 || tab->nb_bits == 64);
  if (old_nb_bits != tab->nb_bits)
  {
    fprintf (stderr, "Warning, computed value of nb_bits (%d) is different "
                     "from the read value of nb_bits (%d).\n", old_nb_bits,
                     tab->nb_bits);
  }

  ASSERT_ALWAYS (tmp_size != 0);
  if (tab->nb_bits == 32)
    ASSERT_ALWAYS (!(tmp_size >> 32));
  tab->size = (index_t) tmp_size;

  // Allocating memory
  renumber_alloc(tab);

  // Begin to read the renumbering table
  i = 0;

  if (tab->add_full_col)
  {
    tab->table[i] = RENUMBER_SPECIAL_VALUE;
    i++;
  }

  for (int k = 0; k < tab->bad_ideals.n; k++)
  {
    bytes_read += get_one_line(tab->file, s);
    parse_one_line_bad_ideals (&tab->bad_ideals, s, k);
    for (int j = 0; j < tab->bad_ideals.nb[k]; j++)
    {
      tab->table[i] = RENUMBER_SPECIAL_VALUE;
      i++;
    }
  }

  // we cached value below 2^MAX_LOG_CACHED
  while (i < tab->size)
  {
    bytes_read += get_one_line(tab->file, s);
    v = parse_one_line(s);
    tab->table[i] = v;

    if ((v >> MAX_LOG_CACHED))
    {
      i++;
      break;
    }
    if (v > prev_v)
      tab->cached[v] = i;

    prev_v = v;
    i++;
  }

  tab->first_not_cached = i;

  for (; i < tab->size; i++)
  {
    bytes_read += get_one_line(tab->file, s);
    tab->table[i] = parse_one_line(s);

    if ((i >> 23) != (report >> 23))
    {
      report = i;
      report_read (t, i, bytes_read, 0);
    }
  }

  report_read (t, i, bytes_read, 1);
  ASSERT_ALWAYS(i == tab->size);

  print_info (stderr, tab);
  fprintf(stderr, "Renumbering struct: nprimes=%"PRid"\n", tab->size);
  fprintf(stderr, "Renumbering struct: first_not_cached=%"PRid"\n",
                                                        tab->first_not_cached);

  fclose_maybe_compressed (tab->file, filename);
}

int renumber_is_bad (int *nb, index_t *first, renumber_t rn, p_r_values_t p,
                     p_r_values_t r, int side)
{
  *first = 0;
  int bad = 0;
  for (int i = 0; i < rn->bad_ideals.n; ++i)
  {
    if (p == rn->bad_ideals.p[i] && r == rn->bad_ideals.r[i]
                                 && side == rn->bad_ideals.side[i])
    {
      *nb = rn->bad_ideals.nb[i];
      bad = 1;
      break;
    }
    else
      *first += rn->bad_ideals.nb[i];
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
      fprintf(renumber_info->file, "%lx\n", r[1][i]);
    for (i = 0; i < k[0]; i++)
      fprintf(renumber_info->file, "%lx\n", r[0][i]);

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
      fprintf(renumber_info->file, "%lx\n", add_value);
      renumber_info->size++;
    }

    for (i = 0; i < k[alg]; i++)
      fprintf(renumber_info->file, "%lx\n", r[alg][i]);

    renumber_info->size += k[alg];
  }
}

#define DEBUG_RENUMB
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
#ifdef DEBUG_RENUMB
    int nstep = 0;
#endif
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
#ifdef DEBUG_RENUMB
      nstep++;
      ASSERT_ALWAYS (nstep < 64);
#endif
    }
  }
  else //p is cached
  {
    i = renumber_info->cached[vp];
#ifdef DEBUG_RENUMB
    ASSERT_ALWAYS (renumber_info->table[i] == vp);
#endif
  }

  /* now i points at the beginning of a decreasing sequence of values of vr */
  if (side != renumber_info->rat)
  {
    while (i < renumber_info->size)
    {
      if (vr < tab[i])
      {
        if (i == renumber_info->size - 1 || vr > tab[i+1])
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

void
renumber_get_p_r_from_index (renumber_t renumber_info, p_r_values_t *p,
                             p_r_values_t * r, int *side, index_t i,
                             cado_poly pol)
{
  index_t j;
  p_r_values_t *tab = renumber_info->table;

  for (j = i; j > 0 && tab[j-1] > tab[j] && tab[j-1] != RENUMBER_SPECIAL_VALUE;)
    j--;

  if (renumber_info->rat == -1)
  {
    *p = (tab[j] - 1) >> 1;
    if (i == j)
    {
      int ret;
      //if we have at least one root on side 1, it is the largest
      if (get_largest_root_mod_p(r, pol->pols[1]->f, pol->pols[1]->degree, *p))
        *side = 1;
      else // else this is the largest on side 0
      {
        ret=get_largest_root_mod_p(r, pol->pols[0]->f, pol->pols[0]->degree, *p);
        ASSERT_ALWAYS (ret > 0);
        *side = 0;
      }
    }
    else
    {
      if (tab[i] <= *p)
      {
        *r = tab[i];
        *side = 0;
      }
      else
      {
        *r = tab[i] - *p - 1;
        *side = 1;
      }
    }
  }
  else
  {
    *p = tab[j] - 1;
    if (*p > (1UL << pol->rat->lpb) && i == j)
    {
      // Case where there is only alg side (p >= lpbr) and we are on the largest
      // root on alg side (i == j)
      int ret = get_largest_root_mod_p(r, pol->alg->f, pol->alg->degree, *p);
      ASSERT_ALWAYS (ret > 0);
      *side = 1 - renumber_info->rat;
    }
    else if (i == j)
    {
      *r = *p + 1; // old convention (r = p+1 in rat side). Has no meaning now.
      *side = renumber_info->rat;
    }
    else
    {
      *r = tab[i];
      *side = 1 - renumber_info->rat;
    }
  }
}

//for DEBUG, should be remove later
void renumber_debug_print_tab (FILE *output, const char *filename,
        cado_poly pol)
{
  renumber_t tab;
  index_t i;
  p_r_values_t p, r;
  int side;

  renumber_init (tab, pol);
  renumber_read_table (tab, filename);

  for (i = 0; i < tab->size; i++)
  {
    if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
    {
      if (i == 0 && tab->add_full_col)
        fprintf (output, "i=0 tab[i]=#   added column\n");
      else
        fprintf (output, "i=%" PRxid " tab[i]=#   above a bad ideals\n", i);
    }
    else
    {
      renumber_get_p_r_from_index (tab, &p, &r, &side, i, pol);
      fprintf (output, "i=%" PRxid " tab[i]=%" PRpr " p=%" PRpr "",
                       i, tab->table[i], p);
      if (side == tab->rat)
        fprintf (output, " rat side\n");
      else if (r == p)
        fprintf (output, " r=%" PRpr " alg side proj\n", r);
      else
        fprintf (output, " r=%" PRpr " alg side\n", r);
    }
  }

  if (tab->bad_ideals.n != 0) {
    fprintf (output, "Bad ideals:\n");
    for (int i = 0; i < tab->bad_ideals.n; ++i) {
      p_r_values_t p = tab->bad_ideals.p[i];
      p_r_values_t r = tab->bad_ideals.r[i];
      int nb = tab->bad_ideals.nb[i];
      fprintf(output, "p=%" PRpr " r=%" PRpr " nb=%d\n", p, r, nb);
    }
  }

  renumber_free(tab);
}
