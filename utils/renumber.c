#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <inttypes.h>
#include <stdio.h>
#include <ctype.h>
#include "portability.h"
#include "utils.h"

//#define RENUMBER_DO_EXPENSIVE_CHECK

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


static uint64_t next_prime_of_powers_of_2[64] = { 0x2, 0x3, 0x5, 0xb, 0x11,
        0x25, 0x43, 0x83, 0x101, 0x209, 0x407, 0x805, 0x1003, 0x2011, 0x401b,
        0x8003, 0x10001, 0x2001d, 0x40003, 0x80015, 0x100007, 0x200011,
        0x40000f, 0x800009, 0x100002b, 0x2000023, 0x400000f, 0x800001d,
        0x10000003, 0x2000000b, 0x40000003, 0x8000000b, 0x10000000f,
        0x200000011, 0x400000019, 0x800000035, 0x100000001f, 0x2000000009,
        0x4000000007, 0x8000000017, 0x1000000000f, 0x2000000001b, 0x4000000000f,
        0x8000000001d, 0x100000000007, 0x20000000003b, 0x40000000000f,
        0x800000000005, 0x1000000000015, 0x2000000000045, 0x4000000000037,
        0x8000000000015, 0x10000000000015, 0x20000000000005, 0x4000000000009f,
        0x80000000000003, 0x100000000000051, 0x200000000000009,
        0x400000000000045, 0x800000000000083, 0x1000000000000021,
        0x200000000000000f, 0x4000000000000087, 0x800000000000001d };

static uint64_t previous_prime_of_powers_of_2[65] = { 0x0, 0x0, 0x3, 0x7, 0xd,
        0x1f, 0x3d, 0x7f, 0xfb, 0x1fd, 0x3fd, 0x7f7, 0xffd, 0x1fff, 0x3ffd,
        0x7fed, 0xfff1, 0x1ffff, 0x3fffb, 0x7ffff, 0xffffd, 0x1ffff7, 0x3ffffd,
        0x7ffff1, 0xfffffd, 0x1ffffd9, 0x3fffffb, 0x7ffffd9, 0xfffffc7,
        0x1ffffffd, 0x3fffffdd, 0x7fffffff, 0xfffffffb, 0x1fffffff7,
        0x3ffffffd7, 0x7ffffffe1, 0xffffffffb, 0x1fffffffe7, 0x3fffffffd3,
        0x7ffffffff9, 0xffffffffa9, 0x1ffffffffeb, 0x3fffffffff5, 0x7ffffffffc7,
        0xfffffffffef, 0x1fffffffffc9, 0x3fffffffffeb, 0x7fffffffff8d,
        0xffffffffffc5, 0x1ffffffffffaf, 0x3ffffffffffe5, 0x7ffffffffff7f,
        0xfffffffffffd1, 0x1fffffffffff91, 0x3fffffffffffdf, 0x7fffffffffffc9,
        0xfffffffffffffb, 0xfffffffffffff3, 0x3ffffffffffffe5, 0xffffffffffffc9,
        0xfffffffffffffa3, 0x1fffffffffffffff, 0x3fffffffffffffc7,
        0x7fffffffffffffe7, 0xffffffffffffffc5 };



/* Skip line beginning with '#'. Return 0 if fgets return NULL.*/
static inline size_t
get_one_line (FILE *f, char *s)
{
  char *rets;
  size_t n;
  do
  {
    rets = fgets(s, RENUMBER_MAXLINE, f);
    if (rets == NULL)
    {
      n = 0;
      break;
    }
    else if (rets[0] != '#')
    {
      n = strnlen(s, RENUMBER_MAXLINE);
      ASSERT_ALWAYS(n != RENUMBER_MAXLINE);
      break;
    }
    // else we skip the line
  } while (1);
  return n;
}

static void
parse_bad_ideals_file (FILE *badfile, renumber_t renum)
{
  renum->bad_ideals.n = 0;
  char s[RENUMBER_MAXLINE];

  while (get_one_line (badfile, s) != 0)
  {
    renum->bad_ideals.n++;
    const char *ptr = s;
    int nb, t;
    for (int count = 0 ; *ptr != '\n'; ptr++)
    {
      if (!(isspace(ptr[0])))
      {
        if (ptr[0] == ':')
          count++;
        else if (count >= 2)
          break;
      }
    }

    for (nb = 0; (t = ugly[(unsigned char) *ptr]) >= 0; ptr++)
      nb = (nb * 10) + t;

    ASSERT_ALWAYS (*ptr == '\n');
    renum->size += nb;
  }
  ASSERT_ALWAYS (feof (badfile));
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
print_info (FILE * f, renumber_t r, int after_reading)
{
  char pre[9] = "# INFO: ";
  fprintf (f, "# Information on renumber struct:\n%ssizeof(p_r_values_t) = %zu\n"
              "%snb_bits = %" PRIu8 "\n%srat = %d %s\n%s#badideals = %d\n"
              "%sadd_full_col = %d\n%slpb0 = %lu\n%slpb1 = %lu\n",
              pre, sizeof(p_r_values_t), pre, r->nb_bits, pre, r->rat,
              (r->rat == -1) ? "(no rational side)" : "", pre, r->bad_ideals.n,
              pre, r->add_full_col, pre, r->lpb[0], pre, r->lpb[1]);

  if (after_reading) /* there is more info to print*/
  {
    fprintf (f, "%ssize = %" PRIu64 "\n"
                "%ssmallest prime not cached = %" PRpr " at index %" PRid "\n"
                "%sbiggest prime below lbp0 = %" PRpr " at index %" PRid "\n"
                "%sbiggest prime below lbp1 = %" PRpr " at index %" PRid "\n",
                pre, r->size, pre, r->smallest_prime_not_cached,
                r->index_smallest_prime_not_cached, pre,
                r->biggest_prime_below_lpb[0],
                r->index_biggest_prime_below_lpb[0], pre,
                r->biggest_prime_below_lpb[1],
                r->index_biggest_prime_below_lpb[1]);
  }
  fflush (stdout);
}

/* sort in decreasing order. Fastest for ~ < 15 values in r[] vs qsort */
inline void
renumber_sort_ul (unsigned long *r, size_t n)
{
  unsigned long rmin;

  if (UNLIKELY (n < 2))
    return;

  if (UNLIKELY (n == 2)) {
    if (r[0] < r[1]) {
      rmin = r[0];
      r[0] = r[1];
      r[1] = rmin;
    }
    return;
  }

  for (size_t i = n; --i;) {
    size_t min = i;
    rmin = r[min];
    for (size_t j = i; j--;) {
      unsigned long rj = r[j];
      if (UNLIKELY (rj < rmin)) {
	min = j;
	rmin = rj;
      }
    }
    if (LIKELY (min != i)) {
      r[min] = r[i];
      r[i] = rmin;
    }
  }
}

/* return zero if no roots mod p, else non-zero */
static int
get_largest_root_mod_p (p_r_values_t *r, mpz_t *pol, int deg, p_r_values_t p)
{
  // if there is a proj root, this the largest (r = p by convention)
  if (mpz_divisible_ui_p (pol[deg], p)) {
    *r = p;
    return 1;
  }

  mpz_poly_t f = {{ .coeff = pol, .deg = deg }};
  unsigned long roots[deg];
  size_t k = (size_t) mpz_poly_roots_ulong (roots, f, p);
  if (k) {
    unsigned long max = roots[--k];
    while (k--)
      if (UNLIKELY (roots[k] > max)) max = roots[k];
    *r = max;
    return 2;
  }

  return 0;
}

static void
renumber_write_first_line (renumber_t renum)
{
  fprintf (renum->file, "%" PRIu8 " %d %d %d %lu %lu\n", renum->nb_bits,
           renum->rat, renum->bad_ideals.n, renum->add_full_col, renum->lpb[0],
           renum->lpb[1]);
}

static void
renumber_read_first_line (renumber_t renum)
{
  int ret;
  ret = fscanf (renum->file, "%" SCNu8 " %d %d %d %lu %lu\n", &(renum->nb_bits),
                &(renum->rat), &(renum->bad_ideals.n), &(renum->add_full_col),
                &(renum->lpb[0]), &(renum->lpb[1]));
  ASSERT_ALWAYS (ret == 6);
  ASSERT_ALWAYS (-1 <= renum->rat && renum->rat <= 1);
  ASSERT_ALWAYS (renum->add_full_col == 0 || renum->add_full_col == 1);
  ASSERT_ALWAYS (renum->nb_bits <= 8 * sizeof(p_r_values_t));
  ASSERT_ALWAYS (renum->nb_bits == 32 || renum->nb_bits == 64);
}

/* Assume v correspond to p + 1 or 2p + 1 (i.e. v > 0) */
static inline p_r_values_t
get_p_from_table_value (renumber_t tab, p_r_values_t v)
{
  if (tab->rat >= 0) /* One alg and one rat side -> p is v-1 */
      return (v - 1);
  else               /* Two alg sides -> p is (v-1)/2 */
      return (v >> 1);
}

/*********************** End internal functions  ******************************/

void
renumber_init_for_reading (renumber_t renumber_info)
{
  memset(renumber_info, 0, sizeof(renumber_t));
  /* Will be set later, by renumber_read_table, with the values of the first
   * line of the renumber file */
}

/* rat contains the rational side (0 or 1) or -1 if two algebraic sides.
 * add_full_col is non-zero if we need to add a column of 1 in the matrix, 0
 * otherwise (for factorization, always 0, for DL 1 if one of the polynomials is
 * not monic). */
void
renumber_init_for_writing (renumber_t renumber_info, int rat, int add_full_col,
                           unsigned long lpb[])
{
  memset(renumber_info, 0, sizeof(renumber_t));

  ASSERT_ALWAYS (-1 <= rat && rat <= 1);
  ASSERT_ALWAYS (add_full_col == 0 || add_full_col == 1);
  ASSERT_ALWAYS (lpb != NULL);
  renumber_info->rat = rat;
  renumber_info->add_full_col = add_full_col;
  renumber_info->lpb[0] = lpb[0];
  renumber_info->lpb[1] = lpb[1];

  int max_nb_bits = MAX(lpb[0], lpb[1]);
  if (renumber_info->rat == -1) /* for two alg side, we need an extra bit. */
    max_nb_bits++;

  if (max_nb_bits <= 32)
    renumber_info->nb_bits = 32;
  else
    renumber_info->nb_bits = 64;
  ASSERT_ALWAYS (renumber_info->nb_bits <= 8 * sizeof(p_r_values_t));
}

void
renumber_clear (renumber_t renumber_info)
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

/* The renumber_t struct _must_ have been initialized before
 * poly = NULL is accepted. It will not print the polynomials on the file */
void
renumber_write_open (renumber_t tab, const char *tablefile, const char *badfile,
                     cado_poly poly)
{
  printf ("# Opening %s to write the renumbering table\n", tablefile);
  fflush (stdout);
  tab->file = fopen_maybe_compressed (tablefile, "w");
  ASSERT_ALWAYS(tab->file != NULL);
  FILE *fbad = NULL;
  if (badfile != NULL)
  {
    printf ("# Opening %s to read the bad ideals\n", badfile);
    fbad = fopen(badfile, "r"); /* never compressed, always small */
    ASSERT_ALWAYS (fbad != NULL);
  }

  tab->size = (tab->add_full_col) ? 1 : 0;

  /* Read bad ideals files */
  if (badfile != NULL)
    parse_bad_ideals_file (fbad, tab); /* update size et bad_ideals.n */

  /* Write the first line */
  renumber_write_first_line (tab);

  /* Print info on stdout (~ what is written on the first line of the file) */
  print_info (stdout, tab, 0);

  /* Write the two polynomials on a line beginning by #, if given */
  if (poly != NULL)
  {
    fprintf (tab->file, "# ");
    mpz_poly_fprintf (tab->file, poly->pols[0]);
    fprintf (tab->file, "# ");
    mpz_poly_fprintf (tab->file, poly->pols[1]);
  }

  /* Write first the bad ideals information at the beginning of file */
  if (badfile != NULL)
  {
    char s[RENUMBER_MAXLINE];
    rewind (fbad);
    while (get_one_line (fbad, s) != 0)
      fputs (s, tab->file);
    ASSERT_ALWAYS (feof (fbad));
  }

  if (fbad != NULL)
    fclose (fbad);
}

void
renumber_write_close (renumber_t tab, const char *tablefile)
{
  fclose_maybe_compressed (tab->file, tablefile);
}


/* The renumber_t struct _must_ have been initialized before */
void
renumber_read_table (renumber_t tab, const char * filename)
{
  char s[RENUMBER_MAXLINE];
  size_t bytes_read = 0, bytes_line;
  stats_data_t infostats;  /* for displaying progress */

  /* open file for reading */
  printf ("# Opening %s to read the renumbering table\n", filename);
  fflush (stdout);
  tab->file = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK (tab->file == NULL, "Cannot open file for reading");

  /* read size of renumbering table */
  renumber_read_first_line (tab);

  /* Allocating memory */
  size_t badideals_pr_size = tab->bad_ideals.n * sizeof (p_r_values_t);
  size_t badideals_int_size = tab->bad_ideals.n * sizeof (int);
  size_t cached_table_size = (2 << MAX_LOG_CACHED) * sizeof (index_t);

    /* Do not know the size yet. Reallocating while reading */
    /* We assume that RENUMBER_DEFAULT_SIZE is enough to hold at least the bad
     * ideals and the added column (if add_full_col is set) */
  uint64_t allocated = RENUMBER_DEFAULT_SIZE;
  size_t default_size = RENUMBER_DEFAULT_SIZE * sizeof (p_r_values_t);

  tab->table           = (p_r_values_t *) malloc (default_size);
  tab->cached          = (index_t *)      malloc (cached_table_size);
  tab->bad_ideals.p    = (p_r_values_t *) malloc (badideals_pr_size);
  tab->bad_ideals.r    = (p_r_values_t *) malloc (badideals_pr_size);
  tab->bad_ideals.side = (int *)          malloc (badideals_int_size);
  tab->bad_ideals.nb   = (int *)          malloc (badideals_int_size);

  ASSERT_ALWAYS (tab->table != NULL);
  ASSERT_ALWAYS (tab->cached != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.p != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.r != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.nb != NULL);
  ASSERT_ALWAYS (tab->bad_ideals.side != NULL);

  memset (tab->cached, 0, cached_table_size);

  if (tab->add_full_col)
  {
    tab->table[0] = RENUMBER_SPECIAL_VALUE;
    tab->size = 1;
  }
  else
    tab->size = 0;

  /* Reading the bad ideals at the top of the renumbering file */
  for (int k = 0; k < tab->bad_ideals.n; k++)
  {
    bytes_read += get_one_line(tab->file, s);
    parse_one_line_bad_ideals (&tab->bad_ideals, s, k);
    for (int j = 0; j < tab->bad_ideals.nb[k]; j++)
    {
      tab->table[tab->size] = RENUMBER_SPECIAL_VALUE;
      tab->size++;
    }
  }

  p_r_values_t prime_cache_limit = next_prime_of_powers_of_2[MAX_LOG_CACHED];
  p_r_values_t expected_biggest_prime_lpb[2];
  expected_biggest_prime_lpb[0] = previous_prime_of_powers_of_2[tab->lpb[0]];
  expected_biggest_prime_lpb[1] = previous_prime_of_powers_of_2[tab->lpb[1]];
  int has_smallest = 0;

  /* Reading the renumbering table */
  stats_init (infostats, stdout, &(tab->size), 24, "Read", "elements", "",
              "elts");

  while ((bytes_line = get_one_line(tab->file, s)) > 0)
  {
    bytes_read += bytes_line;
    if (tab->size >= allocated) /* Not enough space, reallocated tab->table */
    {
      allocated += RENUMBER_DEFAULT_SIZE;
      size_t new_size = allocated * sizeof (p_r_values_t);
      tab->table = (p_r_values_t *) realloc (tab->table, new_size);
      ASSERT_ALWAYS (tab->table != NULL);
    }
    tab->table[tab->size] = parse_one_line(s);

    if (tab->size == 0 || tab->table[tab->size-1] == RENUMBER_SPECIAL_VALUE
                       || tab->table[tab->size] > tab->table[tab->size-1])
    {
      /* We just switch to a new prime in the renumbering table, see if we need
       * to cache it (we cached primes below 2^MAX_LOG_CACHED)
       */
      p_r_values_t p = get_p_from_table_value (tab, tab->table[tab->size]);
      if (p < prime_cache_limit) /* p < 2^MAX_LOG_CACHED */
        tab->cached[p] = tab->size;
      else if (!has_smallest)
      {
        has_smallest = 1;
        tab->index_smallest_prime_not_cached = tab->size;
        tab->smallest_prime_not_cached = p;
      }

      if (p <= expected_biggest_prime_lpb[0])
      {
        tab->index_biggest_prime_below_lpb[0] = tab->size;
        tab->biggest_prime_below_lpb[0] = p;
      }
      if (p <= expected_biggest_prime_lpb[1])
      {
        tab->index_biggest_prime_below_lpb[1] = tab->size;
        tab->biggest_prime_below_lpb[1] = p;
      }
    }
    tab->size++;

    if (stats_test_progress(infostats))
      stats_print_progress (infostats, tab->size, 0, bytes_read, 0);
  }

  if (!has_smallest) /* Every prime is cached. */
  {
    tab->index_smallest_prime_not_cached = tab->size;
    tab->smallest_prime_not_cached =
           next_prime_of_powers_of_2[MAX(tab->lpb[0],tab->lpb[1])];
  }

  stats_print_progress (infostats, tab->size, 0, bytes_read, 1);
  size_t final_size = tab->size * sizeof (p_r_values_t);
  tab->table = (p_r_values_t *) realloc (tab->table, final_size);


  ASSERT_ALWAYS (feof (tab->file));

  print_info (stdout, tab, 1);

  fclose_maybe_compressed (tab->file, filename);
}

int renumber_is_bad (int *nb, index_t *first, renumber_t rn, p_r_values_t p,
                     p_r_values_t r, int side)
{
  *first = (rn->add_full_col) ? 1 : 0;
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

inline size_t
renumber_write_p_2algs (unsigned long p, unsigned long *roots_alg0, size_t nb_roots_alg0,
		       unsigned long *roots_alg1, size_t nb_roots_alg1, char *buffer) {
  size_t size_buffer = 0;
  size_t i;

  renumber_sort_ul(roots_alg0, nb_roots_alg0);
  renumber_sort_ul(roots_alg1, nb_roots_alg1);

  if (LIKELY (nb_roots_alg1)) {
    size_buffer += sprintf (buffer, "%lx\n", (p << 1) + 1); // The largest roots become 2p+1
    for (i = 1; i < nb_roots_alg1; ++i)                      // Add p + 1 on side 1
      size_buffer += sprintf (buffer + size_buffer, "%lx\n", roots_alg1[i] + p + 1);
  }
  else
    *roots_alg0 = (p << 1) + 1;                             // The largest roots become 2p+1

  for (i = 0; i < nb_roots_alg0; ++i)
    size_buffer += sprintf (buffer + size_buffer, "%lx\n", roots_alg0[i]);

  return size_buffer;
}

inline size_t
renumber_write_p_rat_alg (unsigned long p, size_t nb_roots_rat,
			  unsigned long *roots_alg, size_t nb_roots_alg, char *buffer) {
  size_t size_buffer;
  size_t i;

  renumber_sort_ul(roots_alg, nb_roots_alg);

  if (LIKELY(nb_roots_rat)) {
    size_buffer = sprintf (buffer, "%lx\n", p + 1);
  }
  else {
    size_buffer = 0;
    *roots_alg = p + 1; // lpbr < p < lpba, we put p+1 instead of the largest roots
  }

  for (i = 0; i < nb_roots_alg; ++i) {
    size_buffer += sprintf (buffer + size_buffer, "%lx\n", roots_alg[i]);
  }

  return size_buffer;
}

void
renumber_write_p (renumber_t renumber_info, unsigned long p, unsigned long *r[2], int k[2])
{
  size_t size_buffer;
  char buffer[512];

  if (renumber_info->rat == -1)
    size_buffer = renumber_write_p_2algs (p, r[0], (size_t) k[0], r[1], (size_t) k[1], buffer);
  else
    size_buffer = renumber_write_p_rat_alg (p, (size_t) k[renumber_info->rat], r[renumber_info->rat ^ 1], (size_t) k[renumber_info->rat ^ 1], buffer);

  fwrite ((void *) buffer, size_buffer, 1, renumber_info->file);
  renumber_info->size += k[0] + k[1];
}

/* side is 0 if (p,r) corresponds to the left part in the relation,
 * side is 1 if (p,r) corresponds to the right part.
 * If side corresponds to the rational side (if it exists), the value of r is
 * meaningless. */
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

  /**************************************************************************/
  /* Search for i such that
        renumber_info->table[i] = vp
        this is the beginning of a decreasing sequence
  */

  /* For small value of p, the corresponding value of i is cached */
  if (p < renumber_info->smallest_prime_not_cached)
  {
    i = renumber_info->cached[p];
    if (UNLIKELY(tab[i] != vp))
    {
      /* There is a problem, most probably p is not prime. */
      fprintf(stderr, "Fatal error in %s at %s:%d\nError with the cached part of"
                      " the renumbering table\n  p = %" PRpr "\n  vp = %" PRpr
                      "\n  i = cached[p] = %" PRid "\n  tab[i] = %" PRpr "\n",
                      __func__, __FILE__, __LINE__, p, vp, i, tab[i]);
      abort();
    }
  }
  /* p is not cached and below the lpb[side] */
  else if (p <= renumber_info->biggest_prime_below_lpb[side])
  {
#ifdef RENUMBER_DO_EXPENSIVE_CHECK
    int nstep = 0;
#endif
    index_t max = renumber_info->index_biggest_prime_below_lpb[side];
    index_t min = renumber_info->index_smallest_prime_not_cached;

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
#ifdef RENUMBER_DO_EXPENSIVE_CHECK
      {
        /* Stop infinite loops (which indicate bug) */
        nstep++;
        char tmp[256];
        snprintf (tmp, 256, "ntep=%d >= 64 (p=%" PRpr " not prime?)", nstep, p);
        FATAL_ERROR_CHECK ((nstep >= 64), tmp);
      }
#endif
    }
  }
  else /* Error */
  {
    /* prime p is bigger than lpb[side] => Fatal error */
    fprintf(stderr, "Fatal error in %s at %s:%d\nIdeal (p, r, side) = (%" PRpr
                    ", %" PRpr ", %d) is bigger that large prime bound 2^%ld\n",
                    __func__, __FILE__, __LINE__, p, r, side,
                    renumber_info->lpb[side]);
    abort();
  }


  /**************************************************************************/
  /* Now i points at the beginning of a decreasing sequence of values of vr */

  /* Return i in 4 cases:
    first case: an ideal on rational side always corresponds to the first
                element of a sequence
    second case: i is the last index of the table
    third case: the sequence contains only one value
    fourth case: next element of the sequence is too small to correspond to vr
  */
  if (side == renumber_info->rat || i == renumber_info->size - 1
                                 || tab[i] <= tab[i+1]
                                 || vr > tab[i+1])
    return i;
  else /* else we go through the sequence until we find vr */
  {
    while(i != renumber_info->size - 1 && tab[i] > tab[i+1])
    {
      i++;
      if (vr == tab[i])
        return i;
    }
    /* if we arrive here, there is a problem, the ideal was not found in the
       renumbering table */
    fprintf(stderr, "Fatal error in %s at %s:%d\nIdeal (p, r, side) = (%" PRpr
                    ", %" PRpr ", %d) was not found on the renumbering table\n",
                    __func__, __FILE__, __LINE__, p, r, side);
    abort();
  }
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
      if (get_largest_root_mod_p(r, pol->pols[1]->coeff, pol->pols[1]->deg, *p))
        *side = 1;
      else // else this is the largest on side 0
      {
        ret=get_largest_root_mod_p(r, pol->pols[0]->coeff, pol->pols[0]->deg, *p);
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
    unsigned long lpbr = renumber_info->lpb[renumber_info->rat];
    if (*p > (1UL << lpbr) && i == j)
    {
      // Case where there is only alg side (p >= lpbr) and we are on the largest
      // root on alg side (i == j)
      int ret = get_largest_root_mod_p(r, pol->alg->coeff, pol->alg->deg, *p);
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

