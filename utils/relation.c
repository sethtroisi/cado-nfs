#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> /* for PRId64 */
#include <ctype.h> /* for isxdigit */
#include <string.h>
#include <errno.h>

#include "cado_poly.h"
#include "mod_ul.h"
#include "relation.h"
#include "gzip.h"
#include "timing.h"
#include "portability.h"

/*
  Convention for I/O of rels:
    a and b are printed in decimal
    primes are printed in hexadecimal.
*/

// Sometimes, we want our space back...!
void
relation_init (relation_t *rel)
{
    memset(rel, 0, sizeof(relation_t));
}

void
relation_clear (relation_t *rel)
{
    free(rel->rp); rel->rp = NULL;
    free(rel->ap); rel->ap = NULL;
    memset(rel, 0, sizeof(relation_t));
}

/* FIXME: The following interface still strongly relies on the fact that
 * the rational side is [0] and the algebraic side is [1] */
void
relation_add_prime (relation_t *rel, int side, unsigned long p)
{
    int prov[2] = {
        [RATIONAL_SIDE] = rel->nb_rp,
        [ALGEBRAIC_SIDE] = rel->nb_ap, };
    prov[side]++;
    relation_provision_for_primes(rel, prov[0], prov[1]);
    if (side == RATIONAL_SIDE) {
        rel->rp[rel->nb_rp++] = (rat_prime_t) { .p=p, .e=1 };
    } else {
        rel->ap[rel->nb_ap++] = (alg_prime_t) { .p=p, .e=1 };
    }
}

void relation_copy (relation_t *s, relation_t * r)
{
    s->a = r->a;
    s->b = r->b;
    s->nb_rp = r->nb_rp;
    s->nb_ap = r->nb_ap;
    s->nb_rp_alloc = r->nb_rp;
    s->nb_ap_alloc = r->nb_ap;
    s->rp = realloc(s->rp, s->nb_rp_alloc * sizeof(rat_prime_t));
    s->ap = realloc(s->ap, s->nb_ap_alloc * sizeof(alg_prime_t));
    memcpy(s->rp, r->rp, r->nb_rp * sizeof(rat_prime_t));
    memcpy(s->ap, r->ap, r->nb_ap * sizeof(alg_prime_t));
}


void skip_relations_in_file(FILE * f, int n)
{
    char str[RELATION_MAX_BYTES];

    for( ; n-- ; ) {
	char * p = fgets(str, RELATION_MAX_BYTES, f);
        if (!p) {
	    fprintf(stderr, "short read in skip_relations_in_file\n");
            abort();
        }
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Too long string in skip_relations_in_file\n");
	}
    }
}


/* return a/b mod p, and p when gcd(b,p) <> 1: this corresponds to a
   projective root */
unsigned long
findroot(long a, unsigned long b, unsigned long p) {
  int sig, inv;
  unsigned long root;
  modulusul_t m;
  residueul_t r, t, pa, br;

  modul_initmod_ul (m, p);
  modul_init (pa, m);

  if (a >= 0) {
    modul_set_ul (pa, (unsigned long)a, m); /* Does reduction mod p */
    sig = 1;
  } else {
    modul_set_ul (pa, (unsigned long)(-a), m); /* Does reduction mod p */
    sig = -1;
  }

  modul_init (r, m);
  modul_init (t, m);
  modul_init (br, m);
  modul_set_ul (br, b, m); /* Does reduction mod p */
  inv = modul_inv(t, br, m);
  if (inv) {
    modul_mul(r, t, pa, m);
    if (sig == -1)
      modul_neg(r, r, m);
    root = modul_get_ul (r, m);
  } else {
    root = p;
  }
  
  modul_clear (pa, m); /* No-ops. Here for the sake of pedantry */
  modul_clear (r, m);
  modul_clear (t, m);
  modul_clear (br, m);
  modul_clearmod (m);
  return root;
}

// root = p if we don't know the result (p divides leading coeff)
void
computeroots (relation_t *rel)
{
  unsigned long r;
  int i;

  for (i = 0; i < rel->nb_ap; ++i)
    {
      r = findroot (rel->a, rel->b, rel->ap[i].p);
      rel->ap[i].r = r;
  }
}

/* {{{ a few conversion relations. We happen to export these */
static inline void
sswap(char *pd, char *pe)
{
  char c;
  while (pe > pd)
    {
      c = *pe; *pe-- = *pd; *pd++ = c;
    }
}

char *
u64toa16 (char *p, uint64_t m)
{
  char *op;
  static char v[] = "0123456789abcdef";

  op = p;
  do
    {
      *p++ = v[m & 0xf];
      m >>= 4;
    }
  while (m);
  sswap(op, &(p[-1]));
  return (p);
}

char *
u64toa10 (char *p, uint64_t m)
{
  char *op;
  uint64_t n;

  op = p;
  do
    {
      n = m / 10;
      *p++ = (unsigned char) ((m - n * 10) + '0');
      m = n;
    }
  while (m);
  sswap(op, &(p[-1]));
  return (p);
}

char *
d64toa10 (char *p, int64_t m)
{
  if (m < 0)
    {
      *p++ = '-';
      m = -m;
    }
  return (u64toa10 (p, (uint64_t) m));
}

char *
d64toa16 (char *p, int64_t m)
{
  if (m < 0)
    {
      *p++ = '-';
      m = -m;
    }
  return (u64toa16 (p, (uint64_t) m));
}
/*}}}*/

void
fprint_relation (FILE *file, relation_t * rel)
{
  char buf[1<<10], *p, *op;
  size_t lg;
  int i, j;
  
  p = d64toa10(buf, rel->a);
  *p++ = ',';
  p = u64toa10(p, rel->b);
  *p++ = ':';
  for (i = 0; i < rel->nb_rp; ++i)
    if (rel->rp[i].e)
      {
	op = p;
	p = u64toa16 (p, rel->rp[i].p);
	*p++ = ',';
	lg = (p - op);
	for (j = 0; j < rel->rp[i].e - 1; j++)
	  {
	    memcpy (p, op, lg);
	    p += lg;
	  }
      }
  p[-1] = ':';
  for (i = 0; i < rel->nb_ap; ++i)
    if (rel->ap[i].e)
      {
	op = p;
	p = u64toa16 (p, rel->ap[i].p);
	*p++ = ',';
	lg = (p - op);
	for (j = 0; j < rel->ap[i].e - 1; j++)
	  {
	    memcpy (p, op, lg);
	    p += lg;
	  }
      }
  p[-1] = '\n';
  fwrite (buf, sizeof(*buf), p - buf, file);
}

/* same as fprint_relation, but assumes all exponents are > 0 */
void
fprint_relation_raw (FILE *file, relation_t * rel)
{
  char buf[1<<10], *p, *op;
  size_t lg;
  int i, j;
  
  p = d64toa10(buf, rel->a);
  *p++ = ',';
  p = u64toa10(p, rel->b);
  *p++ = ':';
  for (i = 0; i < rel->nb_rp; ++i)
    {
      ASSERT (rel->rp[i].e);
      op = p;
      p = u64toa16 (p, rel->rp[i].p);
      *p++ = ',';
      lg = (p - op);
      for (j = 0; j < rel->rp[i].e - 1; j++)
	{
	  memcpy (p, op, lg);
	  p += lg;
	}
    }
  p[-1] = ':';
  for (i = 0; i < rel->nb_ap; ++i)
    {
      ASSERT (rel->ap[i].e);
      op = p;
      p = u64toa16 (p, rel->ap[i].p);
      *p++ = ',';
      lg = (p - op);
      for (j = 0; j < rel->ap[i].e - 1; j++)
        {
          memcpy (p, op, lg);
          p += lg;
        }
    }
  p[-1] = '\n';
  fwrite (buf, sizeof(*buf), p - buf, file);
}

/* reduces exponents mod 2, and discards primes with even exponent */
/* This assumes that relations have already been sorted */
void
reduce_exponents_mod2 (relation_t *rel)
{
  int i, j;

  if (rel->nb_rp == 0)
    fprintf (stderr, "WARNING: nb_rp = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_rp; i++)
    {
      rel->rp[i].e &= 1; /* reduce exponent mod 2 */
      if (rel->rp[i].e != 0)
        {
          rel->rp[j].p = rel->rp[i].p;
          rel->rp[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
    fprintf(stderr, "WARNING: j_rp=0 in reduce_exponents_mod2. k=%d\n", rel->nb_rp);
  rel->nb_rp = j;

  if(rel->nb_ap == 0)
      fprintf(stderr, "WARNING: nb_ap = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_ap; i++)
    {
      rel->ap[i].e &= 1; /* reduce exponent mod 2 */
      if (rel->ap[i].e != 0)
        {
          rel->ap[j].p = rel->ap[i].p;
          /* rel->ap[j].r = rel->ap[i].r;  */ // in fact useless at this point.
          rel->ap[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
      fprintf(stderr, "WARNING: j_ap = 0 in reduce_exponents_mod2\n");
  rel->nb_ap = j;
}

static int rat_prime_cmp(rat_prime_t * a, rat_prime_t * b)
{
    return (b->p < a->p) - (a->p < b->p);
}

static int alg_prime_cmp(alg_prime_t * a, alg_prime_t * b)
{
    int r = (b->p < a->p) - (a->p < b->p);
    if (r) return r;
    return (b->r < a->r) - (a->r < b->r);
}

typedef int (*sortfunc_t) (const void *, const void *);

void
relation_provision_for_primes (relation_t * rel, int nr, int na)
{
  if (nr > 0 && rel->nb_rp_alloc < nr)
    {
      rel->nb_rp_alloc = nr + nr / 2;
      rel->rp = (rat_prime_t *) realloc (rel->rp,
                                       rel->nb_rp_alloc * sizeof(rat_prime_t));
    }
  if (na > 0 && rel->nb_ap_alloc < na)
    {
      rel->nb_ap_alloc = na + na / 2;
      rel->ap = (alg_prime_t *) realloc (rel->ap,
                                       rel->nb_ap_alloc * sizeof(alg_prime_t));
    }
}

/* assumes all the exponents are initially 1 */
static void
sort_rat_primes (rat_prime_t *rp, int nb_rp)
{
  int i, j;
  unsigned long pi;

  /* insertion sort */
  for (i = 1; i < nb_rp; i++)
    {
      pi = rp[i].p;
      for (j = i - 1; (j >= 0) && (pi < rp[j].p); j--)
        rp[j+1].p = rp[j].p;
      /* the prime pi should go at index j+1 */
      rp[j+1].p = pi;
    }
}

/* assumes all the exponents are 1 */
static void
sort_alg_primes (alg_prime_t *ap, int nb_ap)
{
  int i, j;
  unsigned long pi, ri;

  /* insertion sort: for a given relation (a,b), r is uniquely determined
     by (a,b), thus we only need to compare p */
  for (i = 1; i < nb_ap; i++)
    {
      pi = ap[i].p;
      ri = ap[i].r;
      for (j = i - 1; (j >= 0) && (pi < ap[j].p); j--)
        {
          ap[j+1].p = ap[j].p;
          ap[j+1].r = ap[j].r;
        }
      ap[j+1].p = pi;
      ap[j+1].r = ri;
    }
}

void relation_compress_rat_primes(relation_t * rel)
{
#if 0
  qsort(rel->rp, rel->nb_rp, sizeof(rat_prime_t), (sortfunc_t) &rat_prime_cmp);
#else
  sort_rat_primes (rel->rp, rel->nb_rp);
#endif
    int j = 0;
    for(int i = 0 ; i < rel->nb_rp ; j++) {
        if (i-j) memcpy(rel->rp + j, rel->rp + i, sizeof(rat_prime_t));
        int e = 0;
        for( ; i < rel->nb_rp && rat_prime_cmp(rel->rp+i, rel->rp+j) == 0 ; i++)
            e+=rel->rp[i].e;
        rel->rp[j].e = e;
    }
    rel->nb_rp = j;
}

void relation_compress_alg_primes(relation_t * rel)
{
#if 0
    /* We're considering the list as containing possibly distinct (p,r)
     * pairs, although in reality it cannot happen. I doubt this causes
     * any performance hit though.
     */
  qsort(rel->ap, rel->nb_ap, sizeof(alg_prime_t), (sortfunc_t) &alg_prime_cmp);
#else
  sort_alg_primes (rel->ap, rel->nb_ap);
#endif
    int j = 0;
    for(int i = 0 ; i < rel->nb_ap ; j++) {
        if (i-j) memcpy(rel->ap + j, rel->ap + i, sizeof(alg_prime_t));
        int e = 0;
        for( ; i < rel->nb_ap && alg_prime_cmp(rel->ap+i, rel->ap+j) == 0 ; i++)
            e+=rel->ap[i].e;
        rel->ap[j].e = e;
    }
    rel->nb_ap = j;
}

/* sets a,b. Unless rs->parse_only_ab is set, also fill rs->rel with the
 * (sorted) relation which is read in the input file.
 * Returns the number of characters read in the line, -1 in case of end of
 * file, and 0 in case of an error (when forced_read is non-zero).
 * In case of an error:
 * - if forced_read is zero, fails by an ASSERT
 * - if forced_read is non-zero, return 0.
 *
 * The flag allow_comment allows this function to return also in the case
 * of comments in the relation files.
 */
int relation_stream_get(relation_stream_ptr rs,
        char * supplied_line, size_t supplied_size,
        int forced_read, unsigned int ab_base, int allow_comment)
{
    FILE * f = rs->source;
    char tbuf[RELATION_MAX_BYTES];
    char * line = supplied_line ? supplied_line : tbuf;
    size_t lsize = supplied_line ? supplied_size : sizeof(tbuf);
    int64_t * pa = &rs->rel.a;
    uint64_t * pb = &rs->rel.b;

    static const int ugly[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0,1,2,3,4,5,6,7,8,9,-1, -1, -1, -1, -1, -1,
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

    int c, expected;
    char * p;
    
another_line:
    rs->lnum++;

    p = line;

    *p++ = (c = getc_unlocked(f));
    if (c == EOF) return -1;
    if (c == '\n') goto another_line;
    if (c == '#') {
        for( ; c != EOF && c != '\n' && (size_t) (p-line) < lsize; *p++ = (c=getc_unlocked(f))) ;
        if ((size_t) (p-line) == lsize) {
            fprintf(stderr, "Warning: in %s: encountered very long comment, will have rubbish on next line\n", __func__);
        }
        if (allow_comment) {
            ASSERT_ALWAYS(supplied_line);
            p[-1]='\0'; /* replace newline by '\0' */
            rs->pos += p-line;
            return p-line;
        }

        goto another_line;
    }

    *pa = 0;
    *pb = 0;
    int s = 1;
    int v;
    if (c == '-') { s=-1; *p++ = (c=getc_unlocked(f)); }
    for( ; (v=ugly[(unsigned char) c]) >= 0 ; *p++ = (c=getc_unlocked(f)))
        *pa=*pa*ab_base+v;
    expected = ',';
    if (forced_read && c != expected) {
      for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
      return 0;
    }
    else
      ASSERT_ALWAYS(c == expected);
    *p++ = (c=getc_unlocked(f));
    *pa*=s;
    for( ; (v=ugly[(unsigned char) c]) >= 0 ; *p++ = (c=getc_unlocked(f)))
        *pb=*pb*ab_base+v;
    expected = ':';
    if (forced_read && c != expected) {
      for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
      return 0;
    }
    else
      ASSERT_ALWAYS(c == expected);

    if (!rs->parse_only_ab) {
        /* Do something if we're also interested in primes */
        char * q;
        char * base;
        int n,k;

        base = p;
        n = 1;
        *p++ = (c = getc_unlocked(f));
        for (; c != EOF && c != '\n' && c != ':'; *p++ = (c = getc_unlocked(f)))
            n += c == ',';
        expected = ':';
        if (forced_read && c != expected) {
          for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
          return 0;
        }
        else
          ASSERT_ALWAYS(c == expected);

        relation_provision_for_primes(&rs->rel, n, 0);
        k = 0;
        for (q = base; q != p;) {
            unsigned long pr = 0;
            char * q0 = q;
            for (; (v = ugly[(unsigned char) (c = *q++)]) >= 0;)
                pr = pr * 16 + v;
            if (q>q0+1)
                rs->rel.rp[k++] = (rat_prime_t) { .p = pr,.e = 1};
        }

        rs->rel.nb_rp = k;
        relation_compress_rat_primes(&rs->rel);

        expected = ':';
        if (forced_read && c != expected) {
          for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
          return 0;
        }
        else
          ASSERT_ALWAYS(c == expected);
        ASSERT_ALWAYS(q == p);

        base = p;
        n = 1;
        *p++ = (c = getc_unlocked(f));
        for (; c != EOF && c != '\n' && c != ':'; *p++ = (c = getc_unlocked(f)))
            n += c == ',';
        expected = '\n';
        if (forced_read && c != expected) {
          for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
          return 0;
        }
        else
          ASSERT_ALWAYS(c == expected);

        relation_provision_for_primes(&rs->rel, 0, n);
        k = 0;
        for (q = base; q != p;) {
            unsigned long pr = 0;
            char * q0 = q;
            for (; (v = ugly[(unsigned char) (c = *q++)]) >= 0;)
                pr = pr * 16 + v;
            if (q>q0+1)
                rs->rel.ap[k++] = (alg_prime_t) { .p = pr,.r = -1,.e = 1};
        }

        rs->rel.nb_ap = k;
        relation_compress_alg_primes(&rs->rel);

        expected = '\n';
        if (forced_read && c != expected) {
          for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f))) ;
          return 0;
        }
        else
          ASSERT_ALWAYS(c == expected);
        ASSERT_ALWAYS(q == p);
    }

    /* skip rest of line -- a no-op if we've been told to parse
     * everything. */
    for( ; c != EOF && c != '\n' ; *p++ = (c=getc_unlocked(f)));

    size_t nread =  p-line;
    *p++='\0';

    rs->pos += nread;
    rs->nrels++;

    return nread;
}

/* Same as relation_stream_get, but for a skipped line,
   with supplied_line=NULL and forced_read=0.
   Since this line was already read before, it is necessarily correct. */
int
relation_stream_get_skip (relation_stream_ptr rs)
{
    FILE *f = rs->source;
    int c;
    size_t nread;
    
another_line:
    rs->lnum++;

    c = getc_unlocked (f);
    nread = 1;

    if (c == EOF)
      {
        rs->pos += nread;
        return -1; /* end of file */
      }

    if (c == '#')
      {
        for( ; c != EOF && c != '\n' ; c = getc_unlocked (f), nread++);
        goto another_line;
      }

    /* read entire line */
    for( ; c != EOF && c != '\n' ; c = getc_unlocked (f), nread++);

    rs->pos += nread;
    rs->nrels++;

    return 1;
}

void relation_stream_init(relation_stream_ptr rs)
{
    memset(rs, 0, sizeof(relation_stream));
    rs->t0 = rs->t1 = wct_seconds();
}

void relation_stream_closefile(relation_stream_ptr rs)
{
    ASSERT_ALWAYS(rs->pipe != -1);
    if (rs->source) {
        if (rs->pipe) pclose(rs->source); else fclose(rs->source);
    }
    rs->source = NULL;
    rs->pipe = 0;
}

void relation_stream_clear(relation_stream_ptr rs)
{
    relation_stream_closefile(rs);
    free(rs->rel.rp); rs->rel.rp = NULL;
    free(rs->rel.ap); rs->rel.ap = NULL;
}

void relation_stream_openfile(relation_stream_ptr rs, const char * name)
{
    /* enforce using relation_stream_closefile, which is better practice */
    ASSERT_ALWAYS(rs->source == NULL);
    rs->source = fopen_maybe_compressed2(name, "r", &(rs->pipe), NULL);
    if (rs->source == NULL) {
        fprintf(stderr, "opening %s: %s\n", name, strerror(errno));
        exit(1);
    }
}

void relation_stream_bind(relation_stream_ptr rs, FILE * f)
{
    rs->source = f;
    rs->pipe = -1;
}

void relation_stream_unbind(relation_stream_ptr rs)
{
    rs->source = NULL;
    rs->pipe = 0;
}


int relation_stream_disp_progress_now_p(relation_stream_ptr rs)
{
  static unsigned long change;
  double t;

  if ((rs->nrels >> 18) == (change >> 18))
    return 0;

  t = wct_seconds();
  change = rs->nrels;
  if (t >= rs->t1 + 1)
    {
      rs->dt = t - rs->t0;
      rs->mb_s = rs->dt > 0.01 ? (rs->pos/rs->dt * 1.0e-6) : 0;
      rs->rels_s = rs->dt > 0.01 ? rs->nrels/rs->dt : 0;
      rs->t1 = t;
      return 1;
    }
  return 0;
}

void relation_stream_trigger_disp_progress(relation_stream_ptr rs)
{
    double t = wct_seconds();
    rs->dt = t - rs->t0;
    rs->mb_s = rs->dt > 0.01 ? (rs->pos/rs->dt * 1.0e-6) : 0;
    rs->rels_s = rs->dt > 0.01 ? rs->nrels/rs->dt : 0;
    rs->t1 = t;
}
