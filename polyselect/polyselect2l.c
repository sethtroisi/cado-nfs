/*
  polyselect2l.c is a variant of Paul Zimmermann's polyselect2.c

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-np xxx" denotes the number of special-q's trials for each ad;

  "-lq xxx" denotes the number of small factors (< 251) in the special-q;

  "-maxnorm xxx" only optimize raw polynomials with size <= xxx.
  If the raw polynomial is not good enough, we will still stream
  it to STDERR for further reference.

  Please report bugs to shi.bai AT anu.edu.au.
*/

#include "polyselect2l.h"

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */
#define NEW_ROOTSIEVE
#define MAX_THREADS 16
#define SQ_BATCH_SIZE 10
//#define DEBUG_POLYSELECT2L

#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#endif

/* Read-Only */
uint32_t *Primes = NULL;
unsigned long lenPrimes = 1; // length of Primes[]
int nq = LEN_SPECIAL_Q;
int lq = 1;
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
static int verbose = 0, incr = DEFAULT_INCR, default_MAX_k;
char *out = NULL; /* output file for msieve input (msieve.dat.m) */
cado_poly best_poly, curr_poly;
double best_E = 0.0; /* Murphy's E (the larger the better) */
int seed = 0; /* seed */

/* read-write global variables */
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */
int tot_found = 0; /* total number of polynomials */
int found = 0; /* number of polynomials below maxnorm */
double potential_collisions = 0.0, aver_lognorm = 0.0;
unsigned long collisions = 0;
unsigned long collisions_good = 0;
double total_adminus2;
double best_logmu[11];
double rootsieve_time = 0.0;


/* crt, set r and qqz */
void
crt_sq ( mpz_t qqz,
         mpz_t r,
         unsigned long *q,
         unsigned long *rq )
{
  mpz_t prod, pprod, mod, inv, sum;
  int i;
  unsigned long qq[lq];

  mpz_init_set_ui (prod, 1);
  mpz_init (pprod);
  mpz_init (mod);
  mpz_init (inv);
  mpz_init_set_ui (sum, 0);

  for (i = 0; i < lq; i ++) {
    qq[i] = q[i] * q[i];
    mpz_mul_ui (prod, prod, qq[i]);
  }

  for (i = 0; i < lq; i ++) {
    mpz_divexact_ui (pprod, prod, qq[i]);
    mpz_set_ui (mod, qq[i]);
    mpz_invert (inv, pprod, mod);
    mpz_mul_ui (inv, inv, rq[i]);
    mpz_mul (inv, inv, pprod);
    mpz_add (sum, sum, inv);
  }

  mpz_mod (sum, sum, prod);
  mpz_set (r, sum);
  mpz_set (qqz, prod);

  mpz_clear (prod);
  mpz_clear (pprod);
  mpz_clear (mod);
  mpz_clear (inv);
  mpz_clear (sum);
}


/* print poly info */
void
print_poly_info ( mpz_t *f,
                  unsigned int d,
                  mpz_t g[2] )
{
  unsigned int i, nroots;
  double skew, logmu, alpha;

  gmp_printf ("Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
  for (i = d + 1; i -- != 0; )
    gmp_printf ("c%u: %Zd\n", i, f[i]);

  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  printf ("# lognorm %1.2f, skew %1.2f, alpha %1.2f, E %1.2f,  exp_E %1.2f, %u rroots\n",
          logmu, skew, alpha, logmu + alpha,
          logmu - sqrt (2.0 * exp_rot[d] * log (skew)),
          nroots);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       unsigned long ad, unsigned int d, mpz_t N, unsigned long q,
       mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], qq;
  unsigned int j;
  int cmp;
  double skew, logmu, E;
  /* the expected rotation space is S^5 for degree 6 */

#ifdef DEBUG_POLYSELECT2L
  gmp_printf ("Found match: (%lu,%ld) (%lu,%ld) for ad=%lu, q=%lu, rq=%Zd\n",
              p1, i, p2, i, ad, q, rq);
  gmp_printf ("m0=%Zd\n", m0);
#endif

  mpz_init (l);
  mpz_init (m);
  mpz_init (t);
  mpz_init (k);
  mpz_init (qq);
  mpz_init (adm1);
  mpz_init (mtilde);
  mpz_init (g[0]);
  mpz_init (g[1]);
  mpz_init (gold[0]);
  mpz_init (gold[1]);
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (f == NULL || fold == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in match\n");
    exit (1);
  }
  for (j = 0; j <= d; j++) {
    mpz_init (f[j]);
    mpz_init (fold[j]);
  }
  /* we have l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_mul_ui (l, l, q);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set_si (mtilde, i);
  mpz_set_ui (qq, q);
  mpz_mul_ui (qq, qq, q);
  mpz_mul (mtilde, mtilde, qq);
  mpz_add (mtilde, mtilde, rq);
  mpz_add (mtilde, mtilde, m0);
  /* we want mtilde = d*ad*m + a_{d-1}*l with 0 <= a_{d-1} < d*ad.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_set_ui (m, ad);
  mpz_mul_ui (m, m, d);
  if (mpz_invert (adm1, l, m) == 0)
  {
    fprintf (stderr, "Error in 1/l mod (d*ad)\n");
    exit (1);
  }
  mpz_mul (adm1, adm1, mtilde);
  mpz_mod (adm1, adm1, m);
  mpz_mul (m, adm1, l);
  mpz_sub (m, mtilde, m);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_ui_p (m, d) == 0)
  {
    fprintf (stderr, "Error: m-a_{d-1}*l not divisible by d\n");
    exit (1);
  }
#endif
  mpz_divexact_ui (m, m, d);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_ui_p (m, ad) == 0)
  {
    fprintf (stderr, "Error: (m-a_{d-1}*l)/d not divisible by ad\n");
    exit (1);
  }
#endif
  mpz_divexact_ui (m, m, ad);
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set_ui (f[d], ad);
  mpz_pow_ui (t, m, d);
  mpz_mul_ui (t, t, ad);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif

  mpz_divexact (t, t, l);
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (j = d - 2; j > 0; j--)
  {
#ifdef DEBUG_POLYSELECT2L
    if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
#endif
    mpz_divexact (t, t, l);
    /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
    mpz_pow_ui (mtilde, m, j);
    mpz_fdiv_q (adm1, t, mtilde); /* t -> adm1 * mtilde + t */
    mpz_invert (k, mtilde, l); /* search adm1 + k such that
                                  t = (adm1 + k) * m^j mod l */
    mpz_mul (k, k, t);
    mpz_sub (k, k, adm1);
    mpz_mod (k, k, l);

    mpz_mul_2exp (k, k, 1);
    cmp = mpz_cmp (k, l);
    mpz_div_2exp (k, k, 1);
    if (cmp >= 0)
      mpz_sub (k, k, l);
    mpz_add (adm1, adm1, k);
    mpz_set (f[j], adm1);
    /* subtract adm1*m^j */
    mpz_submul (t, mtilde, adm1);
  }
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif
  mpz_divexact (t, t, l);
  mpz_set (f[0], t);

  /* information on all polynomials */
  pthread_mutex_lock (&lock);
  total_adminus2 += (double) mpz_sizeinbase (f[d-2], 2);
  collisions ++;
  pthread_mutex_unlock (&lock);
  // gmp_printf ("# a_{d-2}=%Zd\n", f[d-2]);

  /* save unoptimized polynomial to fold */
  for (i = d + 1; i -- != 0; )
    mpz_set (fold[i], f[i]);
  mpz_set (gold[1], g[1]);
  mpz_set (gold[0], g[0]);

  /* old lognorm */
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

  /* if the polynomial has norm < "-maxnorm", we optimize it */
  if (logmu <= max_norm)
  {

#ifdef DEBUG_POLYSELECT2L
    /* unoptimized poly */
    printf ("# Raw polynomial:\n");
    print_poly_info (fold, d, gold);
#endif

    /* optimize size */
    optimize (f, d, g, 0, 1);

/* root sieve */
#ifndef NEW_ROOTSIEVE
    unsigned long alim = 2000;
    long jmin, kmin;
#endif
    mpz_neg (m, g[0]);

    rootsieve_time -= seconds_thread ();

#ifdef NEW_ROOTSIEVE
    if (d > 4) {
      ropt_polyselect (f, d, m, g[1], N, MAX_k, 0); // verbose = 2 to see details.
      mpz_neg (g[0], m);
    }
    else {
      unsigned long alim = 2000;
      long jmin, kmin;
      rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
      mpz_neg (g[0], m);
      /* optimize again, but only translation */
      optimize_aux (f, d, g, 0, 0, CIRCULAR);
    }
#else
    rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
    mpz_neg (g[0], m);
    /* optimize again, but only translation */
    optimize_aux (f, d, g, 0, 0, CIRCULAR);
#endif

    rootsieve_time += seconds_thread ();

    skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

    for (i = 10; i > 0 && logmu < best_logmu[i-1]; i--)
      best_logmu[i] = best_logmu[i-1];
    best_logmu[i] = logmu;

    pthread_mutex_lock (&lock);
    collisions_good ++;
    aver_lognorm += logmu;
    tot_found ++;
    pthread_mutex_unlock (&lock);

    /* MurphyE */
    mpz_set (curr_poly->rat->f[0], g[0]);
    mpz_set (curr_poly->rat->f[1], g[1]);
    for (j = d + 1; j -- != 0; )
      mpz_set (curr_poly->alg->f[j], f[j]);
    curr_poly->skew = skew;
    E =  MurphyE (curr_poly, BOUND_F, BOUND_G, AREA, MURPHY_K);

    mpz_neg (m, g[0]);

    pthread_mutex_lock (&lock);
    if (E > best_E)
    {
      best_E = E;
      cado_poly_set (best_poly, curr_poly);
    }
    if (out != NULL) /* msieve output */
    {
      FILE *fp;
      fp = fopen (out, (found == 0) ? "w" : "a");
      if (fp == NULL)
      {
        fprintf (stderr, "Error, cannot open file %s\n", out);
        exit (1);
      }
      fprintf (fp, "0");
      for (j = d + 1; j -- != 0; )
        gmp_fprintf (fp, " %Zd", f[j]);
      mpz_neg (m, g[0]);
      gmp_fprintf (fp, " %Zd %Zd\n", g[1], m);
      fclose (fp);
    }
    found ++;
    pthread_mutex_unlock (&lock);

    /* print optimized polynomial */
#ifdef DEBUG_POLYSELECT2L
    gmp_printf ("# Optimized polynomial:\n");
#endif

    gmp_printf ("n: %Zd\n", N);
    print_poly_info ( f, d, g );
    printf ("# Murphy's E(Bf=%.0f,Bg=%.0f,area=%.2e)=%1.2e (best so far %1.2e)\n",
            BOUND_F, BOUND_G, AREA, E, best_E);
    printf ("\n");
    fflush (stdout);
  }
  else {
    gmp_fprintf (stderr, "# Skip polynomial: %.2f, ad: %lu, l: %Zd, m: %Zd\n", logmu, ad, l, m);
  }

  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (t);
  mpz_clear (k);
  mpz_clear (qq);
  mpz_clear (adm1);
  mpz_clear (mtilde);
  mpz_clear (g[0]);
  mpz_clear (g[1]);
  mpz_clear (gold[0]);
  mpz_clear (gold[1]);
  for (j = 0; j <= d; j++) {
    mpz_clear (f[j]);
    mpz_clear (fold[j]);
  }
  free (f);
  free (fold);
}


/* find collisions between "P" primes */
static inline void
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long i, j, nprimes, p, *rp, nrp;
  long ppl = 0;
  double pc1;
  mpz_t *f, tmp;

  /* init f for roots computation */
  mpz_init_set_ui (tmp, 0);
  f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
  if (f == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }
  for (i = 0; i <= header->d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[header->d], 1);

  rp = (unsigned long*) malloc (header->d * sizeof (unsigned long));
  if (rp == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }

  hash_t H;
  hash_init (H, 2.2 * lenPrimes);

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    ppl = (long) (p*p);

    /* add faked roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    mpz_mod_ui (f[0], header->Ntilde, p);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrp = poly_roots_ulong (rp, f, header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);

    for (j = 0; j < nrp; j++) {
      /* only consider r[j] and r[j] - pp */
      hash_add (H, p, rp[j], header->m0, header->ad, header->d, header->N, 1, tmp);
      hash_add (H, p, rp[j] - ppl, header->m0, header->ad, header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", cputime () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

  /* if the hash table contains n entries, each one smaller than (2P)^2,
     the number of potential collisions is about 1/2n^2/(2P)^2 */
  pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 2], 2.0);
  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

  pthread_mutex_lock (&lock);
  potential_collisions += pc1;
  pthread_mutex_unlock (&lock);
}


/* collision on each special-q */
static inline void
collision_on_each_sq ( header_t header,
                       proots_t R,
                       unsigned long q,
                       mpz_t rqqz,
                       unsigned long *inv_qq ) // inv_qq is of length #Primes
{
  unsigned int nr, j;
  unsigned long nprimes, p, rp, pp;
  long ppl, u;
  mpz_t rppz;

  mpz_init (rppz);

#ifdef DEBUG_POLYSELECT2L
  mpz_t tmp_debug, tmp_debug2;
  mpz_init_set (tmp_debug, header->Ntilde);
  mpz_mod (tmp_debug, tmp_debug, qqz);
  mpz_init_set (tmp_debug2, header->m0);
  mpz_add (tmp_debug2, tmp_debug2, rqqz);
  mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
  mpz_mod (tmp_debug2, tmp_debug2, qqz);
  if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
    fprintf (stderr, "Error: crt root is wrong in collision_on_each_sq\n");
    exit (1);
  }
  mpz_clear (tmp_debug);
  mpz_clear (tmp_debug2);
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  hash_t H;
  hash_init (H,  2.2 * lenPrimes);

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;
    /* set p, p^2, ppl */
    pp = p * p;
    ppl = (long) pp;
    nr = R->nr[nprimes];

    for (j = 0; j < nr; j++)
    {
      rp = R->roots[nprimes][j];
      mpz_set_ui (rppz, rp);
      mpz_sub (rppz, rppz, rqqz);
      mpz_mul_ui (rppz, rppz, inv_qq[nprimes]);
      u = (long) mpz_fdiv_ui (rppz, pp);

#ifdef DEBUG_POLYSELECT2L
      mpz_t Ntmp, m0tmp, tmp, qqz;
      mpz_init_set_ui (qqz, q);
      mpz_mul_ui (qqz, qqz, q);
      mpz_init_set (Ntmp, header->Ntilde);
      mpz_mod_ui (Ntmp, Ntmp, ppl); // tildeN (mod p^2)

      mpz_init_set (m0tmp, header->m0);
      mpz_add (m0tmp, m0tmp, rqqz); // m0 + rq
      mpz_init_set_ui (tmp, u);
      mpz_mul (tmp, tmp, qqz); // i*q^2
      mpz_add (m0tmp, m0tmp, tmp); // m0 + rq + i*q^2

      mpz_pow_ui (m0tmp, m0tmp, header->d);
      mpz_mod_ui (m0tmp, m0tmp, ppl);
      if (mpz_cmp (m0tmp, Ntmp) != 0) {
        fprintf (stderr, "Error: i computation is wrong in collision_on_each_sq\n");

        gmp_printf ("Details: (p=%lu, i(mod p^2)=%ld) for ad=%lu, q=%lu, rq=%Zd, rp=%lu\n",
                    p, u, header->ad, q, rqqz, rp);
        gmp_printf ("m0=%Zd\n", header->m0);
        gmp_printf ("Ntilde=%Zd\n", header->Ntilde);

        exit (1);
      }
      mpz_clear (Ntmp);
      mpz_clear (m0tmp);
      mpz_clear (tmp);
      mpz_clear (qqz);
#endif

      hash_add (H, p, u, header->m0, header->ad, header->d, header->N, q, rqqz);
      hash_add (H, p, u - ppl, header->m0, header->ad, header->d, header->N, q, rqqz);

    }  // next rp

  } // next p

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_each_sq took %dms\n", cputime () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  double pc2 = 0.5 * pow ((double) H->size / (double) Primes[nprimes-1],
                          2.0);

  mpz_clear (rppz);
  hash_clear (H);

  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* collision on batch of special-q's */
static inline void
collision_on_batch_sq ( header_t header,
                        proots_t R,
                        unsigned long *q,
                        mpz_t *qqz,
                        mpz_t *rqqz,
                        unsigned long size )
{
  if (size == 0)
    return;

  unsigned int i;
  unsigned long nprimes, p, pp, tmpul;
  unsigned long **invqq = malloc (size * sizeof (unsigned long *));

  if (invqq) {
    for (i = 0; i < size; i++)
      invqq[i] = malloc (lenPrimes * sizeof (unsigned long));
  }
  else {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

  mpz_t tmp, qprod[size];

  mpz_init (tmp);
  for (i = 0; i < size; i++)
    mpz_init (qprod[i]);

#ifdef DEBUG_POLYSELECT2L  // check roots
  for (i = 0; i < size; i ++) {
    mpz_t tmp_debug, tmp_debug2;
    mpz_init_set (tmp_debug, header->Ntilde);
    mpz_mod (tmp_debug, tmp_debug, qqz[i]);
    mpz_init_set (tmp_debug2, header->m0);
    mpz_add (tmp_debug2, tmp_debug2, rqqz[i]);
    mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
    mpz_mod (tmp_debug2, tmp_debug2, qqz[i]);
    if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
      fprintf (stderr, "Error: crt root is wrong in collision_on_each_sq\n");
      exit (1);
    }
    mpz_clear (tmp_debug);
    mpz_clear (tmp_debug2);
    fprintf (stderr, "i: %lu OK\n", i);
  }
#endif

  // (size -1) multiplications
  mpz_set (qprod[0], qqz[0]);
  for (i = 1; i < size; i ++)
    mpz_mul(qprod[i], qqz[i], qprod[i-1]);

  // int st = cputime();

  /* Step 1: batch inversion */
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    pp = p*p;
    if ((header->d * header->ad) % p == 0)
      continue;

    // the inversion, also reduce qprod[size-1] (mod pp).
    mpz_set (tmp, qprod[size-1]);
    mpz_mod_ui (tmp, tmp, pp);
    tmpul =  mpz_get_ui (tmp);
    tmpul = invert (tmpul, pp);

    //gmp_fprintf (stderr, "primes[%lu]: %lu, %Zd, %lu\n", nprimes, p, qqz, size);

    // for each (q, r) \in a batch
    for (i = size-1; i > 0; i --) {
      mpz_mul_ui(tmp, qprod[i-1], tmpul);
      invqq[i][nprimes] = mpz_fdiv_ui (tmp, pp);
      mpz_mul_ui (tmp, qqz[i], tmpul);
      tmpul = mpz_fdiv_ui (tmp, pp);
    }
    invqq[0][nprimes] = tmpul;

#ifdef DEBUG_POLYSELECT2L // check inversions
    for (i = 0; i < size; i ++) {
      mpz_t tmp_debug, tmp_debug2;
      mpz_init_set_ui (tmp_debug, invqq[i][nprimes]);
      mpz_init (tmp_debug2);
      mpz_mul (tmp_debug2, tmp_debug, qqz[i]);
      unsigned long re = mpz_fdiv_ui (tmp_debug2, pp);
      if (re != 1) {
        fprintf (stderr, "Error, batch inversion is wrong in %s, iter: %d\n", __FUNCTION__, i);
        exit (1);
      }
      mpz_clear (tmp_debug);
      mpz_clear (tmp_debug2);
    }
#endif

  } // next prime

  //fprintf (stderr, "# one batch inversion took %dms\n", cputime () - st);

  /* Step 2: find collisions on q. */
  for (i = 0; i < size; i ++) {

    // int st = cputime();
    collision_on_each_sq ( header,
                           R,
                           q[i],
                           rqqz[i],
                           invqq[i] );

    // printf ("# collision_on_each_sq took %dms\n", cputime () - st);

  }

  for (i = 0; i < size; i++)
    mpz_clear (qprod[i]);
  mpz_clear (tmp);
  for (i = 0; i < size; i++)
    free (invqq[i]);
  free (invqq);
}


/* collision on special-q, call collisio_on_batch_sq */
static inline void
collision_on_sq ( header_t header,
                  proots_t R )
{
  // init special-q roots
  qroots_t SQ_R;
  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  // qroots_print (SQ_R);

  unsigned long K = lq, N = SQ_R->size, tot, i, l;
  unsigned long idx_q[K], q[SQ_BATCH_SIZE];
  mpz_t qqz[SQ_BATCH_SIZE], rqqz[SQ_BATCH_SIZE];

  for (l = 0; l < SQ_BATCH_SIZE; l++) {
    mpz_init (qqz[l]);
    mpz_init (rqqz[l]);
  }

  // less than lq special primes having roots for this ad
  if (N == 0 || N < K)
    return;

  tot =  binom (N, K);
  //fprintf (stderr, "n=%lu, k=%lu, (n,k)=%lu, nq:%d\n", N, K, tot, nq);

  // if tot (n, k) < wanted, use tot as wanted
  if (tot > (unsigned long) nq)
    tot = (unsigned long) nq;

  i = 0;
  while ( i <= (tot-SQ_BATCH_SIZE) ) {

    l = i; // why do I use an extra l here?
    if (l == 0) {

      // enumerate first combination
      first_comb (K, idx_q);
      // print_comb (k, idx_q);
      q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

      for (l = 1; l < SQ_BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
      }
    }
    else {
      for (l = 0; l < SQ_BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
      }
    }

#ifdef DEBUG_POLYSELECT2L
    for (j = 0; j < SQ_BATCH_SIZE; j++)
      gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n", q[j], qqz[j], rqqz[j]);
#endif

    // collision batch
    collision_on_batch_sq ( header,
                            R,
                            q,
                            qqz,
                            rqqz,
                            SQ_BATCH_SIZE );

    i += SQ_BATCH_SIZE;
  }

  // tail batch
  for (l = 0; l < (tot % SQ_BATCH_SIZE); l++) {
    next_comb (N, K, idx_q);
    q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

#ifdef DEBUG_POLYSELECT2L
    gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n", q[l], qqz[l], rqqz[l]);
#endif

  }

  collision_on_batch_sq ( header,
                          R,
                          q,
                          qqz,
                          rqqz,
                          tot % SQ_BATCH_SIZE );

  for (l = 0; l < SQ_BATCH_SIZE; l++) {
    mpz_clear (qqz[l]);
    mpz_clear (rqqz[l]);
  }
  qroots_clear (SQ_R);
}


static void
newAlgo (mpz_t N, unsigned long d, unsigned long ad)
{

  header_t header;
  header_init (header, N, d, ad);

  proots_t R;
  proots_init (R, lenPrimes);

  collision_on_p (header, R);

  collision_on_sq (header, R);

  proots_clear (R, lenPrimes);
  header_clear (header);
}


void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
  return NULL;
}

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */


static void
usage (char *argv)
{
  fprintf (stderr, "Usage: %s [options] P\n", argv);
  fprintf (stderr, "Parameters and options:\n");
  fprintf (stderr, "P            --- degree-1 coefficient of g(x) has\n");
  fprintf (stderr, "                 two prime factors in [P,2P]\n");
  fprintf (stderr, "-v           --- verbose mode\n");
  fprintf (stderr, "-q           --- quiet mode\n");
  fprintf (stderr, "-t nnn       --- use n threads (default 1)\n");
  fprintf (stderr, "-admin nnn   --- start from ad=nnn (default 0)\n");
  fprintf (stderr, "-admax nnn   --- stop at ad=nnn\n");
  fprintf (stderr, "-incr nnn    --- forced factor of ad (default 60)\n");
  fprintf (stderr, "-N nnn       --- input number\n");
  fprintf (stderr, "-degree nnn  --- wanted polynomial degree\n");
  fprintf (stderr, "-nq nnn      --- number of special-q's considered for each coefficient a_d\n");
  fprintf (stderr, "-lq nnn      --- number of factors in the special-q\n");
  fprintf (stderr, "-seed nnn    --- seed for srand()\n");
  fprintf (stderr, "-kmax nnn    --- rotation bound (default %d)\n",
           default_MAX_k);
  fprintf (stderr, "-save xxx    --- save state in file xxx\n");
  fprintf (stderr, "-resume xxx  --- resume state from file xxx\n");
  fprintf (stderr, "-maxnorm xxx --- only optimize polynomials with norm <= xxx\n");
  fprintf (stderr, "-maxtime xxx --- stop the search after xxx seconds\n");
  fprintf (stderr, "-out xxx     --- for msieve-format output\n");
  exit (1);
}


int
main (int argc, char *argv[])
{
  int argc0 = argc;
  char **argv0 = argv, *save = NULL, *resume = NULL;
  double st0 = seconds (), maxtime = DBL_MAX;
  mpz_t N;
  unsigned int d = 0;
  unsigned long P, admin = 0, admax = ULONG_MAX;
  int tries = 0, i, nthreads = 1, st, target_time = TARGET_TIME;
  tab_t *T;
  FILE *fp;
#ifdef MAX_THREADS
  pthread_t tid[MAX_THREADS];
#endif

  /* printf command-line */
  printf ("#");
  for (i = 0; i < argc; i++)
    printf (" %s", argv[i]);
  printf ("\n");
  fflush (stdout);

  mpz_init (N);
  default_MAX_k = MAX_k;
  cado_poly_init (best_poly);
  cado_poly_init (curr_poly);

  /* read params */
  if (argc == 1)
    usage (argv0[0]);

  while (argc >= 2 && argv[1][0] == '-')
  {
    if (strcmp (argv[1], "-t") == 0)
    {
      nthreads = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-maxnorm") == 0)
    {
      max_norm = atof (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-maxtime") == 0)
    {
      maxtime = atof (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-admin") == 0)
    {
      double d;
      d = strtod (argv[2], NULL);
      if (d > (double) ULONG_MAX)
      {
        fprintf (stderr, "Error, too large value of admin\n");
        exit (1);
      }
      admin = (unsigned long) d;
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-admax") == 0)
    {
      double d;
      d = strtod (argv[2], NULL);
      if (d > (double) ULONG_MAX)
      {
        fprintf (stderr, "Error, too large value of admax\n");
        exit (1);
      }
      admax = (unsigned long) d;
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-incr") == 0)
    {
      incr = strtoul (argv[2], NULL, 10);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-N") == 0)
    {
      mpz_set_str (N, argv[2], 10);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-degree") == 0)
    {
      d = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-nq") == 0)
    {
      nq = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-lq") == 0)
    {
      lq = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-seed") == 0)
    {
      seed = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-save") == 0)
    {
      save = argv[2];
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-resume") == 0)
    {
      resume = argv[2];
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-out") == 0)
    {
      out = argv[2];
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 3 && strcmp (argv[1], "-kmax") == 0)
    {
      MAX_k = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (strcmp (argv[1], "-v") == 0)
    {
      verbose ++;
      argv += 1;
      argc -= 1;
    }
    else if (strcmp (argv[1], "-q") == 0)
    {
      verbose = -1;
      argv += 1;
      argc -= 1;
    }
    else
    {
      fprintf (stderr, "Invalid option: %s\n", argv[1]);
      exit (1);
    }
  }

  /* verify params */
  if (mpz_cmp_ui (N, 0) == 0)
  {
    int ret;

    ret = gmp_scanf ("n: %Zd\n", N);
    if (ret != 1)
    {
      fprintf (stderr, "Error, input number N cannot be read\n");
      exit (1);
    }
  }

  if (lq < 1 || nq < 1)
  {
    fprintf (stderr, "Error, number of factors in special-q should >= 1 and/or number of special-q's should >=1\n");
    exit (1);
  }

  if (seed < 1) {
    seed = time(NULL);
    srand(seed);
  }
  else
    srand(seed);

  if (lq == 1)
  {
    nq = LEN_SPECIAL_Q;
  }

  /* set cpoly */
  mpz_set (best_poly->n, N);
  mpz_set (curr_poly->n, N);

  if (argc != 2)
    usage (argv0[0]);

#ifdef MAX_THREADS
  if (nthreads > MAX_THREADS)
  {
    fprintf (stderr, "Error, nthreads should be <= %d\n", MAX_THREADS);
    exit (1);
  }
#endif

  if (mpz_cmp_ui (N, 0) <= 0)
  {
    fprintf (stderr, "Error, missing input number (-N option)\n");
    exit (1);
  }

  if (d == 0)
  {
    fprintf (stderr, "Error, missing degree (-d option)\n");
    exit (1);
  }

  best_poly->alg->degree = d;
  best_poly->rat->degree = 1;
  curr_poly->alg->degree = d;
  curr_poly->rat->degree = 1;

  if (resume != NULL)
  {
    fp = fopen (resume, "r");
    if (fp == NULL)
    {
      fprintf (stderr, "Cannot open resume file %s\n", resume);
      exit (1);
    }
    if (fscanf (fp, "%lu", &admin) != 1)
    {
      fprintf (stderr, "Cannot read ad value from resume file %s\n",
               resume);
      exit (1);
    }
    fclose (fp);
  }

  for (i = 0; i <= 10; i++)
    best_logmu[i] = 999.9;

  /* init primes */
  P = atoi (argv[1]);
  st = cputime ();
  lenPrimes = initPrimes (P, &Primes);

  printf ( "# Info: initializing %lu P primes took %dms, seed=%d\n",
           lenPrimes,
           cputime () - st,
           seed );
  printf ( "# Info: estimated peak memory=%.2fMB (%d threads, batch %d inversions)\n",
           (double) (nthreads * SQ_BATCH_SIZE * lenPrimes * (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           SQ_BATCH_SIZE );

  //printPrimes (Primes, lenPrimes);

  /* init tabs_t for threads */
  T = malloc (nthreads * sizeof (tab_t));
  if (T == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in main\n");
    exit (1);
  }
  for (i = 0; i < nthreads ; i++)
  {
    mpz_init_set (T[i]->N, N);
    T[i]->d = d;
    T[i]->thread = i;
  }

  if (incr <= 0)
  {
    fprintf (stderr, "Error, incr should be positive\n");
    exit (1);
  }

  /* force admin to be divisible by incr */
  if (admin == 0)
    admin = incr;
  admin = ((admin + incr - 1) / incr) * incr; /* incr * ceil (admin/incr) */

  while (admin <= admax && seconds () - st0 <= maxtime)
  {
    for (i = 0; i < nthreads ; i++)
    {
      tries ++;
      if (verbose >= 1)
      {
        gmp_printf ("%d ad=%lu\r", tries, admin);
        fflush (stdout);
      }
      T[i]->ad = admin;
#ifndef MAX_THREADS
      newAlgo (N, d, admin);
#else
      pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
#endif
      admin += incr;
    }
#ifdef MAX_THREADS
    for (i = 0 ; i < nthreads ; i++)
      pthread_join(tid[i], NULL);
#endif

    if (save != NULL)
    {
      fp = fopen (save, "w");
      if (fp == NULL)
      {
        fprintf (stderr, "Cannot open save file %s\n", save);
        exit (1);
      }
      if (fprintf (fp, "%lu\n", admin) < 0)
      {
        fprintf (stderr, "Cannot print ad value to save file %s\n",
                 save);
        exit (1);
      }
      fclose (fp);
    }

    if (cputime () > target_time || verbose > 0)
    {
      printf ("# ad=%lu time=%dms exp. coll.=%1.1f, got %lu (%0.3f/s) with %lu good ones, av. lognorm=%1.2f\n",
              admin,
              cputime (),
              potential_collisions,
              collisions,
              1000.0 * (double) collisions / cputime (),
              collisions_good,
              aver_lognorm / collisions );
      fflush (stdout);
      target_time += TARGET_TIME;
    }
  }

  printf ("# Tried %d ad-value(s), found %d polynomial(s), %d below maxnorm\n",
          tries, tot_found, found);
  printf ("# potential collisions=%1.2e (%1.2e/s)\n",
          potential_collisions, 1000.0 * potential_collisions
          / (double) cputime ());
  printf ("# av. adm2:%1.0f", total_adminus2 / (double) collisions);

  /* print best 10 values of logmu */
  printf ("# best logmu:");
  for (i = 0; i < 10; i++)
    printf (" %1.2f", best_logmu[i]);
  printf ("\n");

  /* print total time (format for cpu_time.sh) */
  printf ("# Total phase took %.2fs\n", seconds () - st0);
  printf ("# Rootsieve took %.2fs\n", rootsieve_time);


  if (best_E == 0.0)
    printf ("No polynomial found, please increase the ad range or decrease P\n");
  else
    print_cadopoly_extra (stdout, best_poly, argc0, argv0, st0, 1 /* raw */);

  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  mpz_clear (N);
  clearPrimes (&Primes);
  cado_poly_clear (best_poly);
  cado_poly_clear (curr_poly);

  return 0;
}
