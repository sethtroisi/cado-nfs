/* TODO: look at http://sourceforge.net/p/msieve/code/HEAD/tree/branches/RDS/cofactorize_siqs.c */

/* tiny MPQS implementation, specially tuned for 64- to 128-bit input */

// #define TRACE -23830
// #define TRACE_P 937

/* if CHECK_OVERFLOW is defined, check that no overflow occurs in the
   additions on unsigned chars (0..255) */
// #define CHECK_OVERFLOW

/* minimum difference between number of relations and factor base size */
#define WANT_EXCESS 11

/* number of small primes we skip (should be >= 1 since we always skip 2) */
#define SKIP 10

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <sys/types.h>
#include <float.h> /* for DBL_MAX */
#include "gmp.h"
#include "timing.h"
#include "macros.h"
#include "mod_ul_default.h"
#include "modredc_ul.h"
#include "gmp_aux.h"

typedef struct {
  unsigned int p;
  unsigned int r;      /* root of x^2 = k*N (mod p) */
  uint64_t invp;       /* 1/p mod 2^64, needed for trial division */
  unsigned char logp;  /* log(p) / log(radix), needed for sieve */
  unsigned long invp2; /* floor(ULONG_MAX/p), needed for trial division */
  unsigned int Mp;     /* M mod p */
} fb_t;

typedef struct {
  unsigned int alloc; /* size of the hash table */
  unsigned int size;  /* number of values in the table */
  unsigned int *p;    /* values of large primes */
  mpz_t *row;         /* bit-vector over the factor base */
  mpz_t *x;           /* value of a*i+b */
  mpz_t *y;           /* value of x^2-N */
} hash_struct;
typedef hash_struct hash_t[1];

typedef struct
{
  mpz_t z;              /* P[0] * P[1] * ... * P[n-1] */
  mpz_t **x;            /* product tree */
  mpz_t *y;             /* smooth parts of x[0][] */
  mpz_t *axb;           /* values of a*i+b */
  unsigned int logm;
  unsigned int m;       /* m = 2^logm */
  unsigned int size;    /* number of elements stored in x[0][] */
} bernstein_struct;
typedef bernstein_struct bernstein_t[1];

void init_tree (bernstein_t, unsigned int, unsigned long *, unsigned int);
void clear_tree (bernstein_t);
void accumulate (bernstein_t, mpz_t, mpz_t);
void get_smooth (bernstein_t);

/* For i < 50, isprime_table[i] == 1 iff i is prime */
static unsigned char isprime_table[] = {
0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0,
0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

static const size_t isprime_table_size =
    sizeof(isprime_table) / sizeof(isprime_table[0]);

static int
jacobi (unsigned long a, unsigned long b)
{
  mpz_t aa, bb;
  int res;

  mpz_init_set_ui (aa, a);
  mpz_init_set_ui (bb, b);
  res = mpz_jacobi (aa, bb);
  mpz_clear (aa);
  mpz_clear (bb);
  return res;
}

/* return b^e mod n, assuming n has at most 32 bits */
static uint64_t
mod_pow_uint64 (uint64_t b, uint64_t e, uint64_t n)
{
  uint64_t r = 1, f = 1, y = b;

  while (f <= e)
    {
      if (e & f)
        r = (r * y) % n;
      f <<= 1;
      y = (y * y) % n;
    }

  return r;
}

/* Uses Tonelli-Shanks, more precisely Algorithm from Table 1 in [1].
   [1] Adleman-Manders-Miller Root Extraction Method Revisited, Zhengjun Cao,
   Qian Sha, Xiao Fan, 2011, http://arxiv.org/abs/1111.4877.
   Solve x^2 = rr (mod p).
*/
static uint64_t
tonelli_shanks (uint64_t rr, uint64_t p)
{
  uint64_t q, s, i, j, l;
  uint64_t aa, hh, delta, dd, bb, zz;

  if (p == 2)
    return rr;

  ASSERT (p <= 4294967295UL);

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);

  zz = 1;
  i = 1;
  do {
    /* zz is equal to i (mod p) */
    zz += 1;
    i++;
  } while (i < isprime_table_size && (!isprime_table[i] || jacobi (zz, p) != -1));
  ASSERT_ALWAYS (i < isprime_table_size);
  aa = mod_pow_uint64 (zz, q, p);
  bb = mod_pow_uint64 (rr, q, p);
  hh = 1;
  for (j = 1; j < s; j++)
    {
      dd = bb;
      for (l = 0; l < s - 1 - j; l++)
        dd = (dd * dd) % p;
      if (dd != 1)
        {
          hh = (hh * aa) % p;
          aa = (aa * aa) % p;
          bb = (bb * aa) % p;
        }
      else
        aa = (aa * aa) % p;
    }
  delta = mod_pow_uint64 (rr, (q + 1) >> 1, p);
  hh = (hh * delta) % p;

  return hh;
}

/* Given k1 containing initially r such that r^2 = N (mod p),
   return k1 and k2 which are the roots of (a*x+b)^2 = N (mod p),
   i.e., k1 = (r-b)/a+M (mod p) and k2 = (-r-b)/a+M (mod p).
   Assume p is odd, and inva = 1/sqrt(a) mod p.
   Ensures k1 <= k2 at the end.
*/
static unsigned long
findroot (unsigned long *k2, unsigned long bmodp, unsigned long p,
          unsigned long k1, unsigned long inva, unsigned long Mp)
{
  /* the two roots are (k1-b)/a and (-k1-b)/a */
  modulus_t pp;
  residueul_t tt, uu, vv;
  modul_initmod_ul (pp, p);
  modul_init (tt, pp);
  modul_init (uu, pp);
  modul_init (vv, pp);
  modul_set_ul (tt, inva, pp);
  modul_set_ul_reduced (vv, k1, pp);
  modul_neg (uu, vv, pp);
  modul_sub_ul (uu, uu, bmodp, pp); /* -r-b */
  modul_mul (uu, uu, tt, pp);       /* (-r-b)/a */
  *k2 = mod_get_ul (uu, pp);
  *k2 += Mp;                        /* (-r-b)/a + M */
  if (*k2 >= p)
    *k2 -= p;
  modul_sub_ul (uu, vv, bmodp, pp); /* r-b */
  modul_mul (uu, uu, tt, pp);       /* (r-b)/a */
  k1 = mod_get_ul (uu, pp);
  k1 += Mp;                         /* (r-b)/a + M */
  if (k1 >= p)
    k1 -= p;
  modul_clear (tt, pp);
  modul_clear (uu, pp);
  modul_clear (vv, pp);
  modul_clearmod (pp);
  if (*k2 < k1)
    {
      unsigned long tmp = *k2;
      *k2 = k1;
      return tmp;
    }
  return k1;
}

#define MAX_PRIMES 2800
static unsigned int Primes[MAX_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, 8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831, 8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973, 10007, 10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099, 10103, 10111, 10133, 10139, 10141, 10151, 10159, 10163, 10169, 10177, 10181, 10193, 10211, 10223, 10243, 10247, 10253, 10259, 10267, 10271, 10273, 10289, 10301, 10303, 10313, 10321, 10331, 10333, 10337, 10343, 10357, 10369, 10391, 10399, 10427, 10429, 10433, 10453, 10457, 10459, 10463, 10477, 10487, 10499, 10501, 10513, 10529, 10531, 10559, 10567, 10589, 10597, 10601, 10607, 10613, 10627, 10631, 10639, 10651, 10657, 10663, 10667, 10687, 10691, 10709, 10711, 10723, 10729, 10733, 10739, 10753, 10771, 10781, 10789, 10799, 10831, 10837, 10847, 10853, 10859, 10861, 10867, 10883, 10889, 10891, 10903, 10909, 10937, 10939, 10949, 10957, 10973, 10979, 10987, 10993, 11003, 11027, 11047, 11057, 11059, 11069, 11071, 11083, 11087, 11093, 11113, 11117, 11119, 11131, 11149, 11159, 11161, 11171, 11173, 11177, 11197, 11213, 11239, 11243, 11251, 11257, 11261, 11273, 11279, 11287, 11299, 11311, 11317, 11321, 11329, 11351, 11353, 11369, 11383, 11393, 11399, 11411, 11423, 11437, 11443, 11447, 11467, 11471, 11483, 11489, 11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657, 11677, 11681, 11689, 11699, 11701, 11717, 11719, 11731, 11743, 11777, 11779, 11783, 11789, 11801, 11807, 11813, 11821, 11827, 11831, 11833, 11839, 11863, 11867, 11887, 11897, 11903, 11909, 11923, 11927, 11933, 11939, 11941, 11953, 11959, 11969, 11971, 11981, 11987, 12007, 12011, 12037, 12041, 12043, 12049, 12071, 12073, 12097, 12101, 12107, 12109, 12113, 12119, 12143, 12149, 12157, 12161, 12163, 12197, 12203, 12211, 12227, 12239, 12241, 12251, 12253, 12263, 12269, 12277, 12281, 12289, 12301, 12323, 12329, 12343, 12347, 12373, 12377, 12379, 12391, 12401, 12409, 12413, 12421, 12433, 12437, 12451, 12457, 12473, 12479, 12487, 12491, 12497, 12503, 12511, 12517, 12527, 12539, 12541, 12547, 12553, 12569, 12577, 12583, 12589, 12601, 12611, 12613, 12619, 12637, 12641, 12647, 12653, 12659, 12671, 12689, 12697, 12703, 12713, 12721, 12739, 12743, 12757, 12763, 12781, 12791, 12799, 12809, 12821, 12823, 12829, 12841, 12853, 12889, 12893, 12899, 12907, 12911, 12917, 12919, 12923, 12941, 12953, 12959, 12967, 12973, 12979, 12983, 13001, 13003, 13007, 13009, 13033, 13037, 13043, 13049, 13063, 13093, 13099, 13103, 13109, 13121, 13127, 13147, 13151, 13159, 13163, 13171, 13177, 13183, 13187, 13217, 13219, 13229, 13241, 13249, 13259, 13267, 13291, 13297, 13309, 13313, 13327, 13331, 13337, 13339, 13367, 13381, 13397, 13399, 13411, 13417, 13421, 13441, 13451, 13457, 13463, 13469, 13477, 13487, 13499, 13513, 13523, 13537, 13553, 13567, 13577, 13591, 13597, 13613, 13619, 13627, 13633, 13649, 13669, 13679, 13681, 13687, 13691, 13693, 13697, 13709, 13711, 13721, 13723, 13729, 13751, 13757, 13759, 13763, 13781, 13789, 13799, 13807, 13829, 13831, 13841, 13859, 13873, 13877, 13879, 13883, 13901, 13903, 13907, 13913, 13921, 13931, 13933, 13963, 13967, 13997, 13999, 14009, 14011, 14029, 14033, 14051, 14057, 14071, 14081, 14083, 14087, 14107, 14143, 14149, 14153, 14159, 14173, 14177, 14197, 14207, 14221, 14243, 14249, 14251, 14281, 14293, 14303, 14321, 14323, 14327, 14341, 14347, 14369, 14387, 14389, 14401, 14407, 14411, 14419, 14423, 14431, 14437, 14447, 14449, 14461, 14479, 14489, 14503, 14519, 14533, 14537, 14543, 14549, 14551, 14557, 14561, 14563, 14591, 14593, 14621, 14627, 14629, 14633, 14639, 14653, 14657, 14669, 14683, 14699, 14713, 14717, 14723, 14731, 14737, 14741, 14747, 14753, 14759, 14767, 14771, 14779, 14783, 14797, 14813, 14821, 14827, 14831, 14843, 14851, 14867, 14869, 14879, 14887, 14891, 14897, 14923, 14929, 14939, 14947, 14951, 14957, 14969, 14983, 15013, 15017, 15031, 15053, 15061, 15073, 15077, 15083, 15091, 15101, 15107, 15121, 15131, 15137, 15139, 15149, 15161, 15173, 15187, 15193, 15199, 15217, 15227, 15233, 15241, 15259, 15263, 15269, 15271, 15277, 15287, 15289, 15299, 15307, 15313, 15319, 15329, 15331, 15349, 15359, 15361, 15373, 15377, 15383, 15391, 15401, 15413, 15427, 15439, 15443, 15451, 15461, 15467, 15473, 15493, 15497, 15511, 15527, 15541, 15551, 15559, 15569, 15581, 15583, 15601, 15607, 15619, 15629, 15641, 15643, 15647, 15649, 15661, 15667, 15671, 15679, 15683, 15727, 15731, 15733, 15737, 15739, 15749, 15761, 15767, 15773, 15787, 15791, 15797, 15803, 15809, 15817, 15823, 15859, 15877, 15881, 15887, 15889, 15901, 15907, 15913, 15919, 15923, 15937, 15959, 15971, 15973, 15991, 16001, 16007, 16033, 16057, 16061, 16063, 16067, 16069, 16073, 16087, 16091, 16097, 16103, 16111, 16127, 16139, 16141, 16183, 16187, 16189, 16193, 16217, 16223, 16229, 16231, 16249, 16253, 16267, 16273, 16301, 16319, 16333, 16339, 16349, 16361, 16363, 16369, 16381, 16411, 16417, 16421, 16427, 16433, 16447, 16451, 16453, 16477, 16481, 16487, 16493, 16519, 16529, 16547, 16553, 16561, 16567, 16573, 16603, 16607, 16619, 16631, 16633, 16649, 16651, 16657, 16661, 16673, 16691, 16693, 16699, 16703, 16729, 16741, 16747, 16759, 16763, 16787, 16811, 16823, 16829, 16831, 16843, 16871, 16879, 16883, 16889, 16901, 16903, 16921, 16927, 16931, 16937, 16943, 16963, 16979, 16981, 16987, 16993, 17011, 17021, 17027, 17029, 17033, 17041, 17047, 17053, 17077, 17093, 17099, 17107, 17117, 17123, 17137, 17159, 17167, 17183, 17189, 17191, 17203, 17207, 17209, 17231, 17239, 17257, 17291, 17293, 17299, 17317, 17321, 17327, 17333, 17341, 17351, 17359, 17377, 17383, 17387, 17389, 17393, 17401, 17417, 17419, 17431, 17443, 17449, 17467, 17471, 17477, 17483, 17489, 17491, 17497, 17509, 17519, 17539, 17551, 17569, 17573, 17579, 17581, 17597, 17599, 17609, 17623, 17627, 17657, 17659, 17669, 17681, 17683, 17707, 17713, 17729, 17737, 17747, 17749, 17761, 17783, 17789, 17791, 17807, 17827, 17837, 17839, 17851, 17863, 17881, 17891, 17903, 17909, 17911, 17921, 17923, 17929, 17939, 17957, 17959, 17971, 17977, 17981, 17987, 17989, 18013, 18041, 18043, 18047, 18049, 18059, 18061, 18077, 18089, 18097, 18119, 18121, 18127, 18131, 18133, 18143, 18149, 18169, 18181, 18191, 18199, 18211, 18217, 18223, 18229, 18233, 18251, 18253, 18257, 18269, 18287, 18289, 18301, 18307, 18311, 18313, 18329, 18341, 18353, 18367, 18371, 18379, 18397, 18401, 18413, 18427, 18433, 18439, 18443, 18451, 18457, 18461, 18481, 18493, 18503, 18517, 18521, 18523, 18539, 18541, 18553, 18583, 18587, 18593, 18617, 18637, 18661, 18671, 18679, 18691, 18701, 18713, 18719, 18731, 18743, 18749, 18757, 18773, 18787, 18793, 18797, 18803, 18839, 18859, 18869, 18899, 18911, 18913, 18917, 18919, 18947, 18959, 18973, 18979, 19001, 19009, 19013, 19031, 19037, 19051, 19069, 19073, 19079, 19081, 19087, 19121, 19139, 19141, 19157, 19163, 19181, 19183, 19207, 19211, 19213, 19219, 19231, 19237, 19249, 19259, 19267, 19273, 19289, 19301, 19309, 19319, 19333, 19373, 19379, 19381, 19387, 19391, 19403, 19417, 19421, 19423, 19427, 19429, 19433, 19441, 19447, 19457, 19463, 19469, 19471, 19477, 19483, 19489, 19501, 19507, 19531, 19541, 19543, 19553, 19559, 19571, 19577, 19583, 19597, 19603, 19609, 19661, 19681, 19687, 19697, 19699, 19709, 19717, 19727, 19739, 19751, 19753, 19759, 19763, 19777, 19793, 19801, 19813, 19819, 19841, 19843, 19853, 19861, 19867, 19889, 19891, 19913, 19919, 19927, 19937, 19949, 19961, 19963, 19973, 19979, 19991, 19993, 19997, 20011, 20021, 20023, 20029, 20047, 20051, 20063, 20071, 20089, 20101, 20107, 20113, 20117, 20123, 20129, 20143, 20147, 20149, 20161, 20173, 20177, 20183, 20201, 20219, 20231, 20233, 20249, 20261, 20269, 20287, 20297, 20323, 20327, 20333, 20341, 20347, 20353, 20357, 20359, 20369, 20389, 20393, 20399, 20407, 20411, 20431, 20441, 20443, 20477, 20479, 20483, 20507, 20509, 20521, 20533, 20543, 20549, 20551, 20563, 20593, 20599, 20611, 20627, 20639, 20641, 20663, 20681, 20693, 20707, 20717, 20719, 20731, 20743, 20747, 20749, 20753, 20759, 20771, 20773, 20789, 20807, 20809, 20849, 20857, 20873, 20879, 20887, 20897, 20899, 20903, 20921, 20929, 20939, 20947, 20959, 20963, 20981, 20983, 21001, 21011, 21013, 21017, 21019, 21023, 21031, 21059, 21061, 21067, 21089, 21101, 21107, 21121, 21139, 21143, 21149, 21157, 21163, 21169, 21179, 21187, 21191, 21193, 21211, 21221, 21227, 21247, 21269, 21277, 21283, 21313, 21317, 21319, 21323, 21341, 21347, 21377, 21379, 21383, 21391, 21397, 21401, 21407, 21419, 21433, 21467, 21481, 21487, 21491, 21493, 21499, 21503, 21517, 21521, 21523, 21529, 21557, 21559, 21563, 21569, 21577, 21587, 21589, 21599, 21601, 21611, 21613, 21617, 21647, 21649, 21661, 21673, 21683, 21701, 21713, 21727, 21737, 21739, 21751, 21757, 21767, 21773, 21787, 21799, 21803, 21817, 21821, 21839, 21841, 21851, 21859, 21863, 21871, 21881, 21893, 21911, 21929, 21937, 21943, 21961, 21977, 21991, 21997, 22003, 22013, 22027, 22031, 22037, 22039, 22051, 22063, 22067, 22073, 22079, 22091, 22093, 22109, 22111, 22123, 22129, 22133, 22147, 22153, 22157, 22159, 22171, 22189, 22193, 22229, 22247, 22259, 22271, 22273, 22277, 22279, 22283, 22291, 22303, 22307, 22343, 22349, 22367, 22369, 22381, 22391, 22397, 22409, 22433, 22441, 22447, 22453, 22469, 22481, 22483, 22501, 22511, 22531, 22541, 22543, 22549, 22567, 22571, 22573, 22613, 22619, 22621, 22637, 22639, 22643, 22651, 22669, 22679, 22691, 22697, 22699, 22709, 22717, 22721, 22727, 22739, 22741, 22751, 22769, 22777, 22783, 22787, 22807, 22811, 22817, 22853, 22859, 22861, 22871, 22877, 22901, 22907, 22921, 22937, 22943, 22961, 22963, 22973, 22993, 23003, 23011, 23017, 23021, 23027, 23029, 23039, 23041, 23053, 23057, 23059, 23063, 23071, 23081, 23087, 23099, 23117, 23131, 23143, 23159, 23167, 23173, 23189, 23197, 23201, 23203, 23209, 23227, 23251, 23269, 23279, 23291, 23293, 23297, 23311, 23321, 23327, 23333, 23339, 23357, 23369, 23371, 23399, 23417, 23431, 23447, 23459, 23473, 23497, 23509, 23531, 23537, 23539, 23549, 23557, 23561, 23563, 23567, 23581, 23593, 23599, 23603, 23609, 23623, 23627, 23629, 23633, 23663, 23669, 23671, 23677, 23687, 23689, 23719, 23741, 23743, 23747, 23753, 23761, 23767, 23773, 23789, 23801, 23813, 23819, 23827, 23831, 23833, 23857, 23869, 23873, 23879, 23887, 23893, 23899, 23909, 23911, 23917, 23929, 23957, 23971, 23977, 23981, 23993, 24001, 24007, 24019, 24023, 24029, 24043, 24049, 24061, 24071, 24077, 24083, 24091, 24097, 24103, 24107, 24109, 24113, 24121, 24133, 24137, 24151, 24169, 24179, 24181, 24197, 24203, 24223, 24229, 24239, 24247, 24251, 24281, 24317, 24329, 24337, 24359, 24371, 24373, 24379, 24391, 24407, 24413, 24419, 24421, 24439, 24443, 24469, 24473, 24481, 24499, 24509, 24517, 24527, 24533, 24547, 24551, 24571, 24593, 24611, 24623, 24631, 24659, 24671, 24677, 24683, 24691, 24697, 24709, 24733, 24749, 24763, 24767, 24781, 24793, 24799, 24809, 24821, 24841, 24847, 24851, 24859, 24877, 24889, 24907, 24917, 24919, 24923, 24943, 24953, 24967, 24971, 24977, 24979, 24989, 25013, 25031, 25033, 25037, 25057, 25073, 25087, 25097, 25111, 25117, 25121, 25127, 25147, 25153, 25163, 25169, 25171, 25183, 25189, 25219, 25229, 25237, 25243, 25247, 25253, 25261, 25301, 25303, 25307, 25309, 25321, 25339, 25343, 25349, 25357, 25367, 25373, 25391};

#define INDEX 25392 /* should be larger than the last prime above */

static int prime_index[INDEX];

static inline void
setbit (mpz_t row, int shift, unsigned long i)
{
  mpz_setbit (row, shift + i);
}

/* res=2: found a new full relation from two partial relations */
void
smooth_stat (int res)
{
  static unsigned long tested = 0, smooth[2] = {0,0}, full2 = 0;

  if (res >= 0)
    {
      tested ++;
      if (res == 1)
        smooth[0] ++;   /* 0 large prime */
      else if (res == 2)
        full2 ++;
      else if (res > 1)
        smooth[1] ++;   /* 1 large prime */
    }
  else
    printf ("Tested %lu norms, %lu+%lu smooth, %lu 1LP-smooth, total %lu (%.0f%%)\n",
            tested, smooth[0], full2, smooth[1], smooth[0] + smooth[1],
            100.0 * (double) (smooth[0] + smooth[1]) / (double) tested);
}

void
hash_init (hash_t H, unsigned long L)
{
  unsigned long i;
  double pi = (double) L / log ((double) L);

  H->alloc = ulong_nextprime ((unsigned long) pi);
  H->size = 0;
  H->p = malloc (H->alloc * sizeof (unsigned int));
  H->row = malloc (H->alloc * sizeof (mpz_t));
  H->x = malloc (H->alloc * sizeof (mpz_t));
  H->y = malloc (H->alloc * sizeof (mpz_t));
  for (i = 0; i < H->alloc; i++)
    {
      H->p[i] = 0;
      mpz_init (H->row[i]);
      mpz_init (H->x[i]);
      mpz_init (H->y[i]);
    }
}

/* Return 0 when large prime p appears for first time,
   return 1 when it already appeared before, and puts then in 'row' the
   bit-vector corresponding to the sum of both relations, and in 'x' the
   product of the two a*i+b values. */
int
hash_insert (hash_t H, unsigned int p, mpz_t row, mpz_t x, mpz_t y)
{
  unsigned int h;

  ASSERT_ALWAYS (2 * H->size <= H->alloc);
  h = p % H->alloc;
  while (H->p[h] != 0 && H->p[h] != p)
    if (++h == H->alloc)
      h = 0;
  if (H->p[h] == 0)
    {
      H->p[h] = p;
      mpz_set (H->row[h], row);
      mpz_set (H->x[h], x);
      mpz_set (H->y[h], y);
      H->size ++;
      return 0;
    }
  else /* found two relations with same large prime! */
    {
      mpz_xor (row, row, H->row[h]);
      mpz_mul (x, x, H->x[h]);
      mpz_mul (y, y, H->y[h]);
      return 1;
    }
}

void
hash_clear (hash_t H)
{
  unsigned long i;

  for (i = 0; i < H->alloc; i++)
    {
      mpz_clear (H->row[i]);
      mpz_clear (H->x[i]);
      mpz_clear (H->y[i]);
    }
  free (H->p);
  free (H->row);
  free (H->x);
  free (H->y);
}

/* Trial divide r (which should be smooth) over the factor base.
   We put relations in column 'shift' and above.
*/
static void
trialdiv (mpz_t r, fb_t *F, long ncol, int shift, mpz_t row)
{
  unsigned int i;
  unsigned long B = F[ncol-1].p, R;

  mpz_set_ui (row, 0);

  if (mpz_sgn (r) < 0)
    {
      setbit (row, shift, 0); /* sign is in column 0 */
      mpz_neg (r, r);
    }

  for (i = 0; i < ncol; i++)
    {
      unsigned long p = F[i].p;

      if (mpz_divisible_ui_p (r, p))
        {
          int e = 0;

          do {
            mpz_divexact_ui (r, r, p);
            e ++;
          } while (mpz_divisible_ui_p (r, p));

          if (e & 1)
            setbit (row, shift, i + 1);

          /* we don't check for cases where r = 1 or is a prime <= B here,
             since they will have very rarely */

          if (mpz_fits_ulong_p (r))
            {
              i = i + 1;
              break;
            }
        }
    }

  /* now r fits into an unsigned long */
  R = mpz_get_ui (r);
  for (; i < ncol; i++)
    {
      unsigned long p = F[i].p;
      unsigned long q = R * F[i].invp;

      /* F[i].invp is 1/p mod 2^64, thus q = R/p mod 2^64:
         the division R/p is exact iff q*p fits in a 64-bit word,
         i.e., when q <= F[i].invp2 = floor(ULONG_MAX/p) [we assume
         unsigned long has 64 bits here] */
      if (q <= F[i].invp2)
        {
          int e = 0;
          do {
            R = q;
            e ++;
            q = R * F[i].invp;
          } while (q <= F[i].invp2);
          if (e & 1)
            setbit (row, shift, i + 1);
          /* now since R has no factor <= p it is either 1 or a prime > p,
             or a product > p^2 of at least two primes */
          if (R <= B && R <= p * p)
            {
              /* R=1 might happen when the largest prime factor of R appears
                 with exponent 2 */
              if (R > 1)
                {
                  i = prime_index[R];
                  ASSERT(R == F[i].p);
                  setbit (row, shift, i + 1);
                }
              return;
            }
        }
    }
  ASSERT(0);
}

static inline void
update (unsigned char *S, unsigned long i, unsigned long p MAYBE_UNUSED,
        unsigned char logp, unsigned int M MAYBE_UNUSED)
{
#ifdef CHECK_OVERFLOW
  unsigned char tmp = S[i] + logp;
  ASSERT_ALWAYS(tmp >= S[i]);
#endif
  S[i] += logp;
#ifdef TRACE
  if (i == M + TRACE)
    printf ("%d: %lu %u\n", TRACE, p, S[i]);
#endif
}

/* put factor in z */
static void
gauss (mpz_t z, mpz_t *Mat, int nrel, int wrel, int ncol, mpz_t *X, mpz_t *Y,
       const mpz_t N0, int verbose)
{
  int i, j, k, k2, j2;
  int shift = wrel; /* shift for the relations */
  mpz_t t;

  /* we put in the right part of the matrix (ncol columns) the matrix we
     want to put in row echelon form, and in the left part (nrel columns)
     the identity matrix: using this trick, the left part will reveal the
     relations at the end (and the right part will be the identity) */
  mpz_init (t);
  j = ncol - 1;
  i = 0;
  for (; j >= 1; )
    {
      /* try to zero all bits in column j, starting from row i */
      for (k = i; k < nrel; k++)
        if (mpz_tstbit (Mat[k], shift + j))
          break;
      if (k == nrel)
        {
          j -= 1;
          continue; /* go to next column */
        }
      if (i < k) /* swap rows i and k */
        mpz_swap (Mat[i], Mat[k]);
      /* we have found a pivot in (i,j) */
      j2 = j - 1;
    restart:
      if (j2 < 0)
        break; /* we still have to reduce column j */
      /* we must find a 2x2 matrix of rank 2, i.e., either [1,0;0,1] or
         [1,0;1,1] or [1,1;0,1] or [1,1;1,0] but not [1,1;1,1] */
      for (k2 = i + 1; k2 < nrel; k2++)
        {
          /* first reduce bit j using pivot (i,j) */
          if (mpz_tstbit (Mat[k2], shift + j))
            mpz_xor (Mat[k2], Mat[k2], Mat[i]);
          if (mpz_tstbit (Mat[k2], shift + j2))
            break;
        }
      if (k2 == nrel)
        {
          j2 --;
          goto restart;
        }
      if (i+1 < k2) /* swap rows i+1 and k2 */
        mpz_swap (Mat[i+1], Mat[k2]);

      mpz_xor (t, Mat[i], Mat[i+1]);
      /* clear bit j2 in row i */
      if (mpz_tstbit (Mat[i], shift + j2))
        mpz_swap (t, Mat[i]);
      /* now we have a pivot in (i,j) and one in (i+1,j2) */
      for (k = i + 2; k < nrel; k++)
        {
          if (mpz_tstbit (Mat[k], shift + j))
            {
              if (mpz_tstbit (Mat[k], shift + j2))
                mpz_xor (Mat[k], Mat[k], t);
              else
                mpz_xor (Mat[k], Mat[k], Mat[i]);
            }
          else if (mpz_tstbit (Mat[k], shift + j2))
            mpz_xor (Mat[k], Mat[k], Mat[i+1]);
          ASSERT(mpz_sizeinbase(Mat[k],2)<=(unsigned int) shift + j);
        }
      i += 2;
      j = j2 - 1;
    }
  /* finish by one row at a time */
  for (; j >= 0; j--)
    {
      for (k = i; k < nrel; k++)
        if (mpz_tstbit (Mat[k], shift + j))
          break;
      if (k == nrel)
        continue;
      if (i < k) /* swap rows i and k */
        mpz_swap (Mat[i], Mat[k]);
      for (k = i + 1; k < nrel; k++)
        {
          if (mpz_tstbit (Mat[k], shift + j))
            mpz_xor (Mat[k], Mat[k], Mat[i]);
          ASSERT(mpz_sizeinbase(Mat[k],2)<=(unsigned int) shift + j);
        }
      i++;
    }
  mpz_clear (t);
  /* all rows from i to nrel-1 are dependencies */
  mpz_t x, y;
  mpz_init (x);
  mpz_init (y);
  if (verbose)
    printf ("Total %d dependencies\n", nrel - i);
  int i0 = i;
  while (i < nrel)
    {
      if (verbose)
        printf ("Trying dependency %d\n", i - i0);
      mpz_set_ui (x, 1);
      mpz_set_ui (y, 1);
      for (j = 0; j < nrel; j++)
        if (mpz_tstbit (Mat[i], j))
          {
            mpz_mul (x, x, X[j]);
            mpz_mul (y, y, Y[j]);
          }
#if 0
      if (mpz_perfect_square_p (y) == 0)
        {
          for (int p = 2; ; p += 1 + (p >= 3))
            {
              mpz_set_ui (x, p);
              int e = mpz_remove (y, y, x);
              if (e & 1)
                {
                  printf ("prime %d appears with odd exponent %d\n",
                          p, e);
                  abort ();
                }
            }
        }
#endif
      ASSERT (mpz_perfect_square_p (y));
      mpz_sqrt (y, y);
      mpz_sub (z, x, y);
      mpz_gcd (z, z, N0);
      if (mpz_cmp_ui (z, 1) > 0 && mpz_cmp (z, N0) < 0)
        {
          if (verbose)
            gmp_printf ("gcd=%Zd\n", z);
          break;
        }
      i++;
    }

  if (verbose && i == nrel)
    printf ("Found no factor\n");

  mpz_clear (x);
  mpz_clear (y);
}

#if 0
/* implements the modified Knuth-Schroeppel function, as in Silverman,
   "The Multiple Polynomial Quadratic Sieve", page 335.
   This is kept for reference only, since it does not seem to give better
   results than best_multiplier().
 */
static int
best_multiplier2 (const mpz_t N, unsigned long ncol, unsigned long K)
{
  unsigned long k, best_k = 1, i, j, p, Np;
  double alpha, best_alpha = DBL_MAX, g;
  mpz_t kN, P;

  mpz_init (kN);
  mpz_init (P);
  for (k = 1; k <= K; k += 2)
    {
      mpz_mul_ui (kN, N, k);
      alpha = 0.5 * log ((double) k);
      for (i = j = 0; i < MAX_PRIMES && j < ncol; i++)
        {
          p = Primes[i];
          if (p == 2)
            {
              Np = mpz_fdiv_ui (kN, 8);
              g = (Np == 1) ? 2 : 0;
            }
          else
            {
              mpz_set_ui (P, p);
              if (k % p == 0)
                g = 1.0 / (double) p;
              else if (mpz_jacobi (kN, P) == 1)
                g = 2.0 / (double) p;
              else
                g = 0.0;
            }
          alpha -= g * log ((double) p);
          j += (g > 0.0);
        }
      if (alpha < best_alpha)
        {
          best_alpha = alpha;
          best_k = k;
        }
    }
  mpz_clear (kN);
  mpz_clear (P);
  return best_k;
}
#endif

/* consider all primes up to B, and all odd multipliers up to K */
static int
best_multiplier (const mpz_t N, unsigned long B, unsigned long K)
{
  unsigned long i, n, *Q, p, q, j, t, Nq, k, best_k = 1;
  double **X, logp, alpha, best_alpha = DBL_MAX;

  for (i = 0; i < MAX_PRIMES && Primes[i] <= B; i++);
  n = i;
  X = (double**) malloc (n * sizeof (double*));
  Q = (unsigned long*) malloc (n * sizeof (unsigned long));

  /* first compute number of roots for all residues mod q = p^k */
  for (i = 0; i < n; i++)
    {
      p = Primes[i];
      for (q = p, k = 1; q * p <= B; q = q * p, k++);
      Q[i] = q;
      X[i] = malloc (q * sizeof(double));
      for (t = 0; t < q; t++)
        X[i][t] = 0;
      unsigned long qmax = q;
      for (q = p, k = 1; q <= qmax; q *= p, k++)
        {
          logp = log ((double) p) / (double) q;
          for (j = 0; j < q; j++)
            {
              t = (j * j) % q;
              /* we have j^2 = t (mod q) */
              for (; t < qmax; t += q)
                X[i][t] += logp;
            }
        }
    }

  /* now check all odd multiplies */
  for (k = 1; k <= K; k += 2)
    {
      alpha = 0.5 * log ((double) k);
      for (i = 0; i < n; i++)
        {
          q = Q[i];
          Nq = mpz_fdiv_ui (N, q);
          Nq = (Nq * k) % q;
          alpha -= X[i][Nq];
        }
      if (alpha < best_alpha)
        {
          best_alpha = alpha;
          best_k = k;
        }
    }

  for (i = 0; i < n; i++)
    free (X[i]);
  free (X);
  free (Q);

  return best_k;
}

/* Put in f a factor of N0 using MPQS.
   Assume N0 is odd. */
void
mpqs_doit (mpz_t f, const mpz_t N0, int verbose)
{
  mpz_t N, *Mat, *X, *Y;
  unsigned char *S, *T, threshold = 0, mask = 128;
  uint64_t *S8, mask8 = 0;
  unsigned long p, k, Nbits, M, i, wrel, L = 0, *P, *Q;
  long j, nrel = 0, lim, ncolw = 0, ncol;
  mpz_t a, b, c, sqrta, r, axb;
  double radix, logradix = 0;
  long st;
  static long init_time = 0, sieve_time = 0, check_time = 0;
  static long gauss_time = 0, total_time = 0;
  unsigned short *W; /* column weight */
  fb_t *F;
  hash_t H;
  bernstein_t Z;

  /* assume N0 is odd */
  ASSERT_ALWAYS (mpz_fdiv_ui (N0, 2) == 1);

  init_time -= milliseconds ();

  Nbits = mpz_sizeinbase (N0, 2);

  /* Set M (half size of sieve array) and size of factor base.
     Note: when M=2^15, setting M=1<<15 statically is faster with gcc. */
  if (Nbits < 56)
    {
      M = 1 << 10;
      ncol = 40;
    }
  else
    {
      M = 1 << (11 + ((Nbits - 56) >> 4));
      ncol = 43 << ((Nbits - 56) >> 4);
    }

  /* in case M fills the L1 cache, reduce it a bit */
  if (M == 32768)
    M -= M / 12;
  /* ensure M is a multiple of 8 */
  M = 8 * (M / 8);

  mpz_init (N);
  k = best_multiplier (N0, Nbits, Nbits);
  mpz_mul_ui (N, N0, k);

  mpz_init (a);
  mpz_init (sqrta);
  mpz_init (b);
  mpz_init (c);
  mpz_init (axb);

  /* compute factor base */
  F = malloc (ncol * sizeof (fb_t));
  P = malloc (ncol * sizeof (unsigned long));
  Q = malloc (ncol * sizeof (unsigned long));
  W = malloc ((ncol + 1) * sizeof (unsigned short));
  /* a prime p can appear in x^2 mod N only if p is a square modulo N */
  mpz_ui_pow_ui (b, 2, 64);
  memset (prime_index, 0, INDEX * sizeof (int));
  for (i = j = 0; j < ncol && i < MAX_PRIMES; i++)
    {
      unsigned long p = Primes[i];
      mpz_set_ui (a, p);
     /* a prime p can be in the factor base only if N is a square mod p */
      if (p == 2 || mpz_jacobi (N, a) != -1)
        {
          F[j].p = P[j] = p;
          prime_index[p] = j;
          for (int k = p - 1; k >= 0 && prime_index[k] == 0; k--)
            prime_index[k] = j;
          mpz_invert (c, a, b); /* c = 1/p mod 2^64 */
          F[j].invp = mpz_get_ui (c);
          F[j].invp2 = ULONG_MAX / p;
          F[j].Mp = M % p;
          j++;
        }
    }
  double lambda = 1.2;
  L = (unsigned long) pow ((double) F[ncol-1].p, lambda);
  if (verbose)
    printf ("Using multiplier k=%lu, M=%lu, ncol=%lu, B=%u, L=%lu\n",
            k, M, ncol, F[ncol-1].p, L);
  ASSERT_ALWAYS(i < MAX_PRIMES);
  ASSERT_ALWAYS(SKIP >= 1);
  lim = ncol;       /* factor base bound */
  wrel = ncol + 1 + WANT_EXCESS; /* +1 is for -1 */
  Mat = (mpz_t*) malloc (wrel * sizeof (mpz_t));
  for (i = 0; i < wrel; i++)
    mpz_init (Mat[i]);
  X = (mpz_t*) malloc (wrel * sizeof (mpz_t));
  Y = (mpz_t*) malloc (wrel * sizeof (mpz_t));
  for (i = 0; i < wrel; i++)
    {
      mpz_init (X[i]);
      mpz_init (Y[i]);
    }
  /* FIXME: the constant 4 here should depend on the number size */
  init_tree (Z, 4, P, lim);

  /* initialize sieve area [-M, M-1] */
  S = malloc (2 * M * sizeof (char));
  T = malloc (2 * M * sizeof (char));
  ASSERT_ALWAYS(((long) S & 7) == 0);

  /* initialize square roots mod p */
  for (j = 0; j < lim; j++)
    {
      unsigned long Np;
      p = F[j].p;
      Np = mpz_fdiv_ui (N, p);
      F[j].r = tonelli_shanks (Np, p);
    }

  /* we want 'a' near sqrt(2*N)/M */
  mpz_mul_ui (a, N, 2);
  mpz_sqrt (a, a);
  mpz_tdiv_q_ui (a, a, M);

  /* we want 'a' to be a square */
  mpz_sqrt (sqrta, a);

  /* we also want sqrt(a) greater than the largest factor base prime */
  if (mpz_cmp_ui (sqrta, F[lim-1].p) <= 0)
    mpz_set_ui (sqrta, F[lim-1].p + 1);

  memset (W, 0, (ncol + 1) * sizeof (unsigned short));

  hash_init (H, L); /* hash table storing relations with large primes */

  init_time += milliseconds ();

  int pols = 0;
  mpz_init (r);
  while (nrel < ncolw + WANT_EXCESS) {
  sieve_time -= milliseconds ();
  do
    mpz_nextprime (sqrta, sqrta);
  while (mpz_jacobi (N, sqrta) != 1);
  /* we want N to be a square mod a: N = b^2 (mod a) */
  unsigned long aui = mpz_get_ui (sqrta);
  mpz_mul (a, sqrta, sqrta);
  ASSERT(mpz_fits_ulong_p (a));
  /* we want b^2-N divisible by a */
  k = tonelli_shanks (mpz_fdiv_ui (N, aui), aui);
  /* we want b = k + sqrta*t: b^2 = k^2 + 2*k*sqrta*t (mod a), thus
     k^2 + 2*k*sqrta*t = N (mod a): t = ((N-k^2)/sqrta)/(2*k) mod sqrta */
  mpz_set_ui (b, k * k);
  mpz_sub (b, N, b);
  ASSERT (mpz_divisible_p (b, sqrta));
  mpz_divexact (b, b, sqrta);
  mpz_set_ui (c, 2 * k);
  mpz_invert (c, c, sqrta);
  mpz_mul (b, b, c);
  mpz_mod (b, b, sqrta);
  mpz_mul (b, b, sqrta);
  mpz_add_ui (b, b, k);
  /* we want b even to ensure k=1 is a root of (a*k+b)^2 = N (mod 2) */
  if (mpz_tstbit (b, 0) == 1)
    mpz_add (b, b, a);
  /* c = (b^2 - N)/a */
  mpz_mul (c, b, b);
  mpz_sub (c, c, N);
  ASSERT (mpz_divisible_p (c, a));
  mpz_divexact (c, c, a);
#ifdef TRACE
  gmp_printf ("a=%Zd b=%Zd\n", a, b);
#endif
  unsigned long aa, bb, sqrtaa;
  ASSERT (mpz_fits_ulong_p (a));
  aa = mpz_get_ui (a);
  ASSERT (mpz_fits_ulong_p (b));
  bb = mpz_get_ui (b);
  sqrtaa = mpz_get_ui (sqrta);
  aa = sqrtaa * sqrtaa;

  /* we now have (a*x+b)^2 - N = a*Q(x) where Q(x) = a*x^2 + 2*b*x + c */

  /* we initialize radix only once, assuming it will not vary very much */
  if (pols++ == 0)
    {
      double Nd = mpz_get_d (N), ad = (double) aa, maxnorm, maxnorm1;
      /* The maximal norm (in absolute value) is a*M^2-n/a for x=-M or M,
         and n/a for x=0. Keep the largest one in case 'a' is not the optimal
         value sqrt(2N)/M. */
      maxnorm1 = Nd / ad;
      maxnorm = ad * (double) M * (double) M - maxnorm1;
      maxnorm = (maxnorm > maxnorm1) ? maxnorm : maxnorm1;
      /* initialize radix: we want radix^255 = maxnorm ~ a*M^2 */
#define MAXNORM 128
      radix = pow (maxnorm, 1.0 / (double) MAXNORM);
      logradix = log (radix);
#ifdef TRACE
      printf ("radix=%f logradix=%f\n", radix, logradix);
#endif

      /* initialize logp values */
#define ROUND round
      for (j = SKIP; j < lim; j++)
        F[j].logp = (char) ROUND (log ((double) F[j].p) / logradix);

      /* take into account the skipped primes */
      double skip_value;
      /* if N = 1 (mod 8), the average power of 2 dividing x^2-N is 2,
         if N = 3 (mod 8), the average power of 2 dividing x^2-N is 1/2,
         if N = 5 (mod 8), the average power of 2 dividing x^2-N is 1,
         if N = 7 (mod 8), the average power of 2 dividing x^2-N is 1/2
         (cf Status Report on Factoring, by Davis, Holdridge and Simmons,
         Eurocrypt'84, page 203) */
      switch (mpz_getlimbn (N, 0) & 7)
        {
        case 1:
          skip_value = 2.0 * log (2.0);
          break;
        case 5:
          skip_value = 1.0 * log (2.0);
          break;
        case 3:
        case 7:
          skip_value = 0.5 * log (2.0);
          break;
        default:
          ASSERT_ALWAYS(0);
        }
#define SKIP_FACTOR 1.5
      skip_value *= SKIP_FACTOR;
      for (j = 1; j < SKIP; j++)
        /* the factor 1.5 increases the probability to find relations for which
           the contribution of skipped primes is less than average, but on the
           other hand it increases the pressure on trial division for relations
           where that contribution is more than average */
        skip_value += SKIP_FACTOR * 2.0 * log ((double) F[j].p)
          / (double) (F[j].p - 1);

      for (i = 0; i <= M; i++)
        {
          /* At location i, the norm is ((a*i+b)^2-N)/a. Since 0 <= b < a,
             we approximate it by ((a*i)^2-N)/a. Since a is not changing
             much, we assume it only depends on i. */
          double x = ad * (double) i;
          /* Note: we don't consider the +b term here, which makes the norm
             approximation not very reliable near roots of (a*x+b)^2 = N */
          x = (x * x - Nd) / ad;
          x = fabs (x);
          x = (log (x) - skip_value) / logradix;
          T[M - i] = (unsigned char) ceil (x < 0 ? 0 : x);
          if (i < M)
            T[M + i] = T[M - i];
        }

      unsigned char Tmax;
      Tmax = (T[0] >  T[M]) ? T[0] : T[M];
      threshold = Tmax - (char) ROUND (lambda * log ((double) F[ncol-1].p) / logradix);
      if (threshold < 128)
        {
          /* we increase Tmax by the difference to 128, so that in
             T[i] = Tmax - T[i] below T[i] is increased by the same
             amount to reach the increased threshold */
          Tmax += 128 - threshold;
          threshold = 128;
        }
      ASSERT_ALWAYS(threshold >= 128);
      mask = 128;
      mask8 = (mask << 8) | mask;
      mask8 = (mask8 << 16) | mask8;
      mask8 = (mask8 << 32) | mask8;
#ifdef TRACE
      printf ("T[%d]=%u threshold=%u\n", TRACE, T[M + TRACE], threshold);
#endif
      /* We want S[i] >= T[i] - xxx, thus S[i] + Tmax >= T[i] + Tmax - xxx
         thus S[i] + (Tmax - T[i]) >= Tmax - xxx. */
      for (i = 0; i < 2*M; i++)
        T[i] = Tmax - T[i];
#ifdef TRACE
      printf ("T[%d]=%u\n", TRACE, T[M +  TRACE]);
#endif
    }

  memcpy (S, T, 2 * M);

  /* sieve */
  modredcul_batch_Q_to_Fp (Q + SKIP, 1, sqrtaa, 0, P + SKIP, lim - SKIP);
  /* skip the small primes whose average contribution is already taken into
     account */
  for (j = SKIP; j < lim; j++)
    {
      unsigned long k2 = -1, i2;

      p = F[j].p;
      k = findroot (&k2, bb % p, p, F[j].r, Q[j] * Q[j], F[j].Mp);
      /* Note: if x^2 = k*N (mod p) has only one root, which can happen only
         when k*N is divisible by p (and then the root is 0) we will count
         twice this root, but this will be very rare, and by not considering
         this special case we speed up the general case */
      for (i = k, i2 = k2; i2 < 2*M; i += p, i2 += p)
        {
          update (S, i, p, F[j].logp, M);
          update (S, i2, p, F[j].logp, M);
        }
      if (i < 2*M)
        update (S, i, p, F[j].logp, M);
    }

  st = milliseconds ();
  sieve_time += st;
  check_time -= st;

#ifdef TRACE
  printf ("%d: S=%d\n", TRACE, S[M + TRACE]);
#endif

  /* find smooth locations */
  /* ncolw <= ncol+1, and we allocated a matrix of ncol+1+WANT_EXCESS rows */
  for (S8 = (uint64_t*) S, i = 0; i < 2*M; S8++)
    if ((S8[0] & mask8) == 0)
      i += 8;
    else
      for (int ii = 0; ii < 8; ii++, i++)
        {
          if (S[i] >= threshold)
            {
              long ii = (long) i - (long) M;
              unsigned int q;

              mpz_mul_si (axb, a, ii);
              mpz_add_ui (axb, axb, bb);
              mpz_add_ui (r, axb, bb);
              mpz_mul_si (r, r, ii);
              mpz_add (r, r, c);
              accumulate (Z, r, axb);
              if (Z->size == Z->m)
                {
                  get_smooth (Z);
                  for (unsigned int i = 0; i < Z->size; i++)
                    {
                    if (mpz_cmp_ui (Z->x[0][i], L) > 0)
                      smooth_stat (0); /* non smooth */
                    else
                      {
                        mpz_set (X[nrel], Z->axb[i]);
                        mpz_mul (Y[nrel], Z->axb[i], Z->axb[i]);
                        mpz_sub (Y[nrel], Y[nrel], N);
                        q = mpz_get_ui (Z->x[0][i]);
                        /* FIXME: only perform trial division for smooth
                           relations */
                        trialdiv (Z->y[i], F, ncol, wrel, Mat[nrel]);

                        smooth_stat (q);
                        if (q > 1)
                          {
                            q = hash_insert (H, q, Mat[nrel], X[nrel], Y[nrel]);
                            if (q == 1)
                              smooth_stat (2);
                          }

                        if (q == 1) /* new relation */
                          {
                            /* matrix M has nrel rows, the left part has wrel
                               columns and is initially the identity matrix,
                               the right part has ncol+1 columns and contains
                               initially the relations */
                            mpz_setbit (Mat[nrel], nrel);
                            for (int i = 0; i <= ncol; i++)
                              if (mpz_tstbit (Mat[nrel], wrel + i))
                                {
                                  ncolw += W[i] == 0;
                                  W[i] ++;
                                }
                            if (++nrel >= ncolw + WANT_EXCESS)
                              goto end_check;
                          }
                      }
                    }
                  Z->size = 0;
                }
            }
        }
  end_check:
  check_time += milliseconds ();
  }
  mpz_clear (r);
  hash_clear (H);
  if (verbose)
    printf ("%ld rels with %d polynomials: %f per poly\n",
            nrel, pols, (double) nrel / (double) pols);

  gauss_time -= milliseconds ();
  gauss (f, Mat, nrel, wrel, ncol + 1, X, Y, N0, verbose);
  st = milliseconds ();
  gauss_time += st;
  total_time = st;

  free (S);
  free (T);
  for (i = 0; i < wrel; i++)
    mpz_clear (Mat[i]);
  free (Mat);
  for (i = 0; i < wrel; i++)
    {
      mpz_clear (X[i]);
      mpz_clear (Y[i]);
    }
  free (X);
  free (Y);
  free (F);
  free (P);
  free (Q);
  free (W);
  mpz_clear (axb);
  mpz_clear (c);
  mpz_clear (b);
  mpz_clear (sqrta);
  mpz_clear (a);
  mpz_clear (N);
  clear_tree (Z);

  if (verbose)
    printf ("Total time: %ldms (init %ld, sieve %ld, check %ld, gauss %ld)\n",
            total_time, init_time, sieve_time, check_time, gauss_time);
}

/******************* Bernstein smooth part algorithm *************************/

/* put in z the product P[0] * P[1] * ... * P[n-1] */
void
compute_P (mpz_t z, unsigned long *P, unsigned int n)
{
  if (n == 1)
    mpz_set_ui (z, P[0]);
  else
    {
      unsigned int k = (n + 1) / 2;
      mpz_t t;

      compute_P (z, P, k);
      mpz_init (t);
      compute_P (t, P + k, n - k);
      mpz_mul (z, z, t);
      mpz_clear (t);
    }
}

/* m = 2^logm is the number of x[] we accumulate */
void
init_tree (bernstein_t T, unsigned int logm, unsigned long *P, unsigned int n)
{
  unsigned int m = 1 << logm, i, j;

  mpz_init (T->z);
  compute_P (T->z, P, n);
  T->logm = logm;
  T->m = m;
  T->size = 0;
  T->x = malloc ((logm + 1) * sizeof (mpz_t*));
  for (i = 0; i <= logm; i++)
    {
      T->x[i] = malloc ((m >> i) * sizeof (mpz_t));
      for (j = 0; j < (m >> i); j++)
        mpz_init (T->x[i][j]);
    }
  T->y = malloc (m * sizeof (mpz_t));
  for (i = 0; i < m; i++)
    mpz_init (T->y[i]);
  T->axb = malloc (m * sizeof (mpz_t));
  for (i = 0; i < m; i++)
    mpz_init (T->axb[i]);
}

void
clear_tree (bernstein_t T)
{
  unsigned int logm = T->logm, m = T->m, i, j;

  mpz_clear (T->z);
  for (i = 0; i <= logm; i++)
    {
      for (j = 0; j < (m >> i); j++)
        mpz_clear (T->x[i][j]);
      free (T->x[i]);
    }
  free (T->x);
  T->size = 0;
  for (i = 0; i < m; i++)
    mpz_clear (T->y[i]);
  free (T->y);
  for (i = 0; i < m; i++)
    mpz_clear (T->axb[i]);
  free (T->axb);
}

/* put in y[0], ..., y[n-1] the smooth parts of x[0], ..., x[n-1] */
void
get_smooth (bernstein_t T)
{
  unsigned int logm = T->logm, m = T->m, i, j, e;

  /* product tree */
  for (i = 1; i <= logm; i++)
    for (j = 0; j < (m >> i); j++)
      mpz_mul (T->x[i][j], T->x[i-1][2*j], T->x[i-1][2*j+1]);

  /* remainder tree */
  mpz_mod (T->x[logm][0], T->z, T->x[logm][0]);
  for (i = logm; i-- > 0;)
    for (j = 0; j < (m >> i); j++)
      mpz_mod ((i > 0) ? T->x[i][j] : T->y[j], T->x[i+1][j/2], T->x[i][j]);

  e = 8; /* we want to get all prime powers up to p^e in the smooth part */

  /* compute modular powers */
  for (i = 0; i < T->size; i++)
    mpz_powm_ui (T->y[i], T->y[i], e, T->x[0][i]);

  /* compute gcd's */
  for (i = 0; i < T->size; i++)
    {
      mpz_gcd (T->y[i], T->y[i], T->x[0][i]);
      if (mpz_sgn (T->x[0][i]) < 0)
        mpz_neg (T->y[i], T->y[i]);
      mpz_divexact (T->x[0][i], T->x[0][i], T->y[i]);
      /* y[i] contains the smooth part, x[0][i] the non-smooth part */
    }
}

/* add a new value of x */
void
accumulate (bernstein_t T, mpz_t x, mpz_t axb)
{
  ASSERT_ALWAYS(T->size < T->m);
  mpz_set (T->x[0][T->size], x);
  mpz_set (T->axb[T->size], axb);
  T->size ++;
}

