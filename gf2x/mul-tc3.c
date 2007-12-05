/********************************************************************
 * Below this line, experimental code
 * (C) 2007 Marco Bodrato <optimaltoom@bodrato.it>
 * This code is released under GPL 2.0 licence.
 * Modified by Paul Zimmermann, April 2007.
 *
 * Reference: http://bodrato.it/papers/#WAIFI2007
 *
 * "Towards Optimal Toom-Cook Multiplication for Univariate and
 * Multivariate Polynomials in Characteristic 2 and 0." by Marco
 * BODRATO; in C.Carlet and B.Sunar, editors, "WAIFI'07 proceedings",
 * LNCS 4547, pp. 119-136. Springer, Madrid, Spain, June 21-22, 2007.
 */

#include <assert.h>
// #define DEBUG

#ifdef DEBUG
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif

#ifndef MUL_TOOM_THRESHOLD
#define MUL_TOOM_THRESHOLD 17
#endif

#if (MUL_TOOM_THRESHOLD < 17)
#error "MUL_TOOM_THRESHOLD should be at least 17"
#endif

/* c <- a + b */
static
void Add (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b, long n)
{
  long i;
  for (i = 0; i < n; i++)
    c[i] = a[i] ^ b[i];
}

/* c <- c + a + b */
static
void Add3 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b, long n)
{
  long i;
  for (i = 0; i < n; i++)
    c[i] ^= a[i] ^ b[i];
}

/* c <- a + x * b, return carry out */
static
_ntl_ulong AddLsh1 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
                    long n)
{
  _ntl_ulong cy = 0UL, t;
  long i;
  for (i = 0; i < n; i++)
    {
      t = a[i] ^ ((b[i] << 1) | cy);
      cy = b[i] >> (NTL_BITS_PER_LONG - 1);
      c[i] = t;
    }
  return cy;
}

/* c <- x * c, return carry out */
static
_ntl_ulong Lsh1 (_ntl_ulong *c, long n, _ntl_ulong cy)
{
  _ntl_ulong t;
  long i;
  for (i = 0; i < n; i++)
    {
      t = (c[i] << 1) | cy;
      cy = c[i] >> (NTL_BITS_PER_LONG - 1);
      c[i] = t;
    }
  return cy;
}

/* c <- a + cy, return carry out (0 for n > 0, cy for n=0) */
static
_ntl_ulong Add1 (_ntl_ulong *c, const _ntl_ulong *a, long n, _ntl_ulong cy)
{
  if (n)
    {
      long i;
      c[0] = a[0] ^ cy;
      for (i = 1; i < n; i++)
        c[i] = a[i];
      return 0;
    }
  else
    return cy;
}

/* c <- ( c + b )/x, return carry */
static
_ntl_ulong Rsh1Add (_ntl_ulong *c, const _ntl_ulong *b, long n)
{
  _ntl_ulong cy = 0, t;

  for (long i = n - 1; i >= 0; i--)
    {
      t = c[i] ^ b[i];
      cy <<= NTL_BITS_PER_LONG - 1;
      c[i] = (t >> 1) | cy;
      cy = t;
    }
  return cy;
}

/* c <- c + (1+x^3) * b, return carry out */
static
_ntl_ulong AddLsh13 (_ntl_ulong *c, const _ntl_ulong *b, long n)
{
  _ntl_ulong cy = 0UL, t;

  for (long i = 0; i < n; i++)
    {
      t = b[i];
      c[i] ^= t ^ (t << 3) ^ cy;
      cy = t >> (NTL_BITS_PER_LONG - 3);
    }
  return cy;
}

/* c <- c + a + b + d */
static
void Add4 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b, 
           const _ntl_ulong *d, long n)
{
  long i;
  for (i = 0; i < n; i++)
    c[i] ^= a[i] ^ b[i] ^ d[i];
}

/* let c = q*(1+x) + X^n*r with X = x^NTL_BITS_PER_LONG and deg(r) < 1
   then c <- q, returns r.
   (Algorithm from Michel Quercia.)
*/
static
_ntl_ulong DivOnePlusX (_ntl_ulong *c, long n)
{
  _ntl_ulong t = 0;
  long i;

  for (i = 0; i < n; i++) {
    t ^= c[i];
    t ^= t << 1;
    t ^= t << 2;
    t ^= t << 4;
    t ^= t << 8;
    t ^= t << 16;
#if (NTL_BITS_PER_LONG == 64)
    t ^= t << 32;
#elif (NTL_BITS_PER_LONG != 32)
#error "NTL_BITS_PER_LONG should be 32 or 64"
#endif
    c[i] = t;
    t >>= (NTL_BITS_PER_LONG - 1);
  }
  return t;
}

#if (defined(DEBUG))
static void
dump (const _ntl_ulong *a, long n)
{
  for (long i = 0; i < n; i++)
    {
      printf ("+%lu*X^%lu", a[i], i);
      if ((i + 1) % 3 == 0)
	printf ("\n");
    }
  printf (":\n");
}
#endif

/*
\\ gp-pari check code.
default(echo, 1);

A = (a2*x^2 + a1*x + a0)*Mod(1,2)
B = (b2*x^2 + b1*x + b0)*Mod(1,2)
C = A * B
c0 = polcoeff(C, 0)
c1 = polcoeff(C, 1)
c2 = polcoeff(C, 2)
c3 = polcoeff(C, 3)
c4 = polcoeff(C, 4)

\\ --- Evaluation phase. 10 add, 4 shift, 5 mul.

W0 = (a2*y^2+a1*y)*Mod(1,2)
W4 = (b2*y^2+b1*y)*Mod(1,2)
W3 = (a2+a1+a0)   *Mod(1,2)
W2 = (b2+b1+b0)   *Mod(1,2)
W1 = W2 * W3                \\ C(1)
W3 = W3 + W0
W2 = W2 + W4
W0 = W0+a0        *Mod(1,2)
W4 = W4+b0        *Mod(1,2)
W3 = W2 * W3                \\ C(y+1)
W2 = W0 * W4                \\ C(y)
W4 = a2 * b2      *Mod(1,2) \\ C(\infty)
W0 = a0 * b0      *Mod(1,2) \\ C(0)

\\ ------ Interpolation phase. 10 add, 2 shift, 2 div. 
W3 = W3 + W2
W3 == ( c4 + (y^2+y+1)*c3 + c2 + c1 )  \\ check

W2 = ( ( W2 + W0 )/y  + W3 + W4*(y^3+1) ) / (y+1)
\\W2 = ( W2 + W0 + W3*y + W4*(y^4+y) ) / (y^2+y)
W2 == ( c2 + c3 )  \\ check

W1 = W1 + W0
W1 == ( c4 + c3 + c2 + c1 )  \\ check

W3 = ( W3 + W1 ) / (y*(y+1))
W3 == ( c3 )  \\ check

W1 = W1 + W4 + W2
W1 == ( c1 )  \\ check

W2 = W2 + W3
W2 == ( c2 )  \\ check

C == W4*x^4+ W3*x^3+ W2*x^2+ W1*x + W0 \\ check

quit;

 */

/*
  c must have space for 2n words.
  stk must have space for max( KarMem(n), 5k+2 + ToomMem(k+1))
  where k = ceil(n/3).
 */
 
static
void Toom3Mul (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
              long n, _ntl_ulong *stk)
{
  long k = (n + 2) / 3; /* ceil(n/3) */
  long r = n - 2 * k;
  _ntl_ulong cy;
  _ntl_ulong *W0 = c;
  _ntl_ulong *W1 = stk;
  _ntl_ulong *W2 = c + 2 * k;
  _ntl_ulong *W3 = stk + 2 * k;
  _ntl_ulong *W4 = c + 4 * k;

/* \\ --- Evaluation phase. 10 add, 4 shift, 5 mul. */

/* W0 = (a2*y^2+a1*y) */
/* W4 = (b2*y^2+b1*y) */
/*    W0 = (a2*y+a1)*y */
/*    W4 = (b2*y+b1)*y */
  cy = AddLsh1 (W0, a + k, a + 2 * k, r);   /* a1 + x a2 */
  cy = Add1 (W0 + r, a + k + r, k - r, cy);
  W0[k] = (cy << 1) ^ Lsh1 (W0, k, 0); /* x a1 + x^2 a2 */
  cy = AddLsh1 (W4 + 2, b + k, b + 2 * k, r);
  cy = Add1 (W4 + 2 + r , b + k + r, k - r, cy);
  W4[2 + k] = (cy << 1) ^ Lsh1 (W4 + 2, k, 0); /* x b1 + x^2 b2 */
  
  /* using W4[2+k] requires that k+3 words are available at W4=c+4k.
     Since c contains 2n=4k+2r words, then W4 contains 2r words, thus
     we need k+3 <= 2r. This is true for n >= 17.
     Also true for n = 9, 12, 14, 15 but timing tests show that 
     Toom3Mul not the fastest routine for such small n. */
     
  ASSERT (k + 3 <= 2 * r);

  /* {c, k+1}: x*a1+x^2*a2, {c+4k, k+1}: x*b1+x^2*b2 */

/* W3 = ((a2+a1)+a0)    */
/* W2 = ((b2+b1)+b0)    */
  Add (c+k+1,     a, a + k,     k);
  Add (c+k+1, c+k+1, a + 2 * k, r); /* a0 + a1 + a2 */
  Add (W2+2,      b, b + k,     k);
  Add (W2+2,   W2+2, b + 2 * k, r); /* b0 + b1 + b2 */
/* W1 = W2 * W3                \\ C(1) */

  /* {c, k+1}: x*a1+x^2*a2, {c+k+1, k}: a0+a1+a2, {c+2k+2,k}: b0+b1+b2,
     {c+4k, k+1}: x*b1+x^2*b2 */

  Toom (W1, c+k+1, W2+2, k, W3); /* W1 : 2*k */

  /* {c, k+1}: x*a1+x^2*a2, {c+k+1, k}: a0+a1+a2, {c+2k+2,k}: b0+b1+b2,
     {c+4k, k+1}: x*b1+x^2*b2, {stk, 2k}: C(1) */

/* W3 = W3 + W0 */
/* W2 = W2 + W4 */
  Add (c+k+1, c+k+1, W0, k);  c[2*k+1]=W0[k];  /* a0 + (x+1)a1 + (x^2+1)a2 */
  Add (W2+2,  W2+2, W4+2, k); W2[k+2]=W4[k+2]; /* b0 + (x+1)b1 + (x^2+1)b2 */
  /* since we use W2[k+2], we need k+3 words in W2, i.e., 3 <= k */
  //  assert (k >= 3);

  /* {c, k+1}: x*a1+x^2*a2, {c+k+1, k+1}: a0+(1+x)*a1+(1+x^2)*a2,
     {c+2k+2,k+1}: b0+(1+x)*b1+(1+x^2)*b2,
     {c+4k, k+1}: x*b1+x^2*b2, {stk, 2k}: C(1) */

/* W0 = W0+a0    */
/* W4 = W4+b0    */
  Add (W0, W0, a, k); /* a0 + (x)a1 + (x^2)a2 */
  Add (W4+2, W4+2, b, k); /* b0 + (x)b1 + (x^2)b2 */

  /* {c, k+1}: a0+x*a1+x^2*a2, {c+k+1, k+1}: a0+(1+x)*a1+(1+x^2)*a2,
     {c+2k+2,k+1}: b0+(1+x)*b1+(1+x^2)*b2,
     {c+4k, k+1}: b0+x*b1+x^2*b2, {stk, 2k}: C(1) */

/* W3 = W2 * W3                \\ C(y+1) */
/* W2 = W0 * W4                \\ C(y) */
  Toom (W3, W2+2, c+k+1, k + 1, stk + 5*k + 2); /* W3 : 2*k+1 */

  /* {c, k+1}: a0+x*a1+x^2*a2, {c+k+1, k+1}: a0+(1+x)*a1+(1+x^2)*a2,
     {c+2k+2,k+1}: b0+(1+x)*b1+(1+x^2)*b2,
     {c+4k, k+1}: b0+x*b1+x^2*b2, {stk, 2k}: C(1), {stk+2k,2k+2}: C(1+x) */

  Toom (W2,   W0,  W4+2, k + 1, stk + 4*k + 3);
  cy = W4[0]; /* contains at most 3 bits */
  ASSERT (cy <= 7);

  /* {c, k+1}: a0+x*a1+x^2*a2, {c+2k, 2k+2}: C(x),
     {c+4k, k+1}: b0+x*b1+x^2*b2, {stk, 2k}: C(1), {stk+2k,2k+2}: C(1+x) */

/* W4 = a2 * b2                \\ C(\infty) */
/* W0 = a0 * b0                \\ C(0) */
  Toom (W0, a, b, k, stk + 4*k + 3); /* W0 : 2*k */

  /* {c, 2k}: C(0), {c+2k, 2k+2}: C(x),
     {c+4k, k+1}: b0+x*b1+x^2*b2, {stk, 2k}: C(1), {stk+2k,2k+2}: C(1+x) */

  Toom (W4, a+2*k, b+2*k,    r, stk + 4*k +3); /* W4 : 2*r */

  /* {c, 2k}: C(0), {c+2k, 2k}+cy: C(x),
     {c+4k, 2r}: C(Inf), {stk, 2k}: C(1), {stk+2k,2k+2}: C(1+x) */

/* \\ ------ Interpolation phase. 10 add, 2 shift, 2 div.  */

/* W3 = W3 + W2: W3 has at most 2k words + 3 bits, W2 has 2k words +
   at most 3 bits (stored in cy). */
  Add (W3, W3, W2, 2 * k);
  W3[2 * k] ^= cy;
  /* now W3 has at most 2k words + 1 bit, which can be non-zero only
     if r = k one most significant bit from a2 and b2 is set. */
  ASSERT (W3[2 * k] <= 1);

  /* {c, 2k}: C(0), {c+2k, 2k}+cy: C(x),
     {c+4k, 2r}: C(Inf), {stk, 2k}: C(1), {stk+2k,2k+2}: C(1+x)+C(x) */

/* W2 = ( ( W2 + W0 )/y  + W3 + W4*(y^3+1) ) / (y+1) */
/* \\W2 = ( W2 + W0 + W3*y + W4*(y^4+y) ) / (y^2+y) */
  /* W2 has 2k words + at most 3 bits (stored in cy), W0 has 2k words */
  Rsh1Add (W2, W0, 2 * k);
  W2[2*k-1] |= cy << (NTL_BITS_PER_LONG - 1);
  /* now W2 has at most 2k words + 3 bits (cy >> 1), but since the final
     results will have 2k words only, we can ignore cy. */
  Add (W2, W2, W3, 2 * k);
  cy = AddLsh13 (W2, W4, 2 * r);
  if ( r != k )
    W2[2 * r] ^= cy;
  /* else ignore the carry, since W2 should have 2k words, taking into
     account the above ignored cy >> 1. */
  DivOnePlusX (W2, 2 * k);

/* W1 = W1 + W0 */
/* W1 == ( c4 + c3 + c2 + c1 )  \\ check */
  Add (W1, W1, W0, 2 * k);

/* W3 = ( W3 + W1 ) / (y*(y+1)) */
/* W3 == ( c3 )  \\ check */
  Rsh1Add (W3, W1, 2*k);
  W3[2 * k - 1] |= W3[2 * k] << (NTL_BITS_PER_LONG - 1);
  DivOnePlusX (W3, 2 * k);

/* W1 = W1 + W4 */
  Add (W1, W1, W4, 2 * r);

  /* perform simultaneously W1 <- W1 + W2, W2 <- W2 + W3,
     and {c + k, 4k} <- {c + k, 4k} + {W1, 4k} */
  Add3 (c + k, W1, W2, k);
  Add4 (W2, W1 + k, W2 + k, W3, k);
  Add3 (W2 + k, W3, W3 + k, k);
  Add (W4, W4, W3 + k, k);

/* C == W4*x^4+ W3*x^3+ W2*x^2+ W1*x + W0 \\ check */

  /* assume 5*k <= 2*n = 4*k + 2*r, i.e., k <= 2*r, which is true for n >= 8 */
  ASSERT(k <= 2 * r);
}

#undef DEBUG
#undef ASSERT
