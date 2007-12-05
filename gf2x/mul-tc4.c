/********************************************************************
 * Below this line, experimental code
 * (C) 2007 Marco Bodrato <optimaltoom@bodrato.it>
 * Modified by Paul Zimmermann, April 2007.
 * This code is released under GPL 2.0 licence.
 *
 * Reference: http://bodrato.it/papers/#WAIFI2007
 *
 * "Towards Optimal Toom-Cook Multiplication for Univariate and
 * Multivariate Polynomials in Characteristic 2 and 0." by Marco
 * BODRATO; in C.Carlet and B.Sunar, editors, "WAIFI'07 proceedings",
 * LNCS 4547, pp. 119-136. Springer, Madrid, Spain, June 21-22, 2007.
 */

#ifndef MUL_TOOM4_THRESHOLD
#define MUL_TOOM4_THRESHOLD 30
#endif

#if (MUL_TOOM4_THRESHOLD < 30)
#error "MUL_TOOM4_THRESHOLD should be at least 30"
#endif

/* let c = q*(1+x^2) + X^n*r with X = x^NTL_BITS_PER_LONG and deg(r) < 2
   then c <- q, returns r.
   (Algorithm from Michel Quercia.)
*/
static _ntl_ulong DivOnePlusX2 (_ntl_ulong *c, long n)
{
  _ntl_ulong t = 0;
  long i;

  for (i = 0; i < n; i++) {
    t ^= c[i];
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
    t >>= (NTL_BITS_PER_LONG - 2);
  }
  return t;
}

/* c <- x * a + x^2 * b, return carry out */
static
_ntl_ulong AddLsh1Lsh2 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
                    long n)
{
  _ntl_ulong cy = 0UL, t;
  long i;
  for (i = 0; i < n; i++)
    {
      t = (a[i]<<1) ^ ((b[i] << 2) | cy);
      cy =  (a[i] >> (NTL_BITS_PER_LONG - 1)) ^ (b[i] >> (NTL_BITS_PER_LONG - 2));
      c[i] = t;
    }
  return cy;
}

/* c <- a + x^2 * b, return carry out */
static
_ntl_ulong AddLsh2 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
                    long n)
{
  _ntl_ulong cy = 0UL, t;
  long i;
  for (i = 0; i < n; i++)
    {
      t = a[i] ^ ((b[i] << 2) | cy);
      cy = b[i] >> (NTL_BITS_PER_LONG - 2);
      c[i] = t;
    }
  return cy;
}

/* c <- a + x^6 * b, return carry out */
static
_ntl_ulong AddLsh6 (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
                    long n)
{
  _ntl_ulong cy = 0UL, t;

  for (long i = 0; i < n; i++)
    {
      t = a[i] ^ ((b[i] << 6) | cy);
      cy = b[i] >> (NTL_BITS_PER_LONG - 6);
      c[i] = t;
    }
  return cy;
}

/* let c = q*(x+x^4) + X^n*r with X = x^NTL_BITS_PER_LONG and deg(r) < 1
   then c <- q, returns r.
*/
static _ntl_ulong DivExactD1 (_ntl_ulong *c, long n)
{
  _ntl_ulong t = 0;
  long i;

  for (i = 0; i < n; i++) {
    t ^=( c[i] >> 1 ) | ( (i+1<n) ? ( c[i+1] << (NTL_BITS_PER_LONG-1) ) : 0);
    t ^= t << 3 ^ t << 6;
    t ^= t << 9 ^ t << 18;
    t ^= t << 27
#if (NTL_BITS_PER_LONG == 64)
                ^ t << 54
#elif (NTL_BITS_PER_LONG != 32)
#error "NTL_BITS_PER_LONG should be 32 or 64"
#endif
                         ;
    c[i] = t;
    t >>= (NTL_BITS_PER_LONG - 3);
  }
  return t;
}

/* let c = q*(x^2+x^4) + X^n*r with X = x^NTL_BITS_PER_LONG and deg(r) < 1
   then c <- q, returns r.
*/
static _ntl_ulong DivExactD2 (_ntl_ulong *c, long n)
{
  /* c <- c/x^2 */
  _ntl_ulong cy = 0, t;

  for (long i = n - 1; i >= 0; i--)
    {
      t = (c[i] >> 2) | (cy << (NTL_BITS_PER_LONG - 2));
      cy = c[i]; /* no need to mask the low 2 bits, since they will disappear
                    with the next cy << (NTL_BITS_PER_LONG - 2) */
      c[i] = t;
    }
  return DivOnePlusX2(c,n);
}

#if 0
/* Same as DivExactD2, but with one pass only. However, does not seem to give
   a significant speedup, thus disabled for now. */
static _ntl_ulong
DivExactD2a (_ntl_ulong *c, long n)
{
  _ntl_ulong t, ci;
  
      t = c[0];
      t ^= t << 2;
      t ^= t << 4;
      t ^= t << 8;
      t ^= t << 16;
#if (NTL_BITS_PER_LONG == 64)
      t ^= t << 32;
#elif (NTL_BITS_PER_LONG != 32)
#error "NTL_BITS_PER_LONG should be 32 or 64"
#endif
      ci = t;
      t >>= (NTL_BITS_PER_LONG - 2);
  for (long i = 1; i < n; i++)
    {
      t ^= c[i];
      t ^= t << 2;
      t ^= t << 4;
      t ^= t << 8;
      t ^= t << 16;
#if (NTL_BITS_PER_LONG == 64)
      t ^= t << 32;
#endif
      /* now t is the result of the division of c[i] by (1+x^2), and 
         t >> (NTL_BITS_PER_LONG - 2) the corresponding carry */
      c[i - 1] = (ci >> 2) | (t << (NTL_BITS_PER_LONG - 2));
      ci = t;
      t >>= (NTL_BITS_PER_LONG - 2);
    }
  c[n - 1] = ci >> 2;
  return t;
}
#endif

/*
\\ gp-pari check code.
\\ (C) 2007 Marco Bodrato <optimaltoom@bodrato.it>
\\ This code is released under GPL 2.0 licence.

U0=u0*Mod(1,2);U1=u1*Mod(1,2);U2=u2*Mod(1,2);U3=u3*Mod(1,2);
V0=v0*Mod(1,2);V1=v1*Mod(1,2);V2=v2*Mod(1,2);V3=v3*Mod(1,2);

U = U3*Y^3 + U2*Y^2*X + U1*Y*X^2 + U0*X^3
V = V3*Y^3 + V2*Y^2*X + V1*Y*X^2 + V0*X^3

\\ P(X,Y): P0=(1,0); P1=(x+1,1); P2=(x,1); P3=(1,1); P4=(1,x); P5=(1,x+1); P6=(0,1)
\\Evaluation phase: 13*2 add, 7*2 shift, 2Smul; 7 mul (n)

W1 = U0 + U1 + U2 + U3    ; W2 = V0 + V1 + V2 + V3
W0 = U1 +(U2 + U3*x)*x    ; W6 = V1 +(V2 + V3*x)*x
W4 = W1 +(W0 + U3*(x+1))*x; W3 = W2 +(W6 + V3*(x+1))*x
W0 = W0*x + U0            ; W6 = W6*x + V0

W5 = W4 * W3              ; W4 = W0 * W6
W3 = W1 * W2

W0 =(U2 +(U1 + U0*x)*x)*x ; W6 =(V2 +(V1 + V0*x)*x)*x
W1 = W1 + W0 + U0*(x^2+x) ; W2 = W2 + W6 + V0*(x^2+x)
W0 = W0 + U3              ; W6 = W6 + V3

W1 = W1 * W2              ; W2 = W0 * W6
W6 = U3 * V3              ; W0 = U0 * V0

\\Interpolation: 22 add, 4 shift, 5 Smul, 4 div (2n)
d1=(x^4+x)*Mod(1,2)	; d1== (x)^1*(x+1)^1*(x^2+x+1)^1 *Mod(1,2)
d2=(x^4+x^2)*Mod(1,2)	; d2== (x)^2*(x+1)^2*(x^2+x+1)^0 *Mod(1,2)

W1 = W1 + W2 + W0*(x^4+x^2+1)
W5 =(W5 + W4 + W1 + W6*(x^4+x^2+1))/d1
W2 = W2 + W6 + W0*(x^6)
W4 = W4 + W2 + W0 + W6*(x^6)
W4 =(W4 + W5*(x^5+x))/d2
W3 = W3 + W0 + W6
W1 = W1 + W3
W2 = W2 +(W1 + W3*x)*x
W3 = W3 + W4 + W5
W1 =(W1 + W3*(x^2+x))/d1
W5 = W5 + W1
W2 =(W2 + W5*(x^2+x))/d2
W4 = W4 + W2

\\Recomposition
W  = W6*Y^6 + W5*Y^5*X + W4*Y^4*X^2+ W3*Y^3*X^3+ W2*Y^2*X^4+ W1*Y*X^5 + W0*X^6
W == U*V

Memory Usage: stk must have space for sp(n), where
sp(n) = 6k+2 + sp(k+1) with k = ceil(n/4).
*/

static
void Toom4Mul (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
              long n, _ntl_ulong *stk)
{
  long k = (n + 3) / 4; /* ceil(n/4) */
  long r = n - 3 * k;
  _ntl_ulong cy,cy1,cy2,cy3,cy4;
  _ntl_ulong *W0 = c;
  _ntl_ulong *W1 = stk;
  _ntl_ulong *W2 = c + 2 * k;
  _ntl_ulong *W3 = stk + 2 * k;
  _ntl_ulong *W4 = c + 4 * k;
  _ntl_ulong *W5 = stk + 4 * k;
  _ntl_ulong *W6 = c + 6 * k;
  _ntl_ulong *newstk = stk + 6 * k + 2;

/* \\Evaluation phase: 13*2 add, 7*2 shift, 2Smul; 7 mul (n) */

/* W1 = U0 + U1 + U2 + U3    ; W2 = V0 + V1 + V2 + V3 */
  Add (W1,     a, a + 3 * k, r);
  Add1(W1+r, a+r, k-r, 0);
  Add3(W1, a + k, a + 2 * k, k); /* U0 + U1 + U2 + U3 */
  Add (W2+2,   b, b + 3 * k, r);
  Add1(W2+r+2,b+r,k-r, 0);
  Add3(W2+2,b+ k, b + 2 * k, k); /* V0 + V1 + V2 + V3 */
/*   Add (W1,     a, a + k,     k); */
/*   Add (W1,    W1, a + 2 * k, k); */
/*   Add (W1,    W1, a + 3 * k, r); /\* U0 + U1 + U2 + U3 *\/ */
/*   Add (W2+2,   b, b + k,     k); */
/*   Add (W2+2,W2+2, b + 2 * k, k); */
/*   Add (W2+2,W2+2, b + 3 * k, r); /\* V0 + V1 + V2 + V3 *\/ */

/* W0 = U1 +(U2 + U3*x)*x    ; W6 = V1 +(V2 + V3*x)*x */
  cy = AddLsh1 (W0, a + 2*k, a + 3 * k, r);   /* U2 + x U3 */
  cy = Add1 (W0 + r, a + 2*k + r, k - r, cy);
  W0[k] = (cy << 1) ^ AddLsh1 (W0, a+k, W0,k); /* U1+x U2 + x^2 U3 */
  cy = AddLsh1 (W6+2, b + 2*k, b + 3 * k, r);   /* V2 + x V3 */
  cy = Add1 (W6+2 + r, b + 2*k + r, k - r, cy);
  W6[k+2] = (cy << 1) ^ AddLsh1 (W6+2, b+k, W6+2,k); /* V1+x V2 + x^2 V3 */
  /* since we use W6[k+2], and we have space for 2r words in W6, we need
     k+3 <= 2*r, which requires n>=30. */
/* W4 = W1 +(W0 + U3*(x+1))*x; W3 = W2 +(W6 + V3*(x+1))*x */
  cy = AddLsh1 (W4+2, W0, a + 3 * k, r);   /* W0 + x U3 */
  cy = Add1 (W4+2+ r, W0 + r, k+1 - r, cy); /* cy == 0 */
  assert (cy == 0);
  Add (W4+2, W4+2, a + 3 * k, r); /* W0 + x U3 + U3 */
  W4[k+2] = (W4[k+2]<<1) ^ AddLsh1 (W4+2, W1, W4+2,k); /* W1+x(W0 +(x+1) U3) */
  cy = AddLsh1 (W3+2, W6+2, b + 3 * k, r);   /* W6 + x V3 */
  cy = Add1 (W3+2+ r, W6+2 + r, k+1 - r, cy); /* cy == 0 */
  Add (W3+2, W3+2, b + 3 * k, r); /* W6 + x V3 + V3 */
  W3[k+2] = (W3[k+2]<<1)^AddLsh1 (W3+2,W2+2,W3+2,k); /* W2+x(W6 + (x+1) V3) */
/* W0 = W0*x + U0            ; W6 = W6*x + V0 */
  W0[k] = (W0[k] << 1) ^ AddLsh1 (W0,  a,W0,k); /* U0+x W0 */
  W6[k+2]=(W6[k+2]<<1) ^ AddLsh1 (W6+2,b,W6+2,k); /* V0+x W6 */
/* W5 = W4 * W3              ; W4 = W0 * W6 */
  Toom (W5, W4+2, W3+2, k+1, newstk); /* W5 : 2*k+1 */

  Toom (W4, W0, W6+2, k+1, newstk); /* W4 : 2*k+1 */
  cy4=W6[0];/* Take care of overlapping byte. */
/* W3 = W1 * W2 */
  Toom (W3, W1, W2+2, k, newstk); /* W3 : 2*k */

/* W0 =(U2 +(U1 + U0*x)*x)*x ; W6 =(V2 +(V1 + V0*x)*x)*x */
  cy = AddLsh1 (W0, a + 1 * k, a + 0 * k, k);   /* U1 + x U0 */
  W0[k] = (cy << 2) ^ AddLsh1Lsh2(W0, a+2*k, W0,k); /* U2+x U1 + x^2 U0 */
  cy = AddLsh1 (W6+2, b + 1 * k, b + 0 * k, k);   /* V1 + x V0 */
  W6[k+2]=(cy << 2) ^ AddLsh1Lsh2(W6+2, b+2*k, W6+2,k); /* V2+x V1 + x^2 V0 */
/* W1 = W1 + W0 + U0*(x^2+x) ; W2 = W2 + W6 + V0*(x^2+x) */
  W1[k] = AddMul1 (W1, W1, a, k, 4+2);
  Add( W0+k+1, W1, W0, k+1);
  W2[k+2] = AddMul1 (W2+2, W2+2, b, k, 4+2);
  Add( W2+2, W2+2, W6+2, k+1);
/* W0 = W0 + U3              ; W6 = W6 + V3 */
  Add (W0,    W0, a + 3 * k, r); /* + U3 */
  Add (W6+2,W6+2, b + 3 * k, r); /* + V3 */
/* W1 = W1 * W2              ; W2 = W0 * W6 */
  cy=W3[0];cy2=W3[1];/* Take care of overlapping byte. */
  Toom (W1, W0+k+1, W2+2, k+1, newstk); /* W1 : 2*k+1 */
  cy1=W3[0];W3[0]=cy;W3[1]=cy2;
  cy=W4[0];cy2=W4[1];/* Take care of overlapping byte. */
  Toom (W2, W0, W6+2, k+1, newstk); /* W2 : 2*k+1 */
  W4[1]=cy2;cy2=W4[0];W4[0]=cy;
/* W6 = U3 * V3              ; W0 = U0 * V0 */
  Toom (W0, a, b, k, newstk); /* W0 : 2*k */
  Toom (W6, a+3*k, b+3*k,    r, newstk); /* W6 : 2*r */
/* \\Interpolation: 22 add, 4 shift, 5 Smul, 4 div (2n) */
/* d1=(x^4+x)*Mod(1,2)	; d1== (x)^1*(x+1)^1*(x^2+x+1)^1 *Mod(1,2) */
/* d2=(x^4+x^2)*Mod(1,2)	; d2== (x)^2*(x+1)^2*(x^2+x+1)^0 *Mod(1,2) */

/* W1 = W1 + W2 + W0*(x^4+x^2+1) */

  Add(W1,W1,W2,2*k);
  cy1 ^= cy2 ^ AddMul1(W1,W1,W0,2*k,16+4+1);

/* W5 =(W5 + W4 + W1 + W6*(x^4+x^2+1))/d1 */
  Add3(W5,W4,W1,2*k); W5[2*k] ^= cy1 ^ cy4;
  W5[2*r] ^= AddMul1(W5,W5,W6,2*r,16+4+1);
  DivExactD1(W5,2*k+1);

/* W2 = W2 + W6 + W0*(x^6) */
  Add(W2,W2,W6,2*r);
  cy2 ^= AddLsh6(W2,W2,W0,2*k);
/* W4 = W4 + W2 + W0 + W6*(x^6) */
  Add3 (W4, W2, W0, 2 * k);
  cy3 = AddLsh6(W4,W4,W6,2*r);
  cy = W6[0]; /* save W6[0]=W4[2k]: we cannot do it before the AddLsh6 call
                 because W6 is used as input */
  W6[0] = cy4 ^ cy2;
  W4[2*r] ^= cy3; /* must come after W6[0] = cy4 in case r=k */
/* W4 =(W4 + W5*(x^5+x))/d2 */
  AddMul1(W4,W4,W5,2*k+1,32+2);
  DivExactD2 (W4, 2 * k + 1);
  W6[0]=cy;
/* W3 = W3 + W0 + W6 */
  Add3(W3,W0,W6,2*r);
  if ( r != k )
    Add(W3+2*r,W3+2*r,W0+2*r,2*(k-r)); /* warning: 2r instead of r */
/* W1 = W1 + W3 */
  Add(W1,W1,W3,2*k);
/* W2 = W2 +(W1 + W3*x)*x */
  cy2^= AddLsh1(W2,W2,W1,2*k) ^ AddLsh2(W2,W2,W3,2*k) ^ (cy1<<1);
/* W3 = W3 + W4 + W5 */
  Add3(W3,W4,W5,2*k);
/* W1 =(W1 + W3*(x^2+x))/d1 */
  cy=W3[0];
  cy1^=AddMul1(W1,W1,W3,2*k,4+2);
  W3[0]=cy1;
  DivExactD1(W1,2*k+1);
  W3[0]=cy;
/* W5 = W5 + W1 */
  Add(W5,W5,W1,2*k);
/* W2 =(W2 + W5*(x^2+x))/d2 */
  cy=W4[0];
  W4[0]=cy2^AddMul1(W2,W2,W5,2*k,4+2);
  DivExactD2 (W2, 2 * k + 1);
  W4[0]=cy;
/* W4 = W4 + W2 */
  Add(W4,W4,W2,2*k);

/* \\Recomposition */
/* W = W6*Y^6 + W5*Y^5 + W4*Y^4+ W3*Y^3+ W2*Y^2+ W1*Y + W0 */
  Add (c + k, c + k, W1, 6 * k);

}
