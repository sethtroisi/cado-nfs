/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02111-1307, USA.
*/

/* Subquadratic gcd over GF(2)[x]. */

// #define BOGONG	// Selects routines that run well on bogong,
			// a 2.2 Ghz AMD Opteron 275 at ANU

#undef STATS		// If defined, statistics will be gathered
			// and printed on stdout


#include <NTL/GF2X.h>

#include "halfgcd.hpp"

/* This threshold controls the internal calls to HalfGCD. */

#define NTL_GF2X_HalfGCD_CROSSOVER (4*NTL_BITS_PER_LONG)

/* this threshold controls the calls in FastGCD, thus is less sensitive */
#define NTL_GF2X_GCD_CROSSOVER (5*NTL_GF2X_HalfGCD_CROSSOVER)

#ifdef STATS
long mulktm = 0x3fffff;		// mask for mul1kt, mul2kt
long addktm = 0xffff;		// mask for Mul1kt, Add1kt
long mul1kt = 0;		// counts calls to mul1 or mul1rpb
long mul2kt = 0;		// counts calls to mul2 or mul2t
long Mul1kt = 0;		// counts calls to mul_1_n
long Mul1kts = 0;		// sum of sizes sb in calls to mul_1_n
long Add1kt = 0;		// counts calls to addmul_1_n
long Add1kts = 0;		// sum of sizes sb in calls to addmul_1_n
#endif

#include "gf2x.h"

using namespace NTL;

class GF2XMatrix {
private:

   GF2XMatrix(const GF2XMatrix&);  // disable
   GF2X elts[2][2];

public:

   GF2XMatrix() { }
   ~GF2XMatrix() { }

   void operator=(const GF2XMatrix&);
   GF2X& operator() (long i, long j) { return elts[i][j]; }
   const GF2X& operator() (long i, long j) const { return elts[i][j]; }
};

void mul_gen(GF2X& c, const GF2X& a, const GF2X& b)
{
   long sa = a.xrep.length();
   long sb = b.xrep.length();

   if (sa <= 0 || sb <= 0) {
      clear(c);
      return;
   }
   c.xrep.SetLength(sa+sb);
   /* gf2x_mul now properly handles possible aliasing a==c and b==c */
   gf2x_mul(c.xrep.elts(), a.xrep.elts(), sa, b.xrep.elts(), sb);
   c.normalize();
}


void
mul (GF2X& U, GF2X& V, const GF2XMatrix& M)
// (U, V)^T = M*(U, V)^T
{
   GF2X RU, RV, R1, R2;
   // RU = U
   // RV = V

   // R1 = M(0,0)
   mul_gen(R1, M(0,0), U);
   // R2 = M(0,1)
   mul_gen(R2, M(0,1), V);
   add(R2, R1, R2);

   // R1 = M(1,0)
   mul_gen(R1, M(1,0), U);
   U = R2; // previous value
   // R2 = M(1,1)
   mul_gen(R2, M(1,1), V);
   add(V, R1, R2);
}

void
mul (GF2XMatrix& A, GF2XMatrix& B, GF2XMatrix& C)
// A = B*C, B and C are destroyed
{
   GF2X B00, B01, B10, B11, C0, C1, T1, T2;

   mul_gen(T1, B(0,0), C(0,0));
   mul_gen(T2, B(0,1), C(1,0));
   add(A(0,0), T1, T2);

   mul_gen(T1, B(1,0), C(0,0));
   mul_gen(T2, B(1,1), C(1,0));
   add(A(1,0), T1, T2);

   C(0,0).kill();
   C(1,0).kill();

   mul_gen(T1, B(0,0), C(0,1));
   mul_gen(T2, B(0,1), C(1,1));
   add(A(0,1), T1, T2);

   mul_gen(T1, B(1,0), C(0,1));
   mul_gen(T2, B(1,1), C(1,1));
   add(A(1,1), T1, T2);

   C(0,1).kill();
   C(1,1).kill();

   B(0,0).kill();
   B(1,0).kill();
   B(0,1).kill();
   B(1,1).kill();
}

void
strassen_mul (GF2XMatrix& C, GF2XMatrix& A, GF2XMatrix& B)
// C = A*B, A and B are destroyed
/* we follow the code from SAGE 1.6, file strassen.pyx */
{
  GF2X S0, T0, S1, T1, S2, T2, S3, T3, P0, P1, P2, P3, P4, P5, P6;

  add (S0, A(1,0), A(1,1));
  add (T0, B(0,1), B(0,0));
  add (S1, S0, A(0,0));
  add (T1, B(1,1), T0);
  add (S2, A(0,0), A(1,0));
  add (T2, B(1,1), B(0,1));
  add (S3, A(0,1), S1);
  add (T3, B(1,0), T1);
  mul_gen (P0, A(0,0), B(0,0));
  mul_gen (P1, A(0,1), B(1,0));
  mul_gen (P2, S0, T0);
  mul_gen (P3, S1, T1);
  mul_gen (P4, S2, T2);
  mul_gen (P5, S3, B(1,1));
  mul_gen (P6, A(1,1), T3);
  add (C(0,0), P0, P1);     /* U0 */
  add (C(0,1), P0, P3);     /* U1 */
  add (C(1,1), C(0,1), P4); /* U2 */
  add (C(1,0), C(1,1), P6); /* U3 */
  add (C(1,1), C(1,1), P2); /* U4 */
  add (C(0,1), C(0,1), P2); /* U5 */
  add (C(0,1), C(0,1), P5); /* U6 */
  A(0,0).kill();
  A(1,0).kill();
  A(0,1).kill();
  A(1,1).kill();
  B(0,0).kill();
  B(1,0).kill();
  B(0,1).kill();
  B(1,1).kill();
}

void
IterHalfGCD (GF2XMatrix& M_out, GF2X& U, GF2X& V, long d_red)
{
   M_out(0,0).SetMaxLength(d_red);
   M_out(0,1).SetMaxLength(d_red);
   M_out(1,0).SetMaxLength(d_red);
   M_out(1,1).SetMaxLength(d_red);

   set(M_out(0,0));   clear(M_out(0,1));
   clear(M_out(1,0)); set(M_out(1,1));

   long goal = deg(U) - d_red;

   if (deg(V) <= goal)
      return;

   GF2X Q, t(INIT_SIZE, d_red);

   while (deg(V) > goal) {
      PlainDivRem(Q, U, U, V);
      swap(U, V);

      mul_gen(t, Q, M_out(1,0));
      add(t, M_out(0,0), t); // add=sub over GF2
      M_out(0,0) = M_out(1,0);
      M_out(1,0) = t;

      mul_gen(t, Q, M_out(1,1));
      add(t, M_out(0,1), t); // add=sub over GF2
      M_out(0,1) = M_out(1,1);
      M_out(1,1) = t;
   }
}

void GF2XMatrix::operator=(const GF2XMatrix& M)
{
   elts[0][0] = M.elts[0][0];
   elts[0][1] = M.elts[0][1];
   elts[1][0] = M.elts[1][0];
   elts[1][1] = M.elts[1][1];
}

// c <- a mod x^n, with n >= 0
void
RightShiftRem (GF2X& c, const GF2X& a, long n)
{
  if (deg(a) < n) // covers a=0 too
    {
      c = a;
      return;
    }

  long sa = a.xrep.length(); // number of words of a
  long wn = n / NTL_BITS_PER_LONG;
  long bn = n - wn * NTL_BITS_PER_LONG; // #bits of ap[wn]
  if (wn >= sa)
    {
      wn = sa;
      bn = 0;
    }

  c.xrep.SetLength (wn + (bn != 0));

   _ntl_ulong *cp = c.xrep.elts();
   const _ntl_ulong *ap = a.xrep.elts();

   long i;

   for (i = 0; i < wn; i++)
     cp[i] = ap[i];
   if (bn)
     cp[wn] = ap[wn] & ((1UL << bn) - 1UL);

   c.normalize();
}

// same as HalfGCD, except U and V are replaced by U' and V'
// if extended is zero, the matrix M_out is not computed, and the inputs are
// fully reduced.
// We want to reduce deg(V) to <= deg(U) - d_red
// Usually if deg(U) = m, then d_red = m/2
void
HalfGCD2 (GF2XMatrix& M_out, GF2X& U, GF2X& V, long d_red, int extended)
{
  long degU = deg(U);

   if (IsZero(V) || deg(V) <= degU - d_red)
     {
       if (extended != 0)
         {
           set(M_out(0,0));   clear(M_out(0,1));
           clear(M_out(1,0)); set(M_out(1,1));
         }
       return;
     }

   if (d_red <= NTL_GF2X_HalfGCD_CROSSOVER)
     {
       IterHalfGCD(M_out, U, V, d_red);
       return;
     }

   long d1 = (d_red + 1) / 2;       /* d1 is about m/4 */
   if (extended == 0) d1 = deg(U) / 4;
   if (d1 < 1) d1 = 1;
   if (d1 >= d_red) d1 = d_red - 1;

   long n = degU - 2 * d1 + 2;      /* n is the ignored part, about m/2 */
   if (n < 0) n = 0;

   GF2X U1, V1;

   if (n != 0)
     {
       RightShiftRem (U1, U, n); /* U1 = U mod x^n: m/2 bits */
       RightShiftRem (V1, V, n); /* V1 = V mod x^n: m/2 bits */
       RightShift (U, U, n);     /* U = U div x^n has about m/2 bits */
       RightShift (V, V, n);     /* V = V div x^n has about m/2 bits */
     }

   GF2XMatrix M1;

   HalfGCD2 (M1, U, V, d1, extended || (n != 0));
   /* the entries of M1 have m/4 bits, and U,V have been reduced to m/4 bits */

   if (n != 0)
     {
       LeftShift (U, U, n);
       LeftShift (V, V, n);
       mul (U1, V1, M1); /* U1,V1:m/2 M1:m/4 cost: 4 M(m/2,m/4) */
       add (U, U, U1);
       add (V, V, V1);
     }

   /* now U, V have 3m/4 bits */

   // FIXME (24 June 2007): shouldn't it be deg(U) instead of deg(V) below?
   long d2 = deg(V) - (degU - d_red); /* should be about m/2 */

   if (IsZero(V) || d2 <= 0)
     {
       if (extended != 0)
         M_out = M1;
       return;
     }

   GF2X Q;
   GF2XMatrix M2;

   DivRem (Q, U, U, V);
   swap (U, V);

   if (extended)
     {
       GF2X t;

       mul_gen (t, Q, M1(1,0));
       add (t, M1(0,0), t);
       swap (M1(0,0), M1(1,0));
       swap (M1(1,0), t);

       mul_gen (t, Q, M1(1,1));
       add (t, M1(0,1), t);
       swap (M1(0,1), M1(1,1));
       swap (M1(1,1), t);
     }
   else
     d2 = deg(U); // fully reduce V

   if (IsZero(V) || deg(V) <= degU - d_red)
     {
       if (extended)
         M_out = M1;
       return;
     }

   n = deg(U) - 2 * d2 + 2; /* should be about m/4 */
   if (n < 0 || extended == 0) n = 0;

   if (n != 0)
     {
       RightShiftRem (U1, U, n); /* U1,V1 have m/4 bits */
       RightShiftRem (V1, V, n);
       RightShift(U, U, n);      /* U,V have m/2 bits */
       RightShift(V, V, n);
     }

   HalfGCD2 (M2, U, V, d2, extended || (n != 0));
   /* now U,V have m/4 bits, like the entries of M2 */

   if (n != 0)
     {
       LeftShift (U, U, n);
       LeftShift (V, V, n);
       mul (U1, V1, M2); /* 4 M(m/4) */
       add (U, U, U1);
       add (V, V, V1);
     }

   if (extended)
     strassen_mul (M_out, M2, M1); /* 8 M(m/4) */
}

/* in the general case, we have deg(u) and deg(v) about the same here */
void FastGCD(GF2X& d, const GF2X& u, const GF2X& v)
{
   GF2X u1, v1;
   GF2XMatrix M1;

   u1 = u;
   v1 = v;

   if (deg(u1) == deg(v1)) {
      if (IsZero(u1)) {
         clear(d);
         return;
      }

      rem(v1, v1, u1);
   }
   else if (deg(u1) < deg(v1)) {
      swap(u1, v1);
   }

   // deg(u1) > deg(v1)

   /* in practice this while loop will be performed only once,
      since HalfGCD fully reduces the inputs */
   while (deg(u1) > NTL_GF2X_GCD_CROSSOVER && !IsZero(v1))
     {
       HalfGCD2 (M1, u1, v1, deg (u1), 0);

       if (!IsZero(v1))
         {
           rem(u1, u1, v1);
           swap(u1, v1);
         }
     }

   // base case.
   GCD(d, u1, v1);
}

