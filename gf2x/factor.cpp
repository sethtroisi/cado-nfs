/* Finds smallest irreducible factor of trinomial over GF(2)

  Copyright 2007 Richard P. Brent, Paul Zimmermann.

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

#define VERSION 1.72

#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for CHAR_BIT */
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZXFactoring.h>
#include <signal.h>
#include <string.h>

NTL_CLIENT

/* Usage: factor [options] r < file [see usage() function for options]
   where file contains lines of the form

      s nnn
   or s u.

   "s nnn" means that x^r+x^s+1 has a factor of degree nnn,
   and this is output as is. In fact anything except "s uxxx" is output as is
   and "s uxxx" is treated like "s u".

   "s u" means that the status of x^r+x^s+1 is unknown; 
   it is checked and "s k px..y" is output if it has a factor of degree k
   (the factor is x..y in hex encoding);
   otherwise "s irreducible" (or "s primitive") is output.
   If the "-z" option is used then "s u" may be output (in this case
   we need to check the trinomial with another program, e.g. irred).

   r and s should satisfy r > s > 0.
*/

#include "HalfGCD.c"

#define BITS_PER_LONG2 (2*NTL_BITS_PER_LONG)

void fastsqr(GF2X& b, GF2X& a, _ntl_ulong r, _ntl_ulong s)	

/* Returns b = a*a mod (x^r + x^s + 1) 

   Here a and b represent polynomials (in GF(2)[x]) of degree < r,
   using the natural representation (as used in NTL and irred) of 
   one bit per coefficient (least significant bit of word 0 is the
   constant term).
   
   The input a is destroyed and must not overlap the output b.
 
   r and s must be odd. In addition, this implementation assumes that
   min(s, r-s) > 2*wlen, where wlen is the wordlength in bits (assumed 
   to be 32 or 64).
    
   The algorithm is as described in Section 4 of Brent, Larvala and 
   Zimmermann, Math. Comp. 72 (2003), 1443--1452 [rpb199], but without
   any overlap of the arrays.
   
   The algorithm has two steps (reduce and interleave). The first step
   reduces the implicit square of polynomial a of degree < r "in place"
   given x^r = x^s + 1 (mod 2), with the result left in scrambled order 
   0 2 4 ... r-3 r-1 1 3 5 ... r-4 r-2 in a.  The second step then 
   "interleaves" the bits to give the correct order in b.
   
   Uses about (2*load + 1.5*store)*r/wlen memory operations, which is
   better than the obvious "square and reduce" algorithm as in 
   NTL's SqrMod. The saving is because operations are not performed
   on bits that are known to be zero.
   
   RPB, 20070318 (based on reducep, reducer and interlvf in irred.c). */

  {
  _ntl_ulong old, nex, mask1, mask2; 	/* "new" is a reserved word */
  _ntl_ulong t, u, v, w, next1, next2;

  long j, alpha, delta, deltaw, deltaq, deltaqc, q1, q4, q5, s1, s2;
  long smax;
  	
#if (NTL_BITS_PER_LONG == 64)

  _ntl_ulong c0, c1, c2, c3, c4, c5;
  			
  long wd = 6;			/* 64-bit words */
  unsigned long wlen = 64;	/* wordlength (bits) = 2^wd */
  long wlenm = 63;		/* wlen - 1 */

  if ((r&s&1) == 0) 
    {
      printf("Error in fastsqr: r = %lu and s = %lu must be odd\n", r, s);
      exit (1);
    }
    
  if (s <= 2*wlen)
    {  
      printf("Error in fastsqr: s = %lu must be at least %lu\n", s, 2*wlen+1);
      exit (1);
    }
  
  if ((r-s) <= 2*wlen)
    {  
      printf("Error in fastsqr: r-s = %lu must be at least %lu\n",
	     r-s, 2*wlen+1);
      exit (1);
    }
  
  smax = (r+wlen-1) >> wd;	/* max size in words = ceil(r/wlen) */

  long sa = a.xrep.length();	/* size in words */

  if (sa < smax)
    a.xrep.SetLength(smax);	/* Always use length smax for a */

  _ntl_ulong *ap = a.xrep.elts();

  for (j = sa; j < smax; j++) 	/* Clear high words of a if necessary */
    ap[j] = 0;  
    
  b.xrep.SetLength(smax+2);	/* Allow 2 words extra space for b */

  _ntl_ulong *bp = b.xrep.elts();

  bp[smax] = 0;			/* Clear the extra two words of b */
  bp[smax+1] = 0;		/* No need to clear words 0..(smax-1) */

  c0 = 0x00000000FFFFFFFFL;	/* Some 64-bit constants */
  c1 = 0x0000FFFF0000FFFFL;
  c2 = 0x00FF00FF00FF00FFL;
  c3 = 0x0F0F0F0F0F0F0F0FL;
  c4 = 0x3333333333333333L;
  c5 = 0x2222222222222222L;

  alpha = r >> 1;		/* alpha = (r-1)/2 */
  delta = (r - s) >> 1;		/* delta = (r - s)/2 */
  deltaw = delta >> wd;		/* deltaw = delta div wlen */
  deltaq = delta & wlenm;	/* deltaq = delta mod wlen */
  deltaqc = wlen - deltaq;
  q1 = (r-1) >> wd;		/* q1 = (r-1) div wlen = smax - 1 */
  q4 = alpha >> wd;		/* q4 = alpha div wlen */
  
  mask1 = (_ntl_ulong)(~0L) >> (wlenm - ((r-1) & wlenm));
  				/* mask1 has wlen-1 - ((r-1) mod wlen)
  				   zero bits in high positions */
  mask2 = (~1L) << (alpha & wlenm);
  				/* mask2 has (alpha mod wlen) + 1 zero bits
  				   in low positions */

  ap[q1] &= mask1;		/* Mask irrelevant high bits of a */
  nex = 0;	
  
  if (deltaq == 0) 
    {				/* Special case, deltaq == 0 (word-aligned) */
    for (j = q1; j > q4; j--)
      ap[j-deltaw] ^= ap[j];
      
    ap[q4-deltaw] ^= ap[q4] & mask2;
    }

  else 
    {				/* Usual case, 0 < deltaq < wlen */
    for (j = q1; j > q4; j--) 
      {
      old = nex;
      nex = ap[j];
      ap[j-deltaw] ^= (nex >> deltaq) | (old << deltaqc);
      }

    /* The last two iterations are special as need to mask some bits */

    old = nex;
    nex = ap[q4] & mask2;
    ap[q4-deltaw]   ^= (nex >> deltaq) | (old << deltaqc);
    ap[q4-deltaw-1] ^= (nex << deltaqc);
    }
  
/* Interleave starts here */

  q5 = (alpha+1) >> wd;		/* q5 = (alpha+1) div wlen */
  s1 = (alpha+1) & wlenm;	/* s1 = (alpha+1) mod wlen */
  s2 = wlenm - s1;		/* In [0, wlen) */
  
  next1 = ap[0];
  old = ap[q5];
  nex = ap[q5+1];

  for (j = 0; j <= q5; j++) {

    next2 = (old >> s1) | ((nex << 1) << s2); /* Beware case s2 == 63 */

    u = next1>>32;			/* High order 32 bits low part of a */
    w = next2>>32;			/* Ditto high part of a */
    t = next1 & c0;			/* Low order 32 bits low part of a */
    v = next2 & c0;			/* Ditto high part of a */

    next1 = ap[j+1];
    old   = nex;
    nex   = ap[j+q5+2];

    u = (u | u<<16) & c1;		 /* Operations on t,u,v,w   */
    t = (t | t<<16) & c1;		 /* can be done in parallel */
    w = (w | w<<16) & c1;
    v = (v | v<<16) & c1;

    u = (u | u<<8) & c2;
    t = (t | t<<8) & c2;
    w = (w | w<<8) & c2;
    v = (v | v<<8) & c2;

    u = (u | u<<4) & c3;
    t = (t | t<<4) & c3;
    w = (w | w<<4) & c3;
    v = (v | v<<4) & c3;

    u = (u | u<<2) & c4;
    t = (t | t<<2) & c4;
    w = (w | w<<2) & c4;
    v = (v | v<<2) & c4;

    u += u & c5;
    t += t & c5;
    w += w & c5;
    v += v & c5;

    bp[2*j+1] = u | (w << 1);
    bp[2*j]   = t | (v << 1);
    }

  bp[smax+1] = 0;				/* Clear high words of b */
  bp[smax] = 0;					/* Note: smax = q1 + 1 */
  bp[q1] &= mask1;				/* Mask high bits */
  b.normalize(); 	     			/* Normalize output */
  return;
  }

#else
				/* This is for 32-bit words */  			
  _ntl_ulong c1, c2, c3, c4, c5;

  long wd = 5;
  long wlen = 32;		/* wordlength (bits) = 2^wd */
  long wlenm = 31;		/* wlen - 1 */

  if ((r&s&1) == 0) 
    {
    printf("Error in fastsqr: r = %d and s = %d must be odd\n", r, s);
    exit (1);
    }
    
  if (s <= 2*wlen)
    {  
    printf("Error in fastsqr: s = %d must be at least %d\n", s, 2*wlen+1);
    exit (1);
    }
  
  if ((r-s) <= 2*wlen)
    {  
    printf("Error in fastsqr: r-s = %d must be at least %d\n", r-s, 2*wlen+1);
    exit (1);
    }
  
  smax = (r+wlen-1) >> wd;	/* max size in words = ceil(r/wlen) */

  long sa = a.xrep.length();	/* size in words */

  if (sa < smax)
    a.xrep.SetLength(smax);	/* Always use length smax for a */

  _ntl_ulong *ap = a.xrep.elts();

  for (j = sa; j < smax; j++) 	/* Clear high words of a if necessary */
    ap[j] = 0;  
    
  b.xrep.SetLength(smax+2);	/* Allow 2 words extra space for b */

  _ntl_ulong *bp = b.xrep.elts();

  bp[smax] = 0;			/* Clear the extra two words of b */
  bp[smax+1] = 0;		/* No need to clear words 0..(smax-1) */

  c1 = 0x0000FFFFL;		/* Some 32-bit constants */
  c2 = 0x00FF00FFL;
  c3 = 0x0F0F0F0FL;
  c4 = 0x33333333L;
  c5 = 0x22222222L;

  alpha = r >> 1;		/* alpha = (r-1)/2 */
  delta = (r - s) >> 1;		/* delta = (r - s)/2 */
  deltaw = delta >> wd;		/* deltaw = delta div wlen */
  deltaq = delta & wlenm;	/* deltaq = delta mod wlen */
  deltaqc = wlen - deltaq;
  q1 = (r-1) >> wd;		/* q1 = (r-1) div wlen = smax - 1 */
  q4 = alpha >> wd;		/* q4 = alpha div wlen */
  
  mask1 = (_ntl_ulong)(~0L) >> (wlenm - ((r-1) & wlenm));
  				/* mask1 has wlen-1 - ((r-1) mod wlen)
  				   zero bits in high positions */
  mask2 = (~1L) << (alpha & wlenm);
  				/* mask2 has (alpha mod wlen) + 1 zero bits
  				   in low positions */

  ap[q1] &= mask1;		/* Mask irrelevant high bits of a */
  nex = 0;	
  
  if (deltaq == 0) 
    {				/* Special case, deltaq == 0 (word-aligned) */
    for (j = q1; j > q4; j--)
      ap[j-deltaw] ^= ap[j];
      
    ap[q4-deltaw] ^= ap[q4] & mask2;
    }

  else 
    {				/* Usual case, 0 < deltaq < wlen */
    for (j = q1; j > q4; j--) 
      {
      old = nex;
      nex = ap[j];
      ap[j-deltaw] ^= (nex >> deltaq) | (old << deltaqc);
      }

    /* The last two iterations are special as need to mask some bits */

    old = nex;
    nex = ap[q4] & mask2;
    ap[q4-deltaw]   ^= (nex >> deltaq) | (old << deltaqc);
    ap[q4-deltaw-1] ^= (nex << deltaqc);
    }
  
/* Interleave starts here */

  q5 = (alpha+1) >> wd;		/* q5 = (alpha+1) div wlen */
  s1 = (alpha+1) & wlenm;	/* s1 = (alpha+1) mod wlen */
  s2 = wlenm - s1;		/* In [0, wlen) */
  
  next1 = ap[0];
  old = ap[q5];
  nex = ap[q5+1];

  for (j = 0; j <= q5; j++) {

    next2 = (old >> s1) | ((nex << 1) << s2); /* Beware case s2 == 31 */

    u = next1>>16;			/* High order 16 bits low part of a */
    w = next2>>16;			/* Ditto high part of a */
    t = next1 & c1;			/* Low order 16 bits low part of a */
    v = next2 & c1;			/* Ditto high part of a */

    next1 = ap[j+1];
    old   = nex;
    nex   = ap[j+q5+2];

    u = (u | u<<8) & c2;		 /* Operations on t,u,v,v   */
    t = (t | t<<8) & c2;		 /* can be done in parallel */
    w = (w | w<<8) & c2;
    v = (v | v<<8) & c2;

    u = (u | u<<4) & c3;
    t = (t | t<<4) & c3;
    w = (w | w<<4) & c3;
    v = (v | v<<4) & c3;

    u = (u | u<<2) & c4;
    t = (t | t<<2) & c4;
    w = (w | w<<2) & c4;
    v = (v | v<<2) & c4;

    u += u & c5;
    t += t & c5;
    w += w & c5;
    v += v & c5;

    bp[2*j+1] = u | (w << 1);
    bp[2*j]   = t | (v << 1);
    }

  bp[smax+1] = 0;				/* Clear high words of b */
  bp[smax] = 0;					/* Note: smax = q1 + 1 */
  bp[q1] &= mask1;				/* Mask high bits */
  b.normalize(); 	     			/* Normalize output */
  return;
  }

#endif

void fastrem (GF2X& a, _ntl_ulong r, _ntl_ulong s)

// Sets a := a mod (x^r + x^s + 1) using a word operations where possible.
// This is much faster than NTL's rem if s > r/2.
// Assumes 0 < s < r and 32-bit or 64-bit words.	RPB 20070321

{
  long d, jw, krb, krw, ksb, ksw; 
  long apr, aps, lastr, lasts, sh1, sh2;
  unsigned long rs, j, wlenr, jb;

#if (NTL_BITS_PER_LONG == 64)
				/* Assume wordlength is 32 or 64 bits */
  				/* This is for 64 bits */
  long wd = 6;
  unsigned long wlen = 64;	/* wordlength (bits) = 2^wd */
  unsigned long wlenm = 63;	/* wlen - 1 */
#else
  long wd = 5;			/* This is for 32 bits */
  unsigned long wlen = 32;	/* wordlength (bits) = 2^wd */
  unsigned long wlenm = 31;	/* wlen - 1 */
#endif

  _ntl_ulong high;
  int slow;

  _ntl_ulong *ap = a.xrep.elts();

  rs = r - s;

  if ((s <= 0) || (rs <= 0)) 
    {
      printf("Error in fastrem: r = %lu and s = %lu must satisfy r > s > 0\n", 
	     r, s);
      exit (1);
    }

  d = deg(a);
  if (d < (long) r)
    return;

  slow = ((s <= 2*wlen) || (rs <= 2*wlen));	// Avoid word ops if slow
  wlenr = wlen + r;
  
  for (j = d; j >= r;)
    {
    jw = j>>wd;			// (j div wlen) = word index
    jb = j&wlenm;		// (j mod wlen) = bit index in word jw
    if (slow || (jb != wlenm) || (j < wlenr))
      {
      if (rep (coeff (a, j)) != 0)		// This is slow but only
        {					// executed a few times
        SetCoeff(a, j, 0);
        SetCoeff(a, j-r,  !(rep (coeff (a, j-r))));
        SetCoeff(a, j-rs, !(rep (coeff (a, j-rs))));
        }
      j--;					// Process next bit
      }
    else					// Process by words
      {						// Loop precomputation
      krb = (j-r)&wlenm;			// bit index in [0..wlen)
      krw = (j-r)>>wd;				// word index > 0
      ksb = (j-rs)&wlenm;			// bit index in [0..wlen)
      ksw = (j-rs)>>wd;				// word index > 0
      sh1 = wlenm - krb;			// first shift amount
      sh2 = wlenm - ksb;			// second shift amount
      
      for (lastr = lasts = 0; krw > 0; )
        {					// Loop to process by words
        high = ap[jw]; 
        apr  = ap[krw];				// In this version we avoid
        aps  = ap[ksw];				// the " ^= " operator in
        ap[jw--]  = 0;				// the inner loop.
        ap[krw--] = (apr^lastr)^(high >> sh1);	// This reduced dependencies.
        ap[ksw--] = (aps^lasts)^(high >> sh2);
	high  = high << 1;
        lastr = high << krb;			// krb == wlenm possible
        lasts = high << ksb;			// ksb == wlenm possible
        }
      ap[krw] ^= lastr;
      ap[ksw] ^= lasts;  
      j = r + krb;				// krb bits remain
      }    
    }
  a.normalize(); 	     			/* Normalize output */
}

int lex_less (GF2X a, GF2X b, int flip)

/* Returns true if 
  (deg(a) < deg(b)) or ((deg(a) = deg(b)) and (a < b in lex ordering))
  except if flip is true when this is changed appropriately. */
  
{
  long da = deg (a);
  long db = deg (b);
  long i;
  if (da < db) return 1;
  if (da > db) return 0;
  if (flip)
    for (i = 0; i <= da; i++)
      {
      if (rep (coeff (a, i)) < rep (coeff (b, i)))
        return 1;
      if (rep (coeff (a, i)) > rep (coeff (b, i)))
        return 0;
      }
  else
    for (i = da; i >= 0; i--)
      {
      if (rep (coeff (a, i)) < rep (coeff (b, i)))
        return 1;
      if (rep (coeff (a, i)) > rep (coeff (b, i)))
        return 0;
      }
  return 0;		/* Equal here, return false */
  
}

void print_hex (char *s, GF2X a, int flip)

/* prints a in hexadecimal (leading bits first unless flip true) */

{
  long d = deg (a);
  long i, j, c;
  char t[17] = "0123456789abcdef";

  c = 0;
  for (i = d; i >= 0; i--)
    {
      j = (flip) ? (d-i) : i;
      c = (2 * c) + rep (coeff (a, j));
      if (i % 4 == 0)
        {
          *s++ = t[c];
          c = 0;
        }
    }
  *s = '\0';
}

void
fastmulmod (GF2X& c, const GF2X& a, const GF2X& b, const GF2XModulus& F,
	         _ntl_ulong r, _ntl_ulong s)
	         
// Returns c = a*b mod (x^r + x^s + 1) 

{
  if ((s <= 0) || (r <= s)) 
    {
      printf("Error in fastmulmod: r = %lu and s = %lu must satisfy r > s > 0\n", 
	     r, s);
      exit (1);
    }

   GF2X t;
   mul_gen(t, a, b);

   if ((r&s&1) && (s > BITS_PER_LONG2) && ((r-s) > BITS_PER_LONG2))
     {
       fastrem(t, r, s);
       c = t;
     }
   else
     rem(c, t, F);        // This is slow if s > r/2  
}

void fastmulbyxmod (GF2X& c, const GF2X& a, _ntl_ulong r, _ntl_ulong s)
	         
// Returns c = a*x mod (x^r + x^s + 1) 

{
  if ((s <= 0) || (r <= s)) 
    {
      printf("Error in fastmulbyxmod: r = %lu and s = %lu must satisfy r > s > 0\n",
	     r, s);
      exit (1);
    }

   MulByX(c, a);
   if (deg(c) >= (long) r)
     {
     SetCoeff(c, r, 0);
     SetCoeff(c, s,  !(rep (coeff (c, s))));
     SetCoeff(c, 0,  !(rep (coeff (c, 0))));
     c.normalize();
     }
}

void copy_h (vec_GF2X& hold, vec_GF2X& h, long m)

/* copy h[0..m-1] into hold[0..m-1] */

{
  for (long i = 0; i < m; i++)
    hold[i] = h[i];
}

void init_h (vec_GF2X& h, long k, long m, const GF2XModulus& F, 
 	_ntl_ulong r, _ntl_ulong s, int usefs)

/* initialize h[j] to sum (x^(s*2^k), 0 <= s < 2^m, hw(s)=j+1), 
   for 0 <= j < m */

{
  long i, j;
  GF2X a, x, htemp;

  /* we use the following Pascal-triangle algorithm:
     if S(m,j) := { 0 <= s < 2^m, hw(s) = j }
     then an element s of S(m,j) is either in S(m-1,j) if s < 2^(m-1),
     or s = 2^(m-1) + t where t is in S(m-1,j-1):
        S(m,j) = S(m-1,j) union [S(m-1,j-1) + 2^(m-1)]
     Or equivalently:
        S(m,j) = {2s, s in S(m-1,j)} union {1+2s, s in S(m-1,j-1)}.
     With S(m,0) = { 0 } and S(m,m) = { 2^m-1 }.
  */
  clear(htemp);
  SetX (x);
  for (j = 1; j <= m; j++)
    {
      h[j-1] = 0;
      for (i = j-1; i >= 0; i--)
	{
	  if (i > 0)
	    {
	     if (usefs)
	       {
	       fastsqr (htemp, h[i-1], r, s);
	       h[i-1] = htemp;
	       }
	     else
	       SqrMod (h[i-1], h[i-1], F);	// Avoid if possible as slow
             fastmulbyxmod (a, h[i-1], r, s);
	    }
	  else
	    a = x;
	  add (h[i], h[i], a);
	}
    }
  for (j = 0; j < m; j++)
    {
    if (usefs)
      {
      for (i = 0; i < k/2; i++)
        {
        fastsqr (htemp, h[j], r, s);
        fastsqr (h[j], htemp, r, s);
        }
      if (k&1)		// last iteration if k odd
        {
        fastsqr (htemp, h[j], r, s);
        h[j] = htemp;	// copy result
        }  
      }
    else  
      for (i = 0; i < k; i++)
        SqrMod (h[j], h[j], F);	// Avoid if possible as slow
    }
}

void
usage ()
{
  fprintf (stderr, "Usage: factor [options] r < file\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "   -b d    - omit hex encoding of factor if degree < d\n");
  fprintf (stderr, "   -f f    - use interval bounds of f^j (default is f=1);\n");
  fprintf (stderr, "             if f=1 then intervals increase linearly in size\n");
  fprintf (stderr, "   -m m    - use double-blocking size of m\n");
  				// default is as for s if r small, or try to
  				// estimate if r large
  fprintf (stderr, "   -p      - assume 2^r-1 is a Mersenne prime so any irreducible\n");
  fprintf (stderr, "             trinomial is primitive\n");
  fprintf (stderr, "   -q q0   - start with q0 blocks of length m (default is q0=4)\n");
  fprintf (stderr, "   -k k0   - use double-blocking strategy only for degree > k0\n");
  fprintf (stderr, "             (default is <= largest k0 such that 2^k0 < r); assume k0 > 0\n");
  fprintf (stderr, "   -skip d - skip degrees <= d (input trinomials are known to have no\n");
  fprintf (stderr, "             factor of degree <= d)\n");
  fprintf (stderr, "   -maxd d - skip degrees > d\n");
  fprintf (stderr, "   -t      - print cumulative times for GCD etc at end\n");
  fprintf (stderr, "   -v      - verbose\n");
  fprintf (stderr, "   -z q1   - output 'u' if number of blocks of length m would be > q1\n");
  fprintf (stderr, "             (default is 'infinity')\n");
  exit (1);
}

/* copied from GMP-ECM, main.c */
static long last_degree = 0;
static long last_s = 0;
static int exit_asap_signalnr = 0; /* Remembers which signal we received */

void 
signal_handler (int sig)
{
  fprintf (stderr, "Interrupted s=%ld after degree %ld\n",
           last_s, last_degree);
  fflush (stderr);
  if (sig == SIGINT || sig == SIGTERM)
    {
      exit_asap_signalnr = sig;
      /* If one of these two signals arrives again, we'll let the default 
         handler take over,  which will usually terminate the process 
         immediately. */
      signal (SIGINT, SIG_DFL);
      signal (SIGTERM, SIG_DFL);
    }
}

void
save_state (long s, long k2)
{
  last_s = s;
  last_degree = k2;
}

int
main (int argc, char *argv[])
{
  long r, rhigh, maxd = LONG_MAX, k = 0, k1, k2, s, ss, i, q, r3; 
  long m = 0, k0 = 0, b0 = 2, q0 = 4, q1 = 0, q2 = 0;
  long skipd = 0, fineDDFtol;
  long tkt = 0;
  int fineDDF, flip, swan;
  char c, *fact;
  double f = 1.0; /* scaling factor, default 1.0 */
  double f0;

  int verbose = 0, verbose_t = 0, primitive = 0, choose_m = 0;
  int usefastsqr;

  double CanZass_t, GCD_t, SqrMul_t, Exp_t, Total_t;
  double SqrMul_kt1, SqrMul_kt2;
  double GCD_kt1, GCD_kt2;

  /* to catch interrupts */
  signal (SIGINT, &signal_handler);
  signal (SIGTERM, &signal_handler);

  /* print whole command line */
  for (i = 0; i < argc; i++)
    fprintf (stderr, "%s ", argv[i]);
  fprintf (stderr, "\n");

  fprintf (stderr, "This is factor version %1.2f\n", VERSION);

  CanZass_t = (double) 0.0;
  GCD_t = (double) 0.0;
  SqrMul_t = (double) 0.0;
  Exp_t = (double) 0.0;
  Total_t = (double) 0.0;
  SqrMul_kt1 = (double) 0.0;
  SqrMul_kt2 = (double) 0.0;
  GCD_kt1 = (double) 0.0;
  GCD_kt2 = (double) 0.0;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-v") == 0)
        {
          verbose ++;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "-t") == 0)
        {
          verbose_t = 1;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "-p") == 0)
        {
          primitive = 1;
          argc --;
          argv ++;
        }
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
        {
          m = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f") == 0)
        {
          f = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-k") == 0)
        {
          k0 = atoi (argv[2]);
	  if (k0 > NTL_BITS_PER_LONG)
	    {
	      fprintf (stderr, "Error, k0 should be <= NTL_BITS_PER_LONG\n");
	      exit (1);
	    }
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-b") == 0)
        {
          b0 = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-q") == 0)
        {
          q0 = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-z") == 0)
        {
          q1 = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-skip") == 0)
        {
          skipd = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-maxd") == 0)
        {
          maxd = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else
        {
          fprintf (stderr, "Error, invalid option: %s\n", argv[1]);
          exit (1);
        }
    }

  if (argc != 2)
    usage ();

  r = atoi (argv[1]);

  if (r < 2)
    {
      fprintf (stderr, "Error, r < 2\n");
      exit (1);
    }

  if ((primitive) && (!ProbPrime(r, 2)))
    {
      fprintf (stderr, 
	       "Error: -p option and r = %lu is not prime\n", r);
      exit (1);
    }

  if (k0 == 0)
    {
      for (; (2 << k0) < r; k0++)
	;
      if (k0 > 20) 
	k0 = (k0+20)/2;		// Empirical optimum, could check at runtime
    }

  if (m == 0)
    {
      choose_m = (r > 10000);	 // Later try to find optimal m
      for ( ; (2 << m) < r; m++) // Reasonable though probably not optimal.
	;
    }			

  if (r < 17) 			// Tiny case, this avoids a problem
    m = 1;			// with e.g. r = 11, s = 1.

  if (r < 1000)			// Small case, this avoids a problem
    q0 = 1;			// with e.g. r = 521, s = 62, failing
    				// to enter fineDDF mode (it was OK with 
    				// V1.29, but not with V1.30 .. V1.32)
  if (m > r / 2)
    {
      fprintf (stderr, "Error, m = %ld greater than r/2\n", m);
      exit (1);
    }

  if (f < (double)1.0)
    {
      fprintf (stderr, "Error, f < 1.0\n");
      exit (1);
    }

  if (q0 < 1)
    {
      fprintf (stderr, "Error, q0 < 1\n");
      exit (1);
    }

  if (q1 == 0)
    q1 = r;				// Effectively infinity
  
  for (r3 = 0; r3*r3*r3 < r; r3++) ;	// r3 = r^(1/3)
  fineDDFtol = 2*r3*r3;			// 2*r^(2/3)

  if (verbose)
    {
      printf ("Test separately degrees up to k0=%ld\n", k0);
      if (skipd > 0)
        printf ("Skip degrees up to %ld\n", skipd);
      printf ("Using block size of m=%ld\n", m);
      printf ("Start with q0=%ld block(s)\n", q0);
      if (q1 < r)
        printf ("Stop when number of blocks is > q1=%ld\n", q1);
      printf ("Use fine DDF above %ld\n", fineDDFtol);
      fflush (stdout);
    }

  fact = (char*) malloc ((2 + (r / 4)) * sizeof (char));

  /* since all degree k factors divide x^(2^k)-x, if x^r+x^s+1 has a degree
     k factor, then so has x^r+x^(s+2^k-1)+1, thus if 2^k-1 <= floor(r/2)-1,
     i.e. 2^(k+1) <= r, then we can save some computations */

  long j, l;

  while (1)
    {
      s = -1;		/* In case no positive integer found */
      if (scanf ("%ld ", &s) == EOF)
        break;

      c = getchar ();

      if (c != 'u') 	/* just copies rest of line */
        {
          if (s != -1)
            printf ("%ld ", s);
	  printf ("%c", c);
	  
          do
            {
              c = getchar ();
	      putchar (c);
            }
          while (c != '\n');
        }
      else /* unknown status */
        {
          int divisible, skip;
	  double ts;

	  ts = GetTime ();
	  skip = 0;
	  tkt++;		/* count trinomials not skipped */
          do
            {
              c = getchar ();	/* discard rest of line */
            }
          while (c != '\n');
          scanf ("\n");

	  GF2X p;
	  SetCoeff(p, r);
	  SetCoeff(p, 0);

	  if (s < 1)
	    {
	    fprintf (stderr, "Error, s < 1\n");
	      exit (1);
	    }
	  if (s >= r)
	    {
	      fprintf (stderr, "Error, s = %ld >= r\n", s);
	      exit (1);
	    }
	    
// We want to flip (i.e. consider reciprocal trinomial x^r + x^(r-s) + 1)
// if s > r/2, because the NTL mod operation is slow in this case.
// However, we also need to flip if s is even, so that we can use 
// the fast squaring algorithm (fastsqr). Thus sometimes we work with
// s > r/2 (but then SqrMod is avoided where possible).  
	    
	  ss = s;		// Save s (in case flipped)  
	  flip = 0;
	  if (2*s > r)
	    {
	    flip = 1;		// Use reciprocal trinomial
	    s = r - s;		// Original s saved in ss
	    }   
	    
	  usefastsqr = 0;
	  if ((r&1) && ((s&1) == 0) && (s > BITS_PER_LONG2) 
	  			    && ((r-s) > BITS_PER_LONG2))
	    {
	    flip = !flip;	// Flip (again if necessary)
	    s = r - s;		// so we can use fastsqr
	    }
	  usefastsqr = ((r&s&1) && (s > BITS_PER_LONG2) 
	  			&& ((r-s) > BITS_PER_LONG2));

	  if (verbose && flip)
	    printf ("Using reciprocal trinomial\n");

	  if (verbose && usefastsqr)
	    printf ("Using fastsqr\n");

	  // If r = +-1 mod 8, 
	  // and s even and s does not divide 2r,
	  // or s odd and (r-s) does not divide 2r,
	  // then x^r + x^s + 1 had an odd number of irreducible factors
	  // by Swan's theorem.

	  swan = (((r&7) == 1) || ((r&7) == 7)) &&
	    	 ((((s&1)==0) && (((2*r)%s) != 0)) || 
	    	  (((s&1)==1) && (((2*r)%(r-s)) != 0)));

	  rhigh = (swan) ? (r/3) : (r/2);

          /* in case the -skip d option is used, and d >= rhigh, stop at d+1 */
          if (skipd >= rhigh)
            {
              rhigh = skipd + 1;
              if (rhigh > r)
                rhigh = r;
            }
          else if (verbose && swan)
	    printf ("By Swan's theorem search at most to degree %ld\n", rhigh);

          if (maxd < rhigh)
            {
              rhigh = maxd;
              if (verbose)
                printf ("Maximal degree checked reduced to %ld\n", rhigh);
            }

	  SetCoeff(p, s, 1);
	  GF2X a, c, htemp, htemp2;
	  GF2XModulus F;
          int Finit = 0; // non zero if F is build
	  clear(htemp);
	  double st;
	  
	  if (choose_m)  			// Try to estimate optimal m
	    {
	    long jt, jtt, mopt;
	    double times, timem;

            if ((Finit == 0) && (!usefastsqr))
              {
                build (F, p);
                Finit = 1;
              }

	    choose_m = 0;			// Only try once
	    SetCoeff (htemp, 0, 1);
	    SetCoeff (htemp, 1, 1);
	    if (usefastsqr)
	      {
	      for (jt = 0; jt < 50; jt++) 	// Warm up cache
	        {				// and make htemp dense
	        fastsqr (htemp2, htemp, r, s);
	        fastsqr (htemp, htemp2, r, s);
	        }
	      times = GetTime ();
	      for (jt = 0; (GetTime() - times) < (double)1.0 ; jt += 100)
	        for (jtt = 0; jtt < 50; jtt++)	// Avoid many GetTime calls
	          {
	          fastsqr (htemp2, htemp, r, s);
	          fastsqr (htemp, htemp2, r, s);
	          }
	      }
	    else
	      {
	      for (jt = 0; jt < 100; jt++) 
	        SqrMod (htemp, htemp, F);	// Warm up cache
	      times = GetTime ();
	      for (jt = 0; (GetTime() - times) < (double)1.0; jt += 100) 
	        for (jtt = 0; jtt < 100; jtt++)
	          SqrMod (htemp, htemp, F);
	      }
	    times = (GetTime () - times)/(double)jt;
	    fastmulmod (htemp, htemp, htemp, F, r, s);
	    timem = GetTime ();
	    for (jt = 0; (GetTime() - timem) < (double)1.0; jt++)
	      fastmulmod (htemp, htemp, htemp, F, r, s);
	    timem = (GetTime () - timem)/(double)jt;
	    // Prefer even m to avoid copying overhead	
	    mopt = 2*(long)((double)0.5*sqrt(timem/times) + (double)0.5);
	    fprintf (stderr, "M %f msec, S %f msec, M/S %f, sqrt(M/S) %f\n", 
	      1000.0* timem, 1000.0*times, timem/times, sqrt(timem/times));
	    fprintf (stderr, "Estimate of best m %ld\n", mopt);
 	    if ((mopt > 1) && (mopt <= 100))
 	      {
 	      m = mopt;
 	      if (verbose) printf ("Using optimal block size of m=%ld\n", m);
 	      }
	    }

	  vec_GF2X h(INIT_SIZE, m);
	  vec_GF2X hold(INIT_SIZE, m); 	/* copy of h for fine DDF */

	  f0 = (double) 1.0;
	  q  = 0;
	  k2 = 1; /* invariant: k2 is last degree checked */
	  fineDDF = 0;  /* 0 as long as GCDs are trivial */
	  for (a = 1; (deg (a) == 0 && k2 < rhigh) || fineDDF != 0;)
	    {
              save_state (ss, k2);

	      if (fineDDF == 0)
		k = k2 + 1; /* next degree to check */
	      else /* fineDDF mode */
		{
		  if (deg (a) > 0) /* smallest factor in current interval */
		    {
		      a = 1; /* reset a to 1 */
		      q2 = q / 2;
		      q = q - q2;
		      copy_h (h, hold, m); /* revert to initial h */
		    }
		  else /* smallest factor in next interval */
		    {
		      k += q * m;
		      q = q2;
		    }
		  k2 = k + q * m - 1;
		}
	      
	      if (k <= k0) /* sieve only one degree at a time */
		{
		  unsigned long K = (1 << k) - 1;
		  
		  if (verbose)
                    printf ("Degree %ld: ", k);

                  if (k <= skipd)
                    {
                      if (verbose)
                        printf ("skip\n");
                    }
                  else
                    {
                      a = 1;
                      SetCoeff (a, K);
                      c = 1;
                      SetCoeff (c, r % K, coeff (c, r % K) == 0);
                      SetCoeff (c, s % K, coeff (c, s % K) == 0);
                      st = GetTime ();
                      GCD (a, a, c);
		      GCD_t += (GetTime () - st);
                      if (verbose)
                        {
                          printf ("gcd took %f sec\n", GetTime () - st);
                          fflush (stdout);
                        }
                    }
		  k2 = k;
		}
	      else /* use double-blocking strategy */
		{
                  if ((Finit == 0) && (!usefastsqr))
                    {
                      build (F, p);
                      Finit = 1;
                    }
		  /* initialize h-table */
		  if (k == k0 + 1)
		    {
		      if (fineDDF == 0)
		        GCD_kt1 += (double) 1.0; // have to do one GCD
		      st = GetTime ();
		      init_h (h, k, m, F, r, s, usefastsqr);
                      Exp_t += (GetTime () - st);
		      if (verbose)
			printf ("Exponentiation took %f\n", GetTime () - st);
		    }
		  
		  if (fineDDF == 0)
		    {
		      if (f > (double) 1.0)  
		 	q = (long) (f0 * (double) q0 + (double) 0.5); // round
		      else
		        q = q + q0;		// linear increase if f = 1
		      f0 = f * f0;
		      k2 = k + q * m - 1;	// q blocks of size m
		    }
		      
		  if (verbose)
                    {
                      printf ("Current q %ld\n", q);
                      fflush (stdout);
                    }

		  /* k2 should not exceed rhigh + m - 1 < r */
		  if (k2 > rhigh) /* should not happen in fine DDF mode */
		    {
		      q = (rhigh - k + m) / m;	// Round up
		      k2 = k + q * m - 1;
		      if (verbose)
			printf ("Reducing q to %ld\n", q);
		    }

		  /* we need to save h whenever q > 1 and the current gcd=1 */
		  if (q > 1 && deg (a) == 0) /* save h for fineDDF */
                    copy_h (hold, h, m);
		    
		  if (verbose)
                    {
                      printf ("Interval %ld..%ld:\n", k, k2);
                      fflush (stdout);
                    }

		  if (q > q1)		/* Output "u" and skip here */
		    {
		    if (verbose)
		      printf ("q > q1 so skipping further test\n");
		    skip = 1;
		    divisible = 0;
		    k2 = r;
		    }
		  else  
		    {
		    st = GetTime ();
		    for (i = k; i + m - 1 <= k2 ; i += m) /* do i .. i+m-1 */
		      {
                    	SqrMul_kt2 += (double)m; // Count m^2 squares as m
		        c = 1;
		        for (j = 0; j < m; j++)
			  {
			    fastmulbyxmod (c, c, r, s);
			    add (c, c, h[j]);
			    if (usefastsqr)
			      {			// Can use fastsqr here
			      for (l = 0; l < m/2; l++)
			        {
			        fastsqr (htemp, h[j], r, s);
			        fastsqr (h[j], htemp, r, s);
			        }
			      if (m&1)		// Last square if m odd
			        {
			        fastsqr (htemp, h[j], r, s);
			        h[j] = htemp;	// Copy result to h[j]
			        }  
			      }  
			    else		// Have to use SqrMod here
			      for (l = 0; l < m; l++)
			        SqrMod (h[j], h[j], F);
			  }
                        if (i + m - 1 <= skipd)
                          {
                            if (verbose >= 2)
                              printf ("   skip degree interval %ld..%ld\n",
                                      i, i + m - 1);
                          }
                        else
                          fastmulmod (a, a, c, F, r, s);
		      }
                    st = GetTime () - st;
                    SqrMul_t += st;

		    if (verbose)
                      {
                        printf ("   squares/products took %f, per term %f\n", 
                                st, st/(double)(q*m));
                        fflush (stdout);
                      }

                    if (k2 <= skipd)
                      {
                        if (verbose)
                          {
                            printf ("   skip gcd\n");
                            fflush (stdout);
                          }
                      }
                    else
                      {
                        st = GetTime ();
                        FastGCD(a, a, p);
			GCD_t += (GetTime () - st);
			GCD_kt2 += (double)1.0;		// Count FastGCD calls
                        if (verbose)
                          {
                            printf ("   gcd took %f\n", GetTime () - st);
			    fflush (stdout);
			  }
			if (deg(a) > 0)
			  {
			    /* check for GCD error */
			    if (deg(a) < k)
			      {
			      printf ("Error: GCD degree %ld < %ld\n", 
			      	deg(a), k);
			      exit (1);
			      }
			    /* exit the loop if we isolated a single factor
			       or q = 1 (we cannot refine any more) 
			       or deg(a) < 10*sqrt(r) so CanZass cheap */
			    if (deg (a) < 2 * k || q == 1 || 
			        deg (a) < fineDDFtol)
			        break;
			    fineDDF = 1; /* enter fineDDF mode */
			    if (verbose)
			      printf ("Entering fine DDF mode\n");
			  }
                      }
		    }
		}
	    }

	  k1 = deg (a);
          /* if k1 = r, the trinomial is irreducible */
          if (rhigh < r)
            divisible = k1 != 0;
          else /* case rhigh = r: either we found a factor of degree < r,
                  or we went up to degree r */
            divisible = k1 < r;
	  if (divisible)
	    {
	      if (k1 >= 2 * k)
		{
		  vec_pair_GF2X_long factors;
		  st = GetTime ();
		  CanZass (factors, a, 0);
		  CanZass_t += GetTime () - st;
		  long i0 = 0;
		  // find lexicographically least factor of smallest degree
		  // (except if flipped)
		  for (i = 1; i < factors.length(); i++)
		    if (lex_less(factors[i].a, factors[i0].a, flip))
		      i0 = i;
		  a = factors[i0].a; // the lex least factor of smallest
		  		     // degree (unless flipped)
		  if (verbose)
		    printf ("CanZass took %f, total degree %ld, small degree %ld\n", 
			    GetTime () - st, k1, deg(a));
		  k1 = deg(a);       // the smallest degree
		  if (k1 < k)
		    {
		    printf ("Error: GCD degree %ld < %ld \n", k1, k);
		    exit (1);
		    }
		}
              if (k1 != 0) /* k1=0 can happen with -skip option */
                print_hex (fact, a, flip);
	    }

          if (!divisible || k1 == r)
            {
              SqrMul_kt1 += (double) (rhigh-k0);	// best possible
              if (skip)
                printf ("%ld u\n", ss);
              else if (primitive)
                printf ("%ld primitive\n", ss);
              else  
	        printf ("%ld irreducible\n", ss);
            }
          else
            {
              if (k1 > k0)
                SqrMul_kt1 += (double) (k1-k0);		// best possible
              if (k1 != 0)
                printf ("%ld %ld", ss, k1);
              else
                printf ("%ld non-irreducible", ss);
	      if (k1 >= b0) 
	        printf (" p%s", fact);
	      printf("\n");  
            }
	  Total_t += (GetTime () - ts);
        }
      fflush (stdout);
    }

  free (fact);

  if (verbose_t && (GCD_kt2 > (double) 0.0) && (Total_t > (double) 0.0))
    {
      printf ("Normal end of run for r = %ld, m = %ld, version %1.2f\n", 
        r, m, VERSION);
      printf ("Total Exp time        %f = %2.2f %%\n", 
	       Exp_t, 100.0*Exp_t/Total_t);
      printf ("Total Sqr/Mul time    %f = %2.2f %%\n", 
	       SqrMul_t, 100.0*SqrMul_t/Total_t);
      printf ("Total GCD time        %f = %2.2f %%\n", 
	       GCD_t, 100.0*GCD_t/Total_t);
      printf ("Total CanZass time    %f = %2.2f %%\n", 
	       CanZass_t, 100.0*CanZass_t/Total_t);
      printf ("Total computation     %f\n", Total_t);
      if (tkt > 0)
        printf ("Time/trinomial        %f for %ld trinomials\n", 
        	 Total_t/tkt, tkt);
      printf ("Efficiency of Sqr/Mul %f\n", SqrMul_kt1/SqrMul_kt2);
      printf ("Efficiency of GCD     %f\n", GCD_kt1/GCD_kt2);
      printf ("Overall efficiency    %f\n",
        (Exp_t + GCD_t*(GCD_kt1/GCD_kt2) + SqrMul_t*(SqrMul_kt1/SqrMul_kt2))/
        (Exp_t + GCD_t + SqrMul_t + CanZass_t));
    }

  return 0;
}

