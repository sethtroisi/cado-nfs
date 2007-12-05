/*
 * Below this line, experimental code
 * (C) 2007 Richard Brent <trinom@rpbrent.com>
 * This code will be released under the GPL 2.0 licence.
 *
 * Based on Marco Bodrato's mul-tc3.c but with full-word aligned 
 * operations to reduce overheads.
 *
 * Reference: http://bodrato.it/papers/#WAIFI2007
 *
 * "Towards Optimal Toom-Cook Multiplication for Univariate and
 * Multivariate Polynomials in Characteristic 2 and 0." by Marco
 * BODRATO; in C.Carlet and B.Sunar, editors, "WAIFI'07 proceedings",
 * LNCS 4547, pp. 119-136. Springer, Madrid, Spain, June 21-22, 2007.
 */

#ifndef MUL_TOOMW_THRESHOLD
#define MUL_TOOMW_THRESHOLD 10
#endif

#if (MUL_TOOMW_THRESHOLD < 8)
#error "MUL_TOOMW_THRESHOLD should be at least 8"
#endif

// #define CHECK_ASSERT

#ifdef CHECK_ASSERT  
#define ASSERT assert
#else
#define ASSERT(x)
#endif

/*
  c must have space for 2n words and should not overlap the inputs.
  stk must have space for sp(n) = toomspace(n) words

    sp(n) = (n lt 8) ? KarMem(7) : 8*(n/3 + 3) + sp(n/3 + 2)
    
  and  KarMem(7) = 19 is the space required by KarMul.

  A simpler bound on the memory required is 5*n + 17 (equality at n = 19).
*/


static void Toom3WMul (_ntl_ulong *c, const _ntl_ulong *a, const _ntl_ulong *b,
              long n, _ntl_ulong *stk)

{
  long k = (n + 2) / 3; 		// size of a0, a1, b0, b1
  long r = n - 2*k;			// size of a2, b2
  long d = (r < k) ? 1 : 0;		// 1 if r < k, 0 otherwise
  long kd = k - d;
   
  const _ntl_ulong *a0 = a;		// Aliases for three parts of a
  const _ntl_ulong *a1 = a + k;
  const _ntl_ulong *a2 = a + 2*k;
  const _ntl_ulong *b0 = b;		// Ditto for b
  const _ntl_ulong *b1 = b + k;
  const _ntl_ulong *b2 = b + 2*k;

  long k2 = 2*(k+2);			// Size of temporary arrays
 
  _ntl_ulong *W0 = c;			// Overlap W0 (size 2*k) with c
  _ntl_ulong *W1 = stk;
  _ntl_ulong *W2 = c + 2*k;		// Overlap W2 with c + 2*k ...
  _ntl_ulong *W3 = W1 + k2;
  _ntl_ulong *W4 = W3 + k2;		// But not W4 as W2 too large
  _ntl_ulong *W5 = W4 + k2;		// W5 is synonymous with W3 in 
  					// Bodrato's mul-tc3.c
  stk += 4*k2;	  			// 4 temporaries of size k2

  long j;
  _ntl_ulong s, u2, v2;
  
// In the comments y = x**w where w = wordlength = NTL_BITS_PER_LONG
// y can be thought of as a w-bit shift operator.
// Bodrato's code corresponds to w = 1, which minimises the size of
// the arguments in the recursive calls, but requires a lot of fiddly
// bit-operations.  By choosing w = 32 or 64 we simplify the coding
// and obtain opportunities for loop optimisation. Both methods have
// the same asymptotic complexity O(n**(ln(5)/ln(3))) = O(n**1.464).

// We try to combine loops as far as possible to reduce overheads and memory
// references. This often means splitting a loop into the "usual" case and
// "special" cases at the start or end, due to different size arrays etc.

// In the comments " + " means addition in GF(2) and " ^ " means
// exponentiation. 

// Evaluation phase				   Size is (max) size in words

// W0 = a1*y + a2*y^2 == A(y) - a0 == A(1+y) - A(1)
// W4 = b1*y + b2*y^2 == B(y) - b0 == B(1+y) - B(1)
// W5 = a0 + a1 + a2 == A(1)
// W2 = b0 + b1 + b2 == B(1)

  W0[0] = W4[0] = 0;
  W0[1] = a1[0];
  W4[1] = b1[0];				// No a2, b2 here
  W5[0] = a0[0] ^ a1[0] ^ (u2 = a2[0]);
  W2[0] = b0[0] ^ b1[0] ^ (v2 = b2[0]);
  for (j = 1; j < r; j++)			// Next r-1 iterations
    {						// This is the usual case
    _ntl_ulong u1, v1;
    W0[j+1] = (u1 = a1[j]) ^ u2;		// Size(a1) = Size(b1) = k
    W4[j+1] = (v1 = b1[j]) ^ v2;
    W5[j] = a0[j] ^ u1 ^ (u2 = a2[j]);		// Size(a2) = Size(b2) = r
    W2[j] = b0[j] ^ v1 ^ (v2 = b2[j]);
    }
  for (; j < k; j++)				// Last iterations for W5, W2
    {	
    W0[j+1] = a1[j];				// Omit a2, b2 here
    W4[j+1] = b1[j];
    W5[j] = a0[j] ^ a1[j];			// Size(W5) := k
    W2[j] = b0[j] ^ b1[j]; 			// Size(W2) := k;
    }
  W0[k+1] = W4[k+1] = 0;			// In case r == k  
  W0[r+1] ^= a2[r-1];				// Size(W0) := kd+2
  W4[r+1] ^= b2[r-1];				// Size(W4) := kd+2

// Recursive calls mixed with further evaluation
// There are 5 recursive calls with sizes at most k+2.
// Thus it is necessary that n > 4 (but we assume that 
// Karatsuba's method or some other method will be used
// for very small n, say n < MUL_TOOMW_THRESHOLD).

// W1 = W2 * W5	== C(1)			

  Toom (W1, W2, W5, k, stk);			// Size(W1) := 2*k

// W5 += W0 == A(1+y)				// Size(W5) < Size(W0)
// W2 += W4 == B(1+y)				// Size(W2) < Size(W4)
// W0 += a0 == A(y)				// Size(W0) > Size(a0)
// W4 += b0 == B(y)				// Size(W4) > Size(b0)

  for (j = 0; j < k; j++)			// First k iterations
    {
    _ntl_ulong u, v;
    W5[j] ^= (u = W0[j]);
    W2[j] ^= (v = W4[j]);
    W0[j]  = u ^ a0[j];
    W4[j]  = v ^ b0[j];
    }
    
  for (; j < kd+2; j++)				// Last 2-d iterations
    {
    W5[j] = W0[j];				// Size(W5) := kd+2
    W2[j] = W4[j];				// Size(W2) := kd+2
    }
    
// W3 = W2 * W5 == C(1+y)

// Output argument in recursive call must differ from inputs.
// That is why we need both W3 and W5.

// ASSERT ((kd+2) <= (n/3 + 2));		// Explaining the space bound
  
  Toom (W3, W2, W5, kd+2, stk);			// Size(W3) := 2*kd + 4
	
// W2 = W0 * W4 == C(y)

  Toom (W2, W0, W4, kd+2, stk);			// Size(W2) := 2*kd + 4

// W0 = a0 * b0 == c0 == C(0/1) == C(0)
  
  Toom (W0, a0, b0, k, stk);			// Size(W0) := 2*k
						// so c[0..(2k-1)] defined
// W4 = a2 * b2 == c4 == C(1/0) == C(infinity)

  Toom (W4, a2, b2, r, stk);			// Size(W4) := 2*r

// Interpolation phase

// W3 += W2 == c1 + c2 + c3*(1 + y + y^2) + c4
// W2 += W0 == C(y) + C(0)

  for (j = 0; j < 2*k; j++)
    {						// First 2*k iterations
    s = W2[j];
    W3[j] ^= s;					// Size(W0) = 2*k
    W2[j]  = s ^ W0[j];				// other sizes 2*kd + 4
    }				

  for (; j < 2*kd+4; j++)
    W3[j] ^= W2[j];				// Last 4 - 2*d iterations

// ASSERT (W2[0] == 0);				// Division should be exact

// W2 = W2/y + W3

  for (j = 0; j < 2*kd + 3; j++)		
    W2[j] = W2[j+1] ^ W3[j];			
  W2[j] = W3[j];				// Size(W2) := 2*kd + 4
  
// W2 = (W2 + W4*(1+y^3))/(1+y) == c2 + c3  

  for (j = 0, s = 0; j < 3; j++)		
    {						
    s ^= W2[j] ^ W4[j];
    W2[j] = s;					// first 3 iterations special
    }
  for (; j < 2*r; j++)
    {
    s ^= W2[j] ^ W4[j] ^ W4[j-3];		// next 2r-3 are usual case
    W2[j] = s;
    }
    						
  for (; j < 2*r+3; j++)
    {
    s ^= W2[j] ^ W4[j-3];			// next 3 are special
    W2[j] = s;
    }

  for (; j < 2*kd+4; j++)
    {						
    s ^= W2[j]; 				// last (k-r-d) == 0 or 1
    W2[j] = s;					// Size(W2) = 2*kd + 3    
    }

// ASSERT (s == 0);				// Division should be exact
  
// W1 += W0 == c1 + c2 + c3 + c4
// W3 += W1 == c3*y*(1+y)

  for (long j = 0; j < 2*k; j++)
    {
    s = W0[j] ^ W1[j];
    W1[j] = s;					// Size(W0) = Size(W1) = 2*k
    W3[j] ^= s;					// Size(W3) = 2*kd + 4 > 2*k 
    }

// ASSERT (W3[0] == 0);				// Next division exact

// W3 = W3/(y + y^2) == c3

  for (j = 0, s = 0; j < 2*kd + 3; j++)
    {
    s ^= W3[j+1];
    W3[j] = s;					
    }
  W3[j] = 0;	   				
// ASSERT (s == 0);				// Division exact  
			  			// Size(W3) := 2*kd + 2

// W1 += W2 + W4 == c1				// Size(W4) == 2*r
// W2 += W3      == c2				//  <= Size(W1) == 2*k
						//  <= Size(W3) == 2*kd + 2
						//  <  Size(W2) == 2*kd + 4
  for (j = 0; j < 2*r; j++)
    {						// Usual case
    s = W2[j];
    W1[j] ^= s ^ W4[j];
    W2[j]  = s ^ W3[j];
    }
  for (; j < 2*k; j++)			
    {						// Next 2*(k-r) iterations
    s = W2[j];
    W1[j] ^= s;					// No W4[j] here
    W2[j]  = s ^ W3[j];
    }
  for (; j < 2*kd + 2; j++)			
    {						// Next 2*(1-d) iterations
    s = W2[j];
    W1[j] = s;					// Extending size of W1
    W2[j] = s ^ W3[j];
    }
  for (; j < 2*kd + 4; j++)			// Last 2 iterations
    W1[j] = W2[j];				// Size(W1) := 2*kd + 4
    						// Size(W2)  = 2*kd + 4
    
// c = W0 + W1*y + W2*y^2 + W3*y^3 + W4*y^4
// We already have
// W0[j] == c[j] for j = 0 .. 2*k-1 because W0 = c, and
// W2[j] == c[j] for j = 2*k .. 2*k+2*kd+3 because W2 = c + 2*k

  for (j = 0; j < 4 - 2*d; j++)			// 4 - 2*d words of W2
    c[j+4*k] ^= W4[j];    			// overlap the W4 region

  for (; j < 2*r; j++)				// Copy rest of W4	
    c[j+4*k] = W4[j];    			// Here c was undefined
      
  for (long j = 0; j < 2*kd + 4; j++)
    c[j+k]   ^= W1[j];
    
// ASSERT (2*kd + 2 + 3*k <= 2*n);		// True if n >= 8 so need
						// MUL_TOOMW_THRESHOLD >= 8  
  for (long j = 0; j < 2*kd + 2; j++)
    c[j+3*k] ^= W3[j];
    
}

#undef CHECK_ASSERT
#undef ASSERT

