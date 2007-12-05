/*
 * Below this line, experimental code
 * (C) 2007 Richard Brent <trinom@rpbrent.com>
 * This code will be released under the GPL 2.0 licence.
 *
 * Based on Sec. 5.2 of Marco Bodrato's paper (reference below)
 * but with full-word aligned operations to reduce overheads.
 *
 * Reference: http://bodrato.it/papers/#WAIFI2007
 *
 * "Towards Optimal Toom-Cook Multiplication for Univariate and
 * Multivariate Polynomials in Characteristic 2 and 0." by Marco
 * BODRATO; in C.Carlet and B.Sunar, editors, "WAIFI'07 proceedings",
 * LNCS 4547, pp. 119-136. Springer, Madrid, Spain, June 21-22, 2007.
 */

#ifndef MUL_TOOMU_THRESHOLD
#define MUL_TOOMU_THRESHOLD 33
#endif

// Need MUL_TOOMU_THRESHOLD >= 11 for internal reasons
// but calls to Toom should have size at least 8 so
// need MUL_TOOMU_THRESHOLD >= 33.

#if (MUL_TOOMU_THRESHOLD < 33)
#error "MUL_TOOMU_THRESHOLD should be at least 33"
#endif

// #define CHECK_ASSERT

#ifdef CHECK_ASSERT  
#define ASSERT assert
#else
#define ASSERT(x)
#endif

/*
  Unbalanced Toom-Cook multiplication, assumes a takes sa words, 
  b takes n = sb = (sa+1)/2 words, 
  returns product c of sa+sb words using five multiplications of
  size (n/2 + O(1)) by (n/2 + O(1)).  See Bodrato, pg. 125, top right.
 
  c should not overlap the inputs.
  
  stk must have space for sp(sa) = toomuspace(sa) words, where 

    sp(sa) = 2*sa + 32 + toomspace(sa/4 + 4)
          >= 4*(2*ceil(n/2) + 3) + toomspace(floor(n/2) + 3)
              
  and toomspace(n) is the maximum space needed for the Toom-Cook routines
  KarMul, Toom3Mul, Toom3wMul, Toom4Mul.
  
  It is assumed that sa >= 33 so n >= 17.
*/


static void Toom3uMul (_ntl_ulong *c, const _ntl_ulong *a, long sa,
              const _ntl_ulong *b, _ntl_ulong *stk)

{
  ASSERT (sa >= MUL_TOOMU_THRESHOLD);	// n should be at least 6 for internal 
                                        // reasons and 17 so calls to Toom
                                        // have size at least 9, so need
                                        // sa >= 33.
                                        
  long n = (sa + 1)/2;			// Assume sb == n == ceil(sa/2)
  long k = (n + 1) / 2; 		// ceil(n/2)
  long d = n&1;				// d = odd(n) = 2k - n, n = 2k - d
  long rb = n - k;			// Size(b1) = k - d
  long ra = sa - 3*k;			// Size(a3) = ra = k - 2d - odd(sa)
  long sc = sa + n;			// Size(c) = 6k - 3d - odd(sa)
   
  const _ntl_ulong *a0 = a;		// Aliases for four parts of a
  const _ntl_ulong *a1 = a0 + k;
  const _ntl_ulong *a2 = a1 + k;
  const _ntl_ulong *a3 = a2 + k;
  const _ntl_ulong *b0 = b;		// Aliases for two parts of b
  const _ntl_ulong *b1 = b0 + k;

  long k2 = 2*(k+3);			// Size of temporary arrays
 
  _ntl_ulong *W0 = c;			// Overlap W0 (size 2*k) with c
  _ntl_ulong *W1 = stk;
  _ntl_ulong *W2 = c + 2*k;		// Overlap W2 with c + 2*k ...
  _ntl_ulong *W3 = W1 + k2;
  _ntl_ulong *W4 = W3 + k2;		// But not W4 as W2 too large
  _ntl_ulong *W5 = W4 + k2;		// W5 is synonymous with W3 in 
  					// Bodrato's paper
  stk += 4*k2;	  			// 4 temporaries of size k2

  _ntl_ulong s, t;
  long j;
  
// In the comments y = x**w where w = wordlength = NTL_BITS_PER_LONG
// y can be thought of as a w-bit shift operator.
// Bodrato's code corresponds to w = 1, which minimises the size of
// the arguments in the Toom calls, but requires a lot of fiddly
// bit-operations.  By choosing w = 32 or 64 we simplify the coding
// and obtain opportunities for loop optimisation. Both methods have
// the same asymptotic complexity.
//
// If the equal-size multiplication is O(n^alpha) then we expect Toom3uMul
// to be worthwhile when alpha > lg(5/2) = 1.3219...
// TC2 has alpha = lg(3)       = 1.58..., 
// TC3 has alpha = lg(5)/lg(3) = 1.46..., 
// TC4 has alpha = lg(7)/2     = 1.40...,
// thus in all these cases Toom3uMul should be worthwhile on average 
// (saving about 5.5% for the case of TC4, and more for other cases).
// However, this analysis does not take O(n) overheads into account
// so it is inaccurate for small n.
// 
// In the comments " + " means addition in GF(2) and " ^ " means
// exponentiation. 

// Evaluation phase				   Size is (max) size in words

// W5 = a3 + a2 + a1 + a0 == A(1)		// Size(W5) := k
// W2 = b1 + b0 == B(1)				// Size(W2) := k

for (j = 0; j < ra; j++)
  {
  W5[j] = a3[j] ^ a2[j] ^ a1[j] ^ a0[j];
  W2[j] = b1[j] ^ b0[j];
  }

for (; j < rb; j++)
  {
  W5[j] = a2[j] ^ a1[j] ^ a0[j];		// No a3[j] here
  W2[j] = b1[j] ^ b0[j];
  }

for (; j < k; j++)
  {
  W5[j] = a2[j] ^ a1[j] ^ a0[j];
  W2[j] = b0[j];				// No b1[j] here
  }

// Calls to Toom mixed with further evaluation. There are 5 calls 
// to Toom with sizes at most k+3-d = n/2 + 3 = (sa+1)/4 + 3.

// W1 = W2 * W5	== C(1)			

  Toom (W1, W2, W5, k, stk);			// Size(W1) := 2*k

// W0 = a3*y^3 + a2*y^2 + a1*y == A(y) - a0	// Size(W0) := k+3-d (at most)

   W0[0] = 0;
   W0[1] = a1[0];
   W0[2] = a2[0] ^ a1[1];
   
   for (j = 0; j < k-3; j++)			// Assumes k > 2, i.e. n > 4
     W0[j+3] = a3[j] ^ a2[j+1] ^ a1[j+2];

   W0[j+3] = a2[j+1] ^ a1[j+2];			// Fix up a3 contribution later
   j++;
   W0[j+3] = a2[j+1];				// No a1[j+2] here, a3 later
   for (j++; j < rb; j++)			// Size(W0) := k+3-d (at most)
     W0[j+3] = 0;				// Need k+3-d <= 2*k
                                                // which is true if n > 4.
   for (j = k-3; j < ra; j++)
     W0[j+3] ^= a3[j];  			// Fix up a3 contribution
     
// W5 += W0 + a3*(y^2 + y)	     		// Size(W5) := k+3-d
// W0 += a0 == A(y)				// Size(W0)  = k+3-d > k

   for (j = k; j < rb+3; j++)			// rb == k-d so rb+3 == k+3-d
     W5[j] = 0;
     
   W5[0] ^= W0[0];
   W0[0] ^= a0[0];   

   t = 0;
   for (j = 1; j < ra; j++)			// Usual case, ra-1 iterations
     {
     _ntl_ulong u;
     u = W0[j];
     W0[j] = u^a0[j];
     s = a3[j-1];
     W5[j] ^= s^(t^u);
     t = s;
     }

  s = a3[j-1];
  W5[j] ^= W0[j]^s^t;
  j++;
  W5[j] ^= s;  

  for (; j < rb+3; j++)
     W5[j] ^= W0[j];
     
   for (j = ra; j < k; j++)
     W0[j] ^= a0[j];     

// Pad W2 to size k+3-d for future Toom call (which requires equal-sized
// inputs).  This is (asymptotically) more efficient that calling AddMul1.
   
   for (j = k; j < rb+3; j++)
     W2[j] = 0;					// Size(W2) := k+3-d

// W2 += b1*y
// W4  = W2 + b1 == B(y)			// Size(W4) := k+3-d
						// but 2 high words zero

   W4[0] = W2[0] ^ b1[0];
   t = b1[0];  
   for (j = 1; j < rb; j++)
     {
     _ntl_ulong s;
     s = W2[j]^t;
     W2[j] = s;
     t = b1[j];
     W4[j] = s^t;
     }
     
   W2[j] ^= t;

   for (; j < rb+3; j++)
     W4[j] = W2[j];         

// W3 = W5 * W2	== C(1+y)			// Size(W3) := 2*(k+3-d)
						// but high 2 words zero
   Toom (W3, W5, W2, rb+3, stk);
   
// W2 = W0 * W4 == C(y)

   Toom (W2, W0, W4, rb+3, stk);		// Size(W2) := 2*(k+3-d)
						// but 2 high words zero
						
// W4 = a3 * b1 == C(infinity)

   if (ra < rb)
     {
     for (j = 0; j < ra; j++)			// W5 := a3 padded to 
       W5[j] = a3[j];   			// size rb > ra
       
     for (; j < rb; j++)
       W5[j] = 0;  				// Size(W5) := rb

     Toom (W4, W5, b1, rb, stk);		// Size(W4) := 2*rb
     }

   else

     Toom (W4, a3, b1, rb, stk);		// Avoid copy if ra == rb

// W0 = a0 * b0 == C(0)

   Toom (W0, a0, b0, k, stk);			// Size(W0) := 2*k;

// Interpolation phase

// W3 += W2 == c1 + c2 + c3*(1 + y + y^2) + c4
// W2 += W0 == C(y) + C(0)

  for (j = 0; j < 2*k; j++)
    {						// First 2*k iterations
    _ntl_ulong s;
    s = W2[j];
    W3[j] ^= s;					// Size(W0) = 2*k
    W2[j]  = s ^ W0[j];				// other sizes 2*rb + 4
    }						// ignoring known zeros				

  for (; j < 2*rb+4; j++)
    W3[j] ^= W2[j];				// Last 4 - 2*d iterations

   ASSERT (W2[0] == 0);				// Division should be exact

// W2 = W2/y + W3

  for (j = 0; j < 2*rb + 3; j++)		
    W2[j] = W2[j+1] ^ W3[j];			
  W2[j] = W3[j];				// Size(W2) := 2*rb + 4
  
// W2 = (W2 + W4*(1+y^3))/(1+y) == c2 + c3  

  for (j = 0, s = 0; j < 3; j++)		
    {						
    s ^= W2[j] ^ W4[j];
    W2[j] = s;					// first 3 iterations special
    }
  for (; j < 2*rb; j++)			
    {
    s ^= W2[j] ^ W4[j] ^ W4[j-3];		// next 2*rb-3 are usual case
    W2[j] = s;
    }
    						
  for (; j < 2*rb+3; j++)
    {
    s ^= W2[j] ^ W4[j-3];			// next 3 are special
    W2[j] = s;
    }

//  W2[j] = 0;					// Size(W2) = 2*rb + 4    
    						// but last word zero
                                                // so Size(W2) := 2*rb + 3
// W1 += W0 == c1 + c2 + c3 + c4
// W3 += W1 == c3*y*(1+y)

  for (j = 0; j < 2*k; j++)
    {
    _ntl_ulong s;
    s = W0[j] ^ W1[j];
    W1[j] = s;					// Size(W0) = Size(W1) = 2*k
    W3[j] ^= s;					// Size(W3) = 2*rb + 4 > 2*k 
    }

  ASSERT (W3[0] == 0);				// Next division exact

// W3 = W3/(y + y^2) == c3

  for (j = 0, s = 0; j < 2*rb + 3; j++)
    {
    s ^= W3[j+1];
    W3[j] = s;					
    }
// W3[j] = 0;	   				
  ASSERT (s == 0);				// Division exact  
			  			// Size(W3) := 2*rb + 2

// W1 += W2 + W4 == c1				// Size(W4) == 2*rb
// W2 += W3      == c2				//  <= Size(W1) == 2*k
						//  <= Size(W3) == 2*rb + 2
						//  <  Size(W2) == 2*rb + 3
  for (j = 0; j < 2*rb; j++)
    {						// Usual case
    _ntl_ulong s;
    s = W2[j];
    W1[j] ^= s ^ W4[j];
    W2[j]  = s ^ W3[j];
    }
    
  for (; j < 2*k; j++)			
    {						// Next 2*d iterations
    _ntl_ulong s;
    s = W2[j];
    W1[j] ^= s;					// No W4[j] here
    W2[j]  = s ^ W3[j];
    }
    
  for (; j < 2*rb + 2; j++)			
    {						// Next 2*(1-d) iterations
    _ntl_ulong s;
    s = W2[j];
    W1[j] = s;					// Extending size of W1
    W2[j] = s ^ W3[j];
    }
  W1[j] = W2[j];				// Size(W1) := 2*rb + 3
    						// Size(W2)  = 2*rb + 3
    
// c = W0 + W1*y + W2*y^2 + W3*y^3 + W4*y^4
// We already have
// W0[j] == c[j] for j = 0 .. 2*k-1 because W0 = c, and
// W2[j] == c[j] for j = 2*k .. 2*k+2*rb+2 because W2 = c + 2*k

  ASSERT (3 - 2*d + 4*k <= sc);

  for (j = 0; j < 3 - 2*d; j++)			// 3 - 2*d words of W2
    c[j+4*k] ^= W4[j];    			// overlap the W4 region

  for (; j < sc - 4*k; j++)			// Copy rest of W4
    c[j+4*k] = W4[j];    			// Here c was undefined
      
  ASSERT (2*rb + 3 + k <= sc);
  
  for (j = 0; j < 2*rb + 3; j++)
    c[j+k]   ^= W1[j];
    
  ASSERT (2*rb + 2 + 3*k <= sc);		// True if n >= 6 so need
                                                // MUL_TOOMU_THRESHOLD >= 6
  for (j = 0; j < 2*rb + 2; j++)
    c[j+3*k] ^= W3[j];
    
}

#undef CHECK_ASSERT
#undef ASSERT
