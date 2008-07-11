
#ifndef _STAGE2_H

#define _STAGE2_H

#define NEXT_D 254
#define NEXT_PASS 255

typedef struct {
  unsigned int B2;
  unsigned int d;          /* The d parameter for stage 2 */
  unsigned int s1;
  unsigned int s2;         /* s2|eulerphi(d), we do stage 2 in s2 passes. 
                              Let S1 + S2 == (Z/dZ)* (mod d), |S2| = s2.
                              We write each prime p, B1 < p <= B2, as
                              p = id - j, j = k1 + k2, k1 in S1, k2 in S2.
                              The V_k1(x+1/x) are precomputed and stored.
                              Each pass processes p = id - k1 - k2 for a 
                              fixed k2, for increasing i and with one table 
                              lookup for the k1 value. */
  int *S1,                 /* Contains eulerphi(d)/s2 elements, for each 
                              k1 in S1 we precompute V_k(x+1/x) */
      *S2;                 /* Contains s2 elements, we do one pass for each */
  unsigned char *pairs;    /* For each element in S2, a list of i,j values 
                              for which to multiply V_{id}(x+1/x) - V_{j}(x+1/x)
                              j = k1 + k2 to accumulator */
} stage2_plan_t;

void stage2_make_plan (stage2_plan_t *, unsigned int, unsigned int, int);
void stage2_clear_plan (stage2_plan_t *);

#endif
