#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#include <stdint.h>
#include "macros.h"
#include "typedefs.h"

#define HC_T uint8_t
#define HCM UMAX(HC_T) /* maximal value of TYPE_HASHCOUNT */


typedef struct {
  index_t p, r;
} ht_t;

typedef struct {
  ht_t *ht; /* p, r */
  HC_T *hc; /* Occurrences of (p,r) */
  p_r_values_t *hr; /* Renumber table. Size = nrels. Type of elements: depend of
	       hm: if >= 2^32, uin64_t, otherwise uint32_t. */
  p_r_values_t  hm; /* Size of ht and hc array. with n = number of relations (in argument),
	       hm = ulong_nextprime ((unsigned long) ((n>>1)*3)) */
} hashtable_t;

#define GET_HASH_P(H,h) (H)->ht[h].p
#define SET_HASH_P(H,h,v) GET_HASH_P(H,h)=(v)
#define GET_HASH_R(H,h) (H)->ht[h].r
#define SET_HASH_R(H,h,v) GET_HASH_R(H,h)=(v)

/* The constants C0 = floor(Pi*10^17) and C1 = floor(exp(1)*10^17) are given
   in Cavallar's thesis. They guarantee that if a, b have at most 53 bits,
   then H(a,b) = C0*a+C1*b is injective, i.e., has no collision. */

#define HC0_64 314159265358979323UL
#define HC1_64 271828182845904523UL
#define HC0_32 3141592653UL
#define HC1_32 2718281828UL
#define HC0 ((sizeof(unsigned long)>4) ? (HC0_64) : (HC0_32))
#define HC1 ((sizeof(unsigned long)>4) ? (HC1_64) : (HC1_32))
#define HK(A,B) (((unsigned long)(A))*HC0+((unsigned long)(B))*HC1)
#define HKM(A,B,M) ((p_r_values_t) (HK(A,B)%((unsigned long) (M))))
#define HASHINSERT(H,P,R,F) \
  hashInsertWithKey(H, (index_t) P, (index_t) R, (p_r_values_t) HKM(P,R,(H)->hm), (unsigned int *) F)

#ifdef  __cplusplus
extern "C" { 
#endif

  extern void hashClear (hashtable_t *);
  extern void hashInit (hashtable_t *, p_r_values_t, unsigned int);
  extern void hashFree (hashtable_t *);
  extern void hashCheck (hashtable_t *);
  extern p_r_values_t hashInsert (hashtable_t *, index_t, index_t, unsigned int *);
  extern p_r_values_t hashInsertWithKey (hashtable_t *, index_t, index_t, p_r_values_t, unsigned int *);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
