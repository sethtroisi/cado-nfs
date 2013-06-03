#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#include <stdint.h>
#include "macros.h"

#define HC_T uint8_t
#define HCM UMAX(HC_T) /* maximal value of TYPE_HASHCOUNT */


typedef struct {
  HT_T p, r;
} ht_t;

typedef struct {
  ht_t *ht; /* p, r */
  HC_T *hc; /* Occurrences of (p,r) */
  HR_T *hr; /* Renumber table. Size = nrels. Type of elements: depend of
	       hm: if >= 2^32, uin64_t, otherwise uint32_t. */
  HR_T  hm; /* Size of ht and hc array. with n = number of relations (in argument),
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
#define HKM(A,B,M) ((HR_T) (HK(A,B)%((unsigned long) (M))))
#define HASHINSERT(H,P,R,F) \
  hashInsertWithKey(H, (HT_T) P, (HT_T) R, (HR_T) HKM(P,R,(H)->hm), (unsigned int *) F)

#ifdef  __cplusplus
extern "C" { 
#endif

  extern void hashClear (hashtable_t *);
  extern void hashInit (hashtable_t *, HR_T, unsigned int);
  extern void hashFree (hashtable_t *);
  extern void hashCheck (hashtable_t *);
  extern HR_T hashInsert (hashtable_t *, HT_T, HT_T, unsigned int *);
  extern HR_T hashInsertWithKey (hashtable_t *, HT_T, HT_T, HR_T, unsigned int *);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
