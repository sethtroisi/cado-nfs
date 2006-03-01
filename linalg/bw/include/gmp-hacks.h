#ifndef GMP_HACKS_H_
#define GMP_HACKS_H_

/* GMP field access macros from gmp-impl.h (which we include anyway IIRC)
 */
#ifndef PTR
#define PTR(x) ((x)->_mp_d)
#endif
#ifndef SIZ
#define SIZ(x) ((x)->_mp_size)
#endif
#ifndef ALLOC
#define ALLOC(x) ((x)->_mp_alloc)
#endif

#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {                                                                  \
    while ((NLIMBS) > 0)                                                \
      {                                                                 \
        if ((DST)[(NLIMBS) - 1] != 0)                                   \
          break;                                                        \
        (NLIMBS)--;                                                     \
      }                                                                 \
  } while (0)
#endif

/* Note : A better construction is provided by gcc : ({ ... }) ; */
	
#define BEGIN_BLOCK	do {
#define END_BLOCK	} while (0)


/* Useful for the lazy boyz */

#define MPZ_GROW_ALLOC(DST, NLIMBS)					\
	BEGIN_BLOCK							\
	if (ALLOC(DST) < (int) (NLIMBS)) {				\
		ALLOC(DST) = (NLIMBS);					\
		PTR(DST)=(mp_limb_t *) realloc(PTR(DST),		\
				(NLIMBS) * sizeof(mp_limb_t));		\
	}								\
	END_BLOCK

#define MPZ_INIT_SET_MPN(DST, SRC, NLIMBS)				\
	BEGIN_BLOCK							\
	ALLOC(DST) = (NLIMBS);						\
	SIZ(DST) = (NLIMBS);						\
	PTR(DST) = (mp_limb_t *)malloc((NLIMBS) * sizeof(mp_limb_t));	\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	
#define MPZ_SET_MPN(DST, SRC, NLIMBS)					\
	BEGIN_BLOCK							\
	MPZ_GROW_ALLOC(DST, NLIMBS);					\
	SIZ(DST) = (NLIMBS);						\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	   
#define MPN_SET_MPZ(DST,NLIMBS,SRC)					\
	BEGIN_BLOCK							\
	memcpy((DST),PTR(SRC),SIZ(SRC) * sizeof(mp_limb_t));		\
	memset((DST)+SIZ(SRC),0,((NLIMBS)-SIZ(SRC)) * sizeof(mp_limb_t));\
	END_BLOCK

/* This one is rather ill-designed, hence the ASSERTS as a caution */
#define MPN_ADD_MPZ(DST,NLIMBS,SRC)					\
	BEGIN_BLOCK							\
	ASSERT(SIZ(SRC) > 0);						\
	ASSERT((NLIMBS) > SIZ(SRC));					\
	mpn_add((DST), (DST), (NLIMBS), PTR(SRC), SIZ(SRC));		\
	END_BLOCK


#endif /* GMP_HACKS_H_ */	
