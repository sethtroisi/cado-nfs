#ifndef ABASE_API_H_
#define ABASE_API_H_

/* This file describes the abase api. abase is here an abstraction layer
 * intended to control the storage of blocks of vectors. Different
 * implementations may store vectors of different block lengths. It
 * should be possible to extend the interface beyond the binary case
 * (earlier versions did).
 *
 * An abase-controlled data is viewed by the user as an almost opaque
 * memory area of type abt *. Displacements in this are must not be
 * computed directly, but rather with the aboffset macro.
 */

/* This file is included by some internal functions. So there's some
 * housekeeping glue which must be discarded when one is only interested
 * in reading the description of functions.
 *
 * Under normal circumstances, with ABASE_BIND undefined, cpp-ing this
 * file gives empty output.
 *
 * (ABASE_BIND could be used as a global ifdef. My original intention was
 * to also allow the possibility of declaring function pointers in a
 * possible virtual base here. However there's no need for it at the
 * moment, and covariant return types would be clumsy at best).
 */

/* types */

/* abobj_t ; this is the type indicating the global information, which is
 * not data dependent. In particular it gives the stride for instance. In
 * all compile-time determined cases, it's goind to be a dummy type
 */
#ifdef  ABASE_BIND
#define abobj_t   ABASE_BIND(obj_t)
#endif

/* abobj_ptr ; on must be able to get a reference on abobj_t by assigning
 * it to an abobj_ptr. This implies that the abobj_t, if not dummy, is
 * best implemented as a [1] array.
 */
/* abobj_srcptr ; same, with const */
#ifdef  ABASE_BIND
#define abobj_ptr   ABASE_BIND(obj_ptr)
#define abobj_srcptr   ABASE_BIND(obj_srcptr)
#endif

/* abt ; this is the main type that abase-handled data areas are made of.
 * Everything which is outside-visible may be declared as abt* or const
 * abt*.
 */
#ifdef  ABASE_BIND
#define abt     ABASE_BIND(base_type)
#endif

/* functions */

/* abobj_init(x) ; initialize the abobj_t variable x */
/* abobj_clear(x) ; clear it */
/* both may be macros */
#ifdef  ABASE_BIND
#define abobj_init(x)   ABASE_BIND(obj_init)(x)
#define abobj_clear(x)  ABASE_BIND(obj_clear)(x)
#endif

/* abobj_init_set(y,x) ; initialize the abobj_t variable y, and set it to
 * have the exact same characteristics as x. Compared to just doing y=x,
 * the difference here is that we have a new independent object, not a
 * reference */
/* may be a macro */
#ifdef  ABASE_BIND
#define abobj_init_set(y,x)   ABASE_BIND(obj_init_set)(y,x)
#endif

/* abobj_set_nbys(x,nbys) ; assert that the corresponding abase will be
 * used to handle nbys values at a time. This function may fail and abort
 * if the required nbys values can not be obtained with the current abase
 * value
 */
#ifdef  ABASE_BIND
#define abobj_set_nbys(x,nbys)   ABASE_BIND(obj_set_nbys)(x,nbys)
#endif


/* abnbits(x) ; the number of bits handled by the corresponding abase.
 * May or may not be a compile-time constant depending on the abase. Note
 * that this corresponds to the number of bits in an abt value, times the
 * repeat count (next macro).
 */
/* may be a macro */
#ifdef  ABASE_BIND
#define abnbits(x)      ABASE_BIND(nbits)(x)
#endif

/* abrepeat(x) ; the number of abt values forming one record. May or may
 * not be a compile-time constant depending on the abase.
 */
/* may be a macro */
#ifdef  ABASE_BIND
#define abrepeat(x)     ABASE_BIND(repeat)(x)
#endif

/* abmax_accumulate(x) ; how many times is it possible to add values
 * while making sure that no overflow occurs. In the binary case, it's
 * defined to infinity.
 */
/* abmax_accumulate_wide(x) ; same, but for values which have been
 * multiplied by a full-length scalar already
 */
#ifdef  ABASE_BIND
#define abmax_accumulate(x)     ABASE_BIND(max_accumulate)(x)
#define abmax_accumulate_wide(x)        ABASE_BIND(max_accumulate_wide)(x)
#endif

/* abinit(x,n) ; this returns a malloc'ed area large enough to hold n
 * records according to the abase x.  */
/* abclear(x,p,n) ; clear it */
/* abinitf(x,n), abclearf(x,p,n) ; same, but allocate on the stack */
#ifdef  ABASE_BIND
#define abinit(x,n)     abase_generic_init(abbytes(x,1),n)
#define abclear(x,p,n)  abase_generic_clear(abbytes(x,1),p,n)
#define abinitf(x,n)    abase_generic_initf(abbytes(x,1),n)
#define abclearf(x,p,n) abase_generic_clearf(abbytes(x,1),p,n)
#endif

/* abzero(x,p,n) ; zero out the n records pointed to by p */
#ifdef  ABASE_BIND
#define abzero(x,p,n)   abase_generic_zero(abbytes(x,1),p,n)
#endif

/* abis_zero(x,p,n) ; test whether the n records pointed to by p are zero */
#ifdef  ABASE_BIND
#define abis_zero(x,p,n)        ABASE_BIND(is_zero)(x,p,n)
#endif

/* abrandom(x,p,n) ; set the n records pointed to by p to random values */
#ifdef  ABASE_BIND
#define abrandom(x,p,n) abase_generic_random(abbytes(x,1),p,n)
#endif

/* abset_ui(x,p,k,v) ; set bit k of the record pointed to by p to the value v */
#ifdef  ABASE_BIND
#define abset_ui(x,p,k,v) ABASE_BIND(set_ui)(x,p,k,v)
#endif

/* aboffset(x,k) ; indicates the value D such that p+D points to the k-th
 * record after the one pointed to by p */
#ifdef  ABASE_BIND
#define aboffset(x,k)   ABASE_BIND(offset)(x,k)
#endif

/* abbytes(x,k) ; the number of bytes occupied by k contiguous records.
 * It's aboffset(x,k) * sizeof(abt)
 */
#ifdef  ABASE_BIND
#define abbytes(x,k)   ABASE_BIND(bytes)(x,k)
#endif

/* abcopy(x,q,p,n) ; copy the n records pointed to by p to the are
 * pointed to by q */
#ifdef  ABASE_BIND
#define abcopy(x,q,p,n) abase_generic_copy(abbytes(x,1),q,p,n)
#endif

/* abadd(x,q,p,n) ; add the record pointed to by p to the one pointed to
 * by q */
#ifdef  ABASE_BIND
#define abadd(x,q,p) ABASE_BIND(add)(x,q,p)
#endif

/* abread(x,f,p,n) ; read n records from a FILE * f to the are pointed to
 * by p */
/* abwrite(x,f,p,n) ; write them */
#ifdef  ABASE_BIND
#define abread(x,f,p,n) abase_generic_read(abbytes(x,1),f,p,n)
#define abwrite(x,f,p,n) abase_generic_write(abbytes(x,1),f,p,n)
#endif

/* abdotprod(x,w,u,v,n) ; computes the dot product of vectors pointed to
 * by u and v. The vector is stored in w as nbits records consisting of
 * abt values.
 *
 * abvdotprod(x,y,w,u,v,n) ; does the same, but with u of type const
 * abvt *, corresponding to the abvobj_t descriptor y.
 * This interface will be revised to accomodate the need to have u and v
 * of different types.
 */
#ifdef  ABASE_BIND
#define abdotprod(x,w,u,v,n) ABASE_BIND(dotprod)(x,w,u,v,n)
#define abvdotprod(x,y,w,u,v,n) ABASE_BIND(vdotprod)(x,y,w,u,v,n)
#endif

/*
 * Things controlled by ABASE_VBIND correspond to something different.
 * It's the readily accessible workalike type which has runtime-specified
 * width. It's dog slow by design, but comes in quite handy for things
 * such as checks and so on.
 */
#ifdef  ABASE_VBIND
#define abvobj_t   ABASE_VBIND(obj_t)
#define abvobj_ptr   ABASE_VBIND(obj_ptr)
#define abvobj_srcptr   ABASE_VBIND(obj_srcptr)
#define abvt     ABASE_VBIND(base_type)
#define abvobj_init(x)   ABASE_VBIND(obj_init)(x)
#define abvobj_clear(x)  ABASE_VBIND(obj_clear)(x)
#define abvobj_init_set(y,x)   ABASE_VBIND(obj_init_set)(y,x)
#define abvobj_set_nbys(x,nbys)   ABASE_VBIND(obj_set_nbys)(x,nbys)
#define abvnbits(x)      ABASE_VBIND(nbits)(x)
#define abvrepeat(x)     ABASE_VBIND(repeat)(x)
#define abvbytes(x,k)   ABASE_VBIND(bytes)(x,k)

#define abvinit(x,n)     abase_generic_init(abvbytes(x,1),n)
#define abvclear(x,p,n)  abase_generic_clear(abvbytes(x,1),p,n)
#define abvinitf(x,n)    abase_generic_initf(abvbytes(x,1),n)
#define abvclearf(x,p,n) abase_generic_clearf(abvbytes(x,1),p,n)
#define abvzero(x,p,n)   abase_generic_zero(abvbytes(x,1),p,n)
#define abvrandom(x,p,n) abase_generic_random(abvbytes(x,1),p,n)
#define abvcopy(x,q,p,n) abase_generic_copy(abvbytes(x,1),q,p,n)
#define abvread(x,f,p,n) abase_generic_read(abvbytes(x,1),f,p,n)

#define abvwrite(x,f,p,n) ABASE_VBIND(write)(x,f,p,n)
#define abvmax_accumulate(x)     ABASE_VBIND(max_accumulate)(x)
#define abvmax_accumulate_wide(x)        ABASE_VBIND(max_accumulate_wide)(x)
#define abvis_zero(x,p,n)        ABASE_VBIND(is_zero)(x,p,n)
#define abvset_ui(x,p,k,v) ABASE_VBIND(set_ui)(x,p,k,v)
#define abvoffset(x,k)   ABASE_VBIND(offset)(x,k)
#define abvadd(x,q,p) ABASE_VBIND(add)(x,q,p)
/* abvdotprod exists already above... */
#endif  /* ABASE_VBIND */

#endif	/* ABASE_API_H_ */

