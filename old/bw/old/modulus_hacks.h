#ifndef MODULUS_HACKS_H_
#define MODULUS_HACKS_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* comment this out to turn off clever modulus hacking */
#define ACTIVATE_HACKS

#if 0
/* XXX Unfortunately, this is *not* true.
 *
 * With tmp={288708538622263, 0}, mod={4294967291},
 * computing tdiv_qr(tmp+1,buf2,0,tmp,2,mod,1) gives a wrong result.
 * (64bits)
 *
 */
/* #define TDIV_QR_OVERLAP_OK */
#endif
#define TDIV_QR_OVERLAP_OK


#ifdef ACTIVATE_HACKS
extern void do_modulus_precomps(void);
#endif /* ACTIVATE_HACKS */
extern void addmul(mp_limb_t *, mp_limb_t *, mp_limb_t *);
extern void inverse(mp_limb_t *, mp_limb_t *);

#ifdef	__cplusplus
}
#endif

#endif /* MODULUS_HACKS_H_ */
