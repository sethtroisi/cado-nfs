#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#define	MODULUS_BITS	19
#define	BITS_TO_WORDS(B,W)	(((B)+(W)-1)/(W))

#ifdef	__x86_64
#define MODULUS_SIZE	BITS_TO_WORDS(MODULUS_BITS, 64)
#else
#define MODULUS_SIZE    BITS_TO_WORDS(MODULUS_BITS, 32)
#endif

#define	CHECKPOINTS	10

#endif	/* CONSTANTS_HPP_ */
