#ifndef LAS_NORMS_H_
#define LAS_NORMS_H_

#include <stdint.h>
#include "las-types.h"
#include "double_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
void init_norms (sieve_info_ptr si, int side);

/* Initialize lognorms for the bucket_region number J. It's a wrapper.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime.
 * The sieve area S must be preallocated with at least (BUCKET_REGION +
 * MEMSET_MIN) space. Indeed, the algorithm that is used might write a
 * bit beyond the meaningful area.
 */
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info_ptr si, unsigned int side, unsigned int smart);

double get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y);

/* This prepares the auxiliary data which is used by
 * init_rat_norms_bucket_region and init_alg_norms_bucket_region
 */
void sieve_info_init_norm_data(sieve_info_ptr si);

void sieve_info_clear_norm_data(sieve_info_ptr si);

void sieve_info_update_norm_data (sieve_info_ptr, int);

int sieve_info_adjust_IJ(sieve_info_ptr si, int nb_threads);

void sieve_info_init_norm_data_sq (sieve_info_ptr si, unsigned long q);

/* To use this LAS_MEMSET, you have to:
 *  1. Define LAS_MEMSET in sieve/las-config.h and have a x86 64 bits.
 *  2. To optimize LAS_MEMSET you have to include this code in your
 *  executable:

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET)
  max_cache = direct_write_vs_stos ();
  min_stos = stos_vs_write128 ();
  fprintf (stderr, "# Las normalisation memset: movaps from 33(0x21) to %zu(0x%zx); rep stosq until %zu(0x%zx); movntps after\n", min_stos, min_stos, max_cache, max_cache);
#endif

*/
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET)
#include <emmintrin.h>

/* Only to avoid a possible warning (in MacOSX). */
#ifdef memset
#undef memset
#endif
#define memset las_memset

/* strategy used for memset of size n:
 *
 * for n < min_stos : movaps
 * for min_stos <= n < max_cache : rep stosq
 * for n >= max_cache : movntps
 */
  static size_t min_stos = 0x8000; /* 32 KB = max size for Intel L0 cache */
  static size_t max_cache = 0x800000; /* 8 MB = average max size for biggest cache */

  /* The fastest memset for x86 64 & SSE2 */
  static inline void *las_memset(void *S, int c, size_t n) {
    uint64_t rc = 0x0101010101010101 * (uint8_t) c;
    if (LIKELY (n > 0x20)) {
      register __m128 mc __asm__ ("xmm7"); /* Way to ask a "legacy" xmm, from xmm0 to xmm7 ? */
      __asm__ __volatile__ ( "movd %[rc], %[mc]\n" : [mc]"=x"(mc) : [rc]"r"((uint32_t) rc));
      void *cS = S;
      *(int64_t *) ((uintptr_t) S + n - 0x10) = rc;
      *(int64_t *) ((uintptr_t) S + n - 0x08) = rc;
      *(int64_t *) ((uintptr_t) S           ) = rc;
      *(int64_t *) ((uintptr_t) S     + 0x08) = rc;
      __asm__ __volatile__ ( "pshufd $0x00, %[mc], %[mc]\n" : [mc]"+x"(mc));
      if (LIKELY (n < min_stos)) {
	uint64_t jmp, offset;
	n += (uintptr_t) S - 0x10;
	n |= 0x0f;
	S = (void *) ((uintptr_t) S | 0x0f);
	n -= (uintptr_t) S;
	offset = (uint8_t) n;
	S = (void *) ((uintptr_t) S + offset - 0x7f);
	offset >>= 2;
	__asm__ __volatile__ ( "leaq 0f(%%rip), %[jmp]\n"
			       "subq %[offset], %[jmp]\n"
			       "shrq $0x08, %[n]\n"
			       "jmpq *%[jmp]\n"
			       
			       ".p2align 4\n 1:\n"
			       "addq $0x100, %[S]\n"
			       "subq $0x01, %[n]\n"
			       "movaps %[mc], -0x80(%[S])\n"
			       "movaps %[mc], -0x70(%[S])\n"
			       "movaps %[mc], -0x60(%[S])\n"
			       "movaps %[mc], -0x50(%[S])\n"
			       "movaps %[mc], -0x40(%[S])\n"
			       "movaps %[mc], -0x30(%[S])\n"
			       "movaps %[mc], -0x20(%[S])\n"
			       "movaps %[mc], -0x10(%[S])\n"
			       "movdqa %[mc],      (%[S])\n"
			       "movaps %[mc],  0x10(%[S])\n"
			       "movaps %[mc],  0x20(%[S])\n"
			       "movaps %[mc],  0x30(%[S])\n"
			       "movaps %[mc],  0x40(%[S])\n"
			       "movaps %[mc],  0x50(%[S])\n"
			       "movaps %[mc],  0x60(%[S])\n"
			       "movaps %[mc],  0x70(%[S])\n"
			       "0:\n"
			       "jnz 1b\n"
			       : [jmp]"=&r"(jmp), [n]"+r"(n), [S]"+R"(S) : [offset]"r"(offset), [mc]"x"(mc) : "cc", "memory");
      } else if (LIKELY (n < max_cache)) {
	n += (uintptr_t) S;
	n &= -0x10;
	_mm_store_ps ((float *) (uintptr_t) (n - 0x40), mc);
	_mm_store_ps ((float *) (uintptr_t) (n - 0x30), mc);
	_mm_store_ps ((float *) (uintptr_t) (n - 0x20), mc);
	_mm_store_ps ((float *) (uintptr_t) (n - 0x10), mc);
	S += 0x40;
	S = (void *)((uintptr_t) S & -0x10);
	_mm_store_ps ((float *) (uintptr_t) (S - 0x30), mc);
	_mm_store_ps ((float *) (uintptr_t) (S - 0x20), mc);
	_mm_store_ps ((float *) (uintptr_t) (S - 0x10), mc);
	_mm_store_ps ((float *) (uintptr_t) (S       ), mc);
	S = (void *)((uintptr_t) S & -0x40);
	n &= -0x40;
	n -= (uintptr_t) S;
	n >>= 3;
	__asm__ __volatile__ ( "cld\n rep\n stosq\n" : "+c"(n), [S]"+D"(S) : [rc]"a"(rc) : "cc", "memory");
      } else {
	uint64_t jmp, offset;
	n += (uintptr_t) S - 0x10;
	n |= 0x0f;
	S = (void *) ((uintptr_t) S | 0x0f);
	n -= (uintptr_t) S;
	offset = (uint8_t) n;
	S = (void *) ((uintptr_t) S + offset - 0x7f);
	offset = (uint64_t) offset >> 2;
	__asm__ __volatile__ ( "leaq 0f(%%rip), %[jmp]\n"
			       "subq %[offset], %[jmp]\n"
			       "shrq $0x08, %[n]\n"
			       "jmpq *%[jmp]\n"
			       
			       ".p2align 4\n 1:\n"
			       "addq $0x100, %[S]\n"
			       "subq $0x01,%[n]\n"
			       "movntps %[mc], -0x80(%[S])\n"
			       "movntps %[mc], -0x70(%[S])\n"
			       "movntps %[mc], -0x60(%[S])\n"
			       "movntps %[mc], -0x50(%[S])\n"
			       "movntps %[mc], -0x40(%[S])\n"
			       "movntps %[mc], -0x30(%[S])\n"
			       "movntps %[mc], -0x20(%[S])\n"
			       "movntps %[mc], -0x10(%[S])\n"
			       "movntdq %[mc],      (%[S])\n"
			       "movntps %[mc],  0x10(%[S])\n"
			       "movntps %[mc],  0x20(%[S])\n"
			       "movntps %[mc],  0x30(%[S])\n"
			       "movntps %[mc],  0x40(%[S])\n"
			       "movntps %[mc],  0x50(%[S])\n"
			       "movntps %[mc],  0x60(%[S])\n"
			       "movntps %[mc],  0x70(%[S])\n"
			       "0:\n"
			       "jnz 1b\n"
			       : [jmp]"=&r"(jmp), [n]"+r"(n), [S]"+R"(S) : [offset]"r"(offset), [mc]"x"(mc) : "cc", "memory");
      }
      return cS;
    } else if (UNLIKELY ((uint8_t) n & 0x30)) {
      *(int64_t *) ((uint8_t *) S     + 0x08) = rc;
      *(int64_t *) ((uint8_t *) S + n - 0x10) = rc;
    between_0x08_0x10:
      *(int64_t *) ((uint8_t *) S           ) = rc;
      *(int64_t *) ((uint8_t *) S + n - 0x08) = rc;
    } else if (UNLIKELY ((uint8_t) n & 0x08)) goto between_0x08_0x10;
    else if (UNLIKELY ((uint8_t) n & 0x04)) {
      *((uint32_t *) ((uint8_t *) S        )) = (uint32_t) rc;
      *((uint32_t *) ((uint8_t *) S + n - 0x04)) = (uint32_t) rc;
    } else if (LIKELY ((uint8_t) n >= 0x01)) {
      *((uint8_t  *) ((uint8_t *) S        )) = (uint8_t ) rc;
      if (LIKELY ((uint8_t) n > 0x01)) {
	*((uint16_t *) ((uint8_t *) S + n - 0x02)) = (uint16_t) rc;
      }
    }
    return S;
  }
  
  extern void tune_las_memset();
#endif
  
/* These functions are internals. Don't use them. Use the wrapper above.
   It's need to declare them here for units & coverage tests.
 */
void init_degree_one_norms_bucket_region_internal     (unsigned char *S, uint32_t J, uint32_t I, double scale, double u0, double u1, double *cexp2);
void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd);
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd, unsigned int nroots, root_ptr roots);
void init_norms_roots_internal (unsigned int degree, double *coeff, double max_abs_root, double precision, unsigned int *nroots, root_ptr roots);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_NORMS_H_ */
