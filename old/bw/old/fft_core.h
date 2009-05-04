#ifndef FFT_CORE_H_
#define FFT_CORE_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*fft_capable_function) (mp_limb_t *, mp_size_t, mp_limb_t *,
				      unsigned int);
extern fft_capable_function direct_fft, inverse_fft;

extern void prepare_fft_engine(unsigned int, int);
extern void cleanup_fft_engine(void);

/* these ones are mostly private */
extern mp_limb_t *l_roots;
extern int *phi_tab;
extern int tabulated_order;

#ifdef __cplusplus
}
#endif

#endif	/* FFT_CORE_H_ */
