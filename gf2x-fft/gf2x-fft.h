#ifndef GF2X_FFT_H_
#define GF2X_FFT_H_

#include "fake_fft.h"
#include "cantor/cantor128.h"
#include "gf2x/gf2x-tfft.h"

#ifdef __cplusplus
#include "fft_adapter.hpp"
DEFINE_FFT_ADAPTER(fake_fft, fake_)
DEFINE_FFT_ADAPTER(c128_fft, c128_)
DEFINE_FFT_ADAPTER(gf2x_tfft, gf2x_tfft_)
#endif

#endif	/* GF2X_FFT_H_ */
