#include "ops_poly.h"

#if (OPS_POLY == OPS_POLY_FFT)
#include "ops_poly_fft.c"
#elif (OPS_POLY == OPS_POLY_SCBK)
#include "ops_poly_schoolbook.c"
#elif (OPS_POLY == OPS_POLY_IFFT)
#include "ops_poly_integerfft.c"
#elif (OPS_POLY == OPS_POLY_WRAP)
#include "ops_poly_wrap.c"
#endif
