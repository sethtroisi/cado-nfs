#ifndef SSE2_H_
#define SSE2_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint64_t v2di __attribute__ ((vector_size (16)));
static inline v2di SHL(v2di x_, unsigned int r_) {
    __asm__("psllq %1,%0" : "+x"(x_) : "K"(r_));
    return x_;
}
static inline v2di SHR(v2di x_, unsigned int r_) {
    __asm__("psrlq %1,%0" : "+x"(x_) : "K"(r_));
    return x_;
}
static inline v2di SHLD(v2di x_, unsigned int r_) {
    __asm__("pslldq %1,%0" : "+x"(x_) : "K"(r_ >> 3));
    return x_;
}
static inline v2di SHRD(v2di x_, unsigned int r_) {
    __asm__("psrldq %1,%0" : "+x"(x_) : "K"(r_ >> 3));
    return x_;
}

#ifdef __cplusplus
}
#endif

#endif	/* SSE2_H_ */
