#include "fake_fft.h"
#include "gf2x.h"
#include "manu.h"
#include "utils.h"

void fake_setup(fake_info_t p, int dF, int dG)
{
    p->d1 = dF;
    p->d2 = dG;
    p->d3 = dF + dG;
    int nc = -1;
    if (p->d1 > nc) nc = p->d1;
    if (p->d2 > nc) nc = p->d2;
    nc += 1;
    p->size = 2 * BITS_TO_WORDS(nc, ULONG_BITS);
}

void fake_dft(const fake_info_t p MAYBE_UNUSED, fake_t dst, unsigned long * src, int n) {
    ASSERT(n <= p->d1 || n <= p->d2);
    int s = BITS_TO_WORDS(n + 1, ULONG_BITS);
    memcpy(dst, src, s * sizeof(ulong));
}
void fake_ift(const fake_info_t p MAYBE_UNUSED, unsigned long * dst, int n, fake_src_t src) {
    ASSERT(n <= p->d3);
    int t = BITS_TO_WORDS(n + 1, ULONG_BITS);
    memcpy(dst, src, t * sizeof(ulong));
}
void fake_compose(const fake_info_t p MAYBE_UNUSED, fake_t dst, fake_src_t s1, fake_src_t s2) {
    int n1 = BITS_TO_WORDS(p->d1 + 1, ULONG_BITS);
    int n2 = BITS_TO_WORDS(p->d2 + 1, ULONG_BITS);
    mul_gf2x(dst, s1, n1, s2, n2);
}
void fake_add(const fake_info_t p MAYBE_UNUSED, fake_t dst, fake_src_t s1, fake_src_t s2) {
    unsigned int i;
    for(i = 0 ; i < p->size ; i++) {
        dst[i] = s1[i] ^ s2[i];
    }
}

/* vim: set sw=4 sta et: */
