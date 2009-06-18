#include "utils.h"

/*  crc code */
#define RED(l,h) do {                                                   \
        /* Compute crc mod x^32 + x^7 + x^6 + x^2 + 1 */                \
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;                              \
        h  = h >> 30 ^ h >> 26 ^ h >> 25;                               \
        /* h is at most 7 bits now. */                                  \
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;                              \
        h = 0;                                                          \
} while (0)

uint32_t crc32(unsigned long * c, int n)
{
    uint32_t v = 0UL;
    for(int i = 0 ; i < n ; ) {
        unsigned long cj = *c++;
        { uint32_t h = v, l = cj; RED(l,h); i++; v = l; if (i == n) break; }
        if (ULONG_BITS == 64) {
            /* wait, there's more ! */
            cj >>= 32;
            { uint32_t h = v, l = cj; RED(l,h); i++; v = l; if (i == n) break; }
        }
    }
    return v;
}
