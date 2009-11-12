#include "utils.h"

#if 0
/* old code */
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
#endif

/* new code */

uint32_t cado_crc_lfsr_turn1(cado_crc_lfsr l, uint32_t c)
{
    static uint32_t const twist = 3486328325;

    uint32_t w = 0;

    w = (c >> (31^l->r)) ^ (c << (31&-l->r));
    l->i--;
    l->r+=11;
    l->i &= 31;

    w ^= l->t[ l->i           ];
    w ^= l->t[(l->i + 22) & 31];
    w ^= l->t[(l->i +  2) & 31];
    w ^= l->t[(l->i +  1) & 31];

    w = w >> 1 ^ (twist & -(w&1));
    l->t[l->i] = w;

    return w;
}

void cado_crc_lfsr_init(cado_crc_lfsr l)
{
    l->i = 0;
    l->r = 0;
    for(int k = 0 ; k < 32 ; k++) l->t[k] = k;
    for(int k = 0 ; k < 96 ; k++) {
        cado_crc_lfsr_turn1(l, 0xdeadbeef / (k+1));
    }
}

void cado_crc_lfsr_clear(cado_crc_lfsr l MAYBE_UNUSED) { }

uint32_t cado_crc_lfsr_turn(cado_crc_lfsr l, void * data, unsigned int count)
{
    uint8_t * ptr = (uint8_t *) data;
    uint32_t w = 0;

    for( ; count-- ; ) {
        w = cado_crc_lfsr_turn1(l, *ptr++);
    }
    return w;
}

uint32_t crc32(void * data, size_t count)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = cado_crc_lfsr_turn(l, data, count);
    cado_crc_lfsr_clear(l);
    return w;
}
