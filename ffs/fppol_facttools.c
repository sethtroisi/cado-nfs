#include "fppol.h"
#include "fppol_internal.h"

#ifdef USE_F2
/*
 [ Seqint(ChangeUniverse(Eltseq((PolynomialRing(GF(2))!Intseq(i,2))^2), Integers()), 2) : i in [0..255] ];
*/
static const uint64_t Sq[256] = {
    0, 1, 4, 5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85, 256,
    257, 260, 261, 272, 273, 276, 277, 320, 321, 324, 325, 336, 337,
    340, 341, 1024, 1025, 1028, 1029, 1040, 1041, 1044, 1045, 1088,
    1089, 1092, 1093, 1104, 1105, 1108, 1109, 1280, 1281, 1284, 1285,
    1296, 1297, 1300, 1301, 1344, 1345, 1348, 1349, 1360, 1361, 1364,
    1365, 4096, 4097, 4100, 4101, 4112, 4113, 4116, 4117, 4160, 4161,
    4164, 4165, 4176, 4177, 4180, 4181, 4352, 4353, 4356, 4357, 4368,
    4369, 4372, 4373, 4416, 4417, 4420, 4421, 4432, 4433, 4436, 4437,
    5120, 5121, 5124, 5125, 5136, 5137, 5140, 5141, 5184, 5185, 5188,
    5189, 5200, 5201, 5204, 5205, 5376, 5377, 5380, 5381, 5392, 5393,
    5396, 5397, 5440, 5441, 5444, 5445, 5456, 5457, 5460, 5461,
    16384, 16385, 16388, 16389, 16400, 16401, 16404, 16405, 16448,
    16449, 16452, 16453, 16464, 16465, 16468, 16469, 16640, 16641,
    16644, 16645, 16656, 16657, 16660, 16661, 16704, 16705, 16708,
    16709, 16720, 16721, 16724, 16725, 17408, 17409, 17412, 17413,
    17424, 17425, 17428, 17429, 17472, 17473, 17476, 17477, 17488,
    17489, 17492, 17493, 17664, 17665, 17668, 17669, 17680, 17681,
    17684, 17685, 17728, 17729, 17732, 17733, 17744, 17745, 17748,
    17749, 20480, 20481, 20484, 20485, 20496, 20497, 20500, 20501,
    20544, 20545, 20548, 20549, 20560, 20561, 20564, 20565, 20736,
    20737, 20740, 20741, 20752, 20753, 20756, 20757, 20800, 20801,
    20804, 20805, 20816, 20817, 20820, 20821, 21504, 21505, 21508,
    21509, 21520, 21521, 21524, 21525, 21568, 21569, 21572, 21573,
    21584, 21585, 21588, 21589, 21760, 21761, 21764, 21765, 21776,
    21777, 21780, 21781, 21824, 21825, 21828, 21829, 21840, 21841,
    21844, 21845 };

#if 0
static void fppol16_sqr(fppol32_ptr r, fppol16_srcptr p)
{
    uint64_t rr, pp = p[0];
    rr = Sq[pp & 255];
    rr |= Sq[(pp>>8) & 255]<<16;
    r[0] = rr;
}
 
static void fppol32_sqr(fppol64_ptr r, fppol32_srcptr p)
{
    uint64_t rr, pp = p[0];
    rr = Sq[pp & 255];
    rr |= Sq[(pp>>8) & 255]<<16;
    rr |= Sq[(pp>>16) & 255]<<32;
    rr |= Sq[(pp>>24) & 255]<<48;
    r[0] = rr;
}
#endif

static void fppol64_sqr(fppol64_ptr rh, fppol64_ptr rl, fppol64_srcptr p)
{
    uint64_t rr, pp = p[0];
    rr = Sq[pp & 255];
    rr |= Sq[(pp>>8) & 255]<<16;
    rr |= Sq[(pp>>16) & 255]<<32;
    rr |= Sq[(pp>>24) & 255]<<48;
    rl[0] = rr;
    rr  = Sq[(pp>>32) & 255];
    rr |= Sq[(pp>>40) & 255]<<16;
    rr |= Sq[(pp>>48) & 255]<<32;
    rr |= Sq[(pp>>56) & 255]<<48;
    rh[0] = rr;
}

static void fppol_sqr(fppol_ptr r, fppol_srcptr p)
{
    if (fppol_is_zero(p)) {
        fppol_set_zero(r);
        return;
    }
    fppol64_t rl, rh;
    __fppol_realloc_lazy(r, (1+(p->deg>>6))<<7);
    for (int k = p->deg>>6; k>=0; --k) {
        fppol64_sqr(rh, rl, p->limbs[k]);
        r->limbs[2*k][0] = rl[0];
        r->limbs[2*k+1][0] = rh[0];
    }
    r->deg = 2*p->deg;
}

#endif



void fppol_qpow(fppol_ptr r, fppol_srcptr p) {
#if defined(USE_F2)
    fppol_sqr(r, p);
#elif defined (USE_F3)
#warning "this is slow code for the cube of a polynomial in charac 3"
    fppol_t t;
    fppol_ptr rr;
    if (r == p) {
        fppol_init(t);
        rr = &t[0];
    } else {
        rr = r;
    }
    fppol_set_zero(rr);
    for(int i = 0; i <= fppol_deg(p); ++i) {
        fp_t c;
        fppol_get_coeff(c, p, i);
        fppol_set_coeff(rr, c, 3*i);
    }
    if (r == p) {
        fppol_set(r, rr);
        fppol_clear(t);
    }
#else
#error "Please implement q-power of polynomials for this field!"
#endif
}

void fppol_derivative(fppol_ptr r, fppol_srcptr p) {
#if defined(USE_F2)
    fppol_div_ti(r, p, 1);
    uint64_t mask = 6148914691236517205U; // 1+t^2+t^4+...
    for (int k = 0; k <= r->deg>>6; ++k)
        r->limbs[k][0] &= mask;
    __fppol_update_degree(r);
#elif defined(USE_F3)
    fppol_div_ti(r, p, 1);
    // 1+t + t^3+t^4 + ...
    uint64_t mask1 = 13176245766935394011U;
    // t + t^4 + t^7 +...
    uint64_t mask2 = 2635249153387078802U;
    for (int k = 0; k <= r->deg>>6; ++k) {
        fppol64_t tmp;
        tmp[0] = r->limbs[k][0] & mask2;
        tmp[1] = r->limbs[k][1] & mask2;
        r->limbs[k][0] &= mask1;
        r->limbs[k][1] &= mask1;
        fppol64_add(r->limbs[k], r->limbs[k], tmp);
    }
    __fppol_update_degree(r);
#else
#error "Please implement derivative of polynomials for this field!"

#endif
}


