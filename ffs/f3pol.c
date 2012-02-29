#include "f3pol.h"

// Look-up table for converting 5 elements of GF(3) to 8 bits.
//
// Given 5 elements of GF(3) in bitsliced representation { p[0], p[1] }, with
// |p[0]| = |p[1]| = 5 bits, and where each element p_i is represented as
// p_i = p[0]_i + 2*p[1]_i, the table is addressed by ((p[1] << 5) | p[0]).
//
// Generated using Magma:
/*
   &cat[ "  " cat &cat L cat "\n" : L in Partition(
     [ &or[ s0[i] eq 1 and s1[i] eq 1 : i in [1..Min(#s0, #s1)] ]
       select "  X," else "   "[1..3-#s] cat s cat ","
       where s  := IntegerToString(x, 10)
       where x  := Seqint(ChangeUniverse(Eltseq(R!s0+2*R!s1), Integers()), 3)
       where s0 := Intseq(p0, 2)
       where s1 := Intseq(p1, 2)
       where R  := PolynomialRing(GF(3))
     : p0 in [0..2^5-1], p1 in [0..2^5-1] ], 16) ];
*/
#define X 255
const uint8_t  __f3_get_ui_conv[] = {
    0,  1,  3,  4,  9, 10, 12, 13, 27, 28, 30, 31, 36, 37, 39, 40,
   81, 82, 84, 85, 90, 91, 93, 94,108,109,111,112,117,118,120,121,
    2,  X,  5,  X, 11,  X, 14,  X, 29,  X, 32,  X, 38,  X, 41,  X,
   83,  X, 86,  X, 92,  X, 95,  X,110,  X,113,  X,119,  X,122,  X,
    6,  7,  X,  X, 15, 16,  X,  X, 33, 34,  X,  X, 42, 43,  X,  X,
   87, 88,  X,  X, 96, 97,  X,  X,114,115,  X,  X,123,124,  X,  X,
    8,  X,  X,  X, 17,  X,  X,  X, 35,  X,  X,  X, 44,  X,  X,  X,
   89,  X,  X,  X, 98,  X,  X,  X,116,  X,  X,  X,125,  X,  X,  X,
   18, 19, 21, 22,  X,  X,  X,  X, 45, 46, 48, 49,  X,  X,  X,  X,
   99,100,102,103,  X,  X,  X,  X,126,127,129,130,  X,  X,  X,  X,
   20,  X, 23,  X,  X,  X,  X,  X, 47,  X, 50,  X,  X,  X,  X,  X,
  101,  X,104,  X,  X,  X,  X,  X,128,  X,131,  X,  X,  X,  X,  X,
   24, 25,  X,  X,  X,  X,  X,  X, 51, 52,  X,  X,  X,  X,  X,  X,
  105,106,  X,  X,  X,  X,  X,  X,132,133,  X,  X,  X,  X,  X,  X,
   26,  X,  X,  X,  X,  X,  X,  X, 53,  X,  X,  X,  X,  X,  X,  X,
  107,  X,  X,  X,  X,  X,  X,  X,134,  X,  X,  X,  X,  X,  X,  X,
   54, 55, 57, 58, 63, 64, 66, 67,  X,  X,  X,  X,  X,  X,  X,  X,
  135,136,138,139,144,145,147,148,  X,  X,  X,  X,  X,  X,  X,  X,
   56,  X, 59,  X, 65,  X, 68,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  137,  X,140,  X,146,  X,149,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   60, 61,  X,  X, 69, 70,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  141,142,  X,  X,150,151,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   62,  X,  X,  X, 71,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  143,  X,  X,  X,152,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   72, 73, 75, 76,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  153,154,156,157,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   74,  X, 77,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  155,  X,158,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   78, 79,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  159,160,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   80,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  161,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  162,163,165,166,171,172,174,175,189,190,192,193,198,199,201,202,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  164,  X,167,  X,173,  X,176,  X,191,  X,194,  X,200,  X,203,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  168,169,  X,  X,177,178,  X,  X,195,196,  X,  X,204,205,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  170,  X,  X,  X,179,  X,  X,  X,197,  X,  X,  X,206,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  180,181,183,184,  X,  X,  X,  X,207,208,210,211,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  182,  X,185,  X,  X,  X,  X,  X,209,  X,212,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  186,187,  X,  X,  X,  X,  X,  X,213,214,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  188,  X,  X,  X,  X,  X,  X,  X,215,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  216,217,219,220,225,226,228,229,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  218,  X,221,  X,227,  X,230,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  222,223,  X,  X,231,232,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  224,  X,  X,  X,233,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  234,235,237,238,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  236,  X,239,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  240,241,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  242,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
};
#undef X


// Look-up table for converting 8 bits to 5 elements of GF(3).
//
// Addressed by x with 0 <= x < 3^5 = 243, the table retuns 10 bits
// ((p[1] << 5) | p[0]), with |p[0]| = |p[1]| = 5 bits, corresponding
// to the bitsliced representation { p[0], p[1] } of 5 elements of GF(3)
// p_i = p[0]_i + 2*p[1]_i.
//
// Generated using Magma:
/*
   &cat[ " " cat &cat L cat "\n" : L in Partition(
     [ " 0x000"[1..6-#s] cat s cat ","
       where s  := IntegerToString(p0 + p1*2^5, 16)
       where p0 := Seqint([ p[i] eq 1 select 1 else 0 : i in [1..#p] ], 2)
       where p1 := Seqint([ p[i] eq 2 select 1 else 0 : i in [1..#p] ], 2)
       where p  := Intseq(x, 3)
     : x in [0..3^5-1] ], 9) ];
*/
const uint16_t __f3_set_ui_conv[] = {
  0x000, 0x001, 0x020, 0x002, 0x003, 0x022, 0x040, 0x041, 0x060,
  0x004, 0x005, 0x024, 0x006, 0x007, 0x026, 0x044, 0x045, 0x064,
  0x080, 0x081, 0x0a0, 0x082, 0x083, 0x0a2, 0x0c0, 0x0c1, 0x0e0,
  0x008, 0x009, 0x028, 0x00a, 0x00b, 0x02a, 0x048, 0x049, 0x068,
  0x00c, 0x00d, 0x02c, 0x00e, 0x00f, 0x02e, 0x04c, 0x04d, 0x06c,
  0x088, 0x089, 0x0a8, 0x08a, 0x08b, 0x0aa, 0x0c8, 0x0c9, 0x0e8,
  0x100, 0x101, 0x120, 0x102, 0x103, 0x122, 0x140, 0x141, 0x160,
  0x104, 0x105, 0x124, 0x106, 0x107, 0x126, 0x144, 0x145, 0x164,
  0x180, 0x181, 0x1a0, 0x182, 0x183, 0x1a2, 0x1c0, 0x1c1, 0x1e0,
  0x010, 0x011, 0x030, 0x012, 0x013, 0x032, 0x050, 0x051, 0x070,
  0x014, 0x015, 0x034, 0x016, 0x017, 0x036, 0x054, 0x055, 0x074,
  0x090, 0x091, 0x0b0, 0x092, 0x093, 0x0b2, 0x0d0, 0x0d1, 0x0f0,
  0x018, 0x019, 0x038, 0x01a, 0x01b, 0x03a, 0x058, 0x059, 0x078,
  0x01c, 0x01d, 0x03c, 0x01e, 0x01f, 0x03e, 0x05c, 0x05d, 0x07c,
  0x098, 0x099, 0x0b8, 0x09a, 0x09b, 0x0ba, 0x0d8, 0x0d9, 0x0f8,
  0x110, 0x111, 0x130, 0x112, 0x113, 0x132, 0x150, 0x151, 0x170,
  0x114, 0x115, 0x134, 0x116, 0x117, 0x136, 0x154, 0x155, 0x174,
  0x190, 0x191, 0x1b0, 0x192, 0x193, 0x1b2, 0x1d0, 0x1d1, 0x1f0,
  0x200, 0x201, 0x220, 0x202, 0x203, 0x222, 0x240, 0x241, 0x260,
  0x204, 0x205, 0x224, 0x206, 0x207, 0x226, 0x244, 0x245, 0x264,
  0x280, 0x281, 0x2a0, 0x282, 0x283, 0x2a2, 0x2c0, 0x2c1, 0x2e0,
  0x208, 0x209, 0x228, 0x20a, 0x20b, 0x22a, 0x248, 0x249, 0x268,
  0x20c, 0x20d, 0x22c, 0x20e, 0x20f, 0x22e, 0x24c, 0x24d, 0x26c,
  0x288, 0x289, 0x2a8, 0x28a, 0x28b, 0x2aa, 0x2c8, 0x2c9, 0x2e8,
  0x300, 0x301, 0x320, 0x302, 0x303, 0x322, 0x340, 0x341, 0x360,
  0x304, 0x305, 0x324, 0x306, 0x307, 0x326, 0x344, 0x345, 0x364,
  0x380, 0x381, 0x3a0, 0x382, 0x383, 0x3a2, 0x3c0, 0x3c1, 0x3e0,
};


// Look-up table for converting 5 elements of GF(3), where the most significant
// one (if any) is 1, to 8 bits.
// Addressed just as __f3_get_ui_conv.
//
// Generated using Magma:
/*
   &cat[ "  " cat &cat L cat "\n" : L in Partition(
     [ (#s0 le #s1 and #s1 gt 0) or 
       &or[ s0[i] eq 1 and s1[i] eq 1 : i in [1..#s1] ]
       select "  X," else "   "[1..3-#s] cat s cat ","
       where s  := IntegerToString(x - (3^(Max(#s0, 1)-1)-1) div 2, 10)
       where x  := Seqint(ChangeUniverse(Eltseq(R!s0+2*R!s1), Integers()), 3)
       where s0 := Intseq(p0, 2)
       where s1 := Intseq(p1, 2)
       where R  := PolynomialRing(GF(3))
     : p0 in [0..2^5-1], p1 in [0..2^4-1] ], 16) ];
*/
#define X 255
const uint8_t  __f3_monic_get_ui_conv[] = {
    0,  1,  2,  3,  5,  6,  8,  9, 14, 15, 17, 18, 23, 24, 26, 27,
   41, 42, 44, 45, 50, 51, 53, 54, 68, 69, 71, 72, 77, 78, 80, 81,
    X,  X,  4,  X,  7,  X, 10,  X, 16,  X, 19,  X, 25,  X, 28,  X,
   43,  X, 46,  X, 52,  X, 55,  X, 70,  X, 73,  X, 79,  X, 82,  X,
    X,  X,  X,  X, 11, 12,  X,  X, 20, 21,  X,  X, 29, 30,  X,  X,
   47, 48,  X,  X, 56, 57,  X,  X, 74, 75,  X,  X, 83, 84,  X,  X,
    X,  X,  X,  X, 13,  X,  X,  X, 22,  X,  X,  X, 31,  X,  X,  X,
   49,  X,  X,  X, 58,  X,  X,  X, 76,  X,  X,  X, 85,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X, 32, 33, 35, 36,  X,  X,  X,  X,
   59, 60, 62, 63,  X,  X,  X,  X, 86, 87, 89, 90,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X, 34,  X, 37,  X,  X,  X,  X,  X,
   61,  X, 64,  X,  X,  X,  X,  X, 88,  X, 91,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X, 38, 39,  X,  X,  X,  X,  X,  X,
   65, 66,  X,  X,  X,  X,  X,  X, 92, 93,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X, 40,  X,  X,  X,  X,  X,  X,  X,
   67,  X,  X,  X,  X,  X,  X,  X, 94,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   95, 96, 98, 99,104,105,107,108,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
   97,  X,100,  X,106,  X,109,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  101,102,  X,  X,110,111,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  103,  X,  X,  X,112,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  113,114,116,117,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  115,  X,118,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  119,120,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
  121,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
};
#undef X


// Look-up table for converting 8 bits to 5 elements of GF(3), where the most
// significant one (if any) is 1.
// Addressed just as __f3_set_ui_conv.
//
// Generated using Magma:
/*
   &cat[ " " cat &cat L cat "\n" : L in Partition(
     [ " 0x000"[1..6-#s] cat s cat ","
       where s  := IntegerToString(p0 + p1*2^5, 16)
       where p0 := Seqint([ p[i] eq 1 select 1 else 0 : i in [1..#p] ], 2)
       where p1 := Seqint([ p[i] eq 2 select 1 else 0 : i in [1..#p] ], 2)
       where p  := Intseq(x + n, 3)
       where n  := Max([0] cat [ n : i in [1..4]
                               | n lt x where n := (3^i-1) div 2 ])
     : x in [0..(3^5-1) div 2] ], [1, 1, 3] cat [9 : i in [1..13]]) ];
*/
const uint16_t __f3_monic_set_ui_conv[] = {
  0x000,
  0x001,
  0x002, 0x003, 0x022,
  0x004, 0x005, 0x024, 0x006, 0x007, 0x026, 0x044, 0x045, 0x064,
  0x008, 0x009, 0x028, 0x00a, 0x00b, 0x02a, 0x048, 0x049, 0x068,
  0x00c, 0x00d, 0x02c, 0x00e, 0x00f, 0x02e, 0x04c, 0x04d, 0x06c,
  0x088, 0x089, 0x0a8, 0x08a, 0x08b, 0x0aa, 0x0c8, 0x0c9, 0x0e8,
  0x010, 0x011, 0x030, 0x012, 0x013, 0x032, 0x050, 0x051, 0x070,
  0x014, 0x015, 0x034, 0x016, 0x017, 0x036, 0x054, 0x055, 0x074,
  0x090, 0x091, 0x0b0, 0x092, 0x093, 0x0b2, 0x0d0, 0x0d1, 0x0f0,
  0x018, 0x019, 0x038, 0x01a, 0x01b, 0x03a, 0x058, 0x059, 0x078,
  0x01c, 0x01d, 0x03c, 0x01e, 0x01f, 0x03e, 0x05c, 0x05d, 0x07c,
  0x098, 0x099, 0x0b8, 0x09a, 0x09b, 0x0ba, 0x0d8, 0x0d9, 0x0f8,
  0x110, 0x111, 0x130, 0x112, 0x113, 0x132, 0x150, 0x151, 0x170,
  0x114, 0x115, 0x134, 0x116, 0x117, 0x136, 0x154, 0x155, 0x174,
  0x190, 0x191, 0x1b0, 0x192, 0x193, 0x1b2, 0x1d0, 0x1d1, 0x1f0,
};
