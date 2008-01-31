#include "globals.hpp"

using namespace globals;

unsigned int m, n;
unsigned int nr, nc;
unsigned int nbys;

field K;

mpz_class modulus;
uint8_t modulus_u8;
uint16_t modulus_u16;
uint32_t modulus_u32;
unsigned long modulus_ulong;

void set_modulus(mpz_class const& n)
{
    modulus = n;
    modulus_u8 = n.get_ui();
    modulus_u16 = n.get_ui();
    modulus_u32 = n.get_ui();
    modulus_ulong = n.get_ui();
}
