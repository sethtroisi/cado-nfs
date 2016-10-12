#include <iostream>
#include <gmp.h>

int main(int argc, char * argv[])
{
    mpz_t a;
    mpz_init_set_ui(a, 0);
    std::cout << a << std::endl;
    mpz_clear(a);
}
