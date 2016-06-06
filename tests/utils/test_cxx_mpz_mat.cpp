#include "cado.h"
#include "macros.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <map>
#include <algorithm>
#include "mpz_mat.h"
#include "tests_common.h"

using namespace std;

int main(int argc, char * argv[])
{
    if (argc > 1) {
        srand(atoi(argv[1]));
    }

    map<unsigned long, cxx_mpz_mat> v;

    mpz_t det, p;
    mpz_init(det);
    mpz_init_set_ui(p, 1009);

    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_mat M;

        mpz_mat_realloc(M, rand() % 16 + 2, rand() % 16 + 2);
        for(unsigned int i = 0 ; i < M->m ; i++) {
            for(unsigned int j = 0 ; j < M->n ; j++) {
                mpz_set_si(mpz_mat_entry(M, i, j), (rand() - (RAND_MAX / 2)));
            }
        }
        cxx_mpz_mat M1 = M, T, M2;
        mpz_gauss_backend_mod(M1, T, p);
        mpz_mat_fprint(stdout, M1);
        unsigned int d = std::min(M1->m, M1->n);
        mpz_mat_realloc(M2, d, d);
        mpz_mat_submat_swap(M2, 0, 0, M1, 0, 0, d, d);
        mpz_mat_determinant_triangular(det, M2);
        v.insert(make_pair(mpz_get_ui(det) % mpz_get_ui(p), M));
    }

    for(decltype(v)::const_iterator it = v.begin() ; it != v.end() ; it++) {
        printf("[det=%ld] ", it->first);
        mpz_mat_fprint(stdout, it->second);
    }

    mpz_clear(det);
    mpz_clear(p);
}
