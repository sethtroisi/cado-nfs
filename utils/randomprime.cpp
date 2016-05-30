#include "cado.h"
#include <stdlib.h>
#include <gmp.h>

#include <vector>
#include <iostream>

#include "utils.h"

/* This binary generates many primes of a given size, using a precomputed
 * table.
 *
 * Note that by design, when asking for (say) 20-bit primes, the set of
 * primes among which we pick does include up to 56 primes above 2^20
 * too.
 */

using namespace std;

/* This stores primes up to any bound below 2^45, using 1.14*pi(x) bytes.
 * (the actual bound for which we care is 42842283925351 ; see
 * http://dx.doi.org/10.1090/S0025-5718-1989-0947470-1 
 *
 * TODO: use simply half-differences for primes below 2^40, that fits.
 *
 */
struct primerange
{
    typedef unsigned long prime_t;
    struct primegroup {
        static const int width = 64-sizeof(prime_t);
        prime_t p;
        uint8_t diff[width];
    };
private:
    vector<primegroup> groups;
public:
    void push_back(primegroup const& G) { groups.push_back(G); }
    prime_t operator[](size_t n) const {
        size_t i0 = n / primegroup::width;
        primegroup const& G(groups[i0]);
        prime_t p = G.p;
        n = n % primegroup::width;
        if (i0 == 0) {
            if (n == 0) return 2;
            if (n == 1) return 3;
            p = 5;
            for(size_t j = 2 ; j < n ; j++) {
                int d = G.diff[j];
                p += 6 * (d/2);
                if (d&1) p += (p % 6) > 3 ? 2 : 4;
            }
        } else {
            for(size_t j = 0 ; j < n ; j++) {
                int d = G.diff[j];
                p += 6 * (d/2);
                if (d&1) p += (p % 6) > 3 ? 2 : 4;
            }
        }
        return p;
    }
    prime_t random(gmp_randstate_t rstate) const {
        size_t i = gmp_urandomm_ui(rstate, groups.size() * primegroup::width);
        return (*this)[i];
    }
};


int main(int argc, char * argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <nbits> [nprimes to generate, defaults to 100]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int nbits = atoi(argv[1]);
    int ngens = 100;
    if (argc >= 3)
        ngens = atoi(argv[2]);

    prime_info pdata;

    prime_info_init(pdata);

    uint64_t p = 2;

    primerange primes;

    for( ; !(p >> nbits) ; ) {
        /* generate a new batch of primes */
        primerange::primegroup G;
        G.p = p;
        for(int i = 0 ; i < G.width ; i++) {
            unsigned long newp = getprime_mt (pdata);
            int d1 = (newp - p) / 6;
            int d0 = ((newp - p) % 6) >> 1;
            G.diff[i] = d1 * 2 + (d0 != 0);
            p = newp;
        }
        primes.push_back(G);
    }
    prime_info_clear(pdata);

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    for(int i = 0 ; i < ngens ; i++) {
        cout << primes.random(rstate) << "\n";
    }

    gmp_randclear(rstate);

}

