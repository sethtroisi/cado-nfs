/* The norm requires that UINT16_MAX be defined only when this is on. It
 * is used in matrix_repr_binary_sliced.hpp */
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <set>
#include "auxfuncs.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "prep_arguments.hpp"
#include "matmul_toy.hpp"
#include <gmp.h>
#include <gmpxx.h>

#include <iostream>
#include <fstream>

/* bw-prep
 *
 * This program is responsible of the following
 *
 * - Choose X vectors
 * - Choose random Y vectors
 */

namespace globals {
    unsigned int nbys = 1;
    mpz_class modulus;
    uint8_t modulus_u8;
    uint16_t modulus_u16;
    uint32_t modulus_u32;
    unsigned long modulus_ulong;
    unsigned int m,n;
    unsigned int nr;
    std::set<uint32_t> zrows;
    std::set<uint32_t> zcols;
    std::vector<std::vector<uint32_t> > xvecs;
}

using namespace std;

/* {{{ myrand() myseed() */
#ifdef  USE_GMP_RANDOM
static gmp_randstate_t random_state;
static mp_limb_t myrand()
{
    return gmp_urandomb_ui(random_state, GMP_LIMB_BITS);
}
static void myseed(unsigned long int x)
{
    gmp_randseed_ui(random_state, x);
}
#else
static mp_limb_t myrand()
{
    return random();
}
static void myseed(unsigned long int x)
{
    srand(x);
}
#endif
/* }}} */

/* x vectors, if equal to single coordinates, must not have their only one
 * at an index which corresponds to a zero row. Because otherwise we will
 * never have a matrix A(X) of full rank.
 *
 * Again if they equal single coordinates, at least one of them must
 * appear outside of the image of transpose(M). We do not know how to
 * compute this image, of course, but if M has zero columns (it should
 * not) then the coordinate in question must be chosen amongst those --
 * at least one of them, that is.
 *
 * The y vectors must not be im the image of M. More specifically, the
 * intersection of the span of y's and Image(M) must be of dimension n.
 * Again, the image is not know, but if zero rows are known then the
 * projection must have full rank.
 */

void setup_x(unsigned int v = 1)
{
    using namespace globals;
    using namespace std;

    bool zcols_still_to_be_hit = !zcols.empty();
    ofstream o;
    must_open(o, files::x);
    xvecs.clear();
    for (unsigned int i = 0; i < m ; i++) {
        unsigned int idx;
        if (v == 1) {
            for(;;) {
                bool testing_zcol=false;
                if (!zcols.empty()) {
                    idx = * zcols.begin();
                    cout << fmt("// Testing zero col %\n")%idx;
                    zcols.erase(zcols.begin());
                     testing_zcol=true;
                } else {
                    idx= ((myrand()) % nr);
                }
                if (zrows.find(idx) != zrows.end()) {
                    cout << fmt("// Avoiding zero row %\n") % idx;
                    continue;
                }

                if (testing_zcol)
                    zcols_still_to_be_hit=false;
                break;
            }
            o << "e" << idx << "\n";
            xvecs.push_back(vector<uint32_t>(1, idx));
        } else {
            set<uint32_t> vs;
            for (unsigned int j = 0 ; j < v ; j++) {
                if (zcols_still_to_be_hit) {
                    idx = * zcols.begin();
                    cout << fmt("// Testing zero col %\n")%idx;
                    zcols.erase(zcols.begin());
                    zcols_still_to_be_hit = false;
                } else {
                    idx= ((myrand()) % nr);
                }
                vs.insert(idx);
            }
            vector<uint32_t> xs;
            xs.insert(xs.end(), vs.begin(), vs.end());
            for(unsigned int j = 0 ; j < xs.size() ; j++) {
                if (j) o << " ";
                o << "e" << xs[j];
            }
            o << "\n";
            xvecs.push_back(xs);
        }
    }
    if (zcols_still_to_be_hit) {
        die("Could not find starting x vectors."
                " This matrix is too weird\n", 1);
    }
}

unsigned int compute_rank(mpz_class * a, unsigned int nrows, unsigned int ncols, mpz_class& p)
{
    unsigned int rank = 0;
    unsigned int i,j,k,l;

    for(i = 0 ; i < nrows ; i++) {
        /* Find the pivot inside the row */
        mpz_class * row = a + i * ncols;
        for(j = 0 ; j < ncols ; j++) {
            if (row[j] == 0) {
                continue;
            }
            /* make this row more canonical */
            mpz_class w;
            mpz_invert(w.get_mpz_t(), row[j].get_mpz_t(), p.get_mpz_t());
            row[j] = 1;
            for(unsigned int k = j + 1 ; k < ncols ; k++) {
                row[k] = (row[k] * w) % p;
            }
            break;
        }
        if (j == ncols)
            continue;
        ASSERT(rank < nrows && rank < ncols);
        rank++;
        // std::cout << fmt("row % is the %-th pivot\n") % i % rank;
        /* Cancel this coeff in all other rows. */
        for(k = i + 1 ; k < nrows ; k++) {
            mpz_class * rowk = a + k * ncols;
            for(l = 0 ; l < ncols ; l++) {
                if (l == j) continue;
                rowk[l] = (rowk[l] - row[l] * rowk[j]) % p;
            }
            rowk[j] = 0;
        }
    }
    return rank;
}

template<typename traits>
void setup_vectors(simple_matmul_toy<traits>& matmul)
{
    using namespace globals;
    using namespace std;


    typedef typename traits::scalar_t scalar_t;
    typedef typename traits::wide_scalar_t wide_scalar_t;
    typedef scalar_t * vector_bunch;

    scalar_t ystrips[n/nbys][nr];


    /* as ridiculous as it seems, we'll use full blown mpzs for storing
     * coordinates.  the point is that we're far from critical here, and
     * want to be able to do modular reductions simply ! */
    vector<mpz_class> ys(nbys, 0);

    /* How many iterates do we check ? */
#define NBITER  2

    unsigned int nx = 1;
    for(unsigned ntri = 0 ; ; ntri++) {
        cout << "// Generating new x,y vector pair\n";
        if (ntri >= nx * 10) {
            cout << "// Getting bored. Trying " << ++nx << " x vectors\n";
        }
        setup_x(nx);

        /* We'll try several times until we're happy. We're going to
         * compute the first iterates online.  */
        mpz_class a1a2[m][NBITER*n];

        for (unsigned int j = 0 ; j < n ; j += nbys) {
            /* Generate this nbys y vectors at a time, as it's going to
             * help us with the online precomputation of the first
             * iterates.  */

            typedef set<unsigned int>::const_iterator suci_t;
            suci_t sc = zrows.begin();
            int kz = 0;

            for(unsigned int c = 0; c < nr; ) {
                unsigned int next_c;
                if (sc == zrows.end()) {
                    next_c = nr;
                } else {
                    next_c = *sc;
                    sc++;
                }

                /* aside zero rows, there's no constraint */
                for( ; c < next_c ; c++) {
                    for(unsigned int y = 0 ; y < nbys ; y++) {
                        ys[y] = myrand() % globals::modulus;
                    }
                    traits::assign(ystrips[j/nbys][c], ys, 0);
                }
                if (c < nr) {
                    /* This row is special, it's zero. We force the
                     * intersection of the y block with the zero rows to be a
                     * set of I_n matrices. */
                    unsigned int y = (kz++ % n) - j;
                    ys.assign(nbys,0);
                    if (y < nbys) {
                        ys[y] = 1;
                        /* otherwise the 1 is for another vertical strip */
                    }
                    traits::assign(ystrips[j/nbys][c], ys, 0);
                }
            }

            /* OK so now ystrips[j/nbys] holds a candidate y vector. Use
             * it to populate some entries of the a1a2 matrix.
             */

            traits::copy(matmul.v, ystrips[j/nbys], nr);

            wide_scalar_t wsum;

            for(int l = 0 ; l < NBITER ; l++) {
                matmul.multiply();
                BUG_ON(xvecs.size() != m);
                for(unsigned int i = 0 ; i < m ; i++) {
                    scalar_t sum;
                    traits::zero(wsum);
                    for(unsigned int k = 0 ; k < xvecs[i].size() ; k++) {
                        traits::addmul(wsum, matmul.w[xvecs[i][k]], 1);
                    }
                    traits::reduce(sum, wsum);

                    traits::assign(ys, sum);
                    for(unsigned int y = 0 ; y < nbys ; y++) {
                        a1a2[i][j+y+l*n] = ys[y];
                    }
                }
                swap(matmul.v, matmul.w);
            }
        }

        /* At this point we should check the rank of a1a2 */
#if 0
        cout << fmt("amat:=Matrix(%,%,[")%m%(NBITER*n);
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < NBITER * n ; j++) {
                if (i||j) cout << ", ";
                cout << a1a2[i][j];
            }
        }
        cout << "]);\n";
        cout << flush;
#endif

        unsigned int r;
        r = compute_rank(& a1a2[0][0], m, NBITER * n, globals::modulus);

        cout << "// Rank: " << r << endl;

        if (r == m)
            break;
    }

    ofstream o;
    must_open(o, files::y);
    for(unsigned int c = 0; c < nr;c++ ) {
        for (unsigned int j = 0 ; j < n/nbys ; j++) {
            if (j) o << " ";
            traits::print(o, ystrips[j][c]);
        }
        o << "\n";
    }
    o.close();
}

int main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(false);
    cerr.tie(&cout);

    std::cout << "// This is bw-prep, version " VERSION << std::endl;

    common_arguments common;
    prep_arguments mine;

    process_arguments(argc, argv, common, mine);

#ifdef  USE_GMP_RANDOM
    gmp_randinit_mt(random_state);
#endif
    if (mine.seed == 0) {
        myseed(time(NULL));
    } else {
        myseed(mine.seed);
    }

    matrix_stats stats;
    stats.need_zrows(&globals::zrows);
    stats.need_zcols(&globals::zcols);

    stats(files::matrix);

    BUG_ON(stats.nr!= stats.nc);

    globals::m = mine.m;
    globals::n = mine.n;
    globals::nr= stats.nr;
    globals::modulus = stats.modulus;
    globals::modulus_u8 = globals::modulus.get_ui();
    globals::modulus_u16 = globals::modulus.get_ui();
    globals::modulus_u32 = globals::modulus.get_ui();
    globals::modulus_ulong = globals::modulus.get_ui();

    {
        ofstream cfg;
        cfg.open(files::params.c_str());
        cfg << fmt("n=%\n") %  globals::n;
        cfg << fmt("m=%\n") %  globals::m;
    }

    if (globals::modulus == 2) {
        if (mine.n % 128 == 0) {
            typedef binary_sse2_traits T;
            globals::nbys = 128;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        } else if (mine.n % 64 == 0) {
            typedef binary_pod_traits<uint64_t> T;
            globals::nbys = 64;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        } else if (mine.n % 32 == 0) {
            typedef binary_pod_traits<uint32_t> T;
            globals::nbys = 32;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        } else if (mine.n % 16 == 0) {
            typedef binary_pod_traits<uint16_t> T;
            globals::nbys = 16;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        } else if (mine.n % 8 == 0) {
            typedef binary_pod_traits<uint8_t> T;
            globals::nbys = 8;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        } else {
            typedef variable_scalar_traits T; 
            globals::nbys = 1;
            simple_matmul_toy<T> mul_code(stats, files::matrix);
            setup_vectors(mul_code);
        }
    } else {
        typedef variable_scalar_traits T; 
        globals::nbys = 1;
        simple_matmul_toy<T> mul_code(stats, files::matrix);
        setup_vectors(mul_code);
    }
}


/* vim:set sw=4 sta et: */
