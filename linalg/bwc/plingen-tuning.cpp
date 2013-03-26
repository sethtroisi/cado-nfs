#include "cado.h"

#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "abase.h"
#include "lingen-polymat.h"
#include "lingen-matpoly.h"
#include "lingen-bigpolymat.h"
#include "plingen.h"
#include "plingen-tuning.h"

#include <vector>
#include <utility>
#include <ostream>
#include <iostream>


using namespace std;

/* {{{ cutoff_array. Should be merged with polymat_cutoff_info, maybe */
/* This is just an array, but compressed as a cut-off list should be */
struct cutoff_array {
    vector<pair<unsigned int, unsigned int> > x;
    cutoff_array(pair<unsigned int, unsigned int> base = make_pair(1,0)) {
        x.push_back(base);
    }
    int push(unsigned int n, unsigned int winner) {
        if (winner == x.back().second) return n - x.back().first;
        x.push_back(make_pair(n, winner));
        return 0;
    }
    /* XXX final_alg is unused for the moment. This is because our
     * interface has in fact very little expressivity. */
    void export_to_cutoff_info(polymat_cutoff_info * dst, unsigned int final_cut, unsigned int final_alg MAYBE_UNUSED) const {
        unsigned int v = 0;
        for( ; v < x.size() && x[v].first < final_cut ; v++);
        dst->table_size = v;
        dst->cut = final_cut;
        dst->subdivide = final_cut;
        dst->table = (unsigned int (*)[2]) realloc(dst->table,
                dst->table_size * sizeof(unsigned int[2]));
        for(v = 0 ; v < dst->table_size ; v++) {
            dst->table[v][0] = x[v].first;
            dst->table[v][1] = x[v].second;
        }
    }
    friend ostream& operator<<(ostream& o, cutoff_array const& A);
};
ostream& operator<<(ostream& o, cutoff_array const& A) {
    o << "{";
    for(auto y : A.x) {
        o << " { " << y.first << ", " << y.second << " },";
    }
    o << " }";
    return o;
}
/* }}} */

struct timer_rusage {
    inline double operator()() const { return seconds(); }
};
struct timer_wct {
    inline double operator()() const { return wct_seconds(); }
};
/* It is important, when doing MPI benches, that we use this algorithm */
struct timer_wct_synchronized {
    inline double operator()() const { 
        double d = wct_seconds();
        MPI_Allreduce(MPI_IN_PLACE, &d, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return d;
    }
};

/* {{{ cutoff_finder and its child small_bench */
template<typename T0 = timer_rusage>
struct cutoff_finder {
    int enough_repeats;
    double minimum_time;
    double enough_time;
    double scale;
    int stable_cutoff_break;
    cutoff_array cutoffs;
    unsigned int win_end;
    int stable;

    cutoff_finder(unsigned int win_start, unsigned int win_end, unsigned int min_length = 1)
        : cutoffs(make_pair(min_length, win_start)), win_end(win_end)
    {
        /* Fill defaults */
        enough_repeats = 1000;
        minimum_time = 0.5;
        enough_time = 1.0;
        scale = 1.1;
        stable_cutoff_break = 4;
        stable = 0;
    }

    inline unsigned int next(unsigned int k) const {
        return MAX(k+1, scale*k);
    }

    inline int done() const {
        return !(stable < stable_cutoff_break);
    }

    void new_winner(unsigned int k, unsigned int w) {
        if (!cutoffs.push(k, w) || w != win_end) {
            stable = 0;
        } else {
            stable++;
        }
    }

    cutoff_array const& result() const { return cutoffs; }
    void export_to_cutoff_info(polymat_cutoff_info* dst, unsigned int final_cut)
    {
        cutoffs.export_to_cutoff_info(dst, final_cut, win_end);
    }
    void export_to_cutoff_info(polymat_cutoff_info* dst)
    {
        ASSERT_ALWAYS(cutoffs.x.back().second == win_end);
        cutoffs.export_to_cutoff_info(dst, cutoffs.x.back().first, win_end);
    }


    template<typename T = T0>
    struct small_bench {
        cutoff_finder const& dad;
        double& tfinal;
        double tt;
        int n;
        explicit small_bench(cutoff_finder const& dad, double& t) :
            dad(dad), tfinal(t) {tt=n=0;}
        inline operator int() {
            if (tt < dad.enough_time) {
                if (tt < dad.minimum_time || n < dad.enough_repeats)
                    return 1;
            }
            tfinal = tt / n;
            return 0;
        }
        inline small_bench& operator++() { n++; return *this; }
        inline void pop() { tt += T()(); }
        inline void push() { tt -= T()(); }
    };

    template<typename T = T0>
    small_bench<T> micro_bench(double& t) const {
        return small_bench<T>(*this, t);
    }

};
/* }}} */

/* This code will first try to bench the basic operations (first
 * local; global are later) of middle product E*pi and product pi*pi.
 * The goal is to see where is the good value for the thresholds:
 * polymat_mul_kara_threshold
 * polymat_mp_kara_threshold
 */

void plingen_tune_mul(abdst_field ab, unsigned int m, unsigned int n)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);
    /* arguments to the ctor:
     * 0 : initially, we know that method 0 wins (basecase)
     * 1 : in the end, we expect method 1 to win (kara)
     * 1 : size makes sense only for >=1
     */
    /* TODO: support multiple values in my cutoff table finder ? */
    cutoff_finder<timer_rusage> finder(0, 1, 1);

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    polymat_cutoff_info_init(always_basecase);
    polymat_cutoff_info_init(improved);

    /* Beware, k is a length, not a degree. Hence length 1 clearly makes
     * no sense */
    for(unsigned int k = 2 ; !finder.done() ; k=finder.next(k)) {
        /* Note: we are benching degree k, but the degree we are
         * interested in for pi is k*m/(m+n) */
        polymat pi, piL, piR;
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, piR, m+n, m+n, k);
        polymat_init(ab, pi, m+n, m+n, 2*k-1);
        polymat_fill_random(ab, piL, k, rstate);
        polymat_fill_random(ab, piR, k, rstate);
        double ttb;
        /* disable kara for a minute */
        polymat_set_mul_kara_cutoff(always_basecase, NULL);
        for(auto x = finder.micro_bench(ttb); x; ++x) {
            x.push();
            polymat_mul(ab, pi, piL, piR);
            x.pop();
        }

        double ttk;
        /* This yields exactly *one* kara recursion step */
        finder.export_to_cutoff_info(improved, k);
        polymat_set_mul_kara_cutoff(improved, NULL);
        for(auto x = finder.micro_bench(ttk); x; ++x) {
            x.push();
            polymat_mul(ab, pi, piL, piR);
            x.pop();
        }

        double ttm;
        /* The matpoly layer is just completetly different -- and gets
         * faster quite early on... */
        for(auto x = finder.micro_bench(ttm); x; ++x) {
            x.push();
            matpoly_mul(ab, (matpoly_ptr) pi, (matpoly_ptr) piL, (matpoly_ptr) piR);
            x.pop();
        }
        printf("%d %1.6f %1.6f %1.6f %1.1f\n", k, ttb, ttk, ttm, ttk/ttb);
        finder.new_winner(k, ttk < ttb); /* < : kara wins: 1 */

        polymat_clear(ab, piL);
        polymat_clear(ab, piR);
        polymat_clear(ab, pi);
    }

    finder.export_to_cutoff_info(improved);
    polymat_set_mul_kara_cutoff(improved, NULL);

    cout << "/* Cutoffs for "<<(m+n)<<"*"<<(m+n)<<"*"<<(m+n)<<" products: */\n";
    cout << "#define MUL_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.result() << endl;
    gmp_randclear(rstate);

    polymat_cutoff_info_clear(always_basecase);
    polymat_cutoff_info_clear(improved);
}/*}}}*/

void plingen_tune_mp(abdst_field ab, unsigned int m, unsigned int n)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);
    /* arguments to the ctor:
     * 0 : initially, we know that method 0 wins (basecase)
     * 1 : in the end, we expect method 1 to win (kara)
     * 1 : size makes sense only for >=1
     */
    cutoff_finder<timer_rusage> finder(0, 1, 1);

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    polymat_cutoff_info_init(always_basecase);
    polymat_cutoff_info_init(improved);

    /* Beware, k is a length, not a degree. Hence length 1 clearly makes
     * no sense */
    for(unsigned int k = 2 ; !finder.done() ; k=finder.next(k)) {
        /* For a 2\ell lingen call, we need a MP which is
         *
         * (2-n/(m+n))\ell times m/(m+n)\ell --> \ell
         *
         * The small argument is m/(m+n)\ell. It's the degree of piL here.
         * The corresponding degree for E, or the interesting fragment of
         * it, is:
         *   (2-n/(m+n))\ell = (2m+n)/(m+n)\ell = (2+n/m)*k
         *
         * while the result gives has length \ell = (1+n/m)*k
         *
         */
        /* Note: we are benching degree k, but the degree we are
         * interested in for ER is k*m/(m+n) */
        polymat ER, E, piL;
        polymat_init(ab, E, m, m+n, 2*k + n*k/m - 1);
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, ER, m, m+n, k + n*k/m);
        polymat_fill_random(ab, E, 2*k + n*k/m - 1, rstate);
        polymat_fill_random(ab, piL, k, rstate);
        double ttb;
        /* disable kara for a minute */
        polymat_set_mp_kara_cutoff(always_basecase, NULL);
        for(auto x = finder.micro_bench(ttb); x; ++x) {
            x.push();
            polymat_mp(ab, ER, E, piL);
            x.pop();
        }

        double ttk;
        /* This yields exactly *one* kara recursion step */
        finder.export_to_cutoff_info(improved, k);
        polymat_set_mp_kara_cutoff(improved, NULL);
        for(auto x = finder.micro_bench(ttk); x; ++x) {
            x.push();
            polymat_mp(ab, ER, E, piL);
            x.pop();
        }
        printf("%d %1.6f %1.6f %1.1f\n", k, ttb, ttk, ttk/ttb);
        finder.new_winner(k, ttk < ttb); /* < : kara wins: 1 */

        polymat_clear(ab, E);
        polymat_clear(ab, piL);
        polymat_clear(ab, ER);
    }

    finder.export_to_cutoff_info(improved);
    polymat_set_mp_kara_cutoff(improved, NULL);

    cout << "/* Cutoffs for "<<(m+n)<<"*"<<(m+n)<<"*"<<(m+n)<<" middle products: */\n";
    cout << "#define MP_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.result() << endl;
    gmp_randclear(rstate);

    polymat_cutoff_info_clear(always_basecase);
    polymat_cutoff_info_clear(improved);
}/*}}}*/

void plingen_tune_bigmul(abdst_field ab, unsigned int m, unsigned int n, unsigned int m1, unsigned int n1, MPI_Comm comm)/*{{{*/
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);
    /* arguments to the ctor:
     * 0 : initially, we know that method 0 wins (basecase)
     * 1 : in the end, we expect method 1 to win (kara)
     * 1 : size makes sense only for >=1
     */
    cutoff_finder<> finder(0, 1, 1);

    bigpolymat model;
    bigpolymat_init_model(model, comm, m1, n1);

    /* Beware, k is a length, not a degree. Hence length 1 clearly makes
     * no sense */
    for(unsigned int k = 2 ; !finder.done() ; k=finder.next(k)) {
        /* Note: we are benching degree k, but the degree we are
         * interested in for pi is k*m/(m+n) */
        polymat pi, piL, piR;
        bigpolymat bpi, bpiL, bpiR;
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, piR, m+n, m+n, k);
        polymat_init(ab, pi, m+n, m+n, 2*k-1);
        bigpolymat_init(ab, bpiL, model, m+n, m+n, k);
        bigpolymat_init(ab, bpiR, model, m+n, m+n, k);
        bigpolymat_init(ab, bpi, model, m+n, m+n, 2*k-1);
        polymat_fill_random(ab, piL, k, rstate);
        polymat_fill_random(ab, piR, k, rstate);
        polymat_fill_random(ab, bigpolymat_my_cell(bpiL), k, rstate);
        polymat_fill_random(ab, bigpolymat_my_cell(bpiR), k, rstate);

        double ttloc;
        if (rank == 0) {
            for(auto x = finder.micro_bench<timer_wct>(ttloc); x; ++x) {
                x.push();
                polymat_mul(ab, pi, piL, piR);
                x.pop();
            }
        }
        MPI_Bcast(&ttloc, 1, MPI_DOUBLE, 0, comm);

        double ttmpi;
        for(auto x = finder.micro_bench<timer_wct_synchronized>(ttmpi); x; ++x) {
            x.push();
            bigpolymat_mul(ab, bpi, bpiL, bpiR);
            x.pop();
        }
        if (rank == 0)
            printf("%d %1.6f %1.6f %1.1f\n", k, ttloc, ttmpi, ttmpi/ttloc);
        finder.new_winner(k, ttmpi < ttloc); /* < : mpi wins: 1 */

        polymat_clear(ab, piL);
        polymat_clear(ab, piR);
        polymat_clear(ab, pi);
        bigpolymat_clear(ab, bpiL);
        bigpolymat_clear(ab, bpiR);
        bigpolymat_clear(ab, bpi);
    }

    bigpolymat_clear_model(model);

    if (rank == 0) {
        cout << "/* Cutoffs for "<<(m+n)<<"*"<<(m+n)<<"*"<<(m+n)<<" MPI products: */\n";
        cout << "#define MUL_MPI_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
            << " " << finder.result() << endl;
    }
    gmp_randclear(rstate);
}/*}}}*/

void plingen_tuning(abdst_field ab, unsigned int m, unsigned int n, MPI_Comm comm, param_list_ptr pl)
{
    mpz_t p;
    int thr[2] = {1,1};
    int mpi[2] = {1,1};
    int rank;

    MPI_Comm_rank(comm, &rank);

    mpz_init(p);
    abfield_characteristic(ab, p);
    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);
    if (rank == 0) {
        printf("# Tuning for m=%d n=%d p=[%zu %u-bits words]"
                " mpi=%dx%d thr=%dx%d\n",
                m, n, mpz_size(p), GMP_LIMB_BITS,
                mpi[0], mpi[1], thr[0], thr[1]);
    }
    mpz_clear(p);

    /* This should normally be reasonably quick, and running it every
     * time can be considered as an option */
    if (rank == 0) {
        plingen_tune_mul(ab, m, n);
        plingen_tune_mp(ab, m, n);
    }
    extern polymat_cutoff_info polymat_mul_kara_cutoff;
    extern polymat_cutoff_info polymat_mp_kara_cutoff;
    bigpolymat_bcast_polymat_cutoff(&polymat_mul_kara_cutoff, 0, comm);
    bigpolymat_bcast_polymat_cutoff(&polymat_mp_kara_cutoff, 0, comm);

    plingen_tune_bigmul(ab, m, n, mpi[0]*thr[0], mpi[1]*thr[1], comm);

#if 0
    int tune_bm_basecase = 1;
    int tune_mp = 1;
    /* record the list of cutoffs, and the method which wins from there
     */

    if (tune_mp) {
        /* Now for benching mp and plain mul */
        unsigned int maxtune = 10000 / (m*n);
        /* Bench the products which would come together with a k-steps
         * basecase algorithm. IOW a 2k, one-level recursive call incurs
         * twice the k-steps basecase, plus once the timings counted here
         * (presently, this has Karatsuba complexity)
         */
        for(unsigned int k = 10 ; k < maxtune ; k+= k/10) {
            polymat E, piL, piR, pi, Er;
            unsigned int sE = k*(m+2*n)/(m+n);
            unsigned int spi = k*m/(m+n);
            polymat_init(E, m, m+n, sE);
            polymat_init(piL, m+n, m+n, spi);
            polymat_init(piR, m+n, m+n, spi);
            polymat_init(pi, m+n, m+n, spi*2);
            polymat_init(Er, m, m+n, sE-spi+1);
            E->size = sE;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            piL->size = spi;
            piR->size = spi;
            for(unsigned int v = 0 ; v < piL->m * piL->n * piL->size ; v++) {
                abrandom(ab, piL->x[v], rstate);
                abrandom(ab, piR->x[v], rstate);
            }
            double ttmp = 0, ttmul = 0;
            ttmp -= seconds();
            polymat_mp(ab, Er, E, piL);
            ttmp += seconds();
            ttmul -= seconds();
            polymat_mul(ab, pi, piL, piR);
            ttmul += seconds();
            double ttmpq = ttmp / (k*k);
            double ttmulq = ttmul / (k*k);
            double ttmpk = ttmp / pow(k, 1.58);
            double ttmulk = ttmul / pow(k, 1.58);
            printf("%u [%.2e+%.2e = %.2e] [%.2e+%.2e = %.2e]\n",
                    k,
                    ttmpq, ttmulq, ttmpq + ttmulq,
                    ttmpk, ttmulk, ttmpk + ttmulk
                    );
            // (seconds()-tt) / (k*k)); // ((sE-spi) * spi) / (m*(m+n)*(m+n)));
            // printf("%zu %.2e\n", E->size, (seconds()-tt) / (k*k)); // (spi * spi) / ((m+n)*(m+n)*(m+n)));
            polymat_clear(E);
            polymat_clear(piL);
            polymat_clear(piR);
            polymat_clear(pi);
            polymat_clear(Er);
        }
    }
    if (tune_bm_basecase) {
        unsigned int maxtune = 10000 / (m * n);
        for(unsigned int k = 10 ; k < maxtune ; k += k/10) {
            unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
            polymat E, pi;
            polymat_init(pi, 0, 0, 0);
            polymat_init(E, m, m+n, maxtune);
            E->size = k;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            double tt = seconds();
            for(unsigned int j = 0 ; j < m + n ; delta[j++]=1);
            bm->t = 1;
            bw_lingen_basecase(bm, pi, E, delta);
            printf("%zu %.2e\n", E->size, (seconds()-tt) / (k * k));
            polymat_clear(pi);
            polymat_clear(E);
            free(delta);
        }
    }


#endif

    return;
}

