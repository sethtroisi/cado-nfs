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
#include <list>
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <assert.h>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>
#include <iomanip>

#include "cado.h"
#undef  ASSERT
#undef  ASSERT_ALWAYS
#include "utils.h"

#define Lmacro(N, m, n) (iceildiv((N)+2*(n),(m))+iceildiv((N)+2*(n),(n))+10)

unsigned int rec_threshold = 0;

#include "auxfuncs.h"   /* for die, coredump_limit */
#include "constants.hpp"
#include "manu.h"
#include "fmt.hpp"
#include "fft_adapter.hpp"
#include "lingen_mat_types.hpp"
#include "files.hpp"

/* output: F0x is the candidate number x. All coordinates of the
 * candidate are grouped, and coefficients come in order (least
 * significant first), separated by a carriage return.
 */

// applies a permutation on the source indices/*{{{*/
template<typename T>
void permute(std::vector<T>& a, unsigned int p[])
{
    std::vector<T> b(a.size());
    for(unsigned int i = 0 ; i < a.size() ; i++) {
        b[p[i]] = a[i];
    }
    a.swap(b);
}/*}}}*/


/* Needed operations on critical path.
 *
 * compute_ctaf: compute the coefficient of degree t of the product of E
 * (m*b matrix) by PI (a b*b matrix), to be stored into small-e (m*b
 * matrix)
 *
 * zeroing out a column of small-e (m*b)
 *
 * do a column gaussian elimination on small-e:
 *   - column after column,
 *     - find the first non-zero bit in the column.
 *     - add this column to other column.
 *
 * Then polynomial operations.
 *
 * => small-e is better stored column-wise.
 * => big-E as well, and PI also.
 *
 * Extra un-critical stuff includes
 * - initialization -> computation of F0
 *                  -> computation of initial E = A * F0
 * - termination -> computation of final F = F0 * PI
 *
 * (while uncritical, the computation of the initial E and the final F
 * may become a problem for debugging situations when the block size is
 * largish).
 *
 * all this stuff has linear complexity.
 */

    /*{{{*/ /* debug code for polmat */
void dbmat(bmat const * x)
{
    using namespace std;
    for(unsigned int i = 0 ; i < x->nrows ; i++) {
        for(unsigned int j = 0 ; j < x->ncols ; j++) {
            std::cout << x->coeff(i,j);
        }
        std::ostringstream w;
        unsigned int j = 0;
        for( ; j < x->ncols ; ) {
            unsigned long z = 0;
            for(unsigned int k = 0 ; j < x->ncols && k < ULONG_BITS ; k++,j++) {
                z <<= 1;
                z |= x->coeff(i,j);
            }
            w << setw(ULONG_BITS/4) << hex << z;
        }
        w.flush();
        std::cout << w.str() << "\n";
    }
}
void dpmat(polmat const * pa, unsigned int k)
{
    bmat x;
    pa->extract_coeff(x, k);
    dbmat(&x);
}/*}}}*/

namespace globals {
    unsigned int nrows;
    unsigned int ncols;
    unsigned int m,n;
    unsigned int t,t0,total_work;
    std::vector<unsigned int> delta;
    std::vector<unsigned int> chance_list;
    std::vector<unsigned int> gone_mad;

    polmat E;
    bmat e0;

    // F0 is exactly the n x n identity matrix, plus the X^(s-exponent)e_{cnum}
    // vectors. Here we store the cnum,exponent pairs.
    std::vector<std::pair<unsigned int, unsigned int> > f0_data;

    double start_time;
}

// To multiply on the right an m x n matrix A by F0, we start by copying
// A into the first n columns. Since we're also dividing out by X^t0, the
// result has to be shifted t0 positions to the right.
// Afterwards, column n+j of the result is column cnum[j] of A, shifted
// exponent[j] positions to the right.
void compute_E_from_A(polmat const &a)/*{{{*/
{
    using namespace globals;
    polmat tmp_E(n, m + n, a.ncoef - t0);
    for(unsigned int j = 0 ; j < n ; j++) {
        tmp_E.import_col_shift(j, a, j, - (int) t0);
    }
    for(unsigned int j = 0 ; j < m ; j++) {
        unsigned int cnum = f0_data[j].first;
        unsigned int exponent = f0_data[j].second;
        tmp_E.import_col_shift(n + j, a, cnum, - (int) exponent);
    }
    E.swap(tmp_E);
}/*}}}*/

// F is in fact F0 * PI.
// To multiply on the *left* by F, we cannot work directly at the column
// level (we could work at the row level, but it does not fit well with
// the way data is organized). So it's merely a matter of adding
// polynomials.
void compute_final_F_from_PI(polmat& F, polmat const& pi)/*{{{*/
{
    std::cout << "Computing final F from PI" << std::endl;
    using namespace globals;
    // We take t0 rows, so that we can do as few shifts as possible
    polmat tmpmat(t0,1,pi.ncoef);
    polmat tmp_F(n, m + n, globals::t0 + pi.ncoef);
    using namespace std;
    for(unsigned int i = 0 ; i < n ; i++) {
        // What contributes to entries in this row ?
        vector<pair<unsigned int, unsigned int> > l;
        set<unsigned int> sexps;
        // l.push_back(make_pair(i,0));
        for(unsigned int j = 0 ; j < m ; j++) {
            if (f0_data[j].first == i) {
                l.push_back(make_pair(n + j, f0_data[j].second));
                sexps.insert(f0_data[j].second);
            }
        }
        vector<unsigned int> exps(sexps.begin(), sexps.end());
        // So the i-th row of f0 has a 1 at position i, and then
        // X^(t0-f.second) at position n+j whenever f.first == i.
        //
        // Now fill in the row.
        for(unsigned int j = 0 ; j < m + n ; j++) {
            if (gone_mad[j])
                continue;
            // We zero out the whole column, it's less trouble.
            tmpmat.zcol(0);
            for(unsigned int k = 0 ; k < l.size() ; k++) {
                tmpmat.addpoly(l[k].second, 0, pi.poly(l[k].first, j));
            }
            tmp_F.addpoly(i,j,pi.poly(i,j));
            for(unsigned int k = 0 ; k < exps.size() ; k++) {
                tmpmat.xmul_poly(exps[k],0,t0-exps[k]);
                tmp_F.addpoly(i,j,tmpmat.poly(exps[k],0));
            }
        }
    }
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (gone_mad[j])
            continue;
        tmp_F.deg(j) = t0 + pi.deg(j);
    }
    F.swap(tmp_F);
}/*}}}*/


void read_data_for_series(polmat& A, int & rc)/*{{{*/
{
    using namespace globals;

    FILE * files[m][n];
    unsigned int nbys[n];

    polmat a(m,n,total_work+1);

    unsigned int i, j;
    unsigned int k;

    for (j = 0; j < n; j++) { nbys[j] = 0; }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n;) {
            unsigned int y = 0;
            char * ptr;
            int d;
            std::string filename = files::a % i % j;
            files[i][j] = fopen(filename.c_str(),"r");
            if (files[i][j] == NULL) {
                die("fopen(%s) : %s", errno, filename.c_str(),
                        strerror(errno));
            }
            printf("Reading file %s", filename.c_str());

            /* NOTE : we drop the first coefficient, because
             * we mean to work with the sequence generated
             * on x and By, so that we obtain a generator
             * afterwards.
             */

            char row[1024];
            unsigned long blah;
            /* We also use this in order to read the number
             * of coefficients per row */
            fgets(row, sizeof(row), files[i][j]);
            ptr = row;
            for(y = 0 ; sscanf(ptr,"%lu%n",&blah,&d) >= 1 ; ptr += d, y++);
            if (y == 0) {
                fprintf(stderr, "\nproblem while reading %s, line = %s\n",
                        filename.c_str(), row);
                abort();
            }
            if (y > 1) {
                printf(" [ %d values per row ]", y);
            }
            printf("\n");
            if (nbys[j] != 0) {
                BUG_ON(y != nbys[j]);
            }
            nbys[j] = y;
            j += y;
        }
    }

    unsigned int read_coeffs = 0;       /* please gcc */

    for (i = 0; i < m; i++) {
        unsigned long * pol[n];
        for (j = 0; j < n; j ++) {
            pol[j] = a.poly(i,j);
        }
        unsigned long mask = 1UL;
        unsigned int shift = 0;
        for (k = 0; k <= total_work; k++) {
            unsigned int rc = 0;
            for (j = 0; j < n; j += nbys[j]) {
                unsigned int l;
                for(l = 0 ; l < nbys[j] ; l++) {
                    unsigned long blah;
                    rc += 1 == fscanf(files[i][j], "%lu", &blah);
                    pol[j + l][shift] |= mask & -blah;
                }
            }
            if (rc == 0) {
                break;
            } else if (rc != n) {
                fprintf(stderr, "A files not in sync ; "
                        "read only %u coeffs in X^%d, row %d\n", rc, k, i);
                break;
            }
            mask <<= 1;
            shift += mask == 0;
            mask += mask == 0;
        }
        if (i == 0) {
            read_coeffs = k;
        } else {
            /* This happens when A files have gone live */
            BUG_ON(k != read_coeffs);
        }
    }
    for (j = 0; j < n; j ++) {
        a.deg(j) = read_coeffs - 1;
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j+= nbys[j]) {
            fclose(files[i][j]);
        }
    }

    printf("Stopped after reading %u coefficients (not counting 1st)\n",
            read_coeffs);

    rc = read_coeffs;
    A.swap(a);
}/*}}}*/

void write_polmat(polmat const& P, const char * fn)/*{{{*/
{
    std::ofstream f(fn);
    if (!f.is_open()) {
        perror(fn);
        exit(1);
    }

    for(unsigned int k = 0 ; k < P.ncoef ; k++) {
        for(unsigned int i = 0 ; i < P.nrows ; i++) {
            for(unsigned int j = 0 ; j < P.ncols ; j++) {
                if (j) f << " ";
                if ((int) k <= P.deg(j)) {
                    f << P.coeff(i,j,k);
                } else {
                    f << 0;
                }
            }
            f << "\n";
        }
        f << "\n";
    }
    printf("Written f to %s\n",fn);
}/*}}}*/

bool recover_f0_data(const char * fn)/*{{{*/
{
    std::ifstream f(fn);

    using namespace globals;
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned int exponent,cnum;
        if (!(f >> cnum >> exponent))
            return false;
        f0_data.push_back(std::make_pair(cnum,exponent));
    }
    std::cout << fmt("recovered init data % on disk\n") % fn;
    return true;
}/*}}}*/

bool write_f0_data(const char * fn)/*{{{*/
{
    std::ofstream f(fn);

    using namespace globals;
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned int cnum = f0_data[i].first;
        unsigned int exponent = f0_data[i].second;
        if (!(f << " " << cnum << " " << exponent))
            return false;
    }
    f << std::endl;
    std::cout << fmt("written init data % to disk\n") % fn;
    return true;
}/*}}}*/

// the computation of F0 says that the column number cnum of the
// coefficient X^exponent of A increases the rank.
//
// These columns are gathered into an invertible matrix at step t0, t0
// being strictly greater than the greatest exponent above. This is so
// that there is no trivial column dependency in F0.
void set_t0_delta_from_F0()/*{{{*/
{
    using namespace globals;
    t0 = f0_data.back().second + 1;     // see above for the +1
    for(unsigned int j = 0 ; j < m + n ; j++) {
        delta[j] = t0;
        // Even irrespective of the value of the exponent, we get t0 for
        // delta.
    }
}/*}}}*/

void write_F0_from_F0_quick()/*{{{*/
{
    using namespace globals;
    unsigned int s = t0;
    /* Now build f */
    polmat F0(n, m + n, s + 1);

    /* First n columns: identity matrix */
    /* rest: X^(s-exponent)cnum's */
    for(unsigned int i=0;i<n;i++) {
        F0.addcoeff(i,i,0,1);
        F0.deg(i) = 0;
    }
    for(unsigned int j=0;j<m;j++) {
        unsigned int cnum = f0_data[j].first;
        unsigned int exponent = f0_data[j].second;
        F0.addcoeff(cnum,n+j,t0-exponent,1);
        F0.deg(n + j) = t0-exponent;
    }
    write_polmat(F0, "F_INIT");
}/*}}}*/

void bw_commit_f(polmat& F) /*{{{*/
{
    using namespace globals;
    using namespace std;

    for (unsigned int j = 0; j < m + n ; j++) {
        std::string filename = files::f % j;
        if (gone_mad[j]) {
            cout << fmt("not writing % -- gone haywire since step %, "
                    "hence % steps ago\n")
                % filename % (t-gone_mad[j]) % gone_mad[j] ;
            continue;
        }
        ofstream f(filename.c_str());
        cout << "Writing " << filename << endl;
        if (!f.is_open()) {
            perror("writing f");
            die("ugh",1);
        }
        for (int k = 0; k <= F.deg(j); k++) {
            for (unsigned int i = 0; i < n; i++) {
                f << F.coeff(i,j,k) << (i == (n-1) ? "\n" : " ");
            }
        }
    }
}
/*}}}*/

void compute_f_init(polmat& A)/*{{{*/
{
    using namespace globals;

    /* build F_INIT */
    printf("Computing t0\n");

    /* For each integer i between 0 and m-1, we have a column, picked
     * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
     * the other ones, has coefficient at row pivots[i] unequal to zero.
     */
    bcol pcols[m];
    unsigned int pivots[m], exponent[m], cnum[m];

    unsigned int r = 0;
    for(unsigned int k=0;r < m && k<A.ncoef;k++) {
        bmat amat;
        A.extract_coeff(amat, k);
        for(unsigned int j=0;r < m && j<n;j++) {
            /* copy column j, coeff k into acol */
            bcol acol;
            amat.extract_col(acol, j);
            /* kill as many coeffs as we can */
            for(unsigned int v=0;v<r;v++) {
                unsigned int u=pivots[v];
                /* the v-th column in the matrix reduced_rank is known to
                 * kill coefficient u (more exactly, to have a -1 as u-th
                 * coefficient, and zeroes for the other coefficients
                 * referenced in the pivots[0] to pivots[v-1] indices).
                 */

                acol.add(pcols[v],acol.coeff(u));
            }
            unsigned int u = acol.ffs();
            if (u == UINT_MAX) {
                printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                        k,j,r);
                if (k * n > m + 40) {
                    printf("The choice of starting vectors was bad. "
                            "Cannot find %u independent cols within A\n",m);
                    exit(1);
                }
                continue;
            }

            /* Bingo, it's a new independent col. */
            pivots[r]=u;
            cnum[r]=j;
            exponent[r]=k;

            f0_data.push_back(std::make_pair(cnum[r], exponent[r]));

            /* TODO: For non-binary stuff, multiply the column so that
             * acol[u] becomes -1
             */
            pcols[r].swap(acol);
            r++;

            if (r == m)
                printf("[X^%d] A, col %d increases rank to %d (head row %d)\n",
                        k,j,r,u);
        }
    }

    t0 = exponent[r-1] + 1;
    printf("Found satisfying init data for t0=%d\n", t0);
                    
    /*
    printf("Init e0 matrix\n");
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            std::cout << A.coeff(i,j,t0);
        }
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << A.coeff(i,cnum[j],exponent[j]);
        }
        std::cout << "\n";
    }
    */


    if (r!=m) {
        printf("This amount of data is insufficient. "
                "Cannot find %u independent cols within A\n",m);
        exit(1);
    }

}/*}}}*/

void print_deltas()/*{{{*/
{
    using namespace globals;
    std::cout << fmt("[t=%[w4]] delta =") % t;
    unsigned int last = UINT_MAX;
    unsigned int nrep = 0;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        unsigned int d = delta[i];
        if (d == last) {
            nrep++;
            continue;
        }
        // Flush the pending repeats
        if (last != UINT_MAX) {
            std::cout << " " << last;
            if (nrep > 1)
                std::cout << "[" << nrep << "]";
        }
        last = d;
        nrep = 1;
    }
    BUG_ON(last == UINT_MAX);
    std::cout << " " << last;
    if (nrep > 1)
        std::cout << "[" << nrep << "]";
    std::cout << "\n";
}/*}}}*/

/*{{{*//* bw_init
 *
 * fill in the structures with the data available on disk. a(X) is
 * fetched this way. f(X) is chosen in order to make [X^0]e(X)
 * nonsingular. Once this condition is satisfied, e(X) is computed. All
 * further computations are done with e(X).
 *
 * a(X) is thrown away as soon as e(X) is computed.
 *
 * At a given point of the algorithm, one can bet that most of the data
 * for e(X) is useless.
 *
 * XXX Therefore it is wise to shrink e(x) aggressively. The new
 * organisation of the data prevents e(x) from being properly swapped
 * out. XXX TODO !!!
 *
 */
static void bw_init(void)
{
    int read_coeffs;
    polmat A;
    using namespace globals;

    printf("Using A(X) div X in order to consider Y as starting point\n");
    total_work--;


    printf("Reading scalar data in polynomial ``a''\n");
    read_data_for_series(A, read_coeffs);


    /* informational *//*{{{ */
    if (read_coeffs < (int) total_work) {
	printf("Since we do not have the full information yet "
	       "(only %d coefficients while %d were needed)"
	       ", we\n"
	       "merely can tell if the system looks like doable\n",
	       read_coeffs, total_work);
    }





    /*}}} */
    /* Data read stage completed. */
    /* TODO. Prepare the FFT engine for handling polynomial
     * multiplications up to ncmax coeffs */
    // unsigned int ncmax = (total_work<<1)+2;
    // ft_order_t::init(ncmax, ops_poly_args);

    delta.assign(m + n, -1);
    gone_mad.assign(m + n, 0);
    chance_list.assign(m + n, 0);

    if (!recover_f0_data(files::f_initq.c_str())) {
        // This is no longer useful
	// give_poly_rank_info(A, read_coeffs - 1);
	compute_f_init(A);
	write_f0_data(files::f_initq.c_str());
    }
    set_t0_delta_from_F0();

    // this is not used. It's only a debugging aid for the helper magma
    // code.
    write_F0_from_F0_quick();

    /*
        if (!read_f_init(F0, files::f_initq.c_str(), n, m + n)) {
            give_poly_rank_info(A, read_coeffs - 1);
            compute_f_init(A);
            write_f0_data(F0, files::f_initq.c_str());
        }
    */
    t = t0;
    printf("t0 = %d\n", t0);


    /*
       for(j=0;j< m + n ;j++) {
       global_delta[j]  = t_counter;
       chance_list[j]   = 0;
       }
     */

    total_work = read_coeffs;
    // A must be understood as A_computed + O(X^total_work)
    // therefore it does not make sense to compute coefficients of E at
    // degrees above total_work-1
    std::cout << "Computing value of E(X)=A(X)F(X) "
        << fmt("(degree %) [ +O(X^%) ]\n") % (total_work-1) % total_work;
    
    // multiply_slow(E, A, F0, t0, total_work-1);
    compute_E_from_A(A);

    printf("Throwing out a(X)\n");
    A.clear();
}
/*}}} */

static void rearrange_ordering(polmat & PI, unsigned int piv[])/*{{{*/
{
    /* Sort the columns. It might seem merely cosmetic and useless to
     * sort w.r.t both the global and local nominal degrees. In fact, it
     * is crucial for the corectness of the computations. (Imagine a
     * 2-step increase, starting with uneven global deltas, and hitting
     * an even situation in the middle. One has to sort out the local
     * deltas to prevent trashing the whole picture).
     *
     * The positional sort, however, *is* cosmetic (makes debugging
     * easier).
     */
    using namespace std;
    using namespace globals;
    typedef pair<pair<unsigned int,unsigned int>, int> corresp_t;
    vector<corresp_t> corresp(m+n);
    for(unsigned int i = 0 ; i < m + n ; i++) {
        corresp[i] = make_pair(make_pair(delta[i],PI.deg(i)), i);
    }
    sort(corresp.begin(), corresp.end(), less<corresp_t>());
    unsigned int p[m+n];
    for(unsigned int i = 0 ; i < m + n ; i++) {
        p[corresp[i].second] = i;
    }
    permute(delta, p);
    permute(chance_list, p);
    permute(gone_mad, p);
    if (piv) {
        for(unsigned int i = 0 ; i < m ; i++) {
            piv[i] = p[piv[i]];
        }
    }
    PI.perm(p);
    e0.perm(p);
}/*}}}*/

static void gauss(unsigned int piv[], polmat& PI)/*{{{*/
{
    /* Do one step of (column) gaussian elimination on e0. The columns
     * are assumed to be in the exact order that corresponds to the
     * ordering of the columns of PI (therefore relative to delta)
     */

    /* Note that we do *NOT* modify E. It seems to be the fastest option,
     * although the other deserves being investigated as well (see
     * old/lingen2.c, there are two versions of the quadratic algorithm).
     */
    unsigned int i,j,k;
    unsigned int rank;

    rank = 0 ;

    using namespace globals;

    /*
    std::cout << "Input matrix\n";
    dbmat(&e0);
    */

    for(j = 0 ; j < e0.ncols ; j++) {
        /* Find the pivot inside the column. */
        i = e0.ffs(j);
        if (i == UINT_MAX)
            continue;
        assert(rank < e0.nrows && rank < e0.ncols);
        // std::cout << fmt("col % is the %-th pivot\n") % j % rank;
        piv[rank++] = j;
        /* Cancel this coeff in all other columns. */
        for(k = j + 1 ; k < e0.ncols ; k++) {
            /* TODO : Over the binary field, this branch avoiding trick
             * could most probably be deleted, I doubt it gains anything. */
            unsigned long c = e0.coeff(i,k);
            /* add c times column j to column k */
            e0.acol(k,j,i,c);
            // E.acol(k,j,c);
            // This one is tempting, but it's a wrong assert.
            // ASSERT(PI.deg(j) <= PI.deg(k));
            // ASSERT(delta[j] <= delta[k]);
            ASSERT(std::make_pair(delta[j],PI.deg(j)) <= std::make_pair(delta[k],PI.deg(k)));
            PI.acol(k,j,c);
        }
        PI.xmul_col(j);
        // E.xmul_col(j);
        delta[j]++;
        if (PI.deg(j) >= (int) PI.ncoef) {
            if (gone_mad[j] == 0) {
                std::cout << fmt("Column % goes haywire from step %\n")%j%t;
            }
            gone_mad[j]++;
            PI.deg(j) = PI.ncoef-1;
        }
    }
    BUG_ON (rank != m);

    /*
    std::cout << "Invertible e0 matrix\n";
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(unsigned int j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */

    /*
    std::cout << "Invertible e0 matrix\n";
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(unsigned int j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */
}/*}}}*/

void read_mat_file_header(const char *name)/*{{{*/
{
    FILE *f = fopen(name, "r");
    if (f == NULL) {
        die("fopen(%s): %s", errno, name, strerror(errno));
    }
    char modulus[512];
    fscanf(f, "// %u ROWS %u COLUMNS ; MODULUS %s",
            &globals::nrows, &globals::ncols, modulus);
    BUG_ON(strcmp(modulus,"2") != 0);
    BUG_ON(globals::nrows > globals::ncols);

    fclose(f);
}/*}}}*/

#if 0/* {{{ */
const char *pi_meta_filename = "pi-%d-%d";

static int retrieve_pi_files(struct t_poly ** p_pi, int t_start)
{
    DIR * pi_dir;
    struct dirent * curr;
    int n_pi_files;
    struct couple {int s; int e;} * pi_files;
    const char *pattern;
    int i;
    struct t_poly * left = NULL, * right = NULL;
    ft_order_t order;
    int o_i;
    struct dft_bb * dft_left, * dft_right, * dft_prod;
    double tt;

    *p_pi=NULL;

    if ((pi_dir=opendir("."))==NULL) {
        perror(".");
        return t_start;
    }

    printf("Scanning directory %s for pi files\n", ".");

    pattern = strrchr(pi_meta_filename,'/');
    if (pattern == NULL) {
        pattern = pi_meta_filename;
    } else {
        pattern++;
    }

    for(n_pi_files=0;(curr=readdir(pi_dir))!=NULL;) {
        int s,e;
        if (sscanf(curr->d_name,pattern,&s,&e)==2) {
            printf("Found %s\n", curr->d_name);
            if (s>e) {
                printf("but that's a stupid one\n");
                continue;
            }
            n_pi_files++;
        }
    }

    if (n_pi_files==0) {
        printf("Found no pi files\n");
        return t_start;
    }

    pi_files=(struct couple *) malloc(n_pi_files*sizeof(struct couple));

    rewinddir(pi_dir);

    for(i=0;i<n_pi_files && (curr=readdir(pi_dir))!=NULL;) {
        if (sscanf(curr->d_name,
                    pattern,
                    &(pi_files[i].s),
                    &(pi_files[i].e))==2)
        {
            if (pi_files[i].s>pi_files[i].e)
                continue;
            i++;
        }
    }
    n_pi_files=i;
    closedir(pi_dir);

    /* The rule is: only look at the best candidate. It's not worth
     * bothering about more subtle cases */
    for(;;) {
        int t_max;
        int best=-1;
        FILE *f;
        char filename[FILENAME_LENGTH];

        printf("Scanning for data starting at t=%d\n",t_start);
        t_max=-1;
        for(i=0;i<n_pi_files;i++) {
            if (pi_files[i].s==t_start && pi_files[i].e != -1) {
                printf("candidate : ");
                printf(pattern,pi_files[i].s,pi_files[i].e);
                printf("\n");
                if (pi_files[i].e>t_max) {
                    t_max=pi_files[i].e;
                    best=i;
                }
            }
        }
        if (t_max==-1) {
            printf("Could not find such data\n");
            break;
        }

        sprintf(filename,pi_meta_filename,t_start,t_max);
        printf("trying %s\n", filename);
        f=fopen(filename,"r");
        if (f==NULL) {
            perror(filename);
            pi_files[best].e=-1;
            continue;
        }
        /* Which degree can we expect for t_start..t_max ?
         */

        unsigned int pideg;
        pideg = iceildiv(m * (t_max - t_start), (m+n));
        pideg += 10;
        if (t_max > total_work) {
            pideg += t_max - total_work;
        }

        right=tp_read(f, pideg);
        fclose(f);

        if (right==NULL) {
            printf("%s : bad or nonexistent data\n",filename);
            pi_files[best].e=-1;
            continue;
        }

        if (left==NULL) {
            left=right;
            right=NULL;
            t_start=t_max;
            continue;
        }

        printf("Beginning multiplication\n");
        *p_pi=tp_comp_alloc(left,right);
        core_if_null(*p_pi,"*p_pi");

        order.set((*p_pi)->degree+1);
        o_i = order;

        dft_left=fft_tp_dft(left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_left,"dft_left");

        dft_right=fft_tp_dft(right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_right,"dft_right");

        dft_prod=fft_bbb_conv(dft_left,dft_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_prod,"dft_prod");

        fft_tp_invdft(*p_pi,dft_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(left);
        tp_free(right);
        dft_bb_free(dft_left);
        dft_bb_free(dft_right);
        dft_bb_free(dft_prod);
        left=*p_pi;
        right=NULL;
        t_start=t_max;
    }
    free(pi_files);

    *p_pi=left;
    return t_start;
}

static void banner_traditional(int t, int deg, double inner, double * last)
{
    *last+=inner;
    if ((t+1)%print_min==0 || t==deg) {
        reclevel_prolog();
        printf("avg=%.1f	"
                "step:	%.2fs	last %d : %.2fs\n",
                ((double)global_sum_delta)/(m+n),inner,
                1+(t%print_min),*last);
        *last=0.0;
    }
}
#endif/*}}}*/

/* There are two choices for the quadratic algorithm. The first one
 * (default) seems to be a bit faster. Once I've got a more reasonable
 * implementation, it'll be time to do a serious comparison.
 *
 * The complexity curve is steeper with the first version.
 */

#if 0/*{{{*/
static void bw_traditional_algo_1(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    bw_mbmat e;
    unsigned int * pivlist;
    double inner_1,inner_2,last;
    int k;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));
    pivlist=(unsigned int *) malloc(m*sizeof(unsigned int));

    mbmat_alloc(e);
    mbmat_zero(e);

    assert(!ec_is_twisted(ec));
    last=0.0;

    for(t=0;t<=ec->degree;t++) {
        compute_ctaf(e,ec,pi,t,t?pivlist:NULL, &inner_1);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        if (check_chance)
            bw_check_chance(e,pi->clist);
        bw_gauss_onestep(e,ec,pi,delta,pivlist,	&inner_2);
        t_counter++;
        banner_traditional(t, ec->degree, inner_1 + inner_2, &last);
    }
    printf("DELTA : ( ");
    for(k=0;k< m + n ;k++) printf("%d ",delta[k]);
    printf(")\n");

    ec_advance(ec,ec->degree+1);	/* cosmetic */
    mbmat_free(e);
    free(pivlist);
    free(perm);
}

static void bw_traditional_algo_2(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    int deg;
    double inner,last;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));

    assert(!ec_is_twisted(ec));

    deg=ec->degree;
    last=0.0;

    for(t=0;t<=deg;t++) {
        if (check_chance)
            bw_check_chance(mbpoly_coeff(ec->p,0),ec->clist);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        ec_apply_perm(ec,perm);
        bw_gauss_onestep(mbpoly_coeff(ec->p,0),ec,pi,delta,NULL,&inner);
        ec_advance(ec,1);
        t_counter++;
        banner_traditional(t, deg, inner, &last);
    }

    free(perm);
}
#endif/*}}}*/

static unsigned long extract_coeff_degree_t(unsigned int t, unsigned long const * a, unsigned int da, unsigned long const * b, unsigned int db)/*{{{*/
{
    unsigned long c = 0;
    /*
    0 <= s <= t
    0 <= s <= na
    0 <= t-s <= nb
    t-nb <= s <= t
    */
    unsigned int high = std::min(t, da);
    unsigned int low = std::max(0u, t - db);
    for(unsigned int s = low ; s <= high ; s++) {
        unsigned int si = s / ULONG_BITS;
        unsigned int ss = s % ULONG_BITS;
        unsigned int ri = (t-s) / ULONG_BITS;
        unsigned int rs = (t-s) % ULONG_BITS;
        c ^= a[si] >> ss & b[ri] >> rs;
    }
    return c & 1UL;
}/*}}}*/


static void extract_coeff_degree_t(unsigned int tstart, unsigned int dt, unsigned int piv[], polmat const& PI)/*{{{*/
{
    using namespace std;
    using namespace globals;
    vector<bool> known(m+n,false);
    if (piv != NULL) {
        for(unsigned int i = 0 ; i < m ; i++)
            known[piv[i]] = true;
    }

    vector<unsigned int> z;

    for(unsigned int j = 0 ; j < m+n ; j++) {
        if (known[j]) continue;
        e0.zcol(j);
        unsigned long some_nonzero = 0;
        for(unsigned int i = 0 ; i < m ; i++) {
            unsigned long c = 0;
            for(unsigned int k = 0 ; k < m + n ; k++) {
                c ^= extract_coeff_degree_t(dt,
                        E.poly(i, k), E.ncoef + 1,
                        PI.poly(k, j), PI.deg(j));
            }
            e0.addcoeff(i,j,c);
            some_nonzero |= c;
        }
        if (!some_nonzero) {
            z.push_back(j);
        }
    }
    vector<unsigned int> ncha(m + n, 0);
    if (z.empty()) {
        chance_list.swap(ncha);
    } else {
        for(unsigned int i = 0 ; i < z.size() ; i++) {
            ++ncha[z[i]];
        }
        std::cout << fmt("Step %, zero cols:") % (tstart + dt);
        for(unsigned int i = 0 ; i < m + n ; i++) {
            if (ncha[i] == 0) {
                chance_list[i] = 0;
            } else {
                unsigned int w = chance_list[i] += ncha[i];
                std::cout << " " << i;
                if (w >= 2) {
                    std::cout << fmt("[%]") %w;
                }
            }
        }
        std::cout << "\n";
    }
}/*}}}*/

/*
 * Rule : in the following, ec can be trashed at will.
 *
 * Compute pi_left, of degree (ec->degree + 1)*(m/(m+n)), such that ec * pi is
 * divisible by X^(ec->degree + 1) (that is, all coefficients up to
 * degree ec->degree in the product are forced to 0).
 *
 */

static void go_quadratic(polmat& pi)/*{{{*/
{
    using namespace globals;

    unsigned int piv[m];
    unsigned int deg = E.ncoef - 1;

    polmat tmp_pi(m + n, m + n, deg * m / (m + n) + 10);
    for(unsigned int i = 0 ; i < m + n ; i++) {
        tmp_pi.addcoeff(i,i,0,1UL);
        tmp_pi.deg(i) = 0;
    }

    rearrange_ordering(tmp_pi, NULL);

    unsigned int tstart = t;

    for (unsigned int dt = 0; dt <= deg ; dt++) {
        double delta;
        delta = seconds() - start_time;

        double percent = (double) dt / (deg + 1);
        percent = percent * percent;

        double estim_final = delta / percent;

        percent *= 100.0;
            
        printf("%5.0f / est %-7.0f (%2.0f%%) ",
                delta, estim_final, percent);

        print_deltas();
	extract_coeff_degree_t(tstart, dt, dt ? piv : NULL, tmp_pi);
        gauss(piv, tmp_pi);
        rearrange_ordering(tmp_pi, piv);
        t++;
        // if (t % 60 < 10 || deg-dt < 30)
        // write_polmat(tmp_pi,std::string(fmt("pi-%-%") % tstart % t).c_str());
    }
    pi.swap(tmp_pi);
    write_polmat(pi,std::string(fmt("pi-%-%") % tstart % t).c_str());
    print_deltas();
}/*}}}*/


static void compute_lingen(polmat& pi);

template<typename fft_type>
static void go_recursive(polmat& pi)
{
    using namespace globals;
    unsigned int deg = E.ncoef - 1;
    unsigned int ldeg = deg / 2 + 1;
    unsigned int rdeg = (deg + 1) / 2;

    /* Repartition of the job:
     *
     *		left		right
     * deg==0	1		0	(never recursive)
     * deg==1	1		1
     * deg==2	2		1
     * deg==n	n/2 + 1		(n+1)/2
     * 
     * The figures are for the number of steps, each one corres-
     * ponding to a m/(m+n) increase of the average degree of pi.
     */
    /* We aim at computing ec * pi / X^ldeg. The degree of this
     * product will be
     *
     * ec->degree + pi->degree - ldeg
     *
     * (We are actually only interested in the low (ec->degree-ldeg)
     * degree part of the product, but the whole thing is required)
     *
     * The expected value of pi->degree is 
     * 	ceil(ldeg*m/(m+n))
     *
     * The probability that pi exceeds this expected degree
     * depends on the base field, but is actually low.
     * However, by the end of the computations, this does
     * happen because the degrees increase unevenly.
     *
     * The DFTs of e and pi can be computed using only the
     * number of points given above, *even if their actual
     * degree is higher*. The FFT routines need to have
     * provision for this.
     *
     * The number of points will then be the smallest power
     * of 2 above deg+ceil(ldeg*m/(m+n))-ldeg+1
     */

    assert(ldeg && rdeg && ldeg + rdeg == deg + 1);

    unsigned int expected_pi_deg = 10 + iceildiv(ldeg*m, (m+n));
    unsigned int kill;
    unsigned int tstart = t;

#ifdef	HAS_CONVOLUTION_SPECIAL
    kill=ldeg;
#else
    kill=0;
#endif

    std::cout << "Recursive call, degree " << deg << std::endl;
    fft_type o(deg, expected_pi_deg, deg + expected_pi_deg - kill, m + n);

    tpolmat<fft_type> E_hat;
    transform(E_hat, E, o, deg);

    E.resize(ldeg - 1 + 1);
    polmat pi_left;
    compute_lingen(pi_left);
    E.clear();

    if (t < t0 + ldeg) {
        printf("Exceptional situation, small generator ; escaping\n");
        pi.swap(pi_left);
        return;
    }

    int pi_l_deg = pi_left.maxdeg();

    printf("deg(pi_l)=%d, bound is %d\n",pi_l_deg,expected_pi_deg);
    if (pi_l_deg >= (int) expected_pi_deg) {/*{{{*/
        printf("Warning : pi grows above its expected degree...\n");
        /*
        printf("order %d , while :\n"
                "deg=%d\n"
                "deg(pi_left)=%d\n"
                "ldeg-1=%d\n"
                "hence, %d is too big\n",
                so_i,
                deg,pi_left->degree,ldeg - 1,
                deg+pi_left->degree-kill + 1);
                */

        int n_exceptional = 0;
        for(unsigned int i = 0 ; i < m + n ; i++) {
            n_exceptional += chance_list[i] * m;
        }
        if (!n_exceptional) {
            die("This should only happen at the end of the computation\n",1);
        }
        pi.swap(pi_left);
        return;
    }/*}}}*/

    tpolmat<fft_type> pi_l_hat;
    transform(pi_l_hat, pi_left, o, pi_l_deg);
    pi_left.clear();

    tpolmat<fft_type> E_middle_hat;

    /* XXX Do special convolutions here */
    compose(E_middle_hat, E_hat, pi_l_hat, o);
    E_hat.clear();
    /* pi_l_hat is used later on ! */

    itransform(E, E_middle_hat, o, deg + pi_l_deg - kill);
    E_middle_hat.clear();

    /* Make sure that the first ldeg-kill coefficients of all entries of
     * E are zero. It's only a matter of verification, so this does not
     * have to be critical. */

    ASSERT(E.valuation() >= ldeg - kill);

    /* This chops off some data */
    E.xdiv_resize(ldeg - kill, rdeg - 1 + 1);

    polmat pi_right;
    compute_lingen(pi_right);
    int pi_r_deg = pi_right.maxdeg();
    E.clear();

    tpolmat<fft_type> pi_r_hat;
    transform(pi_r_hat, pi_right, o, pi_r_deg);
    pi_right.clear();

    tpolmat<fft_type> pi_hat;
    compose(pi_hat, pi_l_hat, pi_r_hat, o);
    pi_l_hat.clear();
    pi_r_hat.clear();

    itransform(pi, pi_hat, o, pi_l_deg + pi_r_deg);

    write_polmat(pi,std::string(fmt("pi-%-%") % tstart % t).c_str());
}

static void compute_lingen(polmat& pi)
{
    /* reads the data in the global thing, E and delta. ;
     * compute the linear generator from this.
     */
    using namespace globals;
    unsigned int deg = E.ncoef - 1;

    if (deg <= rec_threshold) {
        go_quadratic(pi);
    } else if (deg < 3264) {
        /* The bound is such that deg + deg/4 is 64 words or less */
        go_recursive<fake_fft>(pi);
    } else {
        /* Presently, c128 requires input polynomials that are large
         * enough.
         */
        go_recursive<cantor_fft>(pi);
    }
}


void recycle_old_pi(polmat& pi)
{
#if 0/*{{{*/
    if (pi_left!=NULL && !(!check_input && ec->degree<new_t-t_counter)) {
        int dg_kill;
        struct dft_bb * dft_pi_left;
        struct dft_mb * dft_e_left, * dft_e_middle;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (e_left*pi_left)\n");

        /* That's a pity to bring the DFT so far, but we have to
         * check the input is correct */

        dg_kill=new_t-t_counter;

#ifdef	HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            order.set(1+pi_left->degree+ec->degree);
        } else {
            order.set(1+pi_left->degree+ec->degree-dg_kill);
        }
#else
        order.set(1+pi_left->degree+ec->degree);
#endif
        o_i = order;

        dft_e_left=fft_ec_dft(ec,order,&tt);
        printf("DFT(e_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_left,"dft_e_left");

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        ec_untwist(ec);
#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    0,&tt);
        } else {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    dg_kill,&tt);
            ec_advance(ec,dg_kill);
        }
#else
        dft_e_middle=fft_mbb_conv(dft_e_left,dft_pi_left,&tt);
#endif

        printf("CONV(e_left,pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_middle,"dft_e_middle");

#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
            printf("Verifying the product\n");
            check_zero_and_advance(ec, dg_kill);
            printf("Input OK\n");
        } else {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree-dg_kill,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        }
#else
        fft_mb_invdft(ec->p,dft_e_middle,
                ec->degree,&tt);
        printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        printf("Verifying the product\n");
        check_zero_and_advance(ec, dg_kill);
        printf("Input OK\n");
#endif
        tp_act_on_delta(pi_left,global_delta);

        dft_bb_free(dft_pi_left);
        dft_mb_free(dft_e_left);
        dft_mb_free(dft_e_middle);
        t_counter=new_t;

        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
    if (pi_left!=NULL && !check_input && ec->degree<new_t-t_counter) {
        printf("We are not interested in the computation of e(X)\n");
        ec_advance(ec,new_t-t_counter);
        tp_act_on_delta(pi_left,global_delta);
        t_counter=new_t;
        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
#endif/*}}}*/
}

void block_wiedemann(void)
{
    bw_init();


#if 0/*{{{*/
    if (check_pi) {
        if (retrieve_pi_files(&pi_left)) {
            recycle_old_pi(pi_left);
        }
    } else {
        printf("Not reading pi files due to --nopi option\n");
    }
#endif/*}}}*/

    using namespace globals;
    using namespace std;

    cout << fmt("E: % coeffs, t=%\n") % E.ncoef % t;

    { bmat tmp_e0(m,m + n); e0.swap(tmp_e0); }

    for(unsigned int i = 0 ; i < m + n ; i++) {
        E.deg(i) = E.ncoef - 1;
    }

    start_time = seconds();

    // E.resize(deg + 1);
    polmat pi_left;
    compute_lingen(pi_left);

    polmat F;
    compute_final_F_from_PI(F, pi_left);
    bw_commit_f(F);

    printf("// step %d LOOK [", t);
    for (unsigned int j = 0; j < m + n; j++) {
        if (chance_list[j])
            printf(" %d", j);
    }
    printf(" ]\n");

#if 0/*{{{*/
    if (ec->degree>=0) {
        bw_lingen(ec,global_delta,&pi_right);
    } else {
        pi_right=pi_left;
        pi_left=NULL;
    }

    if (pi_left!=NULL) {
        struct dft_bb * dft_pi_left, * dft_pi_right, * dft_pi_prod;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (pi_left*pi_right)\n");
        pi_prod=tp_comp_alloc(pi_left,pi_right);
        core_if_null(pi_prod,"pi_prod");

        order.set(pi_prod->degree+1);
        o_i = order;

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        dft_pi_right=fft_tp_dft(pi_right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_right,"dft_pi_right");

        dft_pi_prod=fft_bbb_conv(dft_pi_left,dft_pi_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_prod,"dft_pi_prod");

        fft_tp_invdft(pi_prod,dft_pi_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(pi_left);
        tp_free(pi_right);
        dft_bb_free(dft_pi_left);
        dft_bb_free(dft_pi_right);
        dft_bb_free(dft_pi_prod);
        pi_left=NULL;
        pi_right=NULL;
        if (check_pi) {
            save_pi(pi_prod,t_start,new_t,t_counter);
        }
        /* new_t is the inner value that has to be discarded */
    } else {
        pi_prod=pi_right;
        /* Don't save here, since it's already been saved by
         * lingen (or it comes from the disk anyway).
         */
    }
#endif/*}}}*/
    // print_chance_list(total_work, chance_list);
}

void usage()
{
    die("Usage: ./bw-lingen-binary --subdir <dir> -t <threshold> <matrix name> <m> <n>\n",1);
}

int main(int argc, char *argv[])
{
    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    argv++, argc--;

    int pop = 0;

    using namespace globals;

    for( ; argc ; ) {
        if (strcmp(argv[0], "--subdir") == 0) {
            if (argc <= 1) usage();
            int rc = chdir(argv[1]);
            if (rc < 0) {
                perror(argv[1]);
                exit(errno);
            }
            argv+=2;
            argc-=2;
            continue;
        }
        if (strcmp(argv[0], "-t") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            rec_threshold = atoi(argv[0]);
            argv++, argc--;
            continue;
        }
#if 0/*{{{*/
        if (strcmp(argv[0], "-p") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            print_min = atoi(argv[0]);
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--enable-complex-field") == 0
                || strcmp(argv[0], "--fermat-prime") == 0)
        {
            ops_poly_args.push_back(argv[0]);
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--no-check-input") == 0) {
            check_input = 0;
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--nopi") == 0) {
            check_pi = 0;
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "-q") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            preferred_quadratic_algorithm = atoi(argv[0]);
            argv++, argc--;
            continue;
            /* changes the quadratic algorithm employed
               {"coppersmith old uneven",0},
               {"new even",1},
             */
        }
#endif/*}}}*/
        if (pop == 0) {
            read_mat_file_header(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        } else if (pop == 1) {
            m = atoi(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        } else if (pop == 2) {
            n = atoi(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        }

        usage();
    }

    if (m == 0 || n == 0 || rec_threshold == 0) {
        usage();
    }

    total_work = Lmacro(ncols, m, n);

    coredump_limit(1);

    block_wiedemann();

    return 0;
}

/* vim: set sw=4 sta et: */
