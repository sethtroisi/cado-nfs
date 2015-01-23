#include "cado.h"

#include <ostream>
#include <sstream>
#include <algorithm>

#include "lingen_mat_types.hpp"
#include "lingen_qcode.h"
#include "utils.h"

/* {{{ utility */
template < typename iterator >
    std::string intlist_to_string(iterator t0, iterator t1)
{
    std::ostringstream cset;
    for (iterator t = t0; t != t1;) {
	iterator u;
	for (u = t; u != t1; u++) {
	    if ((typename std::iterator_traits <
		 iterator >::difference_type) (*u - *t) != (u - t))
		break;
	}
	if (t != t0)
	    cset << ',';
	if (u - t == 1) {
	    cset << *t;
	} else {
	    cset << *t << "-" << u[-1];
	}
	t = u;
    }
    return cset.str();
}

/* }}} */

static unsigned long extract_coeff_degree_t(unsigned int t,
					    unsigned long const *a,
					    unsigned int da,
					    unsigned long const *b,
					    unsigned int db)
{				/*{{{ */
    unsigned long c = 0;
    /*
       0 <= s <= t
       0 <= s <= na
       0 <= t-s <= nb
       t-nb <= s <= t
     */
    unsigned int high = std::min(t, da);
    unsigned int low = std::max(0u, t - db);
    for (unsigned int s = low; s <= high; s++) {
	unsigned int si = s / ULONG_BITS;
	unsigned int ss = s % ULONG_BITS;
	unsigned int ri = (t - s) / ULONG_BITS;
	unsigned int rs = (t - s) % ULONG_BITS;
	c ^= a[si] >> ss & b[ri] >> rs;
    }
    return c & 1UL;
}				/*}}} */


static void extract_coeff_degree_t(lingen_qcode_data_ptr qq, bmat & e0,
				   unsigned int dt,
				   unsigned int piv[], polmat const &PI)
{				/*{{{ */
    using namespace std;
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;

    vector < bool > known(m + n, false);
    if (piv != NULL) {
	for (unsigned int i = 0; i < m; i++)
	    known[piv[i]] = true;
    }

    vector < unsigned int >z;

    /* Note 1: this routine is not optimal for locality. It would be better
       to first loop on i, then k, then j. But then the if (known[j]) continue
       would be in the inner loop. Otherwise first loop on j, then k, then i.
       Note 2: the inner routine extract_coeff_degree_t() requires to reverse
       the coefficients of PI.poly(k, j) [or E.poly(i, k), since the popcount
       is symmetrical]. If we first loop on (i,k) or (j,k), then we could
       precompute the reverse polynomial, and give it to the inner routine,
       which would then just perform a xor of the relevant bit part, and a
       popcount. */
    for (unsigned int j = 0; j < m + n; j++) {
	if (known[j])
	    continue;
	e0.zcol(j);
	unsigned long some_nonzero = 0;
	for (unsigned int i = 0; i < m; i++) {
	    unsigned long c = 0;
	    for (unsigned int k = 0; k < m + n; k++) {
		c ^= extract_coeff_degree_t(dt,
					    qq->iptrs[i * b + k],
					    qq->length - 1, PI.poly(k, j),
					    PI.deg(j));
	    }
	    e0.addcoeff(i, j, c);
	    some_nonzero |= c;
	}
	if (!some_nonzero) {
	    z.push_back(j);
	}
    }
    vector < unsigned int >ncha(m + n, 0);
    if (z.empty()) {
	memset(qq->ch, 0, (m + n) * sizeof(unsigned int));
    } else {
	for (unsigned int i = 0; i < z.size(); i++) {
	    ++ncha[z[i]];
	}

	/* resets the global chance_list counter */
	for (unsigned int i = 0; i < m + n; i++) {
	    if (ncha[i] == 0) {
		qq->ch[i] = 0;
	    } else {
		qq->ch[i]++;
	    }
	}

	printf("%-8u%zucols=0:", qq->t + dt, z.size());

	vector < pair < unsigned int, unsigned int >>zz;
	for (unsigned int i = 0; i < z.size(); i++) {
	    zz.push_back(make_pair(z[i], qq->ch[z[i]]));
	}

	// Now print this out more nicely.
	sort(zz.begin(), zz.end());
	for (; zz.size();) {
	    unsigned int mi = UINT_MAX;
	    for (unsigned int i = 0; i < zz.size(); i++) {
		if (zz[i].second < mi) {
		    mi = zz[i].second;
		}
	    }
	    printf(" [");
	    for (unsigned int i = 0; i < zz.size();) {
		unsigned int j;
		for (j = i; j < zz.size(); j++) {
		    if (zz[j].first - zz[i].first != j - i)
			break;
		}
		if (i)
		    printf(",");
		if (zz[i].first == zz[j - 1].first - 1) {
		    printf("%u,%u", zz[i].first, zz[j - 1].first);
		} else if (zz[i].first < zz[j - 1].first) {
		    printf("%u..%u", zz[i].first, zz[j - 1].first);
		} else {
		    printf("%u", zz[i].first);
		}
		i = j;
	    }
	    printf("]");
	    if (mi > 1)
		printf("*%u", mi);

	    vector < pair < unsigned int, unsigned int >>zz2;
	    for (unsigned int i = 0; i < zz.size(); i++) {
		if (zz[i].second > mi) {
		    zz2.push_back(make_pair(zz[i].first, zz[i].second));
		}
	    }
	    zz.swap(zz2);
	}
	printf("\n");
    }
}				/*}}} */

// applies a permutation on the source indices/*{{{*/
template < typename T > void permute(T * a, unsigned int n, unsigned int p[])
{
    T *b = new T[n];
    for (unsigned int i = 0; i < n ; i++) {
	b[p[i]] = a[i];
    }
    memcpy(a, b, n * sizeof(T));
    delete[]b;
}				/*}}} */

static void rearrange_ordering(lingen_qcode_data_ptr qq, bmat & e0,
			       polmat & PI, unsigned int piv[])
{				/*{{{ */
    unsigned int m = qq->m;
    unsigned int n = qq->b - m;
    /* Sort the columns. It might seem merely cosmetic and useless to
     * sort w.r.t both the global and local nominal degrees. In fact, it
     * is crucial for the correctness of the computations. (Imagine a
     * 2-step increase, starting with uneven global deltas, and hitting
     * an even situation in the middle. One has to sort out the local
     * deltas to prevent trashing the whole picture).
     *
     * The positional sort, however, *is* cosmetic (makes debugging
     * easier).
     */
    using namespace std;
    typedef pair < pair < unsigned int, unsigned int >, int >corresp_t;
    vector < corresp_t > corresp(m + n);
    for (unsigned int i = 0; i < m + n; i++) {
	int pideg = PI.deg(i);
	if (pideg == -1) {	/* overflowed ! */
	    pideg = INT_MAX;
	}
	corresp[i] = make_pair(make_pair(qq->delta[i], pideg), i);
    }
    sort(corresp.begin(), corresp.end(), less < corresp_t > ());
    unsigned int p[m + n];
    for (unsigned int i = 0; i < m + n; i++) {
	p[corresp[i].second] = i;
    }
    permute(qq->delta, m + n, p);
    permute(qq->ch, m + n, p);
    if (piv) {
	for (unsigned int i = 0; i < m; i++) {
	    piv[i] = p[piv[i]];
	}
    }
    PI.perm(p);
    e0.perm(p);
}				/*}}} */

static bool gauss(lingen_qcode_data_ptr qq, bmat & e0, unsigned int dt, unsigned int piv[],
		  polmat & PI)
{				/*{{{ */
    unsigned int m = qq->m;
    unsigned int n = qq->b - m;
    unsigned int i, j, k;
    unsigned int rank;
    rank = 0;

    std::vector < unsigned int >overflowed;

    for (j = 0; j < e0.ncols; j++) {
	/* Find the pivot inside the column. */
	i = e0.ffs(j);
	if (i == UINT_MAX)
	    continue;
	ASSERT(rank < e0.nrows && rank < e0.ncols);
	// std::cout << fmt("col % is the %-th pivot\n") % j % rank;
	piv[rank++] = j;
	/* Cancel this coeff in all other columns. */
	for (k = j + 1; k < e0.ncols; k++) {
	    /* TODO : Over the binary field, this branch avoiding trick
	     * could most probably be deleted, I doubt it gains anything. */
	    unsigned long c = e0.coeff(i, k);
	    /* add c times column j to column k */
	    e0.acol(k, j, i, c);
	    // E.acol(k,j,c);
	    // This one is tempting, but it's a wrong assert.
	    // ASSERT(PI.deg(j) <= PI.deg(k));
	    // ASSERT(delta[j] <= delta[k]);
	    ASSERT(std::make_pair(qq->delta[j], PI.deg(j)) <=
		   std::make_pair(qq->delta[k], PI.deg(k)));
	    PI.acol(k, j, c);
	}
	PI.xmul_col(j);
	// E.xmul_col(j);
	qq->delta[j]++;
	if (PI.deg(j) >= (int) PI.ncoef - 1) {
	    overflowed.push_back(j);
	    /* Then don't hesitate. Set the degree to -1 altoghether.
	     * There's not much meaning remaining in this vector anyway,
	     * so we'd better trash it */
	}
    }
    ASSERT_ALWAYS(rank == m);

    bool finished = !overflowed.empty();

    if (finished) {
        std::string s = intlist_to_string(overflowed.begin(), overflowed.end());

        printf
            ("%-8u** %zu cols (%s) exceed maxdeg=%ld (normal at the end) **\n",
             qq->t + dt, overflowed.size(), s.c_str(), PI.ncoef - 1);
        unsigned ctot = 0;
        for (unsigned int j = 0; j < m + n; j++) {
            ctot += qq->ch[j];
        }
        ASSERT_ALWAYS(ctot > 0);
    }

    return finished;
}				/*}}} */


unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;

    polmat tmp_pi(m + n, m + n, qq->outlength);
    bmat e0(m, m + n);
    unsigned int piv[m];

    for (unsigned int i = 0; i < m + n; i++) {
	tmp_pi.addcoeff(i, i, 0, 1UL);
	tmp_pi.deg(i) = 0;
    }
    for (unsigned int i = 0; i < m; i++) {
        piv[m] = 0;
    }


    rearrange_ordering(qq, e0, tmp_pi, NULL);

    /* qq->t will remain equal to tstart for the duration of the
     * computation */

    bool finished = false;
    unsigned int dt = 0;
    for (dt = 0; !finished && dt < qq->length ; dt++) {
	extract_coeff_degree_t(qq, e0, dt, dt ? piv : NULL, tmp_pi);
	finished = gauss(qq, e0, dt, piv, tmp_pi);
	rearrange_ordering(qq, e0, tmp_pi, piv);
    }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        for(unsigned int j = 0 ; j < m + n ; j++) {
            memcpy(qq->optrs[i * b + j], tmp_pi.poly(i,j), iceildiv(qq->outlength, ULONG_BITS) * sizeof(unsigned long));
        }
    }
    for(unsigned int j = 0 ; j < m + n ; j++) {
        qq->local_delta[j] = tmp_pi.deg(j);
    }
    qq->t += dt;
    return qq->length;
}
