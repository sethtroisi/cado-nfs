#include "cado.h"

#include <ostream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <vector>

#include "lingen_mat_types.hpp"
#include "lingen_qcode.h"
#include "utils.h"

using namespace std;

/* if SWAP_E is defined, put E[i,j] into E[j][i] */
#define SWAP_E

/* if SWAP_PI is defined, put pi[i,j] into pi[j][i] */
#define SWAP_PI

class ulmat_rowmajor {
    unsigned int n;
    unsigned long * x;
    ulmat_rowmajor(ulmat_rowmajor const&) {}
public:
    ulmat_rowmajor(unsigned int m, unsigned int n) : n(n)
    {
        x = new unsigned long[m*n];
        memset(x, 0, m*n*sizeof(unsigned long));
    }
    ~ulmat_rowmajor() { delete[] x; }
    inline unsigned long& operator()(unsigned int i, unsigned int j) {
        return x[i * n + j];
    }
    inline const unsigned long& operator()(unsigned int i, unsigned int j) const {
        return x[i * n + j];
    }
    inline unsigned long * row(unsigned int i) { return x + i * n; }
    inline const unsigned long * row(unsigned int i) const { return x + i * n; }
};

class ulmat_colmajor {
    unsigned int m;
    unsigned long * x;
    ulmat_colmajor(ulmat_colmajor const&) {}
public:
    ulmat_colmajor(unsigned int m, unsigned int n) : m(m)
    {
        x = new unsigned long[m*n];
        memset(x, 0, m*n*sizeof(unsigned long));
    }
    ~ulmat_colmajor() { delete[] x; }
    inline unsigned long& operator()(unsigned int i, unsigned int j) {
        return x[j * m + i];
    }
    inline const unsigned long& operator()(unsigned int i, unsigned int j) const {
        return x[j * m + i];
    }
    inline unsigned long * column(unsigned int j) { return x + j * m; }
    inline const unsigned long * column(unsigned int j) const { return x + j * m; }
};


bool report_spontaneous_zeros(lingen_qcode_data_ptr qq, unsigned int dt, vector<bool>const& is_modified)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;
    vector<pair<unsigned int, unsigned int> > magic;
    bool changed = false;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        if (is_modified[i]) {
            qq->ch[i] = 0;
            changed=true;
        } else {
            qq->ch[i]++;
            magic.push_back(make_pair(i, qq->ch[i]));
        }
    }

    if (magic.empty()) {
        ASSERT_ALWAYS(changed);
        return changed;
    }

    printf("%-8u%zucols=0:", qq->t + dt, magic.size());

    // Now print this out more nicely.
    sort(magic.begin(),magic.end());
    for( ; magic.size() ; ) {
        unsigned int mi = UINT_MAX;
        for(unsigned int i = 0 ; i < magic.size() ; i++) {
            if (magic[i].second < mi) {
                mi = magic[i].second;
            }
        }
        printf(" [");
        for(unsigned int i = 0 ; i < magic.size() ; ) {
            unsigned int j;
            for(j = i; j < magic.size() ; j++) {
                if (magic[j].first-magic[i].first != j-i) break;
            }
            if (i) printf(",");
            if (magic[i].first == magic[j-1].first - 1) {
                printf("%u,%u", magic[i].first, magic[j-1].first);
            } else if (magic[i].first < magic[j-1].first) {
                printf("%u..%u", magic[i].first, magic[j-1].first);
            } else {
                printf("%u", magic[i].first);
            }
            i = j;
        }
        printf("]");
        if (mi > 1)
            printf("*%u",mi);

        vector<pair<unsigned int, unsigned int> > zz2;
        for(unsigned int i = 0 ; i < magic.size() ; i++) {
            if (magic[i].second > mi) {
                zz2.push_back(make_pair(magic[i].first,magic[i].second));
            }
        }
        magic.swap(zz2);
    }
    printf("\n");

    return changed;
}


unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    unsigned int m = qq->m;
    unsigned int b = qq->b;
    unsigned int n = b - m;

#ifdef  SWAP_E
    ulmat_colmajor E(m, b);
#else
    ulmat_rowmajor E(m, b);
#endif

#ifdef  SWAP_PI
    ulmat_colmajor P(b, b);
#else
    ulmat_rowmajor P(b, b);
#endif

    ASSERT_ALWAYS(qq->length <= ULONG_BITS);
    for (unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            E(i, j) = qq->iptrs[i * b + j][0];
        }
    }

    for (unsigned int i = 0; i < b; i++) {
        P(i, i) = 1;
        qq->local_delta[i] = 0;
    }

    /* We need to keep track of the qq->ch array. It checks how many
     * times in a row a given column turns out to be magically zero. 
     *
     * For this, we'll use a local table (is_modified[] below), and we
     * will use it to report which columns happen to be zero, and update
     * the qq->ch array.
     */

    unsigned int e = 0;
    for ( ; e < qq->length; e++) {
	uint64_t mask = (uint64_t) 1 << e;
        vector<bool> is_modified(m + n, false);
#if 0
	unsigned int md = UINT_MAX;
	for (unsigned int j = 0; j < m + n; j++)
	    if (qq->delta[j] > md)
		md = qq->delta[j];
#endif
	for (unsigned int i = 0; i < m; i++) {
	    unsigned int min_degree = UINT_MAX;
            unsigned int pivot = m + n;
	    for (unsigned int j = 0; j < m + n; j++)
		if ((E(i, j) & mask) && (qq->delta[j] < min_degree)) {
		    min_degree = qq->delta[j];
		    pivot = j;
		}
	    ASSERT_ALWAYS(pivot < m + n);
            is_modified[pivot]=1;
	    for (unsigned int k = 0; k < m + n; k++) {
                if (k == pivot) continue;
		if (!(E(i, k) & mask)) continue;
                is_modified[k] = true;
#ifdef SWAP_E
                mpn_xor_n (E.column(k), E.column(k), E.column(pivot), m);
#else
                for (unsigned int l = 0; l < m; l++)
                    E(l, k) ^= E(l, pivot);
#endif
#ifdef  SWAP_PI
                mpn_xor_n (P.column(k), P.column(k), P.column(pivot), m + n);
#else
                for (unsigned int l = 0; l < m + n; l++)
                    P(l, k) ^= P(l, pivot);
#endif
                if (qq->local_delta[pivot] > qq->local_delta[k])
                    qq->local_delta[k] = qq->local_delta[pivot];
	    }
	    for (unsigned int l = 0; l < m; l++)
                E(l, pivot) <<= 1;
	    for (unsigned int l = 0; l < m + n; l++)
                P(l, pivot) <<= 1;
	    qq->delta[pivot] += 1;
            qq->local_delta[pivot]++;
	}
        /* Columns which have not been changed here are those which
         * are magically zero. Count whether this happens often or
         * not.
         *
         * While we're doing this, we can also detect whether our
         * computation stays still. */
        if (!report_spontaneous_zeros(qq, e, is_modified))
            break;
    }

    for (unsigned int i = 0; i < b; i++) {
        for(unsigned int j = 0 ; j < b ; j++) {
            qq->optrs[i * b + j][0] = P(i, j);
        }
    }
    qq->t += e;
    return e;
}
