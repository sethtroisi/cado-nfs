#include <cstdio>
#include "manu.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "balance_arguments.hpp"

#include <boost/cstdint.hpp>
#include <iterator>
#include <vector>
#include <deque>
#include <algorithm>

/* This program does its best to balance for number of cpus that divide
 * nbuckets, where nbuckets is a parameter defaulting to 24 (see
 * balance_arguments.hpp), and can be overridden.
 */
common_arguments common;
balance_arguments mine;

using namespace std;

typedef unsigned int rel;
typedef unsigned int sz;


/* define this to avoid storing the matrix rows in memory. Several input
 * passes are done instead, which makes the program considerably slower
 */
#define	xxxLOWMEM

struct ad_hoc_cmp {
    const vector<vector<rel> >& b;
    unsigned int lim;
    ad_hoc_cmp(const vector<vector<rel> >& b, unsigned int lim)
	: b(b), lim(lim) {}
    bool operator()(const pair<sz, int>& b1,
	    const pair<sz, int>& b2) {
	unsigned int s1 = b[b1.second].size();
	unsigned int s2 = b[b2.second].size();
	if (s1 == s2)	return b1.first > b2.first;
	if (s1 == lim)	return true;
	if (s2 == lim)	return false;
	return b1.first > b2.first;
    }
};

int main(int argc, char * argv[])
{
    ios_base::sync_with_stdio(false);
    cerr.tie(&cout);
    cout.rdbuf()->pubsetbuf(0,0);
    cerr.rdbuf()->pubsetbuf(0,0);

    process_arguments(argc, argv, common, mine);
    unsigned int nr, nc;
    string mstr;
    unsigned int total = 0;
#ifndef	LOWMEM
    vector<matrix_line> lines;
#endif

    vector<pair<sz, rel> > sizes;
    {
	ifstream mtx;

	must_open(mtx, files::matrix);
	get_matrix_header(mtx, nr, nc, mstr);


	cout << "// reading matrix rows" << endl;
#ifndef	LOWMEM
	lines.reserve(nr);
#endif
	istream_iterator<matrix_line> mit(mtx);
	for(rel p = 0 ;
		mit != istream_iterator<matrix_line>() ;
		mit++, p++)
	{
#ifndef	LOWMEM
	    lines.push_back(*mit);
#endif
	    sizes.push_back(make_pair(mit->size(), p));
	    total += mit->size();
	}
	cout << "// ok" << endl;
    }
    BUG_ON(sizes.size() != nr);

    /* Add rows of weight zero so that the matrix is square */
    for(unsigned int i = nr ; i < nc ; i++) {
	sizes.push_back(make_pair(0,i));
#ifndef LOWMEM
	lines.push_back(matrix_line());
#endif
    }
    BUG_ON(sizes.size() != nc);

    /* Now we have the list of row weights, in increasing order */
    sort(sizes.begin(), sizes.end());

    cout << "// dispatching rows into " << mine.nbuckets << " buckets\n";

    vector<vector<rel> > buckets(mine.nbuckets);

    /* This acts as a priority queue indicating which bucket is currently
     * the lightest
     */

    vector<pair<sz, int> > bucket_sizes;
    for(unsigned int i = 0 ; i < mine.nbuckets ; i++) {
	bucket_sizes.push_back(make_pair(0, i));
    }

    ad_hoc_cmp foo(buckets, nc / mine.nbuckets);
    make_heap(bucket_sizes.begin(), bucket_sizes.end(), foo);
    for( ; sizes.size() > nc % mine.nbuckets ; ) {
        /* pick the heaviest row available */
	pair<sz, rel> top = sizes.back();
	sizes.pop_back();
        /* fill the lightest bucket */
	buckets[bucket_sizes[0].second].push_back(top.second);
	bucket_sizes[0].first += top.first;
	pop_heap(bucket_sizes.begin(), bucket_sizes.end(), foo);
	push_heap(bucket_sizes.begin(), bucket_sizes.end(), foo);
    }

    /* Now sort the bucket sizes, and fit the extra (light) rows into the
     * lightest remaining buckets.
     * XXX : Is this really different from continuing the previous loop ?
     */
    sort(bucket_sizes.begin(), bucket_sizes.end());
    for(unsigned int i = 0 ; !sizes.empty() ; i++) {
	pair<sz, rel> top = sizes.back();
	sizes.pop_back();
	buckets[bucket_sizes[i].second].push_back(top.second);
	bucket_sizes[i].first += top.first;
    }

    deque<pair<sz, int> > big_buckets;
    deque<pair<sz, int> > small_buckets;

    for(unsigned int i = 0 ; i < mine.nbuckets ; i++) {
	unsigned int j = bucket_sizes[i].second;
	pair<sz, int> blah(bucket_sizes[i].first, j);

	if (buckets[j].size() == nc / mine.nbuckets) {
	    small_buckets.push_back(blah);
	} else {
	    big_buckets.push_back(blah);
	}
    }
    sort(big_buckets.begin(), big_buckets.end());
    sort(small_buckets.begin(), small_buckets.end());

    double avg = total / (double) mine.nbuckets;
    unsigned int acc = 0;
    /* reorder, keeping track of average weight */
    vector<pair<sz, int> > reorder;
    for(unsigned int i = 0 ; i < mine.nbuckets ; i++) {
	unsigned int s;
	unsigned int r0 = (i * nc) / mine.nbuckets;;
	unsigned int r1 = ((i+1) * nc) / mine.nbuckets;;

	s = r1 - r0;

	deque<pair<sz, int> > * bpool;

	if (s == nc / mine.nbuckets) {
	    bpool = &small_buckets;
	} else {
	    bpool = &big_buckets;
	}
	BUG_ON(bpool->empty());

	if (acc / (double) mine.nbuckets < avg) {
	    /* heaviest */
	    reorder.push_back(bpool->back());
	    acc += bpool->back().first;
	    bpool->pop_back();
	} else {
	    /* lightest */
	    reorder.push_back(bpool->front());
	    acc += bpool->front().first;
	    bpool->pop_front();
	}
    }

    ofstream newmat;
    must_open(newmat, fmt("%.new") % files::matrix);

    /* Pad with zeroes ! */
    put_matrix_header(newmat, nc, nc, mstr);

#ifdef	LOWMEM
    /* read the matrix NBUCKETS times */
    for(unsigned int i = 0 ; i < mine.nbuckets ; i++) {
	unsigned int j = reorder[i].second;
	unsigned int w = reorder[i].first;

	sort(buckets[j].begin(), buckets[j].end());

	cout << fmt("// bucket % weight % nrows %")
	    % j % w % buckets[j].size() << endl;

	ifstream mtx;
	must_open(mtx, files::matrix);

	istream_iterator<matrix_line> mit(mtx);
	matrix_line tmp;
	rel pos = 0;
	if (mstr != string("2")) {
	    unsigned int k;
	    for(k = 0 ; k < buckets[j].size() && buckets[j][k] < nr ; ) {
		ASSERT_ALWAYS(mit != istream_iterator<matrix_line>());
		tmp = *mit++;
		if (pos++ == buckets[j][k]) {
		    newmat << tmp << "\n";
		    k++;
		}
	    }
	    for( ; k < buckets[j].size() ; k++) {
		newmat << "0\n";
	    }
	} else {
	    unsigned int k;
	    for(k = 0 ; k < buckets[j].size() && buckets[j][k] < nr ; ) {
		ASSERT_ALWAYS(mit != istream_iterator<matrix_line>());
		tmp = *mit++;
		if (pos++ == buckets[j][k]) {
		    print_line_without_ones(newmat, tmp) << "\n";
		    k++;
		}
	    }
	    for( ; k < buckets[j].size() ; k++) {
		newmat << "0\n";
	    }
	}
    }
#else	/* ! LOWMEM */
    for(unsigned int i = 0 ; i < mine.nbuckets ; i++) {
	unsigned int j = reorder[i].second;
	unsigned int w = reorder[i].first;

        // This would have the effect of keeping the ordering relatively
        // close to the original ordering in the input file. We disable
        // this for now.
	// sort(buckets[j].begin(), buckets[j].end());

	cout << fmt("// bucket % weight % nrows %")
	    % j % w % buckets[j].size() << endl;
	if (mstr != string("2")) {
	    for(unsigned int k = 0 ; k < buckets[j].size() ; k++) {
		newmat << lines[buckets[j][k]] << "\n";
	    }
	} else {
	    for(unsigned int k = 0 ; k < buckets[j].size() ; k++) {
		print_line_without_ones(newmat, lines[buckets[j][k]]) << "\n";
	    }
	}
    }
#endif

    string oldname = fmt("%.old") % files::matrix;
    string newname = fmt("%.new") % files::matrix;
    string curname = files::matrix;

    if (rename(curname.c_str(), oldname.c_str()) != 0) {
	BUG();
    }

    if (rename(newname.c_str(), curname.c_str()) != 0) {
	BUG();
    }

    return 0;
}
/* vim: set sw=4: */
