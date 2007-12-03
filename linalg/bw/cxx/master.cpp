#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"
#include "gmp-hacks.h"

#include "addmul.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "detect_params.hpp"
#include "master_arguments.hpp"
#include "matrix_header.hpp"
#include "files.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "state.hpp"
#include "ticks.hpp"
#include "elt.hpp"

#include <vector>
#include "matrices.hpp"
#include "matrices_printing.hpp"

using namespace std;
using namespace boost;

typedef unsigned int uint;

master_arguments mine;

namespace globals {
	uint nr;
	uint nb_iter;
	uint m, n, b;
	mpz_class modulus;
}

/* silly helpers */
std::ostream& operator<<(std::ostream& o, const std::vector<uint>& d)
{
	o << "[ ";
	std::copy(d.begin(), d.end(), ostream_iterator<uint>(o, " "));
	o << "]";
	return o;
}
uint max(const std::vector<uint>& d)
{
	ASSERT(!d.empty());
	return *std::max_element(d.begin(), d.end());
}

struct m_dim { inline uint size() const { return globals::m; } };
struct n_dim { inline uint size() const { return globals::n; } };
struct b_dim { inline uint size() const { return globals::m + globals::n; } };
struct v_dim {
	uint v;
	inline uint& ncoeffs() { return v; } 
	inline uint size() const { return v; }
	v_dim() : v() {}
	v_dim(const v_dim& x) : v(x.v) {}
};

typedef cons<m_dim, void> m_;
typedef cons<n_dim, void> n_;
typedef cons<b_dim, n_> bn_;
typedef cons<b_dim, m_> bm_;
typedef cons<n_dim, m_> nm_;
typedef cons<m_dim, m_> mm_;
typedef cons<m_dim, n_> mn_;
typedef cons<v_dim, nm_> vnm_;	// for A
typedef cons<v_dim, bn_> vbn_;	// for F

/* Short note on the void * conversions : 0 is the only constant having
 * automatic conversion to void * ; I know it's dirty...
 */
template<typename Dims>
struct row_ops<row<elt, Dims > > {
	typedef row<elt, Dims > final;
	final& zero() {
		final& me((final&)*this);
		memset(me.base_ptr(), 0, me.total_size() * sizeof(elt));
		return me;
	}
	final& operator=(void * x) {
		ASSERT(x == NULL);
		return zero();
	}
	bool is_zero() const {
		final& me((final&)*this);
		const elt * p = me.base_ptr();
		for(uint x = 0 ; x < me.total_size() ; x++) {
			// XXX elt::operator!=
			if (p[x] != 0)
				return false;
		}
		return true;
	}
	bool operator==(void * x) const {
		ASSERT(x == NULL);
		return is_zero();
	}
	bool operator!=(void * x) const {
		ASSERT(x == NULL);
		return !is_zero();
	}
};

typedef row<elt, m_> m_row;
typedef row<elt, mn_> mn_mat;
typedef row<elt, mm_> mm_mat;
typedef row<elt, bm_> bm_mat;
typedef row<elt, vnm_> vnm_mat;
typedef row<elt, vbn_> vbn_mat;

typedef row<const elt, m_> c_m_row;
typedef row<const elt, mn_> c_mn_mat;
typedef row<const elt, mm_> c_mm_mat;
typedef row<const elt, bm_> c_bm_mat;
typedef row<const elt, vnm_> c_vnm_mat;
typedef row<const elt, vbn_> c_vbn_mat;

uint detect_ncoeffs()
{
	mpz_class t;
	using namespace globals;

	vector<uint> report;
	for(uint i = 0 ; i < m ; i++) {
		for(uint j = 0 ; j < n ; j++) {
			ifstream a((files::a % i % j).c_str());
			if (!a.is_open()) {
				report.push_back(0);
				continue;
			}
			uint v;
			for(v = 0 ; a >> t ; v++);
			report.push_back(v);
		}
	}
	BUG_ON(report.empty());

	return *min_element(report.begin(), report.end());
}

template<typename T, typename E>
void norm(row<E, cons<T, void> > r, mpz_class& m, uint c)
{
	// elt:: unary inversion to size n.
	// elt:: operator* to size n.
	// assign(m, r[c]);
	m = inv(r[c]);
	r[c] = 1;
	mul(r, m, c + 1);
}
template<typename T, typename E>
void norm(row<E, cons<T, void> > r, mpz_class& m)
{
	uint c;
	for(c = 0 ; c < r.size() && r[c] == 0; c++);
	if (c == r.size()) return 0;
	norm(r, m, c);
}
template<typename T, typename E>
void mul(row<E, cons<T, void> > r,
		const mpz_class& mul, uint c = 0)
{
	// XXX elt::operator*= n * n -> n
	for(uint i = c ; i < r.size() ; i++) {
		r[i] *= mul;
		r[i]  = modred(r[i]);
	}
}
template<typename T, typename E1, typename E2>
void submul(row<E1, cons<T, void> > r,
		row<E2, cons<T, void> > p,
		const mpz_class& m,
		uint c = 0)
{
	// XXX elt::operator- n * 2n -> n
	// XXX p * q -> n looks fair as well.
	for(uint i = c ; i < r.size() ; i++) {
		r[i] -= m * p[i];
	}
}

template<typename T, typename E1, typename E2>
void kill(row<E1, cons<T, void> > r,
		row<E2, cons<T, void> > p,
		mpz_class& mul,
		uint c = 0)
{
	ASSERT(p[c] == 1);
	assign(mul, r[c]);
	r[c] = 0;
	submul(r, p, mul, c + 1);
}

/* Initialize f from [Thome2002JSC, section 3.3].
 *
 * Except that the code here is transposed, so that we work with rows
 * instead
 */
uint f_init(vbn_mat& f_series, const vnm_mat& a_series)
{
	/* This is first done by computing the integer s such that the
	 * rows of a[0] ... a[s-1] span the full m-dimensional space
	 */
	uint rank = 0;
	uint s = 0;
	uint j = 0;
	using namespace globals;

	holder<mm_mat> scrap;
	vector<pair<int, pair<int, int> > > pivots(m,
			make_pair(-1, make_pair(0,0)));

	for( ; s < a_series.size() && rank < m ; ) {
		/* Try to insert row j of a[s] */
		scrap[rank] = a_series[s][j];
		
		/* kill with existing rows */
		uint c;
		mpz_class mul;
		for(c = 0 ; c < m ; c++) {
			if (scrap[rank][c] == 0) continue;
			int p = pivots[c].first;
			if (p == -1) {
				norm(scrap[rank], mul, c);
				break;
			}
			kill (scrap[rank], scrap[p], mul, c);
		}
		if (c < m) {
			/* We've found a new pivot */
			cout << fmt("// a_%, row % -> rank++\n") % s % j;
			cout << flush;
			pivots[c] = make_pair(rank, make_pair(s, j));
			rank++;
		}
		if (rank == m)
			break;
		if (++j == m) {
			j = 0;
			s++;
		}
	}
	if (rank < m) {
		cerr << "// singular data\n";
		exit(1);
	}

	s++;

	BUG_ON((s + 1) >= f_series.size());

	for(uint i = 0 ; i <= s ; i++) f_series[i] = 0;

	/* The first n rows are the identity matrix */
	for(uint i = 0 ; i < n ; i++) {
		f_series[0][i][i] = 1;
	}

	/* The others are X^(s-i_k)r_k */
	for(uint i = 0 ; i < m ; i++) {
		int ik = pivots[i].second.first;
		int rk = pivots[i].second.second;
		f_series[s - ik][n + i][rk] = 1;
	}

	return s;
}

/* {{{ proxy class to print slices of matrix polynomials */
template<class Dims>
struct seq {
	row<const elt, cons<v_dim, Dims> > f;
	uint nc;
	seq(row<const elt, cons<v_dim, Dims> > f, uint nc) : f(f), nc(nc) {}
};

template<class Dims>
std::ostream& operator<<(std::ostream& o, const seq<Dims>& aa)
{
	if (test_magma(o)) {
		for(uint v = 0 ; v < aa.nc ; v++) {
			if (v) o << " + ";
			o << "\nX^" << v << " *\n" << set_magma << aa.f[v];
		}
	} else {
		copy(aa.f.begin(), aa.f.begin() + aa.nc - 1,
		std::ostream_iterator<row<const elt, Dims> >(o, "\n\n"));
		// XXX operator<<(ostream&, const elt&)
		o << aa.f[aa.nc - 1];
	}
	return o;
}
/* }}} */

/* Normally, the e matrix on input is filled with something: the result
 * of the transformation undergone during the latest iteration. This
 * amounts necessarily to putting e in (slightly scrambled) row echelon
 * form). It turns out that in this precise case, the m non-zero rows
 * (forming a submatrix of full rank) are precisely the ones we want to
 * keep because the corresponding rows of f have been shifted by X,
 * whereas the rows equal to zero must be recomputed for the next value
 * of t. Hence the test against 0.
 */
std::vector<uint> ctfa(row<elt, bm_> dst, uint t,
		row<elt, vbn_> const & f,
		row<elt, vnm_> const & a,
		const std::vector<uint>& degnom)
{
	vector<uint> c;
	for(uint i = 0 ; i < dst.size() ; i++) {
		if (dst[i] != 0)
			continue;
		for(uint j = 0 ; j < dst[i].size() ; j++) {
			elt::affine<1, 0>::type u, v;
			elt::affine<2, 1>::type res;

			res = 0;
			for(uint s = t - degnom[i] ; s <= t ; s++) {
				for(uint k = 0 ; k < globals::n ; k++) {
					// u = f[t-s][i][k];
					// v = a[s][k][j];
					// XXX elt::operator* to size 2n
					// XXX elt::+= 2n to 2n+1
					// res += u * v;
					res += f[t-s][i][k] * a[s][k][j];
				}
			}
			// XXX elt::operator=(size 2n+1)
			dst[i][j] = res;
		}
		// XXX elt::operator==
		if (dst[i] == 0)
			c.push_back(i);
	}
	return c;
}

void mul_rows(vbn_mat& f, uint j0, uint d, const mpz_class& mul)
{
	mpz_class x;
	for(uint v = 0 ; v <= d ; v++) {
		for(uint i = 0 ; i < globals::n ; i++) {
			assign(x, f[v][i][j0]);
			x *= mul;
			// XXX elt::operator* to size n (including reduction)
			f[v][i][j0] = modred(x);
		}
	}
}

/* Returns x[] such that d[x[0]] ... d[x[nd-1]] is sorted */
vector<uint> order(const vector<uint>& d)
{
	vector<pair<uint, uint> > ordering;
	for(uint i = 0 ; i < d.size() ; i++) {
		ordering.push_back(make_pair(d[i],i));
	}
	sort(ordering.begin(), ordering.end());
	vector<uint> res;
	res.reserve(d.size());
	for(uint i = 0 ; i < d.size() ; i++) {
		res.push_back(ordering[i].second);
	}
	return res;
}

/* This is the naive quadratic version */
void algo1(uint t, vector<uint>& degnom, bm_mat& e,
		 vbn_mat f,
		 c_vnm_mat a)
{
	/* gaussian elimination on the rows */
	vector<uint> o = order(degnom); /* sort the rows by degnom order */
	
	uint i, j, rank = 0;
	
	for(i = 0 ; rank != globals::m && i < globals::b ; i++) {
		for(j = 0 ; j < globals::m ; j++) {
			if (e[o[i]][j] != 0) break;
		}
		if (j == globals::m) continue;
		rank++;

		/* normalize this row */
		mpz_class m;

		norm(e[o[i]], m, j);

		/* adjust in f also */
		for(uint k = 0 ; k <= degnom[o[i]] ; k++) {
			mul(f[k][o[i]], m);
		}

		/* cancel in other rows as well */
		for(uint ii = i + 1 ; ii < globals::b ; ii++) {
			kill(e[o[ii]], e[o[i]], m, j);
			for(uint k = 0 ; k <= degnom[o[i]] ; k++) {
				submul(f[k][o[ii]], f[k][o[i]], m);
			}
		}
	}
	BUG_ON (rank != globals::m);

	for(i = 0 ; i < globals::b ; i++) {
		/* if e[i] is zero, the corresponding row will be
		 * recomputed entirely from the next ctfa step.
		 * Otherwise, it will be kept */

		if (e[i] == 0) continue;

		rank--;		// bookkeeping
		
		/* shift f by X */
		for(uint k = 0 ; k <= degnom[i] ; k++) {
			f[degnom[i]-k+1][i] = f[degnom[i]-k][i];
		}
		f[0][i] = 0;
		degnom[i]++;
	}
	BUG_ON(rank != 0);
}

#if 0
/* Consider ec mod t^l ; update the pi matrix accordingly */
void compute_pi(polymat_t& pi, polymat_t& ec, vector<uint>& degnom, uint l)
{
}
#endif

/* This wrapper saves me a bit of typing */
struct trace {
	double ticks0;
	double wct0;

	trace() : ticks0(oncpu_ticks()), wct0(wallclock_ticks()) {}
	template<typename T>
	bool step(uint t, const T& x, bool force = false) {
		if (!force && (t % 10 != 0))
			return true;
		double delta = oncpu_ticks() - ticks0;
		double wdelta = wallclock_ticks() - wct0;
		double pcpu = (wdelta > 0.1) ? delta / wdelta : 0.0;

		cout << fmt("// step % %[F.2]s %[F.1]%% %\n")
			% t % delta % (100.0 * pcpu) % x; 
		cout << flush;
		return true;
	}
};

int main(int argc, char *argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	common_arguments common;

	process_arguments(argc, argv, common, mine);

	using namespace globals;

	{
		string mstr;
		ifstream mtx;
		must_open(mtx, files::matrix);
		get_matrix_header(mtx, globals::nr, mstr);
		globals::modulus = mpz_class(mstr);
		if (SIZ(globals::modulus.get_mpz_t()) != MODULUS_SIZE) {
			cerr << fmt("ERROR: RECOMPILE WITH"
					" ``#define MODULUS_SIZE %''\n")
				% SIZ(globals::modulus.get_mpz_t());
			exit(1);
		}
	}

	detect_mn(m, n);
	cout << fmt("// detected m = %\n") % m;
	cout << fmt("// detected n = %\n") % n;
	if (m <= 0 || n <= 0) {
		cout << flush;
		BUG();
	}

	b = m + n;

	nb_iter = nr + 2 * n;
	nb_iter = iceildiv(nb_iter, m) + iceildiv(nb_iter, n) + 10;

	uint tmax = detect_ncoeffs();
	
	cout << fmt("// detected tmax = %\n") % tmax;
	
	if (tmax < 2 * n) {
		cerr << "// insufficient amount of data for A[[X]]" << endl;
		exit(1);
	}
	if ((uint) tmax < nb_iter + 1) {
		cerr << "// WARNING: incomplete data" << endl;
	}

	vnm_ a_size;
	vbn_ f_size;

	/* We have tmax coefficients available. Assume that it is enough,
	 * meaning that s and d are such that (cf theorem 3.2 and section
	 * 5 of [Thome2002JSC]) :
	 * tmax >= s + ceil((m+n)/n d)
	 * d <= (tmax - s) * n / (m+n)
	 *
	 * Assume that s <= ceil(m/n) + 2 (safe guess)
	 *
	 * then d <= (tmax - ceil(m/n) - 2) * n / (m+n)
	 * m/n * d <= (tmax - ceil(m/n) - 2) * mn / (n * (m+n))
	 * ceil(m/n * d) <= ceil((tmax - ceil(m/n) - 2) * mn / (n * (m+n)))
	 *
	 * Now the expected typical values for d and s are ceil(N/m) and
	 * ceil(m/n). This means that the output generator has degree
	 * bounded by the expression below for expected_deg_f
	 */
	uint expected_deg_f;
	{
		int exp_s = iceildiv(m, n);
		int bound_s = exp_s + 2;
		int exp_d = iceildiv(nr, m);
		int bound_mnd = iceildiv((tmax - bound_s) * m * n, n * (m + n));

		a_size.ncoeffs() = tmax;
		f_size.ncoeffs() = 1 + bound_s + bound_mnd;
		expected_deg_f = exp_s + iceildiv(m * exp_d, n);
	}

	cout << "// Assuming deg(F(X)) will not exceed "
		<< (f_size.ncoeffs() - 1)
		<< " (expected " << expected_deg_f << ")\n";

	holder<vnm_mat> a(a_size);
	holder<vbn_mat> f(f_size);

	cout << "// Reading data for A(X)\n" << flush;
	for(uint i = 0 ; i < n ; i++) {
		for(uint j = 0 ; j < m ; j++) {
			ifstream ax;
			mpz_class t;
			must_open(ax, files::a % j % i);

			/* NOTE : we drop the first coefficient, because
			 * we mean to work with the sequence generated
			 * on x and By, so that we obtain a generator
			 * afterwards.
			 */
			ax >> t;

			for(uint v = 0 ; v < tmax ; v++) {
				ax >> t;
				a[v][i][j] = t;
			}
		}
	}

	uint t0 = f_init(f, a);
	vector<uint> degnom(b, t0);

	cout << fmt("// t0 = %\n") % t0;
	BUG_ON(t0 > iceildiv(m, n) + 2);

	/* This triggers a complete initialization of E for the first
	 * iteration */
	holder<bm_mat> e;
	e = 0;

	cout << flush;
	cerr << flush;

	trace foo;
	
	for(uint t = t0 ; t < tmax ; t++) {

		/* First, update e */
		vector<uint> chances = ctfa(e, t, f, a, degnom);

		if (chances.size()) {
			cout << fmt("// step %[] LOOK %\n") % t % chances;
			cout << flush;
		}

		/* Make sure any row, even shifted, will fit */
		if (2 + max(degnom) > f.size()) {
			/* if not, we may have hit the point where rows
			 * grow unevenly. Which means we have a
			 * generator. Make sure this is the case, and
			 * bail out.  */
			BUG_ON(chances.empty());
			foo.step(t, degnom);
			break;
		}

		algo1(t, degnom, e, f, a);
		foo.step(t, degnom);
	}
	foo.step(tmax, degnom, true);
	/* No ``LOOK'' message is output for the last step, this is
	 * because we do not have enough info to compute it. Therefore
	 * the info to be relied on is for steps tmax-1 and before... */

	{
		cout << "// Writing F(X) to disk" << endl;
		for(uint i = 0 ; i < b ; i++) {
			ofstream fx;
			must_open(fx, files::f % i);
			for(uint t = 0 ; t <= degnom[i] ; t++) {
				fx << f[t][i] << "\n";
			}
		}
	}

	return 0;
}
