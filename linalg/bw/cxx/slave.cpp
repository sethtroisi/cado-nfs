#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "barrier.h"
#include "manu.h"
#include "gmp-hacks.h"

#include "addmul.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "detect_params.hpp"
#include "slave_arguments.hpp"
#include "files.hpp"
#include "matrix.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "state.hpp"
#include "threads.hpp"
#include "ticks.hpp"
#include "preconditioner.hpp"

// #include "matmul.hpp"

#include <boost/cstdint.hpp>
#include <iterator>

#include "traits_globals.hpp"

/* TODO: maybe remove checks entirely in characteristic two, as they are
 * hardly useful... Most probably something more sensible would be
 * necessary.
 */

slave_arguments mine;

using namespace std;
using namespace boost;

typedef unsigned int uint;

void load_x_vectors(vector < uint32_t > &bw_x) /*{{{ */
{
    for (int i = 0;; i++) {
	ifstream xvec;
	uint32_t x;
	xvec.open((files::x % i).c_str());
	if (!xvec.is_open()) {
	    /* This is enough to guess the value of m */
	    break;
	}
	if (!(xvec >> wanted('e') >> x)) {
	    die("Parse error reading %s\n", 1, (files::x % i).c_str());
	}
	bw_x.push_back(x);
    }
}				/*}}} */

namespace globals {
    uint m, n;
    uint col;
    uint nr, nc;
    uint nb_iter;
    uint cp_lag;
    uint iter0;
    int nbys;
    mpz_class modulus;
    uint8_t modulus_u8;
    uint16_t modulus_u16;
    uint32_t modulus_u32;
    unsigned long modulus_ulong;
     vector < uint32_t > bw_x;
    uint degf;			// mksol only
    barrier_t main_loop_barrier;
    thread_lock_t console_lock = THREAD_LOCK_INITIALIZER;

     vector < mpz_class > dotprod_diff_check;

    /* only for threads */
     vector < std::streampos > mtxfile_pos;
     vector < uint > nb_coeffs;
     vector < void *>thread_class_ptr;

     template < typename traits > struct vdata {
	typedef typename traits::scalar_t scalar_t;
	scalar_t *v0;
	scalar_t *f;
    };

}
/* {{{ I/O : reading tasks */

void fill_check_vector(uint i0, uint i1,
        int32_t * p,
        const std::string & name)
{
    ifstream f;
    must_open(f, name);
    int32_t x;
    istream_iterator < int32_t > fit(f);

    i1 -= i0;
    for (; i0 && fit != istream_iterator < int32_t > (); x = *fit++, i0--);

    for (; i1 && fit != istream_iterator < int32_t > (); *p++ = *fit++, i1--);

    BUG_ON(i1 != 0);
}

/* }}} */

/* {{{ timing helpers */
bool display_decision(int n, int p)
{
    int per[] = {
	1, 2, 5,
	10, 20, 50,
	100, 200, 500,
	1000, 2000, 5000,
	10000, 20000, 50000,
	100000, 200000, 500000,
	1000000, 2000000, 5000000
    };

    if (n <= p)
	return false;

    for (uint i = 0; i < sizeof(per) / sizeof(per[0]); i++) {
	if (n > p + per[i])
	    continue;
	if (n - p == per[i] - (p % per[i]))
	    return true;
    }
    return false;
}

string pdate(double v)
{
    struct tm toto;
    char s[80];
    time_t t = (unsigned long int) v;

    localtime_r(&t, &toto);
    strftime(s, 80, "%b %d %H:%M", &toto);
    /* sprintf(s+strlen(s),".%04d",tv->tv_usec/10000); */
    return string(s);
}

string pdelta(double t)
{
    const char units[] = "mhdw";
    int val[] = { 60, 60, 24, 7 };
    char t_unit = 's';
    for (int i = 0; i < 5 && t > val[i]; i++) {
	t_unit = units[i];
	t /= (double) val[i];
    }
    return fmt("%[F.2]%") % t % t_unit;
}

/* }}} */

template < typename > int addup(int);

/* The Derived type supplies the use_vector part */
template < typename traits, typename Derived >
    // typename traits = typical_scalar_traits<width> >
struct thread:public traits {
    typedef typename traits::scalar_t scalar_t;
    typedef typename traits::wide_scalar_t wide_scalar_t;
    typedef thread < traits, Derived > self;
    int t;
    typename traits::representation::matrix_rowset mat;
    /*
       uint32_t * idx;
       int32_t  * val;
     */
    uint i0;
    uint i1;
    scalar_t *v;
    scalar_t *w;
    wide_scalar_t *scrap;
    int32_t *check_x0;
    int32_t *check_m0;
    double maxwait;
    uint done;
     vector < mpz_class > dot_part;
    /* reference timing */
    uint go_mark;
    double ticks_ref;
    double wct_ref;
     preconditioner < traits > precond;

     thread(void *ptr) {
	using globals::nc;
	using globals::nr;
	using globals::thread_class_ptr;
	using globals::nb_coeffs;
	using globals::mtxfile_pos;

	 globals::vdata < traits > *vptr = (globals::vdata < traits > *)ptr;

	 t = tseqid();

	 BUG_ON(thread_class_ptr.size() <= (uint) t);
	 thread_class_ptr[t] = (void *) this;

	 i0 = (t * nr) / mine.nt;
	 i1 = ((t + 1) * nr) / mine.nt;
	 maxwait = 0.1;

	 v = new scalar_t[nc];
	 w = new scalar_t[nc];
	 scrap = new wide_scalar_t[i1 - i0];

	 mat.alloc(i1 - i0, nb_coeffs[t]);
	 check_x0 = new int32_t[i1 - i0];
	 check_m0 = new int32_t[i1 - i0];

	 dot_part.assign(globals::nbys, mpz_class());

	thread_lock(&globals::console_lock);
	cout << fmt("T% reading %")
            % t % files::matrix << endl;
	thread_unlock(&globals::console_lock);
	{
            ifstream mtx;
	    must_open(mtx, files::matrix);
	    mat.fill(mtx, mtxfile_pos[t], i0, i1);
	}
	thread_lock(&globals::console_lock);
	cout << fmt("T% reading % and %")
            % t % files::x0 % files::m0 << endl;
	thread_unlock(&globals::console_lock);
	fill_check_vector(i0, i1, check_x0, files::x0);
	fill_check_vector(i0, i1, check_m0, files::m0);

	thread_lock(&globals::console_lock);
	cout << fmt("T% done reading input data") % t << endl;
	thread_unlock(&globals::console_lock);
        /* 0 and nr are for the locations affected in the destination
         * vector. As we're only slicing the matrix horizontally so far,
         * it's pointless 
         */
	precond.init(0, nr, files::precond);

	BUG_ON(globals::iter0 & 1);

	// memcpy(v, vptr->v0, nc * sizeof(scalar_t));
	// memset(w, 0, nc * sizeof(scalar_t));
	traits::copy(v, vptr->v0, nc);
	traits::zero(w, nc);

	go_mark = done = globals::iter0;
	ticks_ref = thread_ticks();
	wct_ref = wallclock_ticks();
    }

    ~thread() {
	delete[]v;
	delete[]w;
	delete[]scrap;
	delete[]check_x0;
	delete[]check_m0;
	/* the matrix rowset is deleted by the dtor */
    }

    void display_stats(bool sync = true) {
	if (!display_decision(done, go_mark))
	    return;

	double ticks_diff = thread_ticks() - ticks_ref;
	double wct_diff = wallclock_ticks() - wct_ref;
	double av;
	double m0_estim;
	double tw;
	double rem;
	double pcpu;

	pcpu = ticks_diff / wct_diff;
	av = ticks_diff / (double) (done - go_mark);
	m0_estim = av / (double) globals::nb_coeffs[t];
	tw = globals::nb_iter * av;
	rem = av * (globals::nb_iter - done);
	string eta;

	string e1 = pdate(wct_ref + rem);
	string e2 = pdate(wct_ref + rem / pcpu);
	if (e1 == e2) {
	    eta = e1;
	} else {
	    if (e1.compare(0, 7, e2, 0, 7) == 0) {
		e2.erase(0, 7);
	    }
	    eta = fmt("% .. %") % e1 % e2;
	}

	thread_lock(&globals::console_lock);
	cout << fmt("T% N=% av=% M0=%[.3]s %[F.1]%% tw=% eta=<%>")
	    % t % done % pdelta(av)
	    % m0_estim % (pcpu * 100.0)
	    % pdelta(tw)
	    % eta << endl;
	thread_unlock(&globals::console_lock);
    }

    void one_dotprod(scalar_t &r, const int32_t * sv,
		     const scalar_t * lv) {
	wide_scalar_t tmp;
	zero(tmp);
	int acc = 0;

	for (uint i = i0; i < i1; i++) {
	    addmul(tmp, lv[i], sv[i - i0]);
	    if (++acc == traits::max_accumulate) {
		reduce(tmp, tmp);
		acc = 1;
	    }
	}
        addmul(tmp, r, 1);
	reduce(r, tmp);
    }

    inline void multiply(wide_scalar_t * dst, scalar_t * src) {
	precond(src);
	mat.template mul < traits > (dst, src);
    }

    void note_waited(double twait) {
	if (twait > maxwait) {
	    cerr << fmt("// thread % : waited %[F.2]s")
		% t % (maxwait = twait) << endl;
	}
    }

    /* src is expected to be B^done * y ; compute B^(done+1) y into dst,
     * increment i */
    void flip_flap(scalar_t * dst, scalar_t * src) {
	double twait;
	bool ok;
	using namespace globals;

	ASSERT(dst == ((done & 1) ? v : w));
	ASSERT(src == ((done & 1) ? w : v));

        scalar_t ds, dt;
        traits::zero(ds);
        traits::zero(dt);

	do {
	    multiply(scrap, src);

	    one_dotprod(ds, check_m0, src);
	    reduce(dst, scrap, i0, i1);
            traits::zero(dst + nr, nc - nr);
	    one_dotprod(dt, check_x0, dst);

	    // traits::subtract(dot_part, ds, dt);
	    for (int i = 0; i < nbys; i++) {
		dot_part[i] = traits::get_y(ds,i) - traits::get_y(dt,i);
	    }

	    ok = barrier_wait(&main_loop_barrier, &twait, &addup < self >);
	    note_waited(twait);
	} while (!ok);

	/* read from the other threads the w segments we don't
	 * have */
	for (int x = 0; x < mine.nt; x++) {
	    if (x == t)
		continue;
	    uint xi0 = (x * nr) / mine.nt;
	    uint xi1 = ((x + 1) * nr) / mine.nt;
	    const self *xptr = (const self *)
		globals::thread_class_ptr[x];
	    /*
	       memcpy(dst + xi0,
	       ((done & 1) ? xptr->v : xptr->w) + xi0,
	       (xi1 - xi0) * sizeof(scalar_t));
	     */
	    traits::copy(dst + xi0,
			 ((done & 1) ? xptr->v : xptr->w) + xi0, xi1 - xi0);
	}

	++done;

	((Derived &) * this).use_vector(dst);

	display_stats();
    }

    void loop() {
	using namespace globals;

	/* Check in (x, B^i0 y) */
	((Derived &) * this).use_vector(v);

	for (; done < nb_iter;) {
	    flip_flap(w, v);
	    if (done == nb_iter)
		break;
	    flip_flap(v, w);
	    if (done == nb_iter)
		break;
	}
	/* Note that ideally, we should put this in the dtor, but
	 * this is something I don't like to do.
	 *
	 * It's critical to have it because flip_flap's tail
	 * relies on all threads to be present. */
	barrier_wait(&main_loop_barrier, NULL, NULL);
    }
};

/* do the work that differs when task == slave and task == mksol */
template < typename traits >
    struct slave_thread:public thread < traits, slave_thread < traits > > {
    typedef slave_thread < traits > self;
    typedef thread < traits, self > super;
    typedef typename traits::scalar_t scalar_t;
    typedef typename traits::wide_scalar_t wide_scalar_t;
    ofstream *a;
     vector < uint > px;
    int cp_lag;
     slave_thread(void *x):super(x) {
	int x0 = (super::t * globals::m) / mine.nt;
	int x1 = ((super::t + 1) * globals::m) / mine.nt;
	 a = new ofstream[x1 - x0];
	for (int x = x0; x < x1; x++) {
	    ofstream v;
	    string nm = files::a % x % globals::col;
	     a[x - x0].open(nm.c_str(),
			    ios_base::out | ios_base::ate | ios_base::app);
	     ASSERT_ALWAYS(a[x - x0].is_open());
	     px.push_back(globals::bw_x[x]);
	} cp_lag = globals::cp_lag;
    }
    ~slave_thread() {
	delete[]a;
    }
    void use_vector(scalar_t * dst) {
	// mpz_class z;

	for (uint x = 0; x < px.size(); x++) {
	    // MPZ_SET_MPN(z.get_mpz_t(), dst[px[x]], width);
	    // a[x] << z << "\n";
	    traits::print(a[x], dst[px[x]]) << "\n";
	}

	if (super::done == 200) {
	    /* exceptional case ; only in order to be able to
	     * check early that the system is doable */
	    flush();
	}
#ifndef NDEBUG
	if (super::done < 20) {
	    flush();
	}
#endif

	if (super::done % cp_lag != 0)
	    return;

	uint r = super::done / cp_lag;

	/* There's no point in flushing the starting rev. */
	if (r == 0)
	    return;

	/* writeout w ; this must be done by only
	 * one thread. It's quite silly to use barrier_wait for
	 * this, a simpler concept would clearly do, but that's
	 * what we have at hand for lowest pain... */
	int rc = barrier_wait(&globals::main_loop_barrier, NULL, NULL);
	if (rc == 1) {
	    {
		thread_lock(&globals::console_lock);
		std::string nm0 = files::v % (globals::col) % r;
                if (global::nbys > 1) {
                    std::string nm1 =
                        files::v % (globals::col + globals::nbys - 1) % r;
                    cout << fmt("T% writes % to %")
                        % super::t % nm0 % nm1 << endl;
                } else {
                    cout << fmt("T% writes % to %")
                        % super::t % nm0 << endl;
                }

		thread_unlock(&globals::console_lock);
	    }
	    for (int i = 0; i < globals::nbys; i++) {
		std::string nm;
		nm = files::v % (globals::col + i) % r;
		ofstream v;
		must_open(v, nm);
		for (uint j = 0; j < globals::nc; j++) {
		    v << traits::get_y(dst[j], i) << "\n";
		}
	    }
	}

	flush();
    }
    void flush() {
	/* flush a */
	for (uint x = 0; x < px.size(); x++) {
	    a[x].flush();
	}
    }
};

template < typename traits >
    struct mksol_thread:public thread < traits, mksol_thread < traits > > {
    typedef mksol_thread < traits > self;
    typedef thread < traits, self > super;
    typedef typename traits::scalar_t scalar_t;
    typedef typename traits::wide_scalar_t wide_scalar_t;
    int cp_lag;
    uint bound_deg_f;
    uint degf;
    scalar_t *sum;
    scalar_t *fptr;
    int accumulate_wide;
     mksol_thread(void *ptr):super(ptr) {
	using namespace globals;
	 fptr = ((globals::vdata < traits > *)ptr)->f;

	 sum = new scalar_t[super::i1 - super::i0];
	// memset(sum, 0, (super::i1 - super::i0) * sizeof(scalar_t));
	 traits::zero(sum, super::i1 - super::i0);
	 cp_lag = globals::cp_lag;
	 accumulate_wide = 0;
    } ~mksol_thread() {
	delete[]sum;
    }
    void use_vector(scalar_t * dst) {
	/* We have the full vector available here. However it is
	 * only available for writing */

	using globals::degf;

	for (uint i = super::i0; i < super::i1; i++) {
//                      mp_limb_t t[2 * width + 1];
//                      mp_limb_t c;
//                      mpn_mul_n(t, dst[i], f[degf - super::done], width);
//                      c = mpn_add_n(t, t, sum[i - super::i0], width);
//                      if (c)
//                              c = mpn_add_1(t + width, t + width, width, c);
//                      t[2 * width] = c;
// 
//                      core_ops::assign<width, 2 * width + 1>(sum[i - super::i0], t);
// 
	    traits::addmul_wide(sum[i - super::i0], dst[i],
				fptr[degf - super::done]);
	}
	if (++accumulate_wide >= traits::max_accumulate_wide) {
	    for (uint i = super::i0; i < super::i1; i++) {
		traits::reduce(sum[i - super::i0], sum[i - super::i0]);
	    }
	}

	if (super::done % cp_lag != 0)
	    return;
	/* There's no point in flushing the starting rev. */
	if (super::done == 0)
	    return;

	flush(true);
    }

    void flush(bool priv = false) {

	uint r = super::done / cp_lag;

	if (!priv)
	    r++;

	for (uint i = super::i0; i < super::i1; i++) {
	    traits::reduce(sum[i - super::i0], sum[i - super::i0]);
	}
	accumulate_wide = 1;

	/* There is potential that the last call here, which
	 * comes right before the dtor, is exactly after a
	 * ``normal'' call at a plain checkpoint. In such a case,
	 * the trailing fxy vector will be written out as a bunch
	 * of zeroes, which should not matter */

	/* One thread is responsible for writing to disk the
	 * accumulated sum up to this cp value */
	int rc = barrier_wait(&globals::main_loop_barrier, NULL, NULL);
	if (rc == 1) {
	    std::string nm = files::fxy % mine.scol % globals::col % r;
	    cout << fmt("T% writes %") % super::t % nm << endl;
	    ofstream v;
	    must_open(v, nm);
	    mpz_class z;
	    for (int x = 0; x < mine.nt; x++) {
		uint xi0 = (x * globals::nr) / mine.nt;
		uint xi1 = ((x + 1) * globals::nr) / mine.nt;
		const self *xptr = (const self *)
		    globals::thread_class_ptr[x];
		for (uint xi = xi0; xi < xi1; xi++) {
		    // MPZ_SET_MPN(z.get_mpz_t(),
		    // xptr->sum[xi - xi0],
		    // width);
		    // v << z << "\n";
		    traits::print(v, xptr->sum[xi - xi0]) << "\n";
		}
	    }
	}
	barrier_wait(&globals::main_loop_barrier, NULL, NULL);
	// memset(sum, 0, (super::i1 - super::i0) * sizeof(scalar_t));
	traits::zero(sum, super::i1 - super::i0);
    }
};

template < typename T > int addup(int k)
{
    using namespace globals;

    if (k == 0) {
	for (int i = 0; i < nbys; i++) {
	    dotprod_diff_check[i] = 0;
	}
    }
    const T & me = *(const T *) (thread_class_ptr[tseqid()]);

    for (int i = 0; i < nbys; i++) {
	dotprod_diff_check[i] += me.dot_part[i];
    }

    if (k != mine.nt - 1)
	return 1;

    for (int i = 0; i < nbys; i++) {
	dotprod_diff_check[i] %= globals::modulus;
	if (dotprod_diff_check[i] != 0) {
	    cerr << fmt("// s[%] = %, not 0 !!!\n")
		% i % dotprod_diff_check[i];
	    return 0;
	}
    }
    return 1;
}

/* mksol stuff only */

template < typename traits >
    uint set_f_coeffs(typename traits::scalar_t * f, uint maxdeg)
{
    uint deg = 0;
    // memset(f, 0, (maxdeg + 1) * sizeof(typename traits::scalar_t));
    traits::zero(f, maxdeg + 1);
    ifstream fx;
    must_open(fx, files::f % mine.scol);
    istream_iterator < mpz_class > inp(fx);
    for (deg = 0; !fx.eof(); deg++) {
	vector < mpz_class > line;
	BUG_ON(deg > maxdeg);
	back_insert_iterator < vector < mpz_class > >foo(line);
	/* read all n coefficients, only select the one we're
	 * really interested in.
	 */
	for (uint i = 0; i < globals::n; i++) {
	    *foo++ = *inp++;
	}
	// core_ops::assign<MODULUS_SIZE>(f[deg], line[globals::col]);
	traits::assign(f[deg], line, globals::col);
    }
    deg--;
    for (; deg >= 0 && traits::is_zero(f[deg]); deg--) {
	/* I don't really have time to check this -- shoudl be
	 * okay though. One thing is for sure, it's that we
	 * violate the fact that x^TB^iy = 0. However it's worth
	 * trying, still, to get a kernel vector.
	 */
	cout << "WARNING -- this needs to be double-checked\n";
	cout << fmt("skipping zero scalar at degree % : ") % deg;
	traits::print(cout, f[deg]);
	cout << endl;
    }


    return deg;
}


template < template < typename > class X, typename traits >
    void *thread_program(void *ptr)
{
    cout.flush();
    cerr.flush();

    cout << fmt("// thread % : setting up temporaries\n") % tseqid();
    cout << flush;
    X < traits > blah(ptr);
    barrier_wait(&globals::main_loop_barrier, NULL, NULL);
    cout << fmt("// thread % : go\n") % tseqid() << flush;
    blah.loop();
    blah.flush();
    barrier_wait(&globals::main_loop_barrier, NULL, NULL);

    return NULL;
}

template < typename traits > int program()
{
    using namespace globals;
    typedef typename traits::scalar_t scalar_t;

    int recover;

    nb_iter = nc + 2 * n;
    nb_iter = iceildiv(nb_iter, m) + iceildiv(nb_iter, n) + 10;
    cp_lag = max(nb_iter / CHECKPOINTS, 1u);
    cp_lag += cp_lag & 1;

    barrier_init(&main_loop_barrier, mine.nt, NULL);

    cout << fmt("// checkpointing is done every % iterations\n") % cp_lag;

    globals::vdata < traits > v0_and_f;
    std::vector < void *>targs(mine.nt, (void *) &v0_and_f);

    if (mine.task == "slave") {

	recover = recoverable_iteration(m, col, nbys, cp_lag);
	iter0 = recover * cp_lag;

	if (iter0 >= nb_iter) {
	    cout << "// everything already done for column " << col << "\n";
	    exit(0);
	}

	cout << fmt("// recovering up to checkpoint % == iteration %\n")
	    % recover % iter0;
	cout << fmt("// % total iterations needed, % to go\n")
	    % nb_iter % (nb_iter - iter0);
	cout << flush;

	recover_iteration(m, col, nbys, iter0);

	v0_and_f.v0 = new scalar_t[nc];
	recover_vector < traits > (nc, col, nbys, recover, v0_and_f.v0);

	start_threads(thread_program < slave_thread, traits >, targs);
	delete[]v0_and_f.v0;
    } else if (mine.task == "mksol") {

	/* See the discussion in master.cpp concerning the bound
	 * on the degree of f */
	uint bound_deg_f = iceildiv(m, n) +
	    iceildiv(m * iceildiv(nc, m), n) + 10;

	v0_and_f.f = new scalar_t[bound_deg_f + 1];
	v0_and_f.v0 = new scalar_t[nc];

	degf = set_f_coeffs < traits > (v0_and_f.f, bound_deg_f);

	cout << fmt("// F%[f0w2] has degree %\n") % mine.scol % degf;
	nb_iter = degf;

	/* FIXME ; the 0 here is intended to be changed */
	recover_vector < traits > (nc, col, nbys, 0, v0_and_f.v0);

	start_threads(thread_program < mksol_thread, traits >, targs);

	delete[]v0_and_f.v0;
	delete[]v0_and_f.f;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(false);
    cerr.tie(&cout);
    cout.rdbuf()->pubsetbuf(0, 0);
    cerr.rdbuf()->pubsetbuf(0, 0);

    common_arguments common;

    process_arguments(argc, argv, common, mine);

    string mstr;
    ifstream mtx;

    using namespace globals;

    must_open(mtx, files::matrix);
    get_matrix_header(mtx, nr, nc, mstr);
    BUG_ON_MSG(nr != nc, "Matrix is not square\n");
    /* It should not be too hard to accomodate zero-padded matrices, but:
     * - it causes problem with checks
     * - it's of little use anyway
     */

    globals::modulus = mpz_class(mstr);
    globals::modulus_u8 = globals::modulus.get_ui();
    globals::modulus_u16 = globals::modulus.get_ui();
    globals::modulus_u32 = globals::modulus.get_ui();
    globals::modulus_ulong = globals::modulus.get_ui();

    /*
       if (SIZ(globals::modulus.get_mpz_t()) != MODULUS_SIZE) {
       cerr << "### Not compiled with proper MODULUS_BITS\n"
       << "Using generic code instead\n";
       }
     */

    detect_mn(m, n);
    cout << fmt("// detected m = %\n") % m;
    cout << fmt("// detected n = %\n") % n;

    if (m == 0 || n == 0) {
	cerr << "First run: bw-prep, bw-balance, bw-secure\n";
	exit(1);
    }

    load_x_vectors(bw_x);
    cout << fmt("// loaded % x-vectors\n") % bw_x.size();
    ASSERT_ALWAYS((uint) m == bw_x.size());

    if (mine.e.first == -1 && mine.b.first == 0) {
	mine.e.first = m - 1;
    } else {
	/* Make sure we agree on the value of m */
	BUG_ON((int) m != (mine.e.first - mine.b.first + 1));
    }

    nbys = mine.e.second - mine.b.second + 1;
    col = mine.b.second;

    dotprod_diff_check.assign(nbys, mpz_class());

    /* Configure the thread data area, so that we can pre-fill it */
    configure_threads(mine.nt);

    cout << "// counting coefficients\n" << flush;
    mtxfile_pos.assign(mine.nt, 0);
    nb_coeffs.assign(mine.nt, 0);
    count_matrix_coeffs(mtx, nr, mtxfile_pos, nb_coeffs);

    for (int j = 0; j < mine.nt; j++) {
	uint i0 = (j * nr) / mine.nt;
	uint i1 = ((j + 1) * nr) / mine.nt;
	cout << fmt("// thread % : % rows % coeffs from pos %\n")
	    % j % (i1 - i0)
	    % nb_coeffs[j]
	    % mtxfile_pos[j];
    }

    thread_class_ptr.assign(mine.nt, (void *) NULL);
    mtx.close();

    cout.flush();
    cerr.flush();

    cout << "// modulus: " << globals::modulus << "\n";
    cout << "// modulus: " << MODLIMBS() << " words\n";
    cout << "// modulus: " << MODBITS() << " bits\n";
    cout << "// nbys: " << nbys << "\n";

    TRY_ALL_BINARY_CODES(program);
#if 0
    if (modbits < 8 && nbys == 8) {
	cout << "// Using SSE-2 code\n";
	return program < sse2_8words_traits > ();
    } else if (modbits < 20 && nbys == 2) {
	cout << "// Using SSE-2 code [ 2 words ]\n";
	return program < sse2_2words_traits > ();
    } else if (nbys == 1 && modbits <= (sizeof(unsigned long) * 8 / 2 - 4)) {
	cout << "// Using half-ulong code\n";
	return program < ulong_traits > ();
    } else if (nbys == 1 && modwords == MODULUS_SIZE) {
	cout << "// Using code for " << MODULUS_SIZE << " words\n";
	return program < typical_scalar_traits < MODULUS_SIZE > >();
    } else if (nbys == 1) {
	/* This amounts to at least a 2x slowdown */
	cout << "// Using generic code\n";
	cout << "// SLOW SLOW SLOW SLOW SLOW SLOW SLOW SLOW SLOW\n";
	return program < variable_scalar_traits > ();
    } else {
	cerr << "no available code\n";
	exit(1);
    }
#endif

    return 0;
}

/* vim: set sw=4 sta et: */
