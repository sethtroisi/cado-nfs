#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include "gmp-hacks.h"
#include "constants.hpp"

#include "expressions.hpp"
#include "manu.h"

/* technical aid */
#include <boost/utility/enable_if.hpp>

typedef unsigned int uint;

struct addition;
struct multiplication;
struct subtraction;
struct inversion;
struct conversion;

template<uint w>
class blob : public expressions::value<blob<w> > {
	typedef blob<w> self;

	friend struct multiplication;
	friend struct addition;
	friend struct subtraction;
protected:
	mp_limb_t x[w];
public:
	blob() {};
	blob(const blob<w>& a) { memcpy(x, a.x, sizeof(x)); }
	blob& operator=(mp_limb_t a) {
		x[0] = a;
		memset(x + 1, 0, (w-1) * sizeof(mp_limb_t));
		return *this;
	}
	bool operator==(const blob<w>& a) const {
		return memcmp(x, a.x, sizeof(x)) == 0;
	}
	bool operator==(mp_limb_t a) const {
		if (x[0] != a) return false;
		bool t = true;
		for(uint i = 1 ; i < w && (t = (x[i] == 0)) ; i++);
		return t;
	}
	inline bool operator!=(mp_limb_t a) const
	{ return !operator==(a); }
	inline bool operator!=(const blob<w>& a) const
	{ return !operator==(a); }
	operator mpz_class() const {
		mpz_class r;
		MPZ_SET_MPN(r.get_mpz_t(), x, w);
		return r;
	}
	template<typename T>
	inline blob(const expressions::expr<T>& e) {
		expressions::value<self>::operator=(e);
	}
	template<typename T>
	inline self& operator=(const expressions::expr<T>& e) {
		return expressions::value<self>::operator=(e);
	}

	TOPLEVEL_COMPOUND_BINARY_OP(operator+=, addition);
	TOPLEVEL_COMPOUND_BINARY_OP(operator*=, multiplication);
};

struct elt : public blob<MODULUS_SIZE>,
	       public expressions::value<elt>
{
	typedef elt self;
	typedef blob<MODULUS_SIZE> base0;

	friend struct inversion;
	friend struct conversion;
	friend struct mg_elt;
public:
	template<uint w, uint x> struct affine {
		typedef blob<w * MODULUS_SIZE + x> type;
	};
	typedef affine<2, 0>::type twice;
	static mpz_class modulus;
	/*
	template<typename T>
	inline self& operator=(const expressions::expr<T>& e) {
		return super::operator=(e);
	}
	*/
	elt& operator=(const elt& a) { base0::operator=(a); return *this; }
	elt& operator=(mp_limb_t a) { base0::operator=(a); return *this; }
	template<typename T>
	inline elt(const expressions::expr<T>& e) {
		expressions::value<self>::operator=(e);
	}
	template<typename T>
	inline self& operator=(const expressions::expr<T>& e) {
		return expressions::value<self>::operator=(e);
	}

	TOPLEVEL_COMPOUND_BINARY_OP(operator+=, addition);
	TOPLEVEL_COMPOUND_BINARY_OP(operator*=, multiplication);
};

struct mg_domain {
	mpz_class r2;	/* R^2 mod N */
	mpz_class mn;	/* -1/N mod R */
	mg_domain() {}
	mg_domain(const mpz_class& n) {
		int s = SIZ(n.get_mpz_t());
		ASSERT(s > 0);
		mpz_class r = mpz_class(1UL) << (mp_bits_per_limb * s);
		r2 = (r * r) % n;
		mpz_class foo;
		foo = -n;
		mpz_invert(mn.get_mpz_t(), foo.get_mpz_t(), r.get_mpz_t());
	}
	mg_domain(const mg_domain& d) : r2(d.r2), mn(d.mn) {}
};

struct mg_elt : public elt, public expressions::value<mg_elt> {
	typedef mg_elt self;
	static mg_domain mg;
	TOPLEVEL_COMPOUND_BINARY_OP(operator*=, multiplication);
};

/* Need to define:
 *
 * unary inversion K -> K
 * ``binary'' inversion K, D -> Kmgm, where D is an instance of the mgm
 * domain.
 *
 * conversion K, D -> Kmgm
 */

#define COMMUTE(r, t1, t2)						\
static inline void eval(r& r__, t1 const & t1__, t2 const & t2__) {	\
	eval(r__, t2__, t1__);						\
}

#define	TYPE_IF(c, t) typename boost::enable_if_c<(c), t>::type

/* {{{ multiplication ************************************************ */

namespace expressions {
template<uint v,uint w>	RESOLVE2(multiplication, blob<v>, blob<w>, blob<v+w>);
template<uint w>	RESOLVE2(multiplication, blob<w>, mp_limb_t, blob<w+1>);
template<> RESOLVE2(multiplication, elt, elt, elt::twice);
}

TOPLEVEL_BINARY_OP(operator*, multiplication)
TOPLEVEL_BINARY_OP_COMPLETE(operator*, multiplication, elt)
BINARY_OP_ALLOW(operator*, multiplication, mp_limb_t)
BINARY_OP_ALLOW(operator*, multiplication, int)
BINARY_OP_ALLOW(operator*, multiplication, bool)

struct multiplication {

template<uint w>
static void eval(mpz_class& r, const blob<w>& a, const blob<w>& b) {
	MPZ_GROW_ALLOC(r.get_mpz_t(), 2 * w);
	mpn_mul_n(PTR(r.get_mpz_t()), a.x, b.x, w);
	SIZ(r.get_mpz_t()) = 2*w;
	MPZ_NORMALIZE(r.get_mpz_t());
}

template<uint w>
static void eval(blob<2*w>& r, const blob<w>& a, const blob<w>& b) {
	mpn_mul_n(r.x, a.x, b.x, w);
}

template<uint w>
static void eval(blob<2*w>& r, const elt& a, const elt& b) {
	mpn_mul_n(r.x, a.x, b.x, w);
}

template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const mp_limb_t& b) {
	r.x[w] = mpn_mul_1(r.x, a.x, b, w);
}
template<uint w> COMMUTE(blob<w+1>, mp_limb_t, blob<w>)
};

/* }}} */
/* {{{ addition ************************************************ */

TOPLEVEL_BINARY_OP(operator+, addition)
BINARY_OP_ALLOW(operator+, addition, mp_limb_t)
BINARY_OP_ALLOW(operator+, addition, int)
BINARY_OP_ALLOW(operator+, addition, bool)

struct addition {

template<uint w>
static void eval(mpz_class& r, const blob<w>& a, const blob<w>& b) {
	MPZ_GROW_ALLOC(r.get_mpz_t(), w + 1);
	PTR(r.get_mpz_t())[w] = mpn_add_n(PTR(r.get_mpz_t()), a.x, b.x, w);
	SIZ(r.get_mpz_t()) = w + 1;
	MPZ_NORMALIZE(r.get_mpz_t());
}

template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const blob<w>& b) {
	r.x[w] = mpn_add_n(r.x, a.x, b.x, w);
}

template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const mp_limb_t& b) {
	 r.x[w] = mpn_add_1(r.x, a.x, w, b);
}
template<uint w> COMMUTE(blob<w+1>, mp_limb_t, blob<w>)

/* NOTE NOTE NOTE !
 * this is a wanted optimization. The program must take care of not
 * triggering so many calls of this function that an overflow occurs */
template<uint w>
static void eval(blob<w+1>& r, const blob<w+1>& a, const blob<w>& b) {
	r.x[w] = a.x[w] + mpn_add_n(r.x, a.x, b.x, w);
	/* we're dead if the latter overflows */
}
template<uint w> COMMUTE(blob<w+1>, blob<w>, blob<w+1>)

template<uint v, uint w>
static TYPE_IF(v < w, void)
eval(blob<w+1>& r, const blob<v>& a, const blob<w>& b) {
	mp_limb_t c;
	c = mpn_add_n(r.x, a.x, b.x, v);
	c = mpn_add_1(r.x + v, b.x + v, w-v, c);
	r.x[w] = c;
}

template<uint v, uint w>
static TYPE_IF(w < v, void)
eval(blob<v+1>& r, const blob<v>& a, const blob<w>& b) {
	mp_limb_t c;
	c = mpn_add_n(r.x, a.x, b.x, w);
	c = mpn_add_1(r.x + w, a.x + w, v-w, c);
	r.x[v] = c;
}

};

namespace expressions {
template<uint v, uint w> RESOLVE2(addition, blob<v>, blob<w>, blob<1+((v<w)?w:v)>);
}

/* }}} */
/* {{{ subtraction ************************************************ */

TOPLEVEL_BINARY_OP(operator-, subtraction)
BINARY_OP_ALLOW_2ND(operator-, subtraction, mp_limb_t)
BINARY_OP_ALLOW_2ND(operator-, subtraction, int)
BINARY_OP_ALLOW_2ND(operator-, subtraction, bool)

struct subtraction {

template<uint w>
static void eval(mpz_class& r, const blob<w>& a, const blob<w>& b) {
	MPZ_GROW_ALLOC(r.get_mpz_t(), w + 1);
	PTR(r.get_mpz_t())[w] = mpn_sub_n(PTR(r.get_mpz_t()), a.x, b.x, w);
	SIZ(r.get_mpz_t()) = w + 1;
	MPZ_NORMALIZE(r.get_mpz_t());
}

template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const blob<w>& b) {
	r.x[w] = - mpn_sub_n(r.x, a.x, b.x, w);
}

template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const mp_limb_t& b) {
	 r.x[w] = - mpn_sub_1(r.x, a.x, w, b);
}
/* commutation not provided ; change calling code instead */

/* NOTE NOTE NOTE !
 * this is a wanted optimization. The program must take care of not
 * triggering so many calls of this function that an underflow occurs */
template<uint w>
static void eval(blob<w+1>& r, const blob<w+1>& a, const blob<w>& b) {
	r.x[w] = a.x[w] - mpn_sub_n(r.x, a.x, b.x, w);
	/* we're dead if the latter underflows */
}
template<uint w>
static void eval(blob<w+1>& r, const blob<w>& a, const blob<w+1>& b) {
	r.x[w] = - mpn_sub_n(r.x, a.x, b.x, w) - b.x[w];
}

/* nothing provided when v < w ; change calling code */
template<uint v, uint w>
static TYPE_IF(w < v, void)
eval(blob<v+1>& r, const blob<v>& a, const blob<w>& b) {
	mp_limb_t c;
	c = mpn_sub_n(r.x, a.x, b.x, w);
	c = mpn_sub_1(r.x + w, a.x + w, v-w, c);
	r.x[v] = -c;
}

};

namespace expressions {
template<uint v, uint w> RESOLVE2(subtraction, blob<v>, blob<w>, blob<1+((v<w)?w:v)>);
}

/* }}} */
/* {{{ inversion ************************************************ */

TOPLEVEL_UNARY_OP(inv, inversion);

struct inversion {


static void eval(elt& r, const elt& a) {
	int rc;
	mpz_class zr;
	mpz_class za = a;
	rc = mpz_invert(zr.get_mpz_t(), za.get_mpz_t(), elt::modulus.get_mpz_t());
	BUG_ON(rc == 0);
}

/*
static void eval(mg_elt& r, const mg_elt& a) {
	eval((elt&)r, (const elt&)a);
}
*/

};

namespace expressions {
template<> RESOLVE1(inversion, elt, elt);
template<> RESOLVE1(inversion, mg_elt, mg_elt);
}

/* }}} */
/* {{{ conversion ************************************************ */

TOPLEVEL_BINARY_OP(convert, conversion);

struct conversion {

};

/* }}} */

#include <ostream>

namespace std {
	template<uint w>
	ostream& operator<<(ostream& o, const blob<w> & x) {
		mpz_class z = x;
		return o << z;
	}
}
