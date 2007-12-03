#ifndef SCALAR_HPP_
#define SCALAR_HPP_

#error "this stuff is useless crap"

#include <gmp.h>
#include <gmpxx.h>
#include "gmp-hacks.h"
#include "manu.h"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/utility/enable_if.hpp>
// #include <boost/mpl/math/is_even.hpp>
// #include <boost/mpl/int.hpp>
// #include <boost/mpl/if.hpp>

template<mp_size_t, typename> class scalar;
template<mp_size_t> class vect;
template<mp_size_t> class scalar_pool;
template<mp_size_t> class scalar_simple;
template<mp_size_t> class vect_simple;

// scalar and vect are true pointers ; they do not hold the data, and
// have to be created from a scalar pool.
// However, you manipulate them as you would manipulate actual long
// integers : In particular, assignment copies the pointed-to data, while
// the copy ctor binds a new reference to the same data.
//
// the _simple variants come with their own storage.

// all this is still quite unpolished. In particular, const-wise, it's a
// total mess.

template<mp_size_t width, typename T = mp_limb_t>
class scalar {
	friend class scalar_pool<width>;
	friend class scalar_simple<width>;
	friend class vect<width>;
	struct enabler {};
protected:
	T * data;
	void init(const scalar_pool<width>& pool, int pos) {
		data = pool.data + (pos * width);
	}
	scalar(const scalar_pool<width>& pool, int pos) { init(pool, pos); }
public:
	/* always use the copy ctor from something reasonable */
	scalar() { data = NULL; }
	inline scalar(const scalar<width>& x) : data(x.data) {}
	~scalar() { data = NULL; }
	operator T *() { return data; }
	operator T const *() const { return data; }
	scalar& operator=(mp_limb_t x) {
		memset(data, 0, width * sizeof(mp_limb_t));
		data[0] = x;
		return *this;
	}
	/* There is willingly NO operator= defined */
	scalar& set(const scalar<width>& x) {
		memcpy(data, (const T *) x, width * sizeof(mp_limb_t));
		return *this;
	}
	/*
#define CONDITION(w)						\
                , typename boost::enable_if_c<                  \
			w & 1,					\
                        enabler                                 \
                >::type = enabler()
	scalar& set(const scalar<(width-4)/2>& x CONDITION(width)) {
		mp_size_t xw = (width-4)/2;
		memcpy(data, (const T *) x, xw * sizeof(mp_limb_t));
		memset(data + xw, 0, (width - xw) * sizeof(mp_limb_t));
		return *this;
	}
	*/
	operator mpz_class() const {
		mpz_class r;
		MPZ_SET_MPN(r.get_mpz_t(), data, width);
		return r;
	}
	scalar& operator=(const mpz_class& x) {
		BUG_ON(SIZ(x.get_mpz_t()) < 0 || SIZ(x.get_mpz_t()) >= width);
		MPN_SET_MPZ(data, width, x.get_mpz_t());
		return *this;
	}
	/* Move the pointer ; no check is done */
	scalar& step(int k = 1) {
		data += k * width;
		return *this;
	}
};

template<mp_size_t width>
class vect : public scalar<width> {
public:
	unsigned int size;
protected:
	void init(const scalar_pool<width>& pool, int pos, unsigned int sz)
	{
		scalar<width>::init(pool, pos);
		size = sz;
	}
	vect(const scalar_pool<width>& pool, int pos, unsigned int sz)
		: scalar<width>(pool, pos), size(sz) {}
private:
	/* iterator boilerplate */
	template<class X> class iter : public boost::iterator_facade<
				       iter<X>,
				       scalar<width, X>,
				       boost::random_access_traversal_tag,
				       scalar<width, X>
				       >
	{
#define COND(X, Y)                                              \
                , typename boost::enable_if<                    \
                        boost::is_convertible<Y *, X *>,        \
                        enabler                                 \
                >::type = enabler()
		public:
		typedef long Distance;
		private:
		scalar<width, X> sc;
		struct enabler {};
		typedef typename boost::is_const<X> isc;
		friend class boost::iterator_core_access;
		public:
		iter() {}
		template<class Y>
		iter(scalar<width, Y> const & a COND(X,Y)) : sc(a) {}
		template<class Y>
		iter(iter<Y> const & a COND(X,Y)) : sc(a.sc) {}
		void increment() { sc.step(); }
		void advance(int n) { sc.step(n); }
                template<class Y>
                int distance_to(iter<Y> const & a COND(X,Y)) const
                {
			return ((Y*) a.sc - (X*)sc) / width;
		}
		template<class Y>
		int equal(iter<Y> const & a COND(X,Y)) const
		{ BUG(); return ((Y*)a.sc - (X*)a.sc); }
		scalar<width, X> dereference() const {
			return sc;
		}
	};
public:
	vect() {}
	inline vect(const vect<width>&) {}
	~vect() {}
	using scalar<width>::operator mp_limb_t *;
	scalar<width> operator[](int pos) {
		scalar<width> res;
		res.data = scalar<width>::data + width * pos;
		return res;
	}
	scalar<width, const mp_limb_t> operator[](int pos) const {
		scalar<width> res;
		res.data = scalar<width>::data + width * pos;
		return res;
	}
	/*********/
	typedef iter<mp_limb_t> iterator;
	typedef iter<const mp_limb_t> const_iterator;
	iterator begin() { return (*this)[0]; }
	const_iterator begin() const { return (*this)[0]; }
	iterator end() { return (*this)[size]; }
	const_iterator end() const { return (*this)[size]; }
	/*********/
	vect<width> set_zero(unsigned int i0, unsigned int i1) {
		memset(scalar<width>::data + i0 * width, 0, (i1 - i0) * width);
		return *this;
	}
	inline vect<width> set_zero() { return set_zero(0, size); }
};

template<mp_size_t width>
class scalar_pool {
	friend class scalar<width>;
	friend class vect_simple<width>;
	mp_limb_t * data;
	void init(mp_size_t n) {
		if (data != NULL) delete[] data;
		data = new mp_limb_t[n * width];
	}
	public:
	scalar_pool() { data = NULL; }
	scalar_pool(mp_size_t n) { data = new mp_limb_t[n * width]; }
	~scalar_pool() { delete[] data; }
	/* We do NOT want a copy ctor */
	scalar<width> operator[](int pos) {
		return scalar<width>(*this, pos);
	}
};

template<mp_size_t width>
class scalar_simple : public scalar<width> {
	scalar_pool<width> buf;
public:
	scalar_simple() : buf(1) { ((scalar<width>&)*this).data = buf[0].data; }
	using scalar<width>::operator=;
	using scalar<width>::operator mpz_class;
	using scalar<width>::operator mp_limb_t *;
};

template<mp_size_t width>
class vect_simple : public vect<width> {
	scalar_pool<width> buf;
public:
	vect_simple() {}
	inline vect_simple(const vect_simple<width>&) {}
	~vect_simple() {}
	void init(unsigned int sz) {
		buf.init(sz);
		vect<width>::init(buf, 0, sz);
	}
	vect_simple(unsigned int sz) : buf(sz) {
		vect<width>::init(buf, 0, sz);
	}
};

namespace std {
	template<mp_size_t w>
	inline std::ostream& operator<<(std::ostream& o, const scalar<w>& x)
	{
		return o << (mpz_class) x;
	}
	template<mp_size_t w>
	inline std::istream& operator>>(std::istream& i, scalar<w>& x)
	{
		mpz_class r;
		i >> r;
		x = r;
		return i;
	}
}


#endif	/* SCALAR_HPP_ */
