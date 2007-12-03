#ifndef MATRICES_HPP_
#define MATRICES_HPP_

#include <cstdlib>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <algorithm>
#include <iterator>

/* This type hierarchy mimicks boost::multi_array, but better suits my
 * purposes, and weighs much less */

/* {{{ Quick type list templates */
#if 0
#include <boost/compressed_pair.hpp>
template<typename Car, typename Cdr> struct cons : public boost::compressed_pair<Car, Cdr> {
	typedef boost::compressed_pair<Car, Cdr> super;
	typedef Car head;
	typedef Cdr tail;
	unsigned int size() const { return super::first().size(); }
	unsigned int total_size() const {
		return size() * super::second().total_size();
	}
};
template<typename Car> struct cons<Car, void> : public Car {
	typedef Car head;
	typedef void tail;
	unsigned int size() const { return Car::size(); }
	unsigned int total_size() const { return Car::size(); }
};
#else
/* This is quite an ugly hack. The idea is to disambiguate the bases
 * using different instantiations. */
template<typename T, int WrapLevel = 0> struct wrap : public T {
	enum { wrap_level__ = WrapLevel };
};

template<typename Car, typename Cdr> struct cons :
        public wrap<Car, Cdr::wrap_level__>, public Cdr
{
	enum { wrap_level__ = Cdr::wrap_level__ + 1 };
        typedef Car head;
        typedef Cdr tail;
	typedef cons<Car, Cdr> self;
	// using wrap<Car, wrap_level__ - 1>::size;
	inline unsigned int size() const {
		return wrap<Car, wrap_level__ - 1>::size();
	}
        unsigned int total_size() const {
		return this->size() * Cdr::total_size();
        }
};
template<typename Car> struct cons<Car, void> : public wrap<Car, 0> {
	enum { wrap_level__ = 1 };
        typedef Car head;
        typedef void tail;
	typedef cons<Car, void> self;
	using wrap<Car, wrap_level__ - 1>::size;
        unsigned int total_size() const { return wrap<Car, 0>::size(); }
};
#endif
/* }}} */

/* {{{ convenience enable_if macros */
/* We use private enabler types instead of void*, in order to avoid
 * errors */
#define IF_CONV(X, Y)						\
                , typename boost::enable_if<			\
                        boost::is_convertible<X, Y>,		\
                        enabler					\
                >::type = enabler()

#define	IF_CONV_RET(X, Y, T)					\
	typename boost::enable_if<boost::is_convertible<X, Y>, T>::type
/* }}} */

/* forward-declare our templates */
template<typename T, typename Dims> class row;
template<typename T, typename Dims> class iter;

/* {{{ iterator type ; it is simpler to keep it apart from the row type */
template<typename T, typename Dims> class iter :
	private row<T, Dims>,
	public boost::iterator_facade<
		iter<T, Dims>,	// self
		row<T, Dims>,	// value type
		boost::random_access_traversal_tag,
		row<T, Dims> >	// reference type
{
	struct enabler {};
	template<typename, typename> friend class row;
	template<typename, typename> friend class iter;
	friend class boost::iterator_core_access;
public:
	iter() {}
	typedef row<T, Dims> super;
	typedef iter<T, Dims> self;
	typedef self type;
	typedef super value_type;
	typedef super reference;
private:
	static T *& base_ptr(type& x) { return x.super::base_ptr(); }
	static const T * base_ptr(const type& x) { return x.super::base_ptr(); }
public:
	template<typename U>
	iter(iter<U, Dims> const& a IF_CONV(U*,T*)) : super(a) {}
	
	template<typename U>
	IF_CONV_RET(U, T, self&) operator=(const iter<U, Dims>& a) {
		super::reseat((typename iter<U, Dims>::super) a);
		return *this;
	}
	// enforce the template
	self& operator=(const self& a) { return operator=<T>(a); }
private:
	reference dereference() const { return (super)*this; }
	void advance(int n) { super::it() += n * super::size(); }
	inline void increment() { advance(1); }
	inline void decrement() { advance(-1); }

	template<typename U>
	IF_CONV_RET(U, T, bool) equal(iter<U, Dims> const& a) const {
		return super::it() == a.it();
	}

	template<typename U>
	IF_CONV_RET(U, T, ptrdiff_t) distance_to(iter<U, Dims> const& a) const {
		return (a.it() - super::it()) / (ptrdiff_t) super::size();
	}
};

/* final specialization */
template<typename T> class iter<T, void>
{
	template<typename, typename> friend class row;
	template<typename, typename> friend class iter;
public:
	typedef T * type;
private:
	static T *& base_ptr(type& x) { return x; }
	static const T * base_ptr(const type& x) { return x; }
};
/* }}} */

template<typename T> struct row_ops {};

/* {{{ row type ; does _not_ own the data. It may be regarded as (sort
 * of) a container.
 */
template<typename T, typename Dims>
class row : public Dims, public row_ops<row<T, Dims> >
{
	struct enabler {};
	template<typename, typename> friend class row;
	template<typename, typename> friend class iter;
	typedef iter<T, typename Dims::tail> iter_pre;
	typedef iter<const T, typename Dims::tail> const_iter_pre;
public:
	typedef row<T, Dims> self;
	friend class row_ops<self>;
	typedef typename iter_pre::type
				iterator;
	typedef typename const_iter_pre::type
				const_iterator;
	typedef typename std::iterator_traits<iterator>::value_type
				value_type;
	typedef typename std::iterator_traits<iterator>::reference
				reference;
	typedef typename std::iterator_traits<const_iterator>::reference
				const_reference;
	iterator ptr_;
	inline iterator& it() { return ptr_; }
	inline const iterator& it() const { return ptr_; }
private:
	row(const iterator& p) : ptr_(p) {}
	template<typename U>
	IF_CONV_RET(U*,T*,self&) reseat(row<U, Dims> const& a) {
		it() = a.it();
		return *this;
	}
protected:
	T*& base_ptr() { return iter_pre::base_ptr(it()); }
	const T* base_ptr() const { return iter_pre::base_ptr(it()); }
public:
	row() : ptr_() {}
	template<typename U>
	row(row<U, Dims> const& a IF_CONV(U, T)) :
		Dims((const Dims&) a),
		/* row_ops<self>(a), */
		ptr_(a.it()) {}
	row(const self& a) :
		Dims((const Dims&) a),
		/* row_ops<self>(a), */
		ptr_(a.it()) {}
	using Dims::total_size;

	/* it is not uncommon for the row ops to include an overload of
	 * operator= ; we conveniently avoid disabling this overload by
	 * bringing this definition to the surface of the class */
	using row_ops<self>::operator=;
	template<typename U>
	IF_CONV_RET(U, T, self&) operator=(const row<U, Dims>& a) {
		std::copy(a.begin(), a.end(), begin());
		return *this;
	}
	// enforce the template
	inline self& operator=(const self& a) { return operator=<T>(a); }

	/* random access container facade */
	using Dims::size;
#define	XS(i)	ASSERT((i) < this->size() && (i) >= 0)
	reference operator[](unsigned int i) { XS(i); return *(it() + i); }
	const_reference operator[](unsigned int i) const { XS(i); return *(it() + i); }
	iterator begin() { return it(); }
	const_iterator begin() const { return it(); }
	iterator end() { return it() + size(); }
	const_iterator end() const { return it() + size(); }
	reference front() { return (*this)[0]; }
	const_reference front() const { return (*this)[0]; }
	reference back() { return (*this)[size()-1]; }
	const_reference back() const { return (*this)[size()-1]; }
};
/* }}} */

/* {{{ holder : This is the derivation of the row type that owns the data
 * */
// #include <iostream>
template<typename T> class holder;
template<typename T, typename Dims>
class holder<row<T, Dims> > : public row<T, Dims>
{
	struct enabler;
	public:
	typedef row<T, Dims> super;
	typedef holder<super> self;
	holder(const Dims& d = Dims()) {
		// std::cout << "Allocating " << super::total_size() << " elements" << std::endl;
		(Dims&) *this = d;
		super::base_ptr() = new T[super::total_size()];
	}
	template<typename U>
	IF_CONV_RET(U, T, self&) operator=(const holder<row<U, Dims> >& a) {
		typedef typename holder<row<U, Dims> >::super a_super;
		ASSERT_ALWAYS((const Dims&) *this == (const Dims&) a);
		((super&)*this) = (const a_super&)a;
		return *this;
	}
	template<typename U>
	holder(holder<row<U, Dims> > const& a IF_CONV(U,T)) {
		// std::cout << "Allocating " << super::total_size() << " elements" << std::endl;
		(Dims&) *this = (const Dims&) a;
		super::base_ptr() = new T[super::total_size()];
		operator=(a);
	}
	inline self& operator=(const self& a) { return operator=<T>(a); }
	/* it is not uncommon for the row ops to include an overload of
	 * operator= ; we conveniently avoid disabling this overload by
	 * bringing this definition to the surface of the class */
	using row_ops<super>::operator=;
	~holder() { delete[] super::base_ptr(); }
};
/* }}} */

#endif	/* MATRICES_HPP_ */
