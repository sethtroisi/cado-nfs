#ifndef EXPRESSIONS_HPP_
#define EXPRESSIONS_HPP_

/* {{{ Reference for the expression template header file */

/*
 * This header file provides a framework for allowing simple, yet
 * efficient arithmetic for user-defined types, using operator
 * overloading.
 *
 * The intent is to have the full power of an instruction like:
 *
 * 	my_mul(dst, src1, src2)
 *
 * combined with the ease of use of:
 *
 * 	dst = src1 * src2
 *
 * No copy of src1 or src2 is produced inbetween. When the code is
 * compiled with optimization turned on, the two syntaxes above are
 * expected to yield identical results.
 *
 * Behind the scenes, this is achieved with templates and lazy
 * evaluation. The right-hand side becomes a complete expression tree,
 * which is reduced only once the left-hand side is available. For
 * complex expression trees, the number of temporaries created is exactly
 * the number of inner nodes of your trees.  (You don't have to know much
 * about what happens behind the scenes, but it certainly doesn't harm
 * if you do).
 *
 * In order to allow such operations with your types, you have to:
 *
 * 1) create your types. Say you're building type blah. The type must
 *    meet the following requirements:
 *    -- blah must inherit from expressions::value<blah>. Yup, blah
 *    itself is used as a template parameter.  That's allowed and
 *    extremely useful. Note that expressions::value<blah> has size zero,
 *    so you're not paying anything extra.
 *
 *    -- in the definition of class blah, the following common code (at
 *    some point this will be a macro).
 *
 *       - typedef blah as self.
 *       - include these templates:
 *
 *    template<typename T>
 *    inline elt(const expressions::expr<T>& e) {
 *    	expressions::value<self>::operator=(e);
 *    }
 *    template<typename T>
 *    inline self& operator=(const expressions::expr<T>& e) {
 *    	return expressions::value<self>::operator=(e);
 *    }
 *
 *    [TODO : review this ; the needed part should be trimmed down. It
 *    seems clumsy to require a ctor, and could actually trigger bugs]
 *
 * 2) create function objects for the binary operations. For example, create
 *    ``struct addition'', ``struct multiplication'', and so on.  each of
 *    these function objects should contain several ::eval() functions,
 *    taking (dst, src1, src2) as reference arguments, the last two being
 *    const. The definition and declaration may be separated if you so
 *    wish.
 *
 *    For instance, we may have:
 *
 *    struct multiplication {
 *    template<uint w> static void eval(elt<2*w>&,
 *   				const elt<w>&, const elt<w>&);
 *    };
 *
 *    Your ::eval() functions are the only way of specifying which kinds
 *    of operations are allowed.
 *
 *    Additionally, specializations of the resolve_binary<T, T1, T2> or
 *    resolve_unary<T, T1> must be provided, with T==the function object.
 *    This indicates the preferred type for temporaries of the type T(T1,
 *    T2) (in the ::type member). For instance:
 *
 *    template<uint v, uint w> struct resolve<multiplication, elt<v>, elt<w> > {
 *    		typedef elt<v+w> type;
 *    };
 *
 * 3) add macros for binding your function objects to the corresponding
 *    operators, e.g:
 *
 *      TOPLEVEL_BINARY_OP(operator*, multiplication)
 *      TOPLEVEL_UNARY_OP(inv, inversion)
 *
 *    This concerns only the types you have defined (inheriting from
 *    value<blah>). If you want to allow mixed-type operations, do
 *
 *      BINARY_OP_ALLOW(operator*, multiplication, unsigned long)
 *      UNARY_OP_ALLOW(inv, inversion, unsigned long)
 *
 *    Of course your ``struct multiplication'' object must include the
 *    proper ::eval() prototypes.
 *
 *    there's also BINARY_OP_ALLOW_1ST and BINARY_OP_ALLOW_2ND if you
 *    want to allow the foreign type only as 1st or second argument.
 *
 * 4) (optional) if you want to have +=, *= operators, use within the
 *    definition of class blah the following macros:
 *    class blah {
 *	[...]
 *      TOPLEVEL_COMPOUND_BINARY_OP(operator+=, addition);
 *      TOPLEVEL_COMPOUND_BINARY_OP(operator*=, multiplication);
 *      TOPLEVEL_COMPOUND_UNARY_OP(inverse, multiplication);
 *    };
 *
 *
 * Further notes.
 *
 * there's no reason to limit the overloads to the language operators.
 * Any function name will work.
 *
 * Operators of arbitrary arity are not supported, but it's not a
 * terrible change if needed.
 *
 * If you want things such as x = pow(a, k) % m which really change the
 * evaluation process in depth (by doing a modular multiplication at each
 * step), then you'll have to mess with this header file. It's possible,
 * but not trivial. you'd have to shortcut value::operator=, either via
 * explicit specialization or other kinds of horrible tricks. The problem
 * really is that you're not sure to want that.
 *
 * 
 * The downsides of all this are:
 * - in any case the compilation time for this stuff is higher than for
 *   an approach not having this zero-copy feature
 * - it may happen that your compiler fails to inline out the function
 *   calls which should be omitted, or leave dummy space on the stack
 *   which will never get touched, just because it presumes that some
 *   objects might live in there...
 *
 * }}} */

namespace expressions {

/* CRTP : provide static inheritance. We use this as a common tag */
template<typename Derived> struct expr {
	Derived& top() { return (Derived&)*this; }
	const Derived& top() const { return (const Derived&)*this; }
};

/* The wrapping structure is here only to allow mixing expressions with
 * foreign types. In order to avoid exposing too much of the gory
 * details, we unwrap the type when needed
 */
template<typename T> struct wrap;
template<typename T> struct unwrap {
	typedef T type;
	inline T& operator()(T& x) const { return x; }
	inline const T& operator()(const T& x) const { return x; }
};
template<typename T> struct unwrap<wrap<T> > {
	typedef T type;
	inline T& operator()(wrap<T>& x) const { return x.x_; }
	inline const T& operator()(const wrap<T>& x) const { return x.x_; }
};

/* Types deriving from expr<> represent expression subtrees. They
 * have the capability of flattening the corresponding tree, in order to
 * obtain the top-level evaluated value, in a type which is chosen in
 * accordance with the resolve<> policies, specified as member templates
 * of the operation class.  Access to this functionality is via the
 * flat() member. Note that on immediate values, this is trivial.  Of
 * course all this is useful for creating temporaries.
 */

template<typename, typename, typename> struct resolve_binary;

template<typename Operation, typename T1, typename T2>
struct binary_op : public expr<binary_op<Operation, T1, T2> > {
	const T1& op1;
	const T2& op2;
	binary_op(const T1& op1, const T2& op2) : op1(op1), op2(op2) {}

	typedef typename resolve_binary<Operation,
			typename unwrap<typename T1::otype>::type,
			typename unwrap<typename T2::otype>::type
		>::type otype;
	/* This function is called only when a temporary is needed */
	otype flat() const {
		otype v = *this;
		return v;
	}
};
#define RESOLVE2(t, t1, t2, r) struct resolve_binary<t, t1, t2 > { typedef r type; }

template<typename, typename> struct resolve_unary;
template<typename Operation, typename T1>
struct unary_op : public expr<unary_op<Operation, T1> > {
	const T1& op1;
	unary_op(const T1& op1) : op1(op1) {}

	typedef typename resolve_unary<Operation,
			typename unwrap<typename T1::otype>::type
		>::type otype;
	/* This function is called only when a temporary is needed */
	otype flat() const {
		otype v = *this;
		return v;
	}
			
};
#define RESOLVE1(t, t1, r) struct resolve_unary<t, t1 > { typedef r type; }

/* Note: values appear only as base tags of derived complete types */
/* They also denote the only terminal type which is allowed */
template<typename Derived>
struct value : public expr<value<Derived> >{
	typedef Derived otype;

	Derived& flat() { return (Derived&) *this; }
	const Derived& flat() const { return (const Derived&) *this; }

	template<typename T>
	inline Derived& operator=(const expr<T>& a) {
		return operator=((const T&) a);
	}

	/* FIXME ; this is not documented so far, and in fact it's barely
	 * useful, as we're better off specifying type conversions
	 * explicitly in the derived type
	 */
	template<typename T>
	inline Derived& operator=(const value<T>& a) {
		convert((Derived&)*this, a.flat());
		return flat();
	}
	
	/* Note : these two functions are responsible of our most valued
	 * benefit. The ``silly way'' would consist in tagging the
	 * function above with anything matching the expr<T> pattern.
	 * This would do with lazy evaluation the same thing that is done
	 * with the much simpler scheme of immediate evaluation.  By
	 * separating value<T> from other types, we get the full benefit
	 * of lazy evaluation: the possibility to avoid creation of
	 * temporaries at top level.  */
	template<typename Operation, typename T1>
	inline Derived& operator=(const unary_op<Operation, value<T1> >& a) {
		Operation::eval((Derived&)*this,
				unwrap<typename T1::otype>()(a.op1.flat()));
		return flat();
	}
	template<typename Operation, typename T1, typename T2>
	inline Derived& operator=(const binary_op<Operation, T1, T2>& a) {
		Operation::eval((Derived&)*this,
				unwrap<typename T1::otype>()(a.op1.flat()),
				unwrap<typename T2::otype>()(a.op2.flat()));
		return flat();
	}
};

/* ------------------------------------------------------------------------ */
/* These helper macros are the recommended way of declaring the operators
 * allowed from our framework.
 *
 * See later for an example use ; wrap<T> is defined just below.
 */

/* The _COMPLETE versions are responsible of handing over the control to
 * the correct templates. In some cases, base types are not considered
 * for overload resolution, resulting in missing matches.
 */

#define TOPLEVEL_BINARY_OP(fun_, op)					\
template<typename T, typename U> inline					\
expressions::binary_op<op, T, U> fun_(					\
		const expressions::expr<T>& a,				\
		const expressions::expr<U>& b)				\
{									\
	return expressions::binary_op<op, T, U>(a.top(), b.top());	\
}

#define CONST_RECAST(v, t)						\
	((const expressions::expr<expressions::value< t > >&) (v))
#define RECAST(v, t)							\
	((expressions::expr<expressions::value< t > >&) (v))

#define TOPLEVEL_BINARY_OP_COMPLETE(fun_, op, t)			\
template<typename T> inline						\
expressions::binary_op<op, T, expressions::value<t > > fun_(		\
		const expressions::expr<T>& a,				\
		const t& b)						\
{									\
	return fun_(a, CONST_RECAST(b, t));				\
}									\
template<typename T> inline						\
expressions::binary_op<op, expressions::value<t >, T> fun_(		\
		const t& a,						\
		const expressions::expr<T>& b)				\
{									\
	return fun_(CONST_RECAST(a, t), b);				\
}									\
inline									\
expressions::binary_op<op, expressions::value<t >,			\
			expressions::value<t > > fun_(			\
		const t& a,						\
		const t& b)						\
{									\
	return fun_(CONST_RECAST(a, t), CONST_RECAST(b, t));		\
}

#define TOPLEVEL_UNARY_OP(fun_, op)					\
template<typename T> inline						\
expressions::unary_op<op, T> fun_(const expressions::expr<T>& a)	\
{									\
	return expressions::unary_op<op, T>(a.top());			\
}

#define TOPLEVEL_UNARY_OP_COMPLETE(fun_, op, t)				\
inline						\
expressions::unary_op<op, expressions::value<t > > fun_(const t& a)	\
{									\
	return fun_(CONST_RECAST(a));					\
}

#define BINARY_OP_ALLOW_1ST(fun_, op, type)				\
template<typename T> inline						\
expressions::binary_op<op, expressions::wrap<type>, T> fun_(		\
		const type& a, const expressions::expr<T>& b)		\
{									\
	return expressions::binary_op<op, expressions::wrap<type>, T>(	\
			expressions::wrap<type>(a), b.top());		\
}

#define BINARY_OP_ALLOW_2ND(fun_, op, type)				\
template<typename T> inline						\
expressions::binary_op<op, T, expressions::wrap<type> > fun_(		\
		const expressions::expr<T>& a,				\
		const type& b)						\
{									\
	return expressions::binary_op<op, T, expressions::wrap<type> >(	\
			a.top(), expressions::wrap<type>(b));		\
}

#define BINARY_OP_ALLOW(fun_, op, type)					\
	BINARY_OP_ALLOW_1ST(fun_, op, type)				\
	BINARY_OP_ALLOW_2ND(fun_, op, type)

#define UNARY_OP_ALLOW(fun_, op, type)					\
inline									\
expressions::unary_op<op, expressions::wrap<type> > fun_(const type& b)	\
{									\
	return expressions::unary_op<op, expressions::wrap<type> >(	\
			expressions::wrap<type>(b));			\
}

#define TOPLEVEL_COMPOUND_BINARY_OP(fun_, op)				\
template<typename T> inline						\
self& fun_(const expressions::expr<T>& a)				\
{									\
	typedef expressions::value<self> super;				\
	return operator=(expressions::binary_op<op, super, T>(		\
				super::top(), a.top()));		\
}
#define TOPLEVEL_COMPOUND_UNARY_OP(fun_, op)				\
self& fun_()								\
{									\
	typedef expressions::value<self> super;				\
	return operator=(expressions::unary_op<op, super>(		\
				super::top()));				\
}


/* ------------------------------------------------------------------------ */
/* In order to plug within our framework, the expression types must be
 * declared as inheriting from value<self>, where self is really the type
 * itself (this is the CRTP technique). It is possible, and supported, to
 * use foreign types as well with this technique. For this purpose, we
 * define the following simple template.  Note that this does not solve
 * the problem of assigning to foreign types.
 */
/* interface foreign types with our framework */
template<typename T>
struct wrap : public value< wrap<T> > {
	T x_;
	wrap(const T& x) : x_(x) {}
};

} // namespace expressions

#endif	/* EXPRESSIONS_HPP_ */
