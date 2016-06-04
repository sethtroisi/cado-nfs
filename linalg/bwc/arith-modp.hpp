#ifndef ARITH_MODP_HPP_
#define ARITH_MODP_HPP_

#include <gmp.h>

#include "gmp-hacks.h"
#include "gmp_aux.h"
#include "macros.h"
#include "memory.h"

#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
#include <x86intrin.h>
#endif


#define  xxxDEBUG_INFINITE_LOOPS

namespace arith_modp {
namespace details {
    /* some of this is provided by <type_traits>, but that is post C++98
     */
    template<bool x> struct is_true {};
    template<> struct is_true<true> { typedef int type; };
    template<typename T, typename U> struct is_same { static const bool value = false; };
    template<typename T> struct is_same<T, T> { static const bool value = true; };
    template<typename T> struct make_signed {};
    template<> struct make_signed<unsigned long> { typedef long type; };
    template<> struct make_signed<unsigned long long> { typedef long long type; };
    typedef make_signed<mp_limb_t>::type signed_mp_limb_t;

    template<int n> struct mpn {
        typedef mpn<n> self;
        mp_limb_t x[n];
        mpn() { memset(x, 0, n * sizeof(mp_limb_t)); }
        mpn(mpn const& a) { memcpy(x, a.x, n * sizeof(mp_limb_t)); }
        mpn(mpz_srcptr a) { MPN_SET_MPZ(x, n, a); }
        self& operator=(mpz_srcptr a) { MPN_SET_MPZ(x, n, a); return *this; }
        void zero() { memset(x, 0, n * sizeof(mp_limb_t)); }
        operator mp_limb_t * () { return x; }
        operator const mp_limb_t * () const { return x; }
        static void zero(self * x, int N) {
            memset(x, 0, n * N * sizeof(mp_limb_t));
        }
        static void copy(self * y, const self * x, int N) {
            memcpy(y, x, n * N * sizeof(mp_limb_t));
        }
        bool operator==(self const& a) {
            return memcmp(x, a.x, n * sizeof(mp_limb_t)) == 0;
        }
        /*
        bool operator<(self const& a) {
            return memcmp(x, a.x, n * sizeof(mp_limb_t)) < 0;
        }
        bool operator>(self const& a) {
            return memcmp(x, a.x, n * sizeof(mp_limb_t)) > 0;
        }
        */
    };

    template<int n_, int extra_, typename T>
        struct gfp_base {
            /* This is only example code. The version used is the one which is
             * specialized for extra == 1.
             */
            static const int n = n_;
            static const typename is_true<n_>= extra_>::type extra = extra_;

            typedef mpn<n> elt;
            struct elt_ur: public mpn<n + extra> {
                elt_ur& operator=(elt const& a) {
                    mpn_copyi(mpn<n + extra>::x, a.x, n);
                    memset(mpn<n + extra>::x + n, 0, extra * sizeof(mp_limb_t));
                    return *this;
                }
            };
            struct preinv : public mpn<extra> { int shift; };

            static void propagate_carry(mp_limb_t * dst, mp_limb_t cy) {
                mpn_add_1(dst, dst, extra, cy);
            }
            static void propagate_borrow(mp_limb_t * dst, mp_limb_t cy) {
                mpn_sub_1(dst, dst, extra, cy);
            }
            static bool upperlimbs_are_zero(mp_limb_t * dst) {
                for(int i = extra ; i-- ; )
                    if (dst[i]) return false;
                return true;
            }

            static inline void stream_store(elt * dst, elt const& src) { *dst = src; }
            static inline void add(elt_ur & dst, elt const & src)
            {
                mp_limb_t cy = mpn_add_n(dst, dst, src, n);
                T::propagate_carry(dst + n, cy);
            }

            static inline void sub(elt_ur & dst, elt const & src)
            {
                mp_limb_t cy = mpn_sub_n(dst, dst, src, n);
                T::propagate_borrow(dst + n, cy);
            }

            static inline void add_ur(elt_ur & dst, elt_ur const & src)
            {
                mpn_add_n(dst, dst, src, n + extra);
            }

            static inline void sub_ur(elt_ur & dst, elt_ur const & src)
            {
                mpn_sub_n(dst, dst, src, n + extra);
            }

            static inline void addmul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&)
            {
                mp_limb_t cy = mpn_addmul_1(dst, src, n, x);
                T::propagate_carry(dst + n, cy);
            }
            static inline void submul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&)
            {
                mp_limb_t cy = mpn_submul_1(dst, src, n, x);
                T::propagate_borrow(dst + n, cy);
            }

            /* Preinverse for Barrett reduction. See also the code for reduction,
             * which is further below.
             *
             * We want to apply this to reduce a mod p, with the following
             * constraints.
             *
             *             2^(m-1) <  p  < 2^m
             *        -2^(ell-1)*p <= a  < 2^(ell-1)*p
             *
             * Let I=floor(2^(m+ell)/p). Because of the bound on p, we have 2^ell
             * < I < 2^(ell+1), so that 0<J=I-2^ell<2^ell (which actually fits
             * within ell bits). The preinverse we compute is this J.
             */

            static void compute_preinv(preinv & j, elt const & p)
            {
                mpz_t big;
                mpz_t pz;
                mpz_init_set_ui(big,1);
                mpz_init(pz);
                MPZ_SET_MPN(pz, p, n);
                size_t m = mpz_sizeinbase(pz, 2);
                ASSERT_ALWAYS(m <= (size_t) n * mp_bits_per_limb);
                size_t ell = extra * mp_bits_per_limb;
                mpz_mul_2exp(big, big, m+ell);
                mpz_fdiv_q(big, big, pz);
                ASSERT_ALWAYS(mpz_sizeinbase(big, 2) == (ell + 1));
                mpz_fdiv_r_2exp(big, big, ell);
                MPN_SET_MPZ(j, extra, big);
                j.shift = (mp_bits_per_limb - m) % mp_bits_per_limb;
                mpz_clear(big);
                mpz_clear(pz);
            }

            /* Signed Barrett reduction (extended from Brent-Zimmermann 2010,
             * theorem 2.4)
             */

            /* input: a = a0 + a1 * 2^m, with         2^(m-1) <  p  < 2^m
             *                                   -2^(ell-1)*p <= a  < 2^(ell-1)*p
             *                                              0 <= a0 < 2^m
             * which imply in particular:
             *                                     -2^(ell-1) <= a1 < 2^(ell-1)
             *
             * Case a1 >= 0.
             *
             * Let q0 = floor(a1*I/2^ell) = floor(a1*J/2^ell) + a1.
             * We have 0 <= q0 < 2^ell.
             *
             * Moreover: q0 <= a1*I/2^ell <= a1*2^m/p <= a/p, so that r0=a-p*q0>=0.
             * use p*I >= 2^(m+ell)-p and 2^ell*q0 >= a1*I-2^ell
             *
             * compute 2^ell*p*q0 >= 2^(m+ell)*a1-a1*p-2^ell*p
             *                    >= 2^ell*(a-a0)-p*(a1+2^ell)
             *                    >  2^ell*a - 4*2^ell*p
             *             a-p*q0 <  4p
             * where the third line used a1 < 2^(ell-1) and a0 <= 2^m <= 2*p.
             *
             * Case a1 < 0.
             *
             * We let b1 = a1 + 2^ell, which is the unsigned limb used to
             * represent a1.
             *
             * Let q0 = floor(a1*I/2^ell) = floor(b1*J/2^ell) + b1 - 2^ell - J.
             *
             * Since a1 < 0, we have q0 < 0. With a1 >= -2^(ell-1) and
             * I<2^(ell+1), we obtaib q0 > -2^ell. Therefore q0 is well
             * represented by the machine word
             *  q'0 = q0+2^ell = floor(b1*J/2^ell) + b1 - J
             *
             * We have p*q0 <= p*a1*I/2^ell < p*a1/2^ell*(2^(m+ell)/p-1)
             *              <  a1*2^m - a1/2^ell * p
             *              <  p*q+r-a0-a1/2^ell * p
             *         q-q0 >  (a0-r)/p + a1/2^ell
             *              >  -1.5   since a0>0, r<p, and a1>-2^(ell-1).
             *              >= -1     since q and q0 are integers.
             * So that q-(q0-1) >= 0.
             *
             * Note that because we have -2^ell < q0 < 0, then q0-1 is properly
             * represented by the unsigned machine word 2^ell-1+q0.
             *
             * Also, we have p*q0 >= p*a1*I/2^ell-p
             *                    >= a1*2^m-p
             *                    >= a-a0-p
             *         a-p*(q0-1) <= a0 + 2p < 4p
             *
             * To compute a-p*(q0-1), we actually compute
             * a-p*(q0+2^ell-1)+2^ell*p, which is a submul_ui followed by one
             * addition.
             */

            /* this reduces a in place, and copies the result to r */
            static void reduce(elt & r, elt_ur & a, elt const & p, preinv const & j)
            {
                mp_limb_t tmp[extra + 1];
                if (j.shift) {
                    mpn_lshift(tmp, a + n - 1, extra + 1, j.shift);
                } else {
                    mpn_copyi(tmp + 1, a + n, extra);
                }
                mp_limb_t a1I[2*extra];
                mpn_mul_n(a1I, tmp + 1, j, extra);
                mpn_add_n(a1I + extra, a1I + extra, tmp + 1, extra);
                mp_limb_t * q0 = a1I + extra;
                typename make_signed<mp_limb_t>::type sa1 = (tmp+1)[extra-1];
                if (sa1 < 0) {
                    mpn_sub_n(q0, q0, j, extra);
                    mpn_sub_1(q0, q0, extra, 1);
                    mpn_add_n(a + extra, a + extra, p, n);
                }
                /* emulate a submul_n ; need to do mul first, then sub... */
                mp_limb_t scratch[n + extra];
                mpn_mul(scratch, p, n, q0, extra);
                mpn_sub_n(a, a, scratch, n + extra);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (!upperlimbs_are_zero(a + n) || mpn_cmp(a, p, n) >= 0) {
                    T::sub(a, p);
                    /*
                       {
                       mp_limb_t cy = mpn_sub_n(a, a, p, n);
                       propagate_borrow(a + n, cy);
                       }
                       */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                    spin++;
                    ASSERT_ALWAYS(spin < 4);
#endif
                }
                mpn_copyi(r, a, n);
            }
        };

    template<int n_, int extra_, typename T>
        struct gfp_middle: public gfp_base<n_, extra_, T> {};

    template<int n_, typename T>
        class gfp_middle<n_, 1, T>: public gfp_base<n_, 1, T>
        {
            typedef gfp_base<n_, 1, T> super;
            public:
            static inline void propagate_carry(mp_limb_t * dst, mp_limb_t cy) {
                *dst += cy;
            }
            static inline void propagate_borrow(mp_limb_t * dst, mp_limb_t cy) {
                *dst -= cy;
            }
            static inline bool upperlimbs_are_zero(mp_limb_t * dst) {
                return !dst[0];
            }

            using super::n;
            /* We would prefer write "using typename super::elt" below,
             * however this is only ok with recent c++ */
            typedef typename super::elt elt;
            typedef typename super::elt_ur elt_ur;
            typedef typename super::preinv preinv;

            /* this reduces a in place, and copies the result to r */
            static void reduce(elt & r, elt_ur & a, elt const & p, preinv const & j) {
                mp_limb_t a1 = a[n] << j.shift;
                if (j.shift) {
                    a1 |= a[n-1] >> (mp_bits_per_limb - j.shift);
                }
                signed_mp_limb_t sa1 = a1;
                mp_limb_t tmp[2];
#ifdef  umul_ppmm
                umul_ppmm(tmp[1], tmp[0], a1, j[0]);
#elif defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
                __asm__ ("mulq %3" : "=a"(tmp[0]), "=d" (tmp[1]) : "0" (a1), "rm" (j[0]));
#else
                mpn_mul_n(tmp, &a1, j, 1);
#endif
                mp_limb_t q0 = tmp[1] + a1;
                if (sa1 < 0) {
                    /* see above for the specificities of the negative case */
                    q0 -= j[0] + 1;
                    mpn_add_n(a + 1, a + 1, p, n);
                }
                T::submul_ui(a, p, q0, p, j);
                /*
                   {
                   mp_limb_t cy = mpn_submul_1(a, p, n, q0);
                   super::propagate_borrow(a + n, cy);
                   }
                   */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (a[n] || mpn_cmp(a, p, n) >= 0) {
                    T::sub(a, p);
                    /*
                       {
                       mp_limb_t cy = mpn_sub_n(a, a, p, n);
                       super::propagate_borrow(a + n, cy);
                       }
                       */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                    spin++;
                    ASSERT_ALWAYS(spin < 4);
#endif
                }
                mpn_copyi(r, a, n);
            }
        };

    template<int n, int extra=1>
        struct gfp : public gfp_middle<n, extra, gfp<n, extra> >
    {
    };


    /* Now for some sizes, we see a clear interest in using auxiliary
     * vector types. We call these "fast" types. The general compromise
     * is that we accept types which may be a little wider, but generally
     * allow for better performance. The specs go typically as follows.
     *
     * - conversion to the "fast" types must be done for both operands
     *   (say, source vector as well as destination vector). We don't
     *   intend to go with the same kind of arithmetic that what we do
     *   with elt and elt_ur above, where a "mixed" add/sub function
     *   exists.
     *
     * - "fast" types are amenable to vector instructions
     *
     * - up to some number of additions or subtractions may be performed
     *   on the fast type before reduction.
     *
     * - type may be ambiguous (so comparison entails conversion).
     *
     * We have two natural choices:
     *
     *  - RNS representation
     *  - carry-save (aka nails).
     *
     * The specializations below work with nails. The idea is as follows.
     * For a p spanning three 64-bit words, we spread data into four
     * 48-bit words in an avx register. Then we can accumulate up to 2^16
     * of these at little cost.
     */

    /* the default version is just making no difference, so that we'll
     * use the simple elt / elt_ur mechanism */
    template<typename T> struct fast_type : public T { };


#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    /* Now some specializations */

    /* {{{ gfp<1,1> */
    template<> struct gfp<1, 1> : public gfp_middle<1,1,gfp<1,1> > {
        static inline void add(elt_ur & dst, elt const & src)
        {
            asm("# gfp<1, 1>::add\n"
                "addq %q2, %q0\n"
                "adcq $0x0, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0])
               );
        }

        static inline void sub(elt_ur & dst, elt const & src)
        {
            asm("# gfp<1, 1>::sub\n"
                "subq %q2, %q0\n"
                "sbbq $0x0, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0])
               );
        }

        static inline void addmul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&)
        {
            mp_limb_t foo, bar;
            asm("# gfp<1, 1>::addmul_ui\n"
                "mulq   %[mult]\n"
                "addq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rdx, %[z1]\n"
            : "=a"(foo), "=d"(bar), [z0]"+rm"(dst[0]), [z1]"+rm"(dst[1])
            : "0"(src[0]), [mult]"r1m"(x)
            );
        }
        static inline void submul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&)
        {
            mp_limb_t foo, bar;
            asm("# gfp<1, 1>::submul_ui\n"
                "mulq   %[mult]\n"
                "subq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rdx, %[z1]\n"
            : "=a"(foo), "=d"(bar), [z0]"+rm"(dst[0]), [z1]"+rm"(dst[1])
            : "0"(src[0]), [mult]"r1m"(x)
            );
        }
    };
    /* }}} */

    /* {{{ gfp<2,1> */
    template<> struct gfp<2, 1> : public gfp_middle<2,1,gfp<2,1> > {
        static inline void add(elt_ur & dst, elt const & src) {
            asm("# gfp<2, 1>::add\n"
                "addq %q[s0], %q[d0]\n"
                "adcq %q[s1], %q[d1]\n"
                "adcq $0x0, %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1])
               );
        }

        static inline void sub(elt_ur & dst, elt const & src) {
            asm("# gfp<2, 1>::sub\n"
                "subq %q[s0], %q[d0]\n"
                "sbbq %q[s1], %q[d1]\n"
                "sbbq $0x0, %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1])
               );
        }

        static inline void addmul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&)
        {
            mp_limb_t foo, bar;
            asm("# gfp<2, 1>::addmul_ui\n"
                "mulq   %[mult]\n"
                "addq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "movq   %%rdx, %%rcx\n"
                "movq   %[s1], %%rax\n"
                "mulq   %[mult]\n"
                "addq   %%rcx, %%rax\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rax, %[z1]\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rdx, %[z2]\n"
            : "=&a"(foo), "=&d"(bar),
            [z0]"+rm"(dst[0]),
            [z1]"+rm"(dst[1]),
            [z2]"+rm"(dst[2])
            : [s0]"0"(src[0]), [s1]"rm"(src[1]), [mult]"rm"(x)
            : "rcx"
            );
        }

        static inline void submul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&) {
            mp_limb_t foo, bar;
            asm("# gfp<2, 1>::submul_ui\n"
                "mulq   %[mult]\n"
                "subq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "movq   %%rdx, %%rcx\n"
                "movq   %[s1], %%rax\n"
                "mulq   %[mult]\n"
                "addq   %%rcx, %%rax\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rax, %[z1]\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rdx, %[z2]\n"
            : "=&a"(foo), "=&d"(bar),
            [z0]"+rm"(dst[0]),
            [z1]"+rm"(dst[1]),
            [z2]"+rm"(dst[2])
            : [s0]"0"(src[0]), [s1]"rm"(src[1]), [mult]"rm"(x)
            : "rcx"
            );
        }
    };
    /* }}} */

    /* {{{ macros for assembly for further specializations */
#define FEED_IN_WITH_S0_IN_RAX(in1, r0, r1) \
        /* status: s0 in rax */                                \
        "mulq   %[mult]\n"             /* rdx:rax = s0 * v */  \
        "movq   %%rax, %%" #r0 "\n"    /* lo contrib to d1 */  \
        "movq   " in1 ", %%rax\n"        /* load s1          */  \
        "movq   %%rdx, %%" #r1 "\n"    /* hi contrib to d1 */
#define FEED_IN(in0, in1, r0, r1) \
        "movq   " in0 ", %%rax\n"       \
        FEED_IN_WITH_S0_IN_RAX(in1, r0, r1)
#define INNER_MUL(op, out, in, r0, r1, r2)   \
        /* status: r1:r0 to be added to d_{i+1}:d_i, rax = s_{i+1} */     \
        "xorq   %%" #r2 ", %%" #r2 "\n"                                   \
        "mulq   %[mult]\n"                   /* rdx:rax = s_{i+1} * v */  \
        "" #op "q %%" #r0 ", " out "\n" /* store d_i             */   \
        "adcq   %%rax, %%" #r1 "\n"         /* lo contrib to d_{i+1} */   \
        "adcq   %%rdx, %%" #r2 "\n"         /* hi contrib to d_{i+2} */   \
        "movq   " in ", %%rax\n"       /* load s_{i+2}          */
#define FINISH(op, opc, out0, out1, out2, r0, r1) \
        /* r1:r0 to be added to d_{i+1}:d_i ; rax = s_{i+2} */	\
        "mulq   %[mult]\n"                   			\
        "" #op "q   %%" #r0 ", " out0 "\n"  			\
        "adcq   %%rax, %%" #r1 "\n"				\
        "adcq   $0x0, %%rdx\n"					\
        "" #op "q   %%" #r1 ", " out1 "\n" 			\
        "" #opc "q   %%rdx, " out2 "\n" 
    /* }}} */
    /* {{{ this macro actually exposes the specialization in itself */
#define EXPOSE_SPECIALIZATION(n)					\
    template<> struct gfp<n, 1> : public gfp_middle<n,1,gfp<n,1> > {	\
        using gfp_middle<n,1,gfp<n,1> >::add_ur;   \
        using gfp_middle<n,1,gfp<n,1> >::sub_ur;   \
        static inline void add(elt_ur & dst, elt const & src) {		\
            asm("# gfp<" #n ", 1>::add\n"				\
                    ADDSUB_CODE ## n(add, adc)				\
               );							\
        }								\
        static inline void sub(elt_ur & dst, elt const & src) {		\
            asm("# gfp<" #n ", 1>::sub\n"					\
                    ADDSUB_CODE ## n(sub, sbb)				\
               );							\
        }								\
        static inline void addmul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&) {\
            mp_limb_t foo MAYBE_UNUSED;					\
            asm ("# gfp<" #n ", 1>::addmul_ui\n"				\
                    ADDSUBMUL_CODE ## n(add, adc)			\
            );								\
        }								\
        static inline void submul_ui(elt_ur & dst, elt const & src, mp_limb_t x, elt const&, preinv const&) {\
            mp_limb_t foo MAYBE_UNUSED;					\
            asm("# gfp<" #n ", 1>::submul_ui\n"				\
                    ADDSUBMUL_CODE ## n(sub, sbb)			\
            );								\
        }								\
    }
    /* }}} */

    /* {{{ code for gfp<3, 1> */
#define ADDSUBMUL_CODE3(op, opc)					\
                FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)			\
                INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)		\
                FINISH(op, opc, "%[z1]", "%[z2]", "%[z3]", r9, r10)	\
                : "=&a"(foo),                                           \
                    [z0]"+rm"(dst[0]),				\
                    [z1]"+rm"(dst[1]),				\
                    [z2]"+rm"(dst[2]),				\
                    [z3]"+rm"(dst[3])					\
                :							\
                    [s0]"0"(src[0]),					\
                    [s1]"rm"(src[1]),					\
                    [s2]"rm"(src[2]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "rdx"

#define ADDSUB_CODE3(op, opc)   \
        "" #op  "q %q[s0], %q[d0]\n"					\
        "" #opc "q %q[s1], %q[d1]\n"					\
        "" #opc "q %q[s2], %q[d2]\n"					\
        "" #opc "q $0x0, %q[d3]\n"					\
                :							\
                [d0]"+rm"(dst[0]),					\
                [d1]"+rm"(dst[1]),					\
                [d2]"+rm"(dst[2]),					\
                [d3]"+rm"(dst[3])					\
                :							\
                [s0]"r"(src[0]),					\
                [s1]"r"(src[1]),					\
                [s2]"r"(src[2])

    EXPOSE_SPECIALIZATION(3);
    /* }}} */

    /* {{{ code for gfp<4, 1> */
    /*
#define ADDSUBMUL_CODE4(op, opc)					\
                FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)			\
                INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)		\
                INNER_MUL(op, "%[z1]", "%[s3]", r9, r10, r11)		\
                FINISH(op, opc, "%[z2]", "%[z3]", "%[z4]", r10, r11)	\
                : "=&a"(foo),                                           \
                    [z0]"+rm"(dst[0]),				\
                    [z1]"+rm"(dst[1]),				\
                    [z2]"+rm"(dst[2]),				\
                    [z3]"+rm"(dst[3]),				\
                    [z4]"+rm"(dst[4])					\
                :							\
                    [s0]"0"(src[0]),					\
                    [s1]"rm"(src[1]),					\
                    [s2]"rm"(src[2]),					\
                    [s3]"rm"(src[3]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rdx"

#define xADDSUB_CODE4(op, opc)   \
"" #op  "q %q[s0], (%[z])\n"					\
"" #opc "q %q[s1], 0x8(%[z])\n"					\
"" #opc "q %q[s2], 0x10(%[z])\n"				\
"" #opc "q %q[s3], 0x18(%[z])\n"				\
"" #opc "q $0x0, 0x20(%[z])\n"					\
        :							\
        :							\
            [z]"r"(&dst[0]),				        \
            [s0]"r"(src[0]),					\
            [s1]"r"(src[1]),					\
            [s2]"r"(src[2]),					\
            [s3]"r"(src[3])                                   \
        : "memory"


                */
#define ADDSUB_CODE4(op, opc)   \
        "" #op  "q %q[s0], %q[d0]\n"					\
        "" #opc "q %q[s1], %q[d1]\n"					\
        "" #opc "q %q[s2], %q[d2]\n"					\
        "" #opc "q %q[s3], %q[d3]\n"					\
        "" #opc "q $0x0, %q[d4]\n"					\
                :							\
                [d0]"+rm"(dst[0]),					\
                [d1]"+rm"(dst[1]),					\
                [d2]"+rm"(dst[2]),					\
                [d3]"+rm"(dst[3]),					\
                [d4]"+rm"(dst[4])					\
                :							\
                [s0]"r"(src[0]),					\
                [s1]"r"(src[1]),					\
                [s2]"r"(src[2]),					\
                [s3]"r"(src[3])


#define ADDSUBMUL_CODE4(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                FINISH(op, opc,                                         \
                        "0x10(%[z])", "0x18(%[z])", "0x20(%[z])",       \
                        r10, r11)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

    EXPOSE_SPECIALIZATION(4);
    /* }}} */

    /* {{{ code for gfp<5, 1> */

#define ADDSUBMUL_CODE5(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                FINISH(op, opc,                                         \
                        "0x18(%[z])", "0x20(%[z])", "0x28(%[z])",       \
                        r11, r8)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE5(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q $0x0, 0x28(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(5);
    /* }}} */

    /* {{{ code for gfp<6, 1> */

#define ADDSUBMUL_CODE6(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                FINISH(op, opc,                                         \
                        "0x20(%[z])", "0x28(%[z])", "0x30(%[z])",       \
                        r8, r9)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE6(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q $0x0, 0x30(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(6);
    /* }}} */

    /* {{{ code for gfp<7, 1> */
#define ADDSUBMUL_CODE7(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                INNER_MUL(op, "0x20(%[z])", "0x30(%[s])", r8, r9, r10)	\
                FINISH(op, opc,                                         \
                        "0x28(%[z])", "0x30(%[z])", "0x38(%[z])",       \
                        r9, r10)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE7(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q %q[s6], 0x30(%[z])\n"				\
        "" #opc "q $0x0, 0x38(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5]),                                  \
                    [s6]"r"(src[6])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(7);
    /* }}} */

    /* {{{ code for gfp<8, 1> */
#define ADDSUBMUL_CODE8(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                INNER_MUL(op, "0x20(%[z])", "0x30(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x28(%[z])", "0x38(%[s])", r9, r10, r11)	\
                FINISH(op, opc,                                         \
                        "0x30(%[z])", "0x38(%[z])", "0x40(%[z])",       \
                        r10, r11)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE8(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q %q[s6], 0x30(%[z])\n"				\
        "" #opc "q %q[s7], 0x38(%[z])\n"				\
        "" #opc "q $0x0, 0x40(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5]),                                  \
                    [s6]"r"(src[6]),                                   \
                    [s7]"r"(src[7])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(8);
    /* }}} */
    
    /* further specialization only seem to bring very marginal
     * improvements. */


#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
    template<> struct fast_type<gfp<3, 1> > {
        typedef gfp<3, 1> super;
        struct elt;
        typedef elt elt_ur;
        struct elt {
            typedef elt self;
#ifdef  HAVE_AVX2
            __m256i data[1];
#else
            __m128i data[2];
#endif
            elt() { zero(); }
            elt(elt const& a) {
                data[0] = a.data[0];
#ifndef HAVE_AVX2
                data[1] = a.data[1];
#endif
            }
            /* we do not construct (nor affect) from mpz, because we're not
             * positional */
            void zero() {
#ifdef  HAVE_AVX2
                data[0] = _mm256_setzero_si256();
#else
                data[0] = _mm_setzero_si128();
                data[1] = _mm_setzero_si128();
#endif
            }
            static void zero(elt * x, int N) {
                memset(x, 0, N * sizeof(data));
            }
            static void copy(elt * y, const elt * x, int N) {
                memcpy(y, x, N * sizeof(data));
            }
            bool operator==(elt const& a) {
                return memcmp(data, a.data, sizeof(data)) == 0;
            }
            elt(super::elt const& a) {
                convert(*this, a);
            }

            operator super::elt_ur() const {
                super::elt_ur carries(conv_backend_get_carries(*this));
                super::add(carries, conv_backend_get_main(*this));
                return carries;
            }

            /* same, but we assume carry is zero */
            operator super::elt() const {
                return conv_backend_get_main(*this);
            }
        };

        static inline void stream_store(elt * dst, elt const& src) {
            /* Do we want to stream that or not ? In fact it's slower
             * when streaming... */
#if 0
#ifdef  HAVE_AVX2
            _mm256_stream_si256(dst->data+0, src.data[0]);
#else
            _mm_stream_si128(dst->data+0, src.data[0]);
            _mm_stream_si128(dst->data+1, src.data[1]);
#endif
#else
#ifdef  HAVE_AVX2
            _mm256_storeu_si256(dst->data+0, src.data[0]);
#else
            _mm_storeu_si128(dst->data+0, src.data[0]);
            _mm_storeu_si128(dst->data+1, src.data[1]);
#endif
#endif
        }
        static inline void add(elt & dst, elt const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_add_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_add_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_add_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub(elt & dst, elt const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub_ur(elt_ur & dst, elt_ur const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        /* conversions are done as a combination of blend & shuffle */

        /* case of 192 bits within 256 bits. Three 64-bit words
         * split into four 48-bit words.
         */
        static void convert(elt& dst, const super::elt& a) {
            /* index of 16-bit word in destination, fetched from
             * which index of 16-bit word in the gfp::elt. This is
             * given for the 256-bit registers
             *
             * 0    0
             * 1    1
             * 2    2
             * 3    <empty>
             * 4    3
             * 5    4
             * 6    5
             * 7    <empty>
             * 8    6
             * 9    7
             * 10   8
             * 11   <empty>
             * 12   9
             * 13   10
             * 14   11
             * 15   <empty>
             */
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_shuffle_epi8(
                    _mm256_loadu_si256((__m128i*) a.x),
                    _mm256_setr_epi8( 
                        0,1,2,3,4,5,-1,-1,
                        6,7,8,9,10,11,-1,-1,
                        12,13,14,15,16,17,-1,-1,
                        18,19,20,21,22,23,-1,-1));

#else   /* SSSE3 !! */
            __m128i lo = _mm_loadu_si128((__m128i*) a.x);
            __m128i hi = _mm_loadu_si128((__m128i*) (a.x + 2));
            /* note that 16bit-wide shuffles use an 8-bit immediate,
             * but do not offer the option to selectively insert
             * zeroes. So we're probably better off shuffling bytes.
             */
            dst.data[0] = _mm_shuffle_epi8(lo, _mm_setr_epi8( 
                        0,1,2,3,4,5,-1,-1,
                        6,7,8,9,10,11,-1,-1));
            dst.data[1] = _mm_xor_si128(
                    _mm_shuffle_epi8(lo, _mm_setr_epi8( 
                            12,13,14,15,-1,-1,-1,-1,
                            -1,-1,-1,-1, -1,-1,-1,-1)),
                    _mm_shuffle_epi8(hi, _mm_setr_epi8( 
                            -1,-1,-1,-1, 0, 1,-1,-1,
                            2, 3, 4, 5, 6, 7, -1, -1)));
#endif
        }

        static super::elt conv_backend_get_main(elt const& src) {
            super::elt main;
#ifdef  HAVE_AVX2
            _mm256_storeu_si256((__m256i*) main.x,
                    _mm256_shuffle_epi8(src.data[0],
                        _mm256_setr_epi8(
                            0,1,2,3,4,5,
                            8,9,10,11,12,13,
                            16,17,18,19,20,21,
                            24,25,26,27,28,29,
                            -1,-1,-1,-1,-1,-1,-1,-1)));
#else
            _mm_storeu_si128((__m128i*) main.x,
                    _mm_xor_si128(
                        _mm_shuffle_epi8(src.data[0],
                            _mm_setr_epi8(
                                0,1,2,3,4,5,
                                8,9,10,11,12,13,
                                -1,-1,-1,-1)),
                        _mm_shuffle_epi8(src.data[1],
                            _mm_setr_epi8(
                                -1,-1,-1,-1,-1,-1,
                                -1,-1,-1,-1,-1,-1,
                                0,1,2,3))));
            _mm_storeu_si128((__m128i*) (main.x + 2),
                    _mm_shuffle_epi8(src.data[1],
                        _mm_setr_epi8(
                            4,5,
                            8,9,10,11,12,13,
                            -1,-1,-1,-1,-1,-1,-1,-1)));
#endif
            return main;
        }
        static super::elt_ur conv_backend_get_carries(elt const& src) {
            super::elt_ur carries;
#ifdef  HAVE_AVX2
            _mm256_storeu_si256((__m256i*) carries.x,
                    _mm256_shuffle_epi8(src.data[0],
                        _mm256_setr_epi8(
                            -1,-1,-1,-1,-1,-1,
                            6,7,
                            -1,-1,-1,-1,
                            14,15,
                            -1,-1,-1,-1,
                            22,23,
                            -1,-1,-1,-1,
                            30,31,
                            -1,-1,-1,-1,-1,-1
                            )));
#else
            _mm_storeu_si128((__m128i*) carries.x,
                    _mm_shuffle_epi8(src.data[0],
                        _mm_setr_epi8(
                            -1,-1,-1,-1,-1,-1,
                            6,7,
                            -1,-1,-1,-1,
                            14,15,
                            -1,-1)));
            _mm_storeu_si128((__m128i*) carries.x,
                    _mm_shuffle_epi8(src.data[1],
                        _mm_setr_epi8(
                            -1,-1,
                            6,7,
                            -1,-1,-1,-1,
                            14,15,
                            -1,-1,-1,-1,-1,-1
                            )));
#endif
            return carries;
        }



        /* (add|sub)mul_ui go through convert, do naively and convert
         * back. Yes, it's slightly painful. Here we assume that src
         * has undergone little to no accumulated additions, so that
         * it can basically be converted lossless to a gfp::elt
         */
        static inline void addmul_ui(elt & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            super::elt_ur z(dst);
            super::addmul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }
        static inline void submul_ui(elt_ur & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            super::elt_ur z(dst);
            super::submul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }

        /* we have *TWO* reduction functions here. One which assigns to a
         * standard gfp::elt, and one which assigns to fast_type::elt */
        static void reduce(super::elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            super::elt_ur z(a);
            super::reduce(r, z, p, j);
        }
        static void reduce(elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            reduce(zr, a, p, j);
            r = zr;
        }
    };
#endif  /* defined(HAVE_AVX2) || defined(HAVE_SSSE3) */
#endif
    }

/* expose only what we have in our public interface */
using details::gfp;
using details::fast_type;
}



#endif
