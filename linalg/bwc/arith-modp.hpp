#ifndef ARITH_MODP_HPP_
#define ARITH_MODP_HPP_

#include <gmp.h>

#include "gmp-hacks.h"
#include "macros.h"

#define  DEBUG_INFINITE_LOOPS

namespace arith_modp_details {
    template<bool x> struct is_true {};
    template<> struct is_true<true> { typedef int type; };
    typedef typename std::make_signed<mp_limb_t>::type signed_mp_limb_t;

    template<int n_, int extra_, typename T>
        struct gfp_base {
            /* This is only example code. The version used is the one which is
             * specialized for extra == 1.
             */
            static const int n = n_;
            static const typename is_true<n_>= extra_>::type extra = extra_;

            struct elt { mp_limb_t x[n]; };
            struct elt_ur { mp_limb_t x[n + extra]; };
            struct preinv { mp_limb_t x[extra]; int shift; };

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

            static inline void add(elt_ur & dst, elt const & src)
            {
                mp_limb_t cy = mpn_add_n(dst.x, dst.x, src.x, n);
                T::propagate_carry(dst.x + n, cy);
            }

            static inline void sub(elt_ur & dst, elt const & src)
            {
                mp_limb_t cy = mpn_sub_n(dst.x, dst.x, src.x, n);
                T::propagate_borrow(dst.x + n, cy);
            }

            static inline void addmul(elt_ur & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t cy = mpn_addmul_1(dst.x, src.x, n, x);
                T::propagate_carry(dst.x + n, cy);
            }
            static inline void submul(elt_ur & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t cy = mpn_submul_1(dst.x, src.x, n, x);
                T::propagate_borrow(dst.x + n, cy);
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
                MPZ_SET_MPN(pz, p.x, n);
                size_t m = mpz_sizeinbase(pz, 2);
                ASSERT_ALWAYS(m <= (size_t) n * mp_bits_per_limb);
                size_t ell = extra * mp_bits_per_limb;
                mpz_mul_2exp(big, big, m+ell);
                mpz_fdiv_q(big, big, pz);
                ASSERT_ALWAYS(mpz_sizeinbase(big, 2) == (ell + 1));
                mpz_fdiv_r_2exp(big, big, ell);
                MPN_SET_MPZ(j.x, extra, big);
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
             * a-p*(q0+2^ell-1)+2^ell*p, which is a submul followed by one
             * addition.
             */


            /* this reduces a in place, and copies the result to r */
            static void reduce(elt & r, elt_ur & a, elt const & p, preinv const & j)
            {
                mp_limb_t tmp[extra + 1];
                if (j.shift) {
                    mpn_lshift(tmp, a.x + n - 1, extra + 1, j.shift);
                } else {
                    mpn_copyi(tmp + 1, a.x + n, extra);
                }
                mp_limb_t a1I[2*extra];
                mpn_mul_n(a1I, tmp + 1, j.x, extra);
                mpn_add_n(a1I + extra, a1I + extra, tmp + 1, extra);
                mp_limb_t * q0 = a1I + extra;
                typename std::make_signed<mp_limb_t>::type sa1 = (tmp+1)[extra-1];
                if (sa1 < 0) {
                    mpn_sub_n(q0, q0, j.x, extra);
                    mpn_sub_1(q0, q0, extra, 1);
                    mpn_add_n(a.x + extra, a.x + extra, p.x, n);
                }
                /* emulate a submul_n ; need to do mul first, then sub... */
                mp_limb_t scratch[n + extra];
                mpn_mul(scratch, p.x, n, q0, extra);
                mpn_sub_n(a.x, a.x, scratch, n + extra);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (!upperlimbs_are_zero(a.x + n) || mpn_cmp(a.x, p.x, n) >= 0) {
                    T::sub(a, p);
                    /*
                       {
                       mp_limb_t cy = mpn_sub_n(a.x, a.x, p.x, n);
                       propagate_borrow(a.x + n, cy);
                       }
                       */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                    spin++;
                    ASSERT_ALWAYS(spin < 4);
#endif
                }
                mpn_copyi(r.x, a.x, n);
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
            using typename super::elt;
            using typename super::elt_ur;
            using typename super::preinv;

            /* this reduces a in place, and copies the result to r */
            static void reduce(elt & r, elt_ur & a, elt const & p, preinv const & j) {
                mp_limb_t a1 = a.x[n] << j.shift;
                if (j.shift) {
                    a1 |= a.x[n-1] >> (mp_bits_per_limb - j.shift);
                }
                signed_mp_limb_t sa1 = a1;
                mp_limb_t tmp[2];
#ifdef  umul_ppmm
                umul_ppmm(tmp[1], tmp[0], a1, j.x[0]);
#elif defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
                __asm__ ("mulq %3" : "=a"(tmp[0]), "=d" (tmp[1]) : "0" (a1), "rm" (j.x[0]));
#else
                mpn_mul_n(tmp, &a1, j.x, 1);
#endif
                mp_limb_t q0 = tmp[1] + a1;
                if (sa1 < 0) {
                    /* see above for the specificities of the negative case */
                    q0 -= j.x[0] + 1;
                    mpn_add_n(a.x + 1, a.x + 1, p.x, n);
                }
                T::submul(a, p, q0);
                /*
                   {
                   mp_limb_t cy = mpn_submul_1(a.x, p.x, n, q0);
                   super::propagate_borrow(a.x + n, cy);
                   }
                   */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (a.x[n] || mpn_cmp(a.x, p.x, n) >= 0) {
                    T::sub(a, p);
                    /*
                       {
                       mp_limb_t cy = mpn_sub_n(a.x, a.x, p.x, n);
                       super::propagate_borrow(a.x + n, cy);
                       }
                       */
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                    spin++;
                    ASSERT_ALWAYS(spin < 4);
#endif
                }
                mpn_copyi(r.x, a.x, n);
            }
        };

    template<int n, int extra=1>
        struct gfp : public gfp_middle<n, extra, gfp<n, extra>>
    {
    };


#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    /* Now some specializations */
    template<>
        struct gfp<1, 1> : public gfp_middle<1,1,gfp<1,1> > {
            static inline void add(elt_ur & dst, elt const & src)
            {
                asm("# gfp<1, 1>::add\n"
                    "addq %q2, %q0\n"
                    "adcq $0x0, %q1\n"
                    : "+r"(dst.x[0]), "+r"(dst.x[1])
                    : "rm"(src.x[0])
                   );
            }

            static inline void sub(elt_ur & dst, elt const & src)
            {
                asm("# gfp<1, 1>::sub\n"
                    "subq %q2, %q0\n"
                    "sbbq $0x0, %q1\n"
                    : "+r"(dst.x[0]), "+r"(dst.x[1])
                    : "rm"(src.x[0])
                   );
            }

            static inline void addmul(elt_ur & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t foo, bar;
                asm("# gfp<1, 1>::addmul\n"
                    "mulq    %[mult]\n"
                    "addq    %%rax, %[z0]\n"
                    "adcq    $0, %%rdx\n"
                    "addq    %%rdx, %[z1]\n"
                : "=a"(foo), "=d"(bar), [z0]"+rm"(dst.x[0]), [z1]"+rm"(dst.x[1])
                : "0"(src.x[0]), [mult]"r1m"(x)
                );
            }
            static inline void submul(elt_ur & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t foo, bar;
                asm("# gfp<1, 1>::submul\n"
                    "mulq    %[mult]\n"
                    "subq    %%rax, %[z0]\n"
                    "adcq    $0, %%rdx\n"
                    "subq    %%rdx, %[z1]\n"
                : "=a"(foo), "=d"(bar), [z0]"+rm"(dst.x[0]), [z1]"+rm"(dst.x[1])
                : "0"(src.x[0]), [mult]"r1m"(x)
                );
            }
        };
#endif
}

/* expose only what we have in our public interface */
using arith_modp_details::gfp;



#endif
