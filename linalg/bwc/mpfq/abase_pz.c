/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_pz.h"

#include <inttypes.h>
#include <limits.h>
static int abase_pz_impl_mpi_attr;     /* for MPI functions */
static MPI_Datatype abase_pz_impl_mpi_datatype;
static MPI_Datatype abase_pz_impl_mpi_datatype_ur;
static MPI_Op abase_pz_impl_mpi_addition_op;
static MPI_Op abase_pz_impl_mpi_addition_op_ur;
static int abase_pz_impl_mpi_use_count;   /* several stacked init()/clear() pairs are supported */
/* Active handler: simd_pz */
/* Automatically generated code  */
/* Active handler: pz */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
   fieldtype=prime,
   tag=pz,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
     ],
    member_templates_restrict={
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
     u64k1=[ u64k1, u64k2, ],
     u64k2=[ u64k1, u64k2, ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=abase_vbase,
    global_prefix=abase_,
    name=abase_vbase,
    substitutions=[
     [ (?^:abase_pz_elt \*), void *, ],
     [ (?^:abase_pz_src_elt\b), const void *, ],
     [ (?^:abase_pz_elt\b), void *, ],
     [ (?^:abase_pz_dst_elt\b), void *, ],
     [ (?^:abase_pz_elt_ur \*), void *, ],
     [ (?^:abase_pz_src_elt_ur\b), const void *, ],
     [ (?^:abase_pz_elt_ur\b), void *, ],
     [ (?^:abase_pz_dst_elt_ur\b), void *, ],
     [ (?^:abase_pz_vec \*), void *, ],
     [ (?^:abase_pz_src_vec\b), const void *, ],
     [ (?^:abase_pz_vec\b), void *, ],
     [ (?^:abase_pz_dst_vec\b), void *, ],
     [ (?^:abase_pz_vec_ur \*), void *, ],
     [ (?^:abase_pz_src_vec_ur\b), const void *, ],
     [ (?^:abase_pz_vec_ur\b), void *, ],
     [ (?^:abase_pz_dst_vec_ur\b), void *, ],
     [ (?^:abase_pz_poly \*), void *, ],
     [ (?^:abase_pz_src_poly\b), const void *, ],
     [ (?^:abase_pz_poly\b), void *, ],
     [ (?^:abase_pz_dst_poly\b), void *, ],
     ],
    },
   vtag=pz,
   w=64,
   } */


/* Functions operating on the field structure */
/* *pz::code_for_field_characteristic */
void abase_pz_field_characteristic(abase_pz_dst_field k, mpz_t z)
{
        int i;
        int n = k->kl;
        mpz_set_ui(z, k->p[n - 1]);
        for (i = n - 2; i >= 0; --i) {
        mpz_mul_2exp(z, z, 64);
        mpz_add_ui(z, z, k->p[i]);
        }
}

/* *pz::code_for_field_init */
void abase_pz_field_init(abase_pz_dst_field k)
{
        k->p = NULL;
        k->bigmul_p = NULL;
        k->io_base = 10;
        mpz_init(k->factor);
        k->ts_info.e = 0;
}

/* *pz::code_for_field_clear */
void abase_pz_field_clear(abase_pz_dst_field k)
{
        free(k->p);
        free(k->bigmul_p);
        if (k->ts_info.e > 0) {
        free(k->ts_info.hh);
        free(k->ts_info.z);
        }
        mpz_clear(k->factor);
}

/* *pz::code_for_field_specify */
void abase_pz_field_specify(abase_pz_dst_field k, unsigned long dummy MAYBE_UNUSED, void * vp)
{
        mpz_ptr p = (mpz_ptr) vp;
        k->kl = mpz_size(p);
        k->p = (mp_limb_t *) malloc(k->kl * sizeof(mp_limb_t));
        if (!k->p)
        MALLOC_FAILED();
        for (unsigned int i = 0; i < k->kl; ++i)
        k->p[i] = mpz_getlimbn(p, i);
        k->url = 1 + 2 * k->kl;
        k->url_margin = LONG_MAX;
        // precompute bigmul_p = largest multiple of p that fits in an elt_ur,
        //   p*Floor( (2^(nn*w)-1)/p )
        if (k->bigmul_p == NULL)
        k->bigmul_p = (mp_limb_t *) malloc(k->url * sizeof(mp_limb_t));
        if (!k->bigmul_p)
        MALLOC_FAILED();
        abase_pz_elt_ur big;
        abase_pz_elt_ur_init(k, &big);
        mp_limb_t q[k->url - k->kl + 1], r[k->kl], tmp[k->url + 1];
    
        for (unsigned int i = 0; i < k->url; ++i)
        big[i] = ~0UL;
        mpn_tdiv_qr(q, r, 0, big, k->url, k->p, k->kl);
        mpn_mul(tmp, q, k->url - k->kl + 1, k->p, k->kl);
        for (unsigned int i = 0; i < k->url; ++i)
        (k->bigmul_p)[i] = tmp[i];
        abase_pz_elt_ur_clear(k, &big);
        assert(tmp[k->url] == 0UL);
}


/* Element allocation functions */
/* *pz::code_for_init */
void abase_pz_init(abase_pz_dst_field k, abase_pz_elt * x)
{
        *x = (abase_pz_elt) malloc(k->kl * sizeof(mp_limb_t));
        if (!(*x))
        MALLOC_FAILED();
}

/* *pz::code_for_clear */
void abase_pz_clear(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_elt * x)
{
        free(*x);
}


/* Elementary assignment functions */
/* *pz::code_for_set_mpn */
void abase_pz_set_mpn(abase_pz_dst_field k, abase_pz_dst_elt z, mp_limb_t * x, size_t n)
{
        if (n < k->kl) {
        memcpy(z, x, n * sizeof(mp_limb_t));
        memset(z + n, 0, (k->kl - n) * sizeof(mp_limb_t));
        } else {
        mp_limb_t *tmp;
        tmp = (mp_limb_t *) malloc((n + 1 - k->kl) * sizeof(mp_limb_t));
        if (!tmp)
            MALLOC_FAILED();
        mpn_tdiv_qr(tmp, z, 0, x, n, k->p, k->kl);
        free(tmp);
        }
}

/* *pz::code_for_set_mpz */
void abase_pz_set_mpz(abase_pz_dst_field k, abase_pz_dst_elt z, mpz_t x)
{
        if (mpz_sgn(x) < 0) {
        abase_pz_set_mpn(k, z, x->_mp_d, -x->_mp_size);
        abase_pz_neg(k, z, z);
        } else
        abase_pz_set_mpn(k, z, x->_mp_d, x->_mp_size);
}


/* Assignment of random values */
/* *pz::code_for_random */
void abase_pz_random(abase_pz_dst_field k, abase_pz_dst_elt z, gmp_randstate_t state)
{
        mpz_t zz;
        mpz_init(zz);
        mpz_urandomb(zz, state, k->kl * GMP_LIMB_BITS);
        abase_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
}

/* *pz::code_for_random2 */
void abase_pz_random2(abase_pz_dst_field k, abase_pz_dst_elt z, gmp_randstate_t state)
{
        mpz_t zz;
        mpz_init(zz);
        mpz_rrandomb(zz, state, k->kl * GMP_LIMB_BITS);
        abase_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
}


/* Arithmetic operations on elements */
/* *pz::code_for_add */
void abase_pz_add(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, abase_pz_src_elt y)
{
        mp_limb_t cy;
        cy = mpn_add_n(z, x, y, k->kl);
        if (cy || mpn_cmp(z, k->p, k->kl) >= 0)
        mpn_sub_n(z, z, k->p, k->kl);
}

/* *pz::code_for_sub */
void abase_pz_sub(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, abase_pz_src_elt y)
{
        mp_limb_t cy;
        cy = mpn_sub_n(z, x, y, k->kl);
        if (cy)
        mpn_add_n(z, z, k->p, k->kl);
}

/* *pz::code_for_neg */
void abase_pz_neg(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x)
{
        if (abase_pz_is_zero(k, x))
        abase_pz_set_zero(k, z);
        else
        mpn_sub_n(z, k->p, x, k->kl);
}

/* *pz::code_for_mul */
void abase_pz_mul(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, abase_pz_src_elt y)
{
        abase_pz_elt_ur zz;
        abase_pz_elt_ur_init(k, &zz);
        abase_pz_mul_ur(k, zz, x, y);
        abase_pz_reduce(k, z, zz);
        abase_pz_elt_ur_clear(k, &zz);
}

/* *pz::code_for_sqr */
void abase_pz_sqr(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x)
{
        abase_pz_elt_ur zz;
        abase_pz_elt_ur_init(k, &zz);
        abase_pz_sqr_ur(k, zz, x);
        abase_pz_reduce(k, z, zz);
        abase_pz_elt_ur_clear(k, &zz);
}

/* *pz::code_for_is_sqr */
int abase_pz_is_sqr(abase_pz_dst_field k, abase_pz_src_elt x)
{
        mpz_t a, p;
        mpz_init(a);
        mpz_init(p);
        abase_pz_field_characteristic(k, p);
        abase_pz_get_mpz(k, a, x);
        int ret = mpz_jacobi(a, p);
        mpz_clear(a);
        mpz_clear(p);
        return (ret >= 0);
}

/* missing sqrt */
/* *pz::code_for_pow */
void abase_pz_pow(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, unsigned long * n, size_t nl)
{
        mpz_t p, xx, nn;
        mpz_init(p);
        mpz_init(xx);
        abase_pz_field_characteristic(k, p);
        abase_pz_get_mpz(k, xx, x);
        nn->_mp_d = n;
        nn->_mp_size = nl;
        nn->_mp_alloc = nl;
        mpz_powm(xx, xx, nn, p);
        abase_pz_set_mpz(k, z, xx);
        mpz_clear(p);
        mpz_clear(xx);
}

/* *pz::code_for_add_ui */
void abase_pz_add_ui(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, unsigned long y)
{
        mp_limb_t cy;
        cy = mpn_add_1(z, x, k->kl, y);
        if (cy || mpn_cmp(z, k->p, k->kl) >= 0)
        mpn_sub_n(z, z, k->p, k->kl);
}

/* *pz::code_for_sub_ui */
void abase_pz_sub_ui(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, unsigned long y)
{
        mp_limb_t cy;
        cy = mpn_sub_1(z, x, k->kl, y);
        if (cy)
        mpn_add_n(z, z, k->p, k->kl);
}

/* *pz::code_for_mul_ui */
void abase_pz_mul_ui(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x, unsigned long y)
{
        mpz_t xx, p;
        mpz_init(xx);
        mpz_init(p);
        abase_pz_field_characteristic(k, p);
        abase_pz_get_mpz(k, xx, x);
        mpz_mul_ui(xx, xx, y);
        mpz_tdiv_r(xx, xx, p);
        abase_pz_set_mpz(k, z, xx);
        mpz_clear(xx);
        mpz_clear(p);
}

/* *pz::code_for_inv */
int abase_pz_inv(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x)
{
        mpz_t xx, p;
        mpz_init(xx);
        mpz_init(p);
        abase_pz_field_characteristic(k, p);
        abase_pz_get_mpz(k, xx, x);
        int ret = mpz_invert(xx, xx, p);
        abase_pz_set_mpz(k, z, xx);
        mpz_clear(xx);
        mpz_clear(p);
        if (!ret) {
        if (abase_pz_is_zero(k, x))
            return 0;
        else {
            fprintf(stderr, "Not implemented\n");
            exit(1);
        }
        }
        return 1;
}


/* Operations involving unreduced elements */
/* *pz::code_for_elt_ur_init */
void abase_pz_elt_ur_init(abase_pz_dst_field k, abase_pz_elt_ur * x)
{
        *x = (abase_pz_elt_ur) malloc(k->url * sizeof(mp_limb_t));
        if (!(*x))
        MALLOC_FAILED();
}

/* *pz::code_for_elt_ur_clear */
void abase_pz_elt_ur_clear(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_elt_ur * x)
{
        free(*x);
}

/* *pz::code_for_elt_ur_add */
void abase_pz_elt_ur_add(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt_ur x, abase_pz_src_elt_ur y)
{
        mpn_add_n(z, x, y, k->url);
}

/* *pz::code_for_elt_ur_neg */
void abase_pz_elt_ur_neg(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt_ur x)
{
        mpn_neg(z, x, k->url);
}

/* *pz::code_for_elt_ur_sub */
void abase_pz_elt_ur_sub(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt_ur x, abase_pz_src_elt_ur y)
{
        mpn_sub_n(z, x, y, k->url);
}

/* *pz::code_for_mul_ur */
void abase_pz_mul_ur(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt x, abase_pz_src_elt y)
{
        mpn_mul_n(z, x, y, k->kl);
        z[k->url - 1] = 0;
}

/* *pz::code_for_sqr_ur */
void abase_pz_sqr_ur(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt x)
{
        mpn_sqr(z, x, k->kl);
        z[k->url - 1] = 0;
}

/* *pz::code_for_reduce */
void abase_pz_reduce(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_dst_elt_ur x)
{
    //    int neg = 0;
        if (x[k->url - 1] >> (GMP_LIMB_BITS - 1)) {	// negative input
    //        neg = 1;
    //        abase_pz_elt_ur_neg(k, x, x);
        mpn_add_n(x, x, k->bigmul_p, k->url);
        }
        mp_limb_t *tmp;
        tmp = (mp_limb_t *) malloc((k->url + 1) * sizeof(mp_limb_t));
        if (!tmp)
        MALLOC_FAILED();
        mpn_tdiv_qr(tmp, z, 0, x, k->url, k->p, k->kl);
    //    if (neg)
    //        abase_pz_neg(k, z, z);
        free(tmp);
}

/* *pz::code_for_normalize */
void abase_pz_normalize(abase_pz_dst_field k, abase_pz_dst_elt z)
{
        mp_limb_t tmp;
        mpn_tdiv_qr(&tmp, z, 0, z, k->kl, k->p, k->kl);
}


/* Comparison functions */

/* Input/output functions */
/* *pz::code_for_asprint */
void abase_pz_asprint(abase_pz_dst_field k, char * * pstr, abase_pz_src_elt x)
{
        // deal with 0
        if (abase_pz_is_zero(k, x)) {
        *pstr = (char *) malloc(2);
        *pstr[0] = '0';
        *pstr[1] = '\0';
        return;
        }
        // allocate enough room for base 2 conversion.
        *pstr = (char *) malloc(((k->kl) * 64 + 1));
        if (*pstr == NULL)
        MALLOC_FAILED();
        mp_limb_t *tmp;
        tmp = (mp_limb_t *) malloc((k->kl) * sizeof(mp_limb_t));
        if (!tmp)
        MALLOC_FAILED();
        int tl = k->kl - 1;
        while (x[tl] == 0)		// x is non-zero, no need to test tl>0
        tl--;
        for (int i = 0; i <= tl; ++i)
        tmp[i] = x[i];
        int n = mpn_get_str((unsigned char *) (*pstr), k->io_base, tmp, tl + 1);
        free(tmp);
        for (int i = 0; i < n; ++i)
        (*pstr)[i] += '0';
        (*pstr)[n] = '\0';
        // remove leading 0s
        int shift = 0;
        while (((*pstr)[shift] == '0') && ((*pstr)[shift + 1] != '\0'))
        shift++;
        if (shift > 0) {
        int i = 0;
        while ((*pstr)[i + shift] != '\0') {
            (*pstr)[i] = (*pstr)[i + shift];
            i++;
        }
        (*pstr)[i] = '\0';
        }
}

/* *pz::code_for_fprint */
void abase_pz_fprint(abase_pz_dst_field k, FILE * file, abase_pz_src_elt x)
{
        char *str;
        abase_pz_asprint(k, &str, x);
        fprintf(file, "%s", str);
        free(str);
}

/* *pz::code_for_sscan */
int abase_pz_sscan(abase_pz_dst_field k, abase_pz_dst_elt z, const char * str)
{
        mpz_t zz;
        mpz_init(zz);
        if (gmp_sscanf(str, "%Zd", zz) != 1) {
        mpz_clear(zz);
        return 0;
        }
        abase_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
        return 1;
}

/* *pz::code_for_fscan */
int abase_pz_fscan(abase_pz_dst_field k, FILE * file, abase_pz_dst_elt x)
{
        mpz_t zz;
        mpz_init(zz);
        if (gmp_fscanf(file, "%Zd", zz) != 1) {
        mpz_clear(zz);
        return 0;
        }
        abase_pz_set_mpz(k, x, zz);
        mpz_clear(zz);
        return 1;
}


/* Vector functions */
/* *pz::code_for_vec_init */
void abase_pz_vec_init(abase_pz_dst_field k, abase_pz_vec * v, unsigned int n)
{
        *v = (abase_pz_vec) malloc((k->kl) * n * sizeof(mp_limb_t));
        if (!(*v))
        MALLOC_FAILED();
}

/* *pz::code_for_vec_reinit */
void abase_pz_vec_reinit(abase_pz_dst_field k, abase_pz_vec * v, unsigned int n MAYBE_UNUSED, unsigned int m)
{
        *v = (abase_pz_vec) realloc(*v, (k->kl) * m * sizeof(mp_limb_t));
        if (!(*v))
        MALLOC_FAILED();
}

/* *pz::code_for_vec_clear */
void abase_pz_vec_clear(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_vec * v, unsigned int n MAYBE_UNUSED)
{
        free(*v);
}

/* *pz::code_for_vec_set */
void abase_pz_vec_set(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec v, unsigned int n)
{
        if (v != w)
        memmove(w, v, n * k->kl * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_set_partial */
void abase_pz_vec_set_partial(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec v, unsigned int bw, unsigned int bv, unsigned int l)
{
        if ((v + bv * k->kl * sizeof(mp_limb_t)) !=
        (w + bw * k->kl * sizeof(mp_limb_t)))
        memmove((w + bw * k->kl * sizeof(mp_limb_t)),
            (v + bv * k->kl * sizeof(mp_limb_t)),
            l * k->kl * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_set_zero */
void abase_pz_vec_set_zero(abase_pz_dst_field k, abase_pz_dst_vec w, unsigned int n)
{
        memset(w, 0, n * k->kl * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_setcoef */
void abase_pz_vec_setcoef(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_elt x, unsigned int i)
{
        abase_pz_set(k, w + i * k->kl, x);
}

/* *pz::code_for_vec_setcoef_ui */
void abase_pz_vec_setcoef_ui(abase_pz_dst_field k, abase_pz_dst_vec w, unsigned long x0, unsigned int i)
{
        abase_pz_set_ui(k, w + i * k->kl, x0);
}

/* *pz::code_for_vec_getcoef */
void abase_pz_vec_getcoef(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_vec v, unsigned int i)
{
        abase_pz_set(k, z, v + i * k->kl);
}

/* *pz::code_for_vec_add */
void abase_pz_vec_add(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec u, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_add(k, w + i * k->kl, u + i * k->kl, v + i * k->kl);
}

/* *pz::code_for_vec_neg */
void abase_pz_vec_neg(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_neg(k, w + i * k->kl, v + i * k->kl);
}

/* *pz::code_for_vec_rev */
void abase_pz_vec_rev(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec v, unsigned int n)
{
        unsigned int nn = n >> 1;
        abase_pz_elt tmp;
        abase_pz_init(k, &tmp);
        for (unsigned int i = 0; i < nn; ++i) {
        abase_pz_set(k, tmp, v + i * k->kl);
        abase_pz_set(k, w + i * k->kl, v + (n - 1 - i) * k->kl);
        abase_pz_set(k, w + (n - 1 - i) * k->kl, tmp);
        }
        if (n & 1)
        abase_pz_set(k, w + nn * k->kl, v + nn * k->kl);
        abase_pz_clear(k, &tmp);
}

/* *pz::code_for_vec_sub */
void abase_pz_vec_sub(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec u, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_sub(k, w + i * k->kl, u + i * k->kl, v + i * k->kl);
}

/* *pz::code_for_vec_scal_mul */
void abase_pz_vec_scal_mul(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec v, abase_pz_src_elt x, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_mul(k, w + i * k->kl, v + i * k->kl, x);
}

/* *pz::code_for_vec_conv */
void abase_pz_vec_conv(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_src_vec u, unsigned int n, abase_pz_src_vec v, unsigned int m)
{
        abase_pz_vec_ur tmp;
        abase_pz_vec_ur_init(k, &tmp, m + n - 1);
        abase_pz_vec_conv_ur(k, tmp, u, n, v, m);
        abase_pz_vec_reduce(k, w, tmp, m + n - 1);
        abase_pz_vec_ur_clear(k, &tmp, m + n - 1);
}

/* *pz::code_for_vec_random */
void abase_pz_vec_random(abase_pz_dst_field k, abase_pz_dst_vec w, unsigned int n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_random(k, w + i * k->kl, state);
}

/* *pz::code_for_vec_random2 */
void abase_pz_vec_random2(abase_pz_dst_field k, abase_pz_dst_vec w, unsigned int n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_random2(k, w + i * k->kl, state);
}

/* *pz::code_for_vec_cmp */
int abase_pz_vec_cmp(abase_pz_dst_field k, abase_pz_src_vec u, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i) {
        int ret = abase_pz_cmp(k, u + i * k->kl, v + i * k->kl);
        if (ret != 0)
            return ret;
        }
        return 0;
}

/* *pz::code_for_vec_is_zero */
int abase_pz_vec_is_zero(abase_pz_dst_field k, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i) {
        if (!abase_pz_is_zero(k, v + i * k->kl))
            return 0;
        }
        return 1;
}

/* *pz::code_for_vec_asprint */
void abase_pz_vec_asprint(abase_pz_dst_field k, char * * pstr, abase_pz_src_vec v, unsigned int n)
{
        if (n == 0) {
        *pstr = (char *) malloc(4 * sizeof(char));
        sprintf(*pstr, "[ ]");
        return;
        }
        int alloc = 100;
        int len = 0;
        *pstr = (char *) malloc(alloc * sizeof(char));
        char *str = *pstr;
        *str++ = '[';
        *str++ = ' ';
        len = 2;
        for (unsigned int i = 0; i < n; ++i) {
        if (i) {
            (*pstr)[len++] = ',';
            (*pstr)[len++] = ' ';
        }
        char *tmp;
        abase_pz_asprint(k, &tmp, v + i * k->kl);
        int ltmp = strlen(tmp);
        if (len + ltmp + 4 > alloc) {
            alloc = len + ltmp + 100;
            *pstr = (char *) realloc(*pstr, alloc * sizeof(char));
        }
        strncpy(*pstr + len, tmp, ltmp + 4);
        len += ltmp;
        free(tmp);
        }
        (*pstr)[len++] = ' ';
        (*pstr)[len++] = ']';
        (*pstr)[len] = '\0';
}

/* *pz::code_for_vec_fprint */
void abase_pz_vec_fprint(abase_pz_dst_field k, FILE * file, abase_pz_src_vec v, unsigned int n)
{
        char *str;
        abase_pz_vec_asprint(k, &str, v, n);
        fprintf(file, "%s", str);
        free(str);
}

/* *pz::code_for_vec_print */
void abase_pz_vec_print(abase_pz_dst_field k, abase_pz_src_vec v, unsigned int n)
{
        abase_pz_vec_fprint(k, stdout, v, n);
}

/* *pz::code_for_vec_sscan */
int abase_pz_vec_sscan(abase_pz_dst_field k, abase_pz_vec * v, unsigned int * n, const char * str)
{
        // start with a clean vector
        abase_pz_vec_reinit(k, v, *n, 0);
        *n = 0;
        while (isspace((int) (unsigned char) str[0]))
        str++;
        if (str[0] != '[')
        return 0;
        str++;
        if (str[0] != ' ')
        return 0;
        str++;
        if (str[0] == ']') {
        return 1;
        }
        unsigned int i = 0;
        for (;;) {
        if (*n < i + 1) {
            abase_pz_vec_reinit(k, v, *n, i + 1);
            *n = i + 1;
        }
        int ret;
        ret = abase_pz_sscan(k, (*v) + i * k->kl, str);
        if (!ret) {
            return 0;
        }
        i++;
        while (isdigit((int) (unsigned char) str[0]))
            str++;
        while (isspace((int) (unsigned char) str[0]))
            str++;
        if (str[0] == ']')
            break;
        if (str[0] != ',')
            return 0;
        str++;
        while (isspace((int) (unsigned char) str[0]))
            str++;
        }
        return 1;
}

/* *pz::code_for_vec_fscan */
int abase_pz_vec_fscan(abase_pz_dst_field k, FILE * file, abase_pz_vec * v, unsigned int * n)
{
        char *tmp;
        int c;
        int allocated, len = 0;
        allocated = 100;
        tmp = (char *) malloc(allocated * sizeof(char));
        if (!tmp)
        MALLOC_FAILED();
        for (;;) {
        c = fgetc(file);
        if (c == EOF)
            return 0;
        if (len == allocated) {
            allocated += 100;
            tmp = (char *) realloc(tmp, allocated * sizeof(char));
        }
        tmp[len] = c;
        len++;
        if (c == ']')
            break;
        }
        if (len == allocated) {
        allocated += 1;
        tmp = (char *) realloc(tmp, allocated * sizeof(char));
        }
        tmp[len] = '\0';
        int ret = abase_pz_vec_sscan(k, v, n, tmp);
        free(tmp);
        return ret;
}

/* missing vec_scan */
/* *pz::code_for_vec_ur_init */
void abase_pz_vec_ur_init(abase_pz_dst_field k, abase_pz_vec_ur * v, unsigned int n)
{
        *v = (abase_pz_vec) malloc((k->url) * n * sizeof(mp_limb_t));
        if (!(*v))
        MALLOC_FAILED();
}

/* *pz::code_for_vec_ur_set_zero */
void abase_pz_vec_ur_set_zero(abase_pz_dst_field k, abase_pz_dst_vec_ur w, unsigned int n)
{
        memset(w, 0, n * k->url * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_set_vec */
void abase_pz_vec_ur_set_vec(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_elt_ur_set_elt(k, w + i * k->url, v + i * k->kl);
}

/* *pz::code_for_vec_ur_reinit */
void abase_pz_vec_ur_reinit(abase_pz_dst_field k, abase_pz_vec_ur * v, unsigned int n MAYBE_UNUSED, unsigned int m)
{
        *v = (abase_pz_vec) realloc(*v, (k->url) * m * sizeof(mp_limb_t));
        if (!(*v))
        MALLOC_FAILED();
}

/* *pz::code_for_vec_ur_clear */
void abase_pz_vec_ur_clear(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_vec_ur * v, unsigned int n MAYBE_UNUSED)
{
        free(*v);
}

/* *pz::code_for_vec_ur_set */
void abase_pz_vec_ur_set(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec_ur v, unsigned int n)
{
        if (v != w)
        memmove(w, v, n * k->url * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_setcoef */
void abase_pz_vec_ur_setcoef(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_elt_ur x, unsigned int i)
{
        abase_pz_elt_ur_set(k, w + i * k->url, x);
}

/* *pz::code_for_vec_ur_getcoef */
void abase_pz_vec_ur_getcoef(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_vec_ur v, unsigned int i)
{
        abase_pz_elt_ur_set(k, z, v + i * k->url);
}

/* *pz::code_for_vec_ur_add */
void abase_pz_vec_ur_add(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec_ur u, abase_pz_src_vec_ur v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_elt_ur_add(k, w + i * k->url, u + i * k->url,
                      v + i * k->url);
}

/* *pz::code_for_vec_ur_sub */
void abase_pz_vec_ur_sub(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec_ur u, abase_pz_src_vec_ur v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_elt_ur_sub(k, w + i * k->url, u + i * k->url,
                      v + i * k->url);
}

/* *pz::code_for_vec_ur_neg */
void abase_pz_vec_ur_neg(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec_ur v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_elt_ur_neg(k, w + i * k->url, v + i * k->url);
}

/* *pz::code_for_vec_ur_rev */
void abase_pz_vec_ur_rev(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec_ur v, unsigned int n)
{
        unsigned int nn = n >> 1;
        abase_pz_elt_ur tmp;
        abase_pz_elt_ur_init(k, &tmp);
        for (unsigned int i = 0; i < nn; ++i) {
        abase_pz_elt_ur_set(k, tmp, v + i * k->url);
        abase_pz_elt_ur_set(k, w + i * k->url, v + (n - 1 - i) * k->url);
        abase_pz_elt_ur_set(k, w + (n - 1 - i) * k->url, tmp);
        }
        if (n & 1)
        abase_pz_elt_ur_set(k, w + nn * k->url, v + nn * k->url);
        abase_pz_elt_ur_clear(k, &tmp);
}

/* *pz::code_for_vec_scal_mul_ur */
void abase_pz_vec_scal_mul_ur(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec v, abase_pz_src_elt x, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_mul_ur(k, w + i * k->url, v + i * k->kl, x);
}

static void abase_pz_vec_conv_ur_ks(abase_pz_field, abase_pz_dst_vec_ur, abase_pz_src_vec, unsigned int, abase_pz_src_vec, unsigned int);
static /* *pz::code_for_vec_conv_ur */
void abase_pz_vec_conv_ur_ks(abase_pz_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec u, unsigned int n, abase_pz_src_vec v, unsigned int m)
{
        // compute base as a power 2^GMP_NUMB_BITS
        // This is the least number of words that can accomodate
        //     log_2( (p-1)^2 * min(n,m) )
        mpz_t p;
        mpz_init(p);
        abase_pz_field_characteristic(k, p);
        mpz_sub_ui(p, p, 1);
        mpz_mul(p, p, p);
        mpz_mul_ui(p, p, MIN(m, n));
    
        long nbits = mpz_sizeinbase(p, 2);
        unsigned long nwords = 1 + ((nbits - 1) / GMP_NUMB_BITS);
        nbits = GMP_NUMB_BITS * nwords;
        mpz_clear(p);
        assert(k->url >= nwords);
    
        // Create big integers
        mp_limb_t *U, *V;
        U = (mp_limb_t *) malloc(n * nwords * sizeof(mp_limb_t));
        V = (mp_limb_t *) malloc(m * nwords * sizeof(mp_limb_t));
        if (!U || !V)
        MALLOC_FAILED();
        memset(U, 0, n * nwords * sizeof(mp_limb_t));
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_get_mpn(k, U + i * nwords, u + i * k->kl);
        memset(V, 0, m * nwords * sizeof(mp_limb_t));
        for (unsigned int i = 0; i < m; ++i)
        abase_pz_get_mpn(k, V + i * nwords, v + i * k->kl);
        // Mul !
        mp_limb_t *W;
        W = (mp_limb_t *) malloc((n + m) * nwords * sizeof(mp_limb_t));
        if (!W)
        MALLOC_FAILED();
        if (n>=m)
            mpn_mul(W, U, n * nwords, V, m * nwords);
        else
            mpn_mul(W, V, m * nwords, U, n * nwords);
        // Put coeffs in w
        memset(w, 0, (n + m - 1) * k->url * sizeof(mp_limb_t));
        for (unsigned int i = 0; i < n + m - 1; ++i)
        memcpy(w + i * k->url, W + i * nwords, nwords * sizeof(mp_limb_t));
        free(U);
        free(V);
        free(W);
}

/* *pz::code_for_vec_conv_ur */
void abase_pz_vec_conv_ur(abase_pz_dst_field k, abase_pz_dst_vec_ur w, abase_pz_src_vec u, unsigned int n, abase_pz_src_vec v, unsigned int m)
{
        abase_pz_vec_conv_ur_ks(k, w, u, n, v, m);
}

/* *pz::code_for_vec_reduce */
void abase_pz_vec_reduce(abase_pz_dst_field k, abase_pz_dst_vec w, abase_pz_dst_vec_ur v, unsigned int n)
{
        for (unsigned int i = 0; i < n; ++i)
        abase_pz_reduce(k, w + i * k->kl, v + i * k->url);
}

/* *pz::code_for_vec_elt_stride */
ptrdiff_t abase_pz_vec_elt_stride(abase_pz_dst_field k, int n)
{
        return n * k->kl * sizeof(mp_limb_t);
}

/* *pz::code_for_vec_ur_elt_stride */
ptrdiff_t abase_pz_vec_ur_elt_stride(abase_pz_dst_field k, int n)
{
        return n * k->url * sizeof(mp_limb_t);
}


/* Polynomial functions */
/* *Mpfq::defaults::poly::code_for_poly_setmonic, pz */
void abase_pz_poly_setmonic(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly q, abase_pz_src_poly p)
{
    long degp = abase_pz_poly_deg(K, p);
    if (degp == -1) {
        q->size = 0;
        return;
    }
    if (degp == 0) {
        abase_pz_elt aux;
        abase_pz_init(K, &aux);
        abase_pz_set_ui(K, aux, 1);
        abase_pz_poly_setcoef(K, q, aux, 0);
        abase_pz_clear(K, &aux);
        q->size = 1;
        return;
    }
    abase_pz_elt lc;
    abase_pz_init(K, &lc);
    abase_pz_poly_getcoef(K, lc, p, degp);
    abase_pz_inv(K, lc, lc);
    abase_pz_poly_setcoef_ui(K, q, 1, degp);
    abase_pz_vec_scal_mul(K, q->c, p->c, lc, degp);
    q->size = degp+1;
    abase_pz_clear(K, &lc);
}

/* *Mpfq::defaults::poly::code_for_poly_divmod, pz */
void abase_pz_poly_divmod(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly q, abase_pz_dst_poly r, abase_pz_src_poly a, abase_pz_src_poly b)
{
    if (b->size == 0) {
        fprintf(stderr, "Error: division by 0\n");
        exit(1);
    }
    if (a->size == 0) {
        q->size = 0; r->size = 0;
        return;
    }
    int dega = abase_pz_poly_deg(K, a);
    if (dega<0) {
        q->size = 0; r->size = 0;
        return;
    }
    // Compute deg b and inverse of leading coef
    int degb = abase_pz_poly_deg(K, b);
    if (degb<0) {
        fprintf(stderr, "Error: division by 0\n");
        exit(1);
    }
    if (degb > dega) {
        q->size=0;
        abase_pz_poly_set(K, r, a);
        return;
    }
    int bmonic;
    abase_pz_elt ilb;
    abase_pz_init(K, &ilb);
    abase_pz_elt temp;
    abase_pz_init(K, &temp);
    abase_pz_poly_getcoef(K, temp, b, degb);
    if (abase_pz_cmp_ui(K, temp, 1) == 0) {
        abase_pz_set_ui(K, ilb, 1);
        bmonic = 1;
    } else {
        abase_pz_inv(K, ilb, temp);
        bmonic = 0;
    }
    
    abase_pz_poly qq, rr;
    abase_pz_poly_init(K, qq, dega-degb+1);
    abase_pz_poly_init(K, rr, dega);
    
    abase_pz_poly_set(K, rr, a);
    abase_pz_elt aux, aux2;
    
    abase_pz_init(K, &aux);
    abase_pz_init(K, &aux2);
    
    int i;
    int j;
    for (i = dega; i >= (int)degb; --i) {
        abase_pz_poly_getcoef(K, aux, rr, i);
        if (!bmonic) 
            abase_pz_mul(K, aux, aux, ilb);
        abase_pz_poly_setcoef(K, qq, aux, i-degb);
        for (j = i-1; j >= (int)(i - degb); --j) {
            abase_pz_poly_getcoef(K, temp, b, j-i+degb);
            abase_pz_mul(K, aux2, aux, temp);
            abase_pz_poly_getcoef(K, temp, rr, j);
    
            abase_pz_sub(K, temp, temp, aux2);
            abase_pz_poly_setcoef(K, rr, temp, j);
        }
    }    
    
    rr->size = degb;
    int degr = abase_pz_poly_deg(K, rr);
    rr->size = degr+1;
    
    if (q != NULL) 
        abase_pz_poly_set(K, q, qq);
    if (r != NULL)
        abase_pz_poly_set(K, r, rr);
    abase_pz_clear(K, &temp);
    abase_pz_clear(K, &aux);
    abase_pz_clear(K, &aux2);
    abase_pz_poly_clear(K, rr);
    abase_pz_poly_clear(K, qq);
}

static void abase_pz_poly_preinv(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, unsigned int);
static /* *Mpfq::defaults::poly::code_for_poly_precomp_mod, pz */
void abase_pz_poly_preinv(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly q, abase_pz_src_poly p, unsigned int n)
{
    // Compute the inverse of p(x) modulo x^n
    // Newton iteration: x_{n+1} = x_n + x_n(1 - a*x_n)
    // Requires p(0) = 1
    // Assume p != q (no alias)
    abase_pz_elt temp;
    abase_pz_init(K, &temp);
    abase_pz_poly_getcoef(K, temp, p, 0);//Should be in the assert
    assert( abase_pz_cmp_ui(K, temp, 1) == 0);
    assert (p != q);
    int m;
    if (n <= 2) {
        abase_pz_poly_setcoef_ui(K, q, 1, 0);
        q->size = 1;
        m = 1;
        if (n == 1)
            return;
    } else {
        // n >= 3: recursive call at prec m = ceil(n/2)
        m = 1 + ((n-1)/2);
        abase_pz_poly_preinv(K, q, p, m);
    }
    // enlarge q if necessary
    if (q->alloc < n) {
        abase_pz_vec_reinit(K, &(q->c), q->alloc, n);
        q->alloc = n;
    }
    // refine value
    abase_pz_vec tmp;
    abase_pz_vec_init(K, &tmp, m+n-1);
    
    abase_pz_vec_conv(K, tmp, p->c, MIN(n, p->size), q->c, m);
    int nn = MIN(n, MIN(n, p->size) + m -1);
    abase_pz_vec_neg(K, tmp, tmp, nn);
    abase_pz_vec_getcoef(K, temp, tmp, 0);
    abase_pz_add_ui(K, temp, temp, 1);
    abase_pz_vec_setcoef(K, tmp, temp, 0);
    abase_pz_vec_conv(K, tmp, q->c, m, tmp, nn);
    abase_pz_vec_set(K, q->c + m, tmp + m, n-m);
    q->size = n;
    
    abase_pz_clear(K, &temp);
    abase_pz_vec_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::poly::code_for_poly_precomp_mod, pz */
void abase_pz_poly_precomp_mod(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly q, abase_pz_src_poly p)
{
    assert(p != q);
    int N = abase_pz_poly_deg(K, p);
    abase_pz_poly rp;
    abase_pz_poly_init(K, rp, N+1);
    abase_pz_vec_rev(K, rp->c, p->c, N+1);
    rp->size = N+1;
    abase_pz_poly_preinv(K, q, rp, N);
    abase_pz_poly_clear(K, rp);
}

/* *Mpfq::defaults::poly::code_for_poly_mod_pre, pz */
void abase_pz_poly_mod_pre(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly r, abase_pz_src_poly q, abase_pz_src_poly p, abase_pz_src_poly irp)
{
    int N = abase_pz_poly_deg(K, p);
    int degq = abase_pz_poly_deg(K, q);
    if (degq < N) {
        abase_pz_poly_set(K, r, q);
        return;
    }
    int m = degq - N;
    assert (degq <= 2*N-2);
    abase_pz_poly revq;
    abase_pz_poly_init(K, revq, MAX(degq+1, m+1));
    abase_pz_vec_rev(K, revq->c, q->c, degq+1);
    revq->size = q->size;
    abase_pz_poly_mul(K, revq, revq, irp);
    abase_pz_vec_rev(K, revq->c, revq->c, m+1);
    revq->size = m+1;
    
    abase_pz_poly_mul(K, revq, revq, p);
    abase_pz_poly_sub(K, r, q, revq);
    r->size = abase_pz_poly_deg(K, r)+1;
    abase_pz_poly_clear(K, revq);
}

/* missing poly_scan */

/* Functions related to SIMD operation */
/* *simd_pz::code_for_dotprod */
void abase_pz_dotprod(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_vec xw, abase_pz_src_vec xu1, abase_pz_src_vec xu0, unsigned int n)
{
        abase_pz_elt_ur s,t;
        abase_pz_elt_ur_init(K, &s);
        abase_pz_elt_ur_init(K, &t);
        abase_pz_elt_ur_set_zero(K, s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_pz_mul_ur(K, t, abase_pz_vec_coeff_ptr_const(K, xu0, i), abase_pz_vec_coeff_ptr_const(K, xu1, i));
            abase_pz_elt_ur_add(K, s, s, t);
        }
        abase_pz_reduce(K, abase_pz_vec_coeff_ptr(K, xw, 0), s);
        abase_pz_elt_ur_clear(K, &s);
        abase_pz_elt_ur_clear(K, &t);
}


/* Member templates related to SIMD operation */

/* MPI interface */
static void abase_pz_mpi_op_inner(void *, void *, int *, MPI_Datatype *);
static /* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_pz_mpi_op_inner(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    abase_pz_dst_field K;
    MPI_Type_get_attr(*datatype, abase_pz_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    abase_pz_vec_add(K, inoutvec, inoutvec, invec, *len);
}

static void abase_pz_mpi_op_inner_ur(void *, void *, int *, MPI_Datatype *);
static /* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_pz_mpi_op_inner_ur(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    abase_pz_dst_field K;
    MPI_Type_get_attr(*datatype, abase_pz_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    abase_pz_vec_ur_add(K, inoutvec, inoutvec, invec, *len);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_pz_mpi_ops_init(abase_pz_dst_field K MAYBE_UNUSED)
{
        if (abase_pz_impl_mpi_use_count++) return;
    MPI_Type_create_keyval(MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN, &abase_pz_impl_mpi_attr, NULL);
    MPI_Type_contiguous(abase_pz_vec_elt_stride(K, 1), MPI_BYTE, &abase_pz_impl_mpi_datatype);
    MPI_Type_commit(&abase_pz_impl_mpi_datatype);
    MPI_Type_contiguous(abase_pz_vec_ur_elt_stride(K, 1), MPI_BYTE, &abase_pz_impl_mpi_datatype_ur);
    MPI_Type_commit(&abase_pz_impl_mpi_datatype_ur);
    MPI_Type_set_attr(abase_pz_impl_mpi_datatype, abase_pz_impl_mpi_attr, K);
    MPI_Type_set_attr(abase_pz_impl_mpi_datatype_ur, abase_pz_impl_mpi_attr, K);
    /* 1 here indicates that our operation is always taken to be
     * commutative */
    MPI_Op_create(&abase_pz_mpi_op_inner, 1, &abase_pz_impl_mpi_addition_op);
    MPI_Op_create(&abase_pz_mpi_op_inner_ur, 1, &abase_pz_impl_mpi_addition_op_ur);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype */
MPI_Datatype abase_pz_mpi_datatype(abase_pz_dst_field K MAYBE_UNUSED)
{
    return abase_pz_impl_mpi_datatype;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype_ur */
MPI_Datatype abase_pz_mpi_datatype_ur(abase_pz_dst_field K MAYBE_UNUSED)
{
    return abase_pz_impl_mpi_datatype_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op */
MPI_Op abase_pz_mpi_addition_op(abase_pz_dst_field K MAYBE_UNUSED)
{
    return abase_pz_impl_mpi_addition_op;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op_ur */
MPI_Op abase_pz_mpi_addition_op_ur(abase_pz_dst_field K MAYBE_UNUSED)
{
    return abase_pz_impl_mpi_addition_op_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_clear */
void abase_pz_mpi_ops_clear(abase_pz_dst_field K MAYBE_UNUSED)
{
        if (--abase_pz_impl_mpi_use_count) return;
    MPI_Op_free(&abase_pz_impl_mpi_addition_op);
    MPI_Op_free(&abase_pz_impl_mpi_addition_op_ur);
    MPI_Type_delete_attr(abase_pz_impl_mpi_datatype, abase_pz_impl_mpi_attr);
    MPI_Type_delete_attr(abase_pz_impl_mpi_datatype_ur, abase_pz_impl_mpi_attr);
    MPI_Type_free(&abase_pz_impl_mpi_datatype);
    MPI_Type_free(&abase_pz_impl_mpi_datatype_ur);
    MPI_Type_free_keyval(&abase_pz_impl_mpi_attr);
}


/* Object-oriented interface */
static const char * abase_pz_wrapper_impl_name();
static const char * abase_pz_wrapper_impl_name()
{
    return abase_pz_impl_name();
}

static unsigned long abase_pz_wrapper_impl_max_characteristic_bits();
static unsigned long abase_pz_wrapper_impl_max_characteristic_bits()
{
    return abase_pz_impl_max_characteristic_bits();
}

static unsigned long abase_pz_wrapper_impl_max_degree();
static unsigned long abase_pz_wrapper_impl_max_degree()
{
    return abase_pz_impl_max_degree();
}

static void abase_pz_wrapper_field_characteristic(abase_vbase_ptr, mpz_t);
static void abase_pz_wrapper_field_characteristic(abase_vbase_ptr vbase MAYBE_UNUSED, mpz_t z MAYBE_UNUSED)
{
    abase_pz_field_characteristic(vbase->obj, z);
}

static int abase_pz_wrapper_field_degree(abase_vbase_ptr);
static int abase_pz_wrapper_field_degree(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_field_degree(vbase->obj);
}

static void abase_pz_wrapper_field_init(abase_vbase_ptr);
static void abase_pz_wrapper_field_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_field_init(vbase->obj);
}

static void abase_pz_wrapper_field_clear(abase_vbase_ptr);
static void abase_pz_wrapper_field_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_field_clear(vbase->obj);
}

static void abase_pz_wrapper_field_specify(abase_vbase_ptr, unsigned long, void *);
static void abase_pz_wrapper_field_specify(abase_vbase_ptr vbase MAYBE_UNUSED, unsigned long dummy MAYBE_UNUSED, void * vp MAYBE_UNUSED)
{
    abase_pz_field_specify(vbase->obj, dummy, vp);
}

static void abase_pz_wrapper_field_setopt(abase_vbase_ptr, unsigned long, void *);
static void abase_pz_wrapper_field_setopt(abase_vbase_ptr vbase MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, void * y MAYBE_UNUSED)
{
    abase_pz_field_setopt(vbase->obj, x, y);
}

static void abase_pz_wrapper_init(abase_vbase_ptr, abase_pz_elt *);
static void abase_pz_wrapper_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_elt * x MAYBE_UNUSED)
{
    abase_pz_init(vbase->obj, x);
}

static void abase_pz_wrapper_clear(abase_vbase_ptr, abase_pz_elt *);
static void abase_pz_wrapper_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_elt * x MAYBE_UNUSED)
{
    abase_pz_clear(vbase->obj, x);
}

static void abase_pz_wrapper_set(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt);
static void abase_pz_wrapper_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_set(vbase->obj, z, x);
}

static void abase_pz_wrapper_set_ui(abase_vbase_ptr, abase_pz_dst_elt, unsigned long);
static void abase_pz_wrapper_set_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED)
{
    abase_pz_set_ui(vbase->obj, z, x0);
}

static void abase_pz_wrapper_set_zero(abase_vbase_ptr, abase_pz_dst_elt);
static void abase_pz_wrapper_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED)
{
    abase_pz_set_zero(vbase->obj, z);
}

static unsigned long abase_pz_wrapper_get_ui(abase_vbase_ptr, abase_pz_src_elt);
static unsigned long abase_pz_wrapper_get_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    return abase_pz_get_ui(vbase->obj, x);
}

static void abase_pz_wrapper_set_mpn(abase_vbase_ptr, abase_pz_dst_elt, mp_limb_t *, size_t);
static void abase_pz_wrapper_set_mpn(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, mp_limb_t * x MAYBE_UNUSED, size_t n MAYBE_UNUSED)
{
    abase_pz_set_mpn(vbase->obj, z, x, n);
}

static void abase_pz_wrapper_set_mpz(abase_vbase_ptr, abase_pz_dst_elt, mpz_t);
static void abase_pz_wrapper_set_mpz(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, mpz_t x MAYBE_UNUSED)
{
    abase_pz_set_mpz(vbase->obj, z, x);
}

static void abase_pz_wrapper_get_mpn(abase_vbase_ptr, mp_limb_t *, abase_pz_src_elt);
static void abase_pz_wrapper_get_mpn(abase_vbase_ptr vbase MAYBE_UNUSED, mp_limb_t * z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_get_mpn(vbase->obj, z, x);
}

static void abase_pz_wrapper_get_mpz(abase_vbase_ptr, mpz_t, abase_pz_src_elt);
static void abase_pz_wrapper_get_mpz(abase_vbase_ptr vbase MAYBE_UNUSED, mpz_t z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_get_mpz(vbase->obj, z, x);
}

static void abase_pz_wrapper_random(abase_vbase_ptr, abase_pz_dst_elt, gmp_randstate_t);
static void abase_pz_wrapper_random(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_random(vbase->obj, z, state);
}

static void abase_pz_wrapper_random2(abase_vbase_ptr, abase_pz_dst_elt, gmp_randstate_t);
static void abase_pz_wrapper_random2(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_random2(vbase->obj, z, state);
}

static void abase_pz_wrapper_add(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
static void abase_pz_wrapper_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    abase_pz_add(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_sub(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
static void abase_pz_wrapper_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    abase_pz_sub(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_neg(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt);
static void abase_pz_wrapper_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_neg(vbase->obj, z, x);
}

static void abase_pz_wrapper_mul(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
static void abase_pz_wrapper_mul(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    abase_pz_mul(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_sqr(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt);
static void abase_pz_wrapper_sqr(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_sqr(vbase->obj, z, x);
}

static int abase_pz_wrapper_is_sqr(abase_vbase_ptr, abase_pz_src_elt);
static int abase_pz_wrapper_is_sqr(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    return abase_pz_is_sqr(vbase->obj, x);
}

static void abase_pz_wrapper_pow(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, unsigned long *, size_t);
static void abase_pz_wrapper_pow(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED, size_t nl MAYBE_UNUSED)
{
    abase_pz_pow(vbase->obj, z, x, n, nl);
}

static void abase_pz_wrapper_frobenius(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt);
static void abase_pz_wrapper_frobenius(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    abase_pz_frobenius(vbase->obj, x, y);
}

static void abase_pz_wrapper_add_ui(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
static void abase_pz_wrapper_add_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    abase_pz_add_ui(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_sub_ui(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
static void abase_pz_wrapper_sub_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    abase_pz_sub_ui(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_mul_ui(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
static void abase_pz_wrapper_mul_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    abase_pz_mul_ui(vbase->obj, z, x, y);
}

static int abase_pz_wrapper_inv(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_elt);
static int abase_pz_wrapper_inv(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    return abase_pz_inv(vbase->obj, z, x);
}

static void abase_pz_wrapper_elt_ur_init(abase_vbase_ptr, abase_pz_elt_ur *);
static void abase_pz_wrapper_elt_ur_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_elt_ur * x MAYBE_UNUSED)
{
    abase_pz_elt_ur_init(vbase->obj, x);
}

static void abase_pz_wrapper_elt_ur_clear(abase_vbase_ptr, abase_pz_elt_ur *);
static void abase_pz_wrapper_elt_ur_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_elt_ur * x MAYBE_UNUSED)
{
    abase_pz_elt_ur_clear(vbase->obj, x);
}

static void abase_pz_wrapper_elt_ur_set(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt_ur);
static void abase_pz_wrapper_elt_ur_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt_ur x MAYBE_UNUSED)
{
    abase_pz_elt_ur_set(vbase->obj, z, x);
}

static void abase_pz_wrapper_elt_ur_set_elt(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt);
static void abase_pz_wrapper_elt_ur_set_elt(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_elt_ur_set_elt(vbase->obj, z, x);
}

static void abase_pz_wrapper_elt_ur_set_zero(abase_vbase_ptr, abase_pz_dst_elt_ur);
static void abase_pz_wrapper_elt_ur_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED)
{
    abase_pz_elt_ur_set_zero(vbase->obj, z);
}

static void abase_pz_wrapper_elt_ur_set_ui(abase_vbase_ptr, abase_pz_dst_elt_ur, unsigned long);
static void abase_pz_wrapper_elt_ur_set_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED)
{
    abase_pz_elt_ur_set_ui(vbase->obj, z, x0);
}

static void abase_pz_wrapper_elt_ur_add(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt_ur, abase_pz_src_elt_ur);
static void abase_pz_wrapper_elt_ur_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt_ur x MAYBE_UNUSED, abase_pz_src_elt_ur y MAYBE_UNUSED)
{
    abase_pz_elt_ur_add(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_elt_ur_neg(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt_ur);
static void abase_pz_wrapper_elt_ur_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt_ur x MAYBE_UNUSED)
{
    abase_pz_elt_ur_neg(vbase->obj, z, x);
}

static void abase_pz_wrapper_elt_ur_sub(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt_ur, abase_pz_src_elt_ur);
static void abase_pz_wrapper_elt_ur_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt_ur x MAYBE_UNUSED, abase_pz_src_elt_ur y MAYBE_UNUSED)
{
    abase_pz_elt_ur_sub(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_mul_ur(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt, abase_pz_src_elt);
static void abase_pz_wrapper_mul_ur(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    abase_pz_mul_ur(vbase->obj, z, x, y);
}

static void abase_pz_wrapper_sqr_ur(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt);
static void abase_pz_wrapper_sqr_ur(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_sqr_ur(vbase->obj, z, x);
}

static void abase_pz_wrapper_reduce(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_dst_elt_ur);
static void abase_pz_wrapper_reduce(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_dst_elt_ur x MAYBE_UNUSED)
{
    abase_pz_reduce(vbase->obj, z, x);
}

static void abase_pz_wrapper_normalize(abase_vbase_ptr, abase_pz_dst_elt);
static void abase_pz_wrapper_normalize(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED)
{
    abase_pz_normalize(vbase->obj, z);
}

static void abase_pz_wrapper_addmul_si_ur(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_elt, long);
static void abase_pz_wrapper_addmul_si_ur(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur w MAYBE_UNUSED, abase_pz_src_elt u MAYBE_UNUSED, long v MAYBE_UNUSED)
{
    abase_pz_addmul_si_ur(vbase->obj, w, u, v);
}

static int abase_pz_wrapper_cmp(abase_vbase_ptr, abase_pz_src_elt, abase_pz_src_elt);
static int abase_pz_wrapper_cmp(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, abase_pz_src_elt y MAYBE_UNUSED)
{
    return abase_pz_cmp(vbase->obj, x, y);
}

static int abase_pz_wrapper_cmp_ui(abase_vbase_ptr, abase_pz_src_elt, unsigned long);
static int abase_pz_wrapper_cmp_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned long y0 MAYBE_UNUSED)
{
    return abase_pz_cmp_ui(vbase->obj, x, y0);
}

static int abase_pz_wrapper_is_zero(abase_vbase_ptr, abase_pz_src_elt);
static int abase_pz_wrapper_is_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    return abase_pz_is_zero(vbase->obj, x);
}

static void abase_pz_wrapper_asprint(abase_vbase_ptr, char * *, abase_pz_src_elt);
static void abase_pz_wrapper_asprint(abase_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_asprint(vbase->obj, pstr, x);
}

static void abase_pz_wrapper_fprint(abase_vbase_ptr, FILE *, abase_pz_src_elt);
static void abase_pz_wrapper_fprint(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_fprint(vbase->obj, file, x);
}

static void abase_pz_wrapper_print(abase_vbase_ptr, abase_pz_src_elt);
static void abase_pz_wrapper_print(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_print(vbase->obj, x);
}

static int abase_pz_wrapper_sscan(abase_vbase_ptr, abase_pz_dst_elt, const char *);
static int abase_pz_wrapper_sscan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return abase_pz_sscan(vbase->obj, z, str);
}

static int abase_pz_wrapper_fscan(abase_vbase_ptr, FILE *, abase_pz_dst_elt);
static int abase_pz_wrapper_fscan(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_dst_elt x MAYBE_UNUSED)
{
    return abase_pz_fscan(vbase->obj, file, x);
}

static int abase_pz_wrapper_scan(abase_vbase_ptr, abase_pz_dst_elt);
static int abase_pz_wrapper_scan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt x MAYBE_UNUSED)
{
    return abase_pz_scan(vbase->obj, x);
}

static void abase_pz_wrapper_vec_init(abase_vbase_ptr, abase_pz_vec *, unsigned int);
static void abase_pz_wrapper_vec_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_init(vbase->obj, v, n);
}

static void abase_pz_wrapper_vec_reinit(abase_vbase_ptr, abase_pz_vec *, unsigned int, unsigned int);
static void abase_pz_wrapper_vec_reinit(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_pz_vec_reinit(vbase->obj, v, n, m);
}

static void abase_pz_wrapper_vec_clear(abase_vbase_ptr, abase_pz_vec *, unsigned int);
static void abase_pz_wrapper_vec_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_clear(vbase->obj, v, n);
}

static void abase_pz_wrapper_vec_set(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_set(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_set_partial(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, unsigned int, unsigned int, unsigned int);
static void abase_pz_wrapper_vec_set_partial(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int bw MAYBE_UNUSED, unsigned int bv MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
{
    abase_pz_vec_set_partial(vbase->obj, w, v, bw, bv, l);
}

static void abase_pz_wrapper_vec_set_zero(abase_vbase_ptr, abase_pz_dst_vec, unsigned int);
static void abase_pz_wrapper_vec_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_set_zero(vbase->obj, w, n);
}

static void abase_pz_wrapper_vec_setcoef(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_elt, unsigned int);
static void abase_pz_wrapper_vec_setcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_vec_setcoef(vbase->obj, w, x, i);
}

static void abase_pz_wrapper_vec_setcoef_ui(abase_vbase_ptr, abase_pz_dst_vec, unsigned long, unsigned int);
static void abase_pz_wrapper_vec_setcoef_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_vec_setcoef_ui(vbase->obj, w, x0, i);
}

static void abase_pz_wrapper_vec_getcoef(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_getcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt z MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_vec_getcoef(vbase->obj, z, v, i);
}

static void abase_pz_wrapper_vec_add(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec u MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_add(vbase->obj, w, u, v, n);
}

static void abase_pz_wrapper_vec_neg(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_neg(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_rev(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_rev(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_rev(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_sub(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec u MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_sub(vbase->obj, w, u, v, n);
}

static void abase_pz_wrapper_vec_scal_mul(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_elt, unsigned int);
static void abase_pz_wrapper_vec_scal_mul(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_scal_mul(vbase->obj, w, v, x, n);
}

static void abase_pz_wrapper_vec_conv(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, unsigned int, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_conv(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_pz_vec_conv(vbase->obj, w, u, n, v, m);
}

static void abase_pz_wrapper_vec_random(abase_vbase_ptr, abase_pz_dst_vec, unsigned int, gmp_randstate_t);
static void abase_pz_wrapper_vec_random(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_vec_random(vbase->obj, w, n, state);
}

static void abase_pz_wrapper_vec_random2(abase_vbase_ptr, abase_pz_dst_vec, unsigned int, gmp_randstate_t);
static void abase_pz_wrapper_vec_random2(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_vec_random2(vbase->obj, w, n, state);
}

static int abase_pz_wrapper_vec_cmp(abase_vbase_ptr, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
static int abase_pz_wrapper_vec_cmp(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec u MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return abase_pz_vec_cmp(vbase->obj, u, v, n);
}

static int abase_pz_wrapper_vec_is_zero(abase_vbase_ptr, abase_pz_src_vec, unsigned int);
static int abase_pz_wrapper_vec_is_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return abase_pz_vec_is_zero(vbase->obj, v, n);
}

static abase_pz_dst_vec abase_pz_wrapper_vec_subvec(abase_vbase_ptr, abase_pz_dst_vec, int);
static abase_pz_dst_vec abase_pz_wrapper_vec_subvec(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_subvec(vbase->obj, v, i);
}

static abase_pz_src_vec abase_pz_wrapper_vec_subvec_const(abase_vbase_ptr, abase_pz_src_vec, int);
static abase_pz_src_vec abase_pz_wrapper_vec_subvec_const(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_subvec_const(vbase->obj, v, i);
}

static abase_pz_dst_elt abase_pz_wrapper_vec_coeff_ptr(abase_vbase_ptr, abase_pz_dst_vec, int);
static abase_pz_dst_elt abase_pz_wrapper_vec_coeff_ptr(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_coeff_ptr(vbase->obj, v, i);
}

static abase_pz_src_elt abase_pz_wrapper_vec_coeff_ptr_const(abase_vbase_ptr, abase_pz_src_vec, int);
static abase_pz_src_elt abase_pz_wrapper_vec_coeff_ptr_const(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_coeff_ptr_const(vbase->obj, v, i);
}

static void abase_pz_wrapper_vec_asprint(abase_vbase_ptr, char * *, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_asprint(abase_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_asprint(vbase->obj, pstr, v, n);
}

static void abase_pz_wrapper_vec_fprint(abase_vbase_ptr, FILE *, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_fprint(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_fprint(vbase->obj, file, v, n);
}

static void abase_pz_wrapper_vec_print(abase_vbase_ptr, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_print(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_print(vbase->obj, v, n);
}

static int abase_pz_wrapper_vec_sscan(abase_vbase_ptr, abase_pz_vec *, unsigned int *, const char *);
static int abase_pz_wrapper_vec_sscan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec * v MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return abase_pz_vec_sscan(vbase->obj, v, n, str);
}

static int abase_pz_wrapper_vec_fscan(abase_vbase_ptr, FILE *, abase_pz_vec *, unsigned int *);
static int abase_pz_wrapper_vec_fscan(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_vec * v MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED)
{
    return abase_pz_vec_fscan(vbase->obj, file, v, n);
}

static void abase_pz_wrapper_vec_ur_init(abase_vbase_ptr, abase_pz_vec_ur *, unsigned int);
static void abase_pz_wrapper_vec_ur_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_init(vbase->obj, v, n);
}

static void abase_pz_wrapper_vec_ur_set_zero(abase_vbase_ptr, abase_pz_dst_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_set_zero(vbase->obj, w, n);
}

static void abase_pz_wrapper_vec_ur_set_vec(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_ur_set_vec(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_set_vec(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_ur_reinit(abase_vbase_ptr, abase_pz_vec_ur *, unsigned int, unsigned int);
static void abase_pz_wrapper_vec_ur_reinit(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_pz_vec_ur_reinit(vbase->obj, v, n, m);
}

static void abase_pz_wrapper_vec_ur_clear(abase_vbase_ptr, abase_pz_vec_ur *, unsigned int);
static void abase_pz_wrapper_vec_ur_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_clear(vbase->obj, v, n);
}

static void abase_pz_wrapper_vec_ur_set(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_set(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_ur_setcoef(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_elt_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_setcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_elt_ur x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_vec_ur_setcoef(vbase->obj, w, x, i);
}

static void abase_pz_wrapper_vec_ur_getcoef(abase_vbase_ptr, abase_pz_dst_elt_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_getcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt_ur z MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_vec_ur_getcoef(vbase->obj, z, v, i);
}

static void abase_pz_wrapper_vec_ur_add(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec_ur u MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_add(vbase->obj, w, u, v, n);
}

static void abase_pz_wrapper_vec_ur_sub(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec_ur u MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_sub(vbase->obj, w, u, v, n);
}

static void abase_pz_wrapper_vec_ur_neg(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_neg(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_ur_rev(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_ur_rev(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_ur_rev(vbase->obj, w, v, n);
}

static void abase_pz_wrapper_vec_scal_mul_ur(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec, abase_pz_src_elt, unsigned int);
static void abase_pz_wrapper_vec_scal_mul_ur(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_scal_mul_ur(vbase->obj, w, v, x, n);
}

static void abase_pz_wrapper_vec_conv_ur(abase_vbase_ptr, abase_pz_dst_vec_ur, abase_pz_src_vec, unsigned int, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_vec_conv_ur(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur w MAYBE_UNUSED, abase_pz_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, abase_pz_src_vec v MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_pz_vec_conv_ur(vbase->obj, w, u, n, v, m);
}

static void abase_pz_wrapper_vec_reduce(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_dst_vec_ur, unsigned int);
static void abase_pz_wrapper_vec_reduce(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec w MAYBE_UNUSED, abase_pz_dst_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_vec_reduce(vbase->obj, w, v, n);
}

static abase_pz_dst_vec_ur abase_pz_wrapper_vec_ur_subvec(abase_vbase_ptr, abase_pz_dst_vec_ur, int);
static abase_pz_dst_vec_ur abase_pz_wrapper_vec_ur_subvec(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_ur_subvec(vbase->obj, v, i);
}

static abase_pz_src_vec_ur abase_pz_wrapper_vec_ur_subvec_const(abase_vbase_ptr, abase_pz_src_vec_ur, int);
static abase_pz_src_vec_ur abase_pz_wrapper_vec_ur_subvec_const(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_ur_subvec_const(vbase->obj, v, i);
}

static abase_pz_dst_elt abase_pz_wrapper_vec_ur_coeff_ptr(abase_vbase_ptr, abase_pz_dst_vec_ur, int);
static abase_pz_dst_elt abase_pz_wrapper_vec_ur_coeff_ptr(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_ur_coeff_ptr(vbase->obj, v, i);
}

static abase_pz_src_elt abase_pz_wrapper_vec_ur_coeff_ptr_const(abase_vbase_ptr, abase_pz_src_vec_ur, int);
static abase_pz_src_elt abase_pz_wrapper_vec_ur_coeff_ptr_const(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return abase_pz_vec_ur_coeff_ptr_const(vbase->obj, v, i);
}

static ptrdiff_t abase_pz_wrapper_vec_elt_stride(abase_vbase_ptr, int);
static ptrdiff_t abase_pz_wrapper_vec_elt_stride(abase_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return abase_pz_vec_elt_stride(vbase->obj, n);
}

static ptrdiff_t abase_pz_wrapper_vec_ur_elt_stride(abase_vbase_ptr, int);
static ptrdiff_t abase_pz_wrapper_vec_ur_elt_stride(abase_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return abase_pz_vec_ur_elt_stride(vbase->obj, n);
}

static void abase_pz_wrapper_poly_init(abase_vbase_ptr, abase_pz_poly, unsigned int);
static void abase_pz_wrapper_poly_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_poly p MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_poly_init(vbase->obj, p, n);
}

static void abase_pz_wrapper_poly_clear(abase_vbase_ptr, abase_pz_poly);
static void abase_pz_wrapper_poly_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_poly p MAYBE_UNUSED)
{
    abase_pz_poly_clear(vbase->obj, p);
}

static void abase_pz_wrapper_poly_set(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED)
{
    abase_pz_poly_set(vbase->obj, w, u);
}

static void abase_pz_wrapper_poly_setmonic(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_setmonic(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly q MAYBE_UNUSED, abase_pz_src_poly p MAYBE_UNUSED)
{
    abase_pz_poly_setmonic(vbase->obj, q, p);
}

static void abase_pz_wrapper_poly_setcoef(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_elt, unsigned int);
static void abase_pz_wrapper_poly_setcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_poly_setcoef(vbase->obj, w, x, i);
}

static void abase_pz_wrapper_poly_setcoef_ui(abase_vbase_ptr, abase_pz_dst_poly, unsigned long, unsigned int);
static void abase_pz_wrapper_poly_setcoef_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_poly_setcoef_ui(vbase->obj, w, x, i);
}

static void abase_pz_wrapper_poly_getcoef(abase_vbase_ptr, abase_pz_dst_elt, abase_pz_src_poly, unsigned int);
static void abase_pz_wrapper_poly_getcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt x MAYBE_UNUSED, abase_pz_src_poly w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_pz_poly_getcoef(vbase->obj, x, w, i);
}

static int abase_pz_wrapper_poly_deg(abase_vbase_ptr, abase_pz_src_poly);
static int abase_pz_wrapper_poly_deg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_poly w MAYBE_UNUSED)
{
    return abase_pz_poly_deg(vbase->obj, w);
}

static void abase_pz_wrapper_poly_add(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, abase_pz_src_poly v MAYBE_UNUSED)
{
    abase_pz_poly_add(vbase->obj, w, u, v);
}

static void abase_pz_wrapper_poly_sub(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, abase_pz_src_poly v MAYBE_UNUSED)
{
    abase_pz_poly_sub(vbase->obj, w, u, v);
}

static void abase_pz_wrapper_poly_add_ui(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, unsigned long);
static void abase_pz_wrapper_poly_add_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    abase_pz_poly_add_ui(vbase->obj, w, u, x);
}

static void abase_pz_wrapper_poly_sub_ui(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, unsigned long);
static void abase_pz_wrapper_poly_sub_ui(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    abase_pz_poly_sub_ui(vbase->obj, w, u, x);
}

static void abase_pz_wrapper_poly_neg(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED)
{
    abase_pz_poly_neg(vbase->obj, w, u);
}

static void abase_pz_wrapper_poly_scal_mul(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_elt);
static void abase_pz_wrapper_poly_scal_mul(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, abase_pz_src_elt x MAYBE_UNUSED)
{
    abase_pz_poly_scal_mul(vbase->obj, w, u, x);
}

static void abase_pz_wrapper_poly_mul(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_mul(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, abase_pz_src_poly v MAYBE_UNUSED)
{
    abase_pz_poly_mul(vbase->obj, w, u, v);
}

static void abase_pz_wrapper_poly_divmod(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_divmod(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly q MAYBE_UNUSED, abase_pz_dst_poly r MAYBE_UNUSED, abase_pz_src_poly a MAYBE_UNUSED, abase_pz_src_poly b MAYBE_UNUSED)
{
    abase_pz_poly_divmod(vbase->obj, q, r, a, b);
}

static void abase_pz_wrapper_poly_precomp_mod(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_precomp_mod(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly q MAYBE_UNUSED, abase_pz_src_poly p MAYBE_UNUSED)
{
    abase_pz_poly_precomp_mod(vbase->obj, q, p);
}

static void abase_pz_wrapper_poly_mod_pre(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_mod_pre(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly r MAYBE_UNUSED, abase_pz_src_poly q MAYBE_UNUSED, abase_pz_src_poly p MAYBE_UNUSED, abase_pz_src_poly irp MAYBE_UNUSED)
{
    abase_pz_poly_mod_pre(vbase->obj, r, q, p, irp);
}

static void abase_pz_wrapper_poly_gcd(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_gcd(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly g MAYBE_UNUSED, abase_pz_src_poly a0 MAYBE_UNUSED, abase_pz_src_poly b0 MAYBE_UNUSED)
{
    abase_pz_poly_gcd(vbase->obj, g, a0, b0);
}

static void abase_pz_wrapper_poly_xgcd(abase_vbase_ptr, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static void abase_pz_wrapper_poly_xgcd(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly g MAYBE_UNUSED, abase_pz_dst_poly u0 MAYBE_UNUSED, abase_pz_dst_poly v0 MAYBE_UNUSED, abase_pz_src_poly a0 MAYBE_UNUSED, abase_pz_src_poly b0 MAYBE_UNUSED)
{
    abase_pz_poly_xgcd(vbase->obj, g, u0, v0, a0, b0);
}

static void abase_pz_wrapper_poly_random(abase_vbase_ptr, abase_pz_dst_poly, unsigned int, gmp_randstate_t);
static void abase_pz_wrapper_poly_random(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_poly_random(vbase->obj, w, n, state);
}

static void abase_pz_wrapper_poly_random2(abase_vbase_ptr, abase_pz_dst_poly, unsigned int, gmp_randstate_t);
static void abase_pz_wrapper_poly_random2(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    abase_pz_poly_random2(vbase->obj, w, n, state);
}

static int abase_pz_wrapper_poly_cmp(abase_vbase_ptr, abase_pz_src_poly, abase_pz_src_poly);
static int abase_pz_wrapper_poly_cmp(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_poly u MAYBE_UNUSED, abase_pz_src_poly v MAYBE_UNUSED)
{
    return abase_pz_poly_cmp(vbase->obj, u, v);
}

static void abase_pz_wrapper_poly_asprint(abase_vbase_ptr, char * *, abase_pz_src_poly);
static void abase_pz_wrapper_poly_asprint(abase_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, abase_pz_src_poly w MAYBE_UNUSED)
{
    abase_pz_poly_asprint(vbase->obj, pstr, w);
}

static void abase_pz_wrapper_poly_fprint(abase_vbase_ptr, FILE *, abase_pz_src_poly);
static void abase_pz_wrapper_poly_fprint(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_src_poly w MAYBE_UNUSED)
{
    abase_pz_poly_fprint(vbase->obj, file, w);
}

static void abase_pz_wrapper_poly_print(abase_vbase_ptr, abase_pz_src_poly);
static void abase_pz_wrapper_poly_print(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_src_poly w MAYBE_UNUSED)
{
    abase_pz_poly_print(vbase->obj, w);
}

static int abase_pz_wrapper_poly_sscan(abase_vbase_ptr, abase_pz_dst_poly, const char *);
static int abase_pz_wrapper_poly_sscan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return abase_pz_poly_sscan(vbase->obj, w, str);
}

static int abase_pz_wrapper_poly_fscan(abase_vbase_ptr, FILE *, abase_pz_dst_poly);
static int abase_pz_wrapper_poly_fscan(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_pz_dst_poly w MAYBE_UNUSED)
{
    return abase_pz_poly_fscan(vbase->obj, file, w);
}

static int abase_pz_wrapper_groupsize(abase_vbase_ptr);
static int abase_pz_wrapper_groupsize(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_groupsize(vbase->obj);
}

static int abase_pz_wrapper_offset(abase_vbase_ptr, int);
static int abase_pz_wrapper_offset(abase_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return abase_pz_offset(vbase->obj, n);
}

static int abase_pz_wrapper_stride(abase_vbase_ptr);
static int abase_pz_wrapper_stride(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_stride(vbase->obj);
}

static void abase_pz_wrapper_set_ui_at(abase_vbase_ptr, abase_pz_dst_elt, int, unsigned long);
static void abase_pz_wrapper_set_ui_at(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_pz_set_ui_at(vbase->obj, p, k, v);
}

static void abase_pz_wrapper_set_ui_all(abase_vbase_ptr, abase_pz_dst_elt, unsigned long);
static void abase_pz_wrapper_set_ui_all(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt p MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_pz_set_ui_all(vbase->obj, p, v);
}

static void abase_pz_wrapper_elt_ur_set_ui_at(abase_vbase_ptr, abase_pz_dst_elt, int, unsigned long);
static void abase_pz_wrapper_elt_ur_set_ui_at(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_pz_elt_ur_set_ui_at(vbase->obj, p, k, v);
}

static void abase_pz_wrapper_elt_ur_set_ui_all(abase_vbase_ptr, abase_pz_dst_elt, unsigned long);
static void abase_pz_wrapper_elt_ur_set_ui_all(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_elt p MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_pz_elt_ur_set_ui_all(vbase->obj, p, v);
}

static void abase_pz_wrapper_dotprod(abase_vbase_ptr, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
static void abase_pz_wrapper_dotprod(abase_vbase_ptr vbase MAYBE_UNUSED, abase_pz_dst_vec xw MAYBE_UNUSED, abase_pz_src_vec xu1 MAYBE_UNUSED, abase_pz_src_vec xu0 MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_pz_dotprod(vbase->obj, xw, xu1, xu0, n);
}

static void abase_pz_wrapper_mpi_ops_init(abase_vbase_ptr);
static void abase_pz_wrapper_mpi_ops_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_mpi_ops_init(vbase->obj);
}

static MPI_Datatype abase_pz_wrapper_mpi_datatype(abase_vbase_ptr);
static MPI_Datatype abase_pz_wrapper_mpi_datatype(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_mpi_datatype(vbase->obj);
}

static MPI_Datatype abase_pz_wrapper_mpi_datatype_ur(abase_vbase_ptr);
static MPI_Datatype abase_pz_wrapper_mpi_datatype_ur(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_mpi_datatype_ur(vbase->obj);
}

static MPI_Op abase_pz_wrapper_mpi_addition_op(abase_vbase_ptr);
static MPI_Op abase_pz_wrapper_mpi_addition_op(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_mpi_addition_op(vbase->obj);
}

static MPI_Op abase_pz_wrapper_mpi_addition_op_ur(abase_vbase_ptr);
static MPI_Op abase_pz_wrapper_mpi_addition_op_ur(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_pz_mpi_addition_op_ur(vbase->obj);
}

static void abase_pz_wrapper_mpi_ops_clear(abase_vbase_ptr);
static void abase_pz_wrapper_mpi_ops_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_mpi_ops_clear(vbase->obj);
}

static void abase_pz_wrapper_oo_field_init(abase_vbase_ptr);
static void abase_pz_wrapper_oo_field_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_oo_field_init(vbase);
}

static void abase_pz_wrapper_oo_field_clear(abase_vbase_ptr);
static void abase_pz_wrapper_oo_field_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_pz_oo_field_clear(vbase);
}

void abase_pz_oo_field_init(abase_vbase_ptr vbase)
{
    memset(vbase, 0, sizeof(struct abase_vbase_s));
    vbase->obj = malloc(sizeof(abase_pz_field));
    abase_pz_field_init((abase_pz_dst_field) vbase->obj);
    vbase->impl_name = (const char * (*) ()) abase_pz_wrapper_impl_name;
    vbase->impl_max_characteristic_bits = (unsigned long (*) ()) abase_pz_wrapper_impl_max_characteristic_bits;
    vbase->impl_max_degree = (unsigned long (*) ()) abase_pz_wrapper_impl_max_degree;
    vbase->field_characteristic = (void (*) (abase_vbase_ptr, mpz_t)) abase_pz_wrapper_field_characteristic;
    vbase->field_degree = (int (*) (abase_vbase_ptr)) abase_pz_wrapper_field_degree;
    vbase->field_init = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_field_init;
    vbase->field_clear = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_field_clear;
    vbase->field_specify = (void (*) (abase_vbase_ptr, unsigned long, void *)) abase_pz_wrapper_field_specify;
    vbase->field_setopt = (void (*) (abase_vbase_ptr, unsigned long, void *)) abase_pz_wrapper_field_setopt;
    vbase->init = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_init;
    vbase->clear = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_clear;
    vbase->set = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_set;
    vbase->set_ui = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_pz_wrapper_set_ui;
    vbase->set_zero = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_set_zero;
    vbase->get_ui = (unsigned long (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_get_ui;
    vbase->set_mpn = (void (*) (abase_vbase_ptr, void *, mp_limb_t *, size_t)) abase_pz_wrapper_set_mpn;
    vbase->set_mpz = (void (*) (abase_vbase_ptr, void *, mpz_t)) abase_pz_wrapper_set_mpz;
    vbase->get_mpn = (void (*) (abase_vbase_ptr, mp_limb_t *, const void *)) abase_pz_wrapper_get_mpn;
    vbase->get_mpz = (void (*) (abase_vbase_ptr, mpz_t, const void *)) abase_pz_wrapper_get_mpz;
    vbase->random = (void (*) (abase_vbase_ptr, void *, gmp_randstate_t)) abase_pz_wrapper_random;
    vbase->random2 = (void (*) (abase_vbase_ptr, void *, gmp_randstate_t)) abase_pz_wrapper_random2;
    vbase->add = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_add;
    vbase->sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_sub;
    vbase->neg = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_neg;
    vbase->mul = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_mul;
    vbase->sqr = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_sqr;
    vbase->is_sqr = (int (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_is_sqr;
    /* missing sqrt */
    vbase->pow = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long *, size_t)) abase_pz_wrapper_pow;
    vbase->frobenius = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_frobenius;
    vbase->add_ui = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long)) abase_pz_wrapper_add_ui;
    vbase->sub_ui = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long)) abase_pz_wrapper_sub_ui;
    vbase->mul_ui = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long)) abase_pz_wrapper_mul_ui;
    vbase->inv = (int (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_inv;
    vbase->elt_ur_init = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_elt_ur_init;
    vbase->elt_ur_clear = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_elt_ur_clear;
    vbase->elt_ur_set = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_elt_ur_set;
    vbase->elt_ur_set_elt = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_elt_ur_set_elt;
    vbase->elt_ur_set_zero = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_elt_ur_set_zero;
    vbase->elt_ur_set_ui = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_pz_wrapper_elt_ur_set_ui;
    vbase->elt_ur_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_elt_ur_add;
    vbase->elt_ur_neg = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_elt_ur_neg;
    vbase->elt_ur_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_elt_ur_sub;
    vbase->mul_ur = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_mul_ur;
    vbase->sqr_ur = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_sqr_ur;
    vbase->reduce = (void (*) (abase_vbase_ptr, void *, void *)) abase_pz_wrapper_reduce;
    vbase->normalize = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_normalize;
    vbase->addmul_si_ur = (void (*) (abase_vbase_ptr, void *, const void *, long)) abase_pz_wrapper_addmul_si_ur;
    vbase->cmp = (int (*) (abase_vbase_ptr, const void *, const void *)) abase_pz_wrapper_cmp;
    vbase->cmp_ui = (int (*) (abase_vbase_ptr, const void *, unsigned long)) abase_pz_wrapper_cmp_ui;
    vbase->is_zero = (int (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_is_zero;
    vbase->asprint = (void (*) (abase_vbase_ptr, char * *, const void *)) abase_pz_wrapper_asprint;
    vbase->fprint = (void (*) (abase_vbase_ptr, FILE *, const void *)) abase_pz_wrapper_fprint;
    vbase->print = (void (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_print;
    vbase->sscan = (int (*) (abase_vbase_ptr, void *, const char *)) abase_pz_wrapper_sscan;
    vbase->fscan = (int (*) (abase_vbase_ptr, FILE *, void *)) abase_pz_wrapper_fscan;
    vbase->scan = (int (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_scan;
    vbase->vec_init = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_init;
    vbase->vec_reinit = (void (*) (abase_vbase_ptr, void *, unsigned int, unsigned int)) abase_pz_wrapper_vec_reinit;
    vbase->vec_clear = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_clear;
    vbase->vec_set = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_set;
    vbase->vec_set_partial = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int, unsigned int, unsigned int)) abase_pz_wrapper_vec_set_partial;
    vbase->vec_set_zero = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_set_zero;
    vbase->vec_setcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_setcoef;
    vbase->vec_setcoef_ui = (void (*) (abase_vbase_ptr, void *, unsigned long, unsigned int)) abase_pz_wrapper_vec_setcoef_ui;
    vbase->vec_getcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_getcoef;
    vbase->vec_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_add;
    vbase->vec_neg = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_neg;
    vbase->vec_rev = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_rev;
    vbase->vec_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_sub;
    vbase->vec_scal_mul = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_scal_mul;
    vbase->vec_conv = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int)) abase_pz_wrapper_vec_conv;
    vbase->vec_random = (void (*) (abase_vbase_ptr, void *, unsigned int, gmp_randstate_t)) abase_pz_wrapper_vec_random;
    vbase->vec_random2 = (void (*) (abase_vbase_ptr, void *, unsigned int, gmp_randstate_t)) abase_pz_wrapper_vec_random2;
    vbase->vec_cmp = (int (*) (abase_vbase_ptr, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_cmp;
    vbase->vec_is_zero = (int (*) (abase_vbase_ptr, const void *, unsigned int)) abase_pz_wrapper_vec_is_zero;
    vbase->vec_subvec = (void * (*) (abase_vbase_ptr, void *, int)) abase_pz_wrapper_vec_subvec;
    vbase->vec_subvec_const = (const void * (*) (abase_vbase_ptr, const void *, int)) abase_pz_wrapper_vec_subvec_const;
    vbase->vec_coeff_ptr = (void * (*) (abase_vbase_ptr, void *, int)) abase_pz_wrapper_vec_coeff_ptr;
    vbase->vec_coeff_ptr_const = (const void * (*) (abase_vbase_ptr, const void *, int)) abase_pz_wrapper_vec_coeff_ptr_const;
    vbase->vec_asprint = (void (*) (abase_vbase_ptr, char * *, const void *, unsigned int)) abase_pz_wrapper_vec_asprint;
    vbase->vec_fprint = (void (*) (abase_vbase_ptr, FILE *, const void *, unsigned int)) abase_pz_wrapper_vec_fprint;
    vbase->vec_print = (void (*) (abase_vbase_ptr, const void *, unsigned int)) abase_pz_wrapper_vec_print;
    vbase->vec_sscan = (int (*) (abase_vbase_ptr, void *, unsigned int *, const char *)) abase_pz_wrapper_vec_sscan;
    vbase->vec_fscan = (int (*) (abase_vbase_ptr, FILE *, void *, unsigned int *)) abase_pz_wrapper_vec_fscan;
    /* missing vec_scan */
    vbase->vec_ur_init = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_ur_init;
    vbase->vec_ur_set_zero = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_ur_set_zero;
    vbase->vec_ur_set_vec = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_set_vec;
    vbase->vec_ur_reinit = (void (*) (abase_vbase_ptr, void *, unsigned int, unsigned int)) abase_pz_wrapper_vec_ur_reinit;
    vbase->vec_ur_clear = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_vec_ur_clear;
    vbase->vec_ur_set = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_set;
    vbase->vec_ur_setcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_setcoef;
    vbase->vec_ur_getcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_getcoef;
    vbase->vec_ur_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_add;
    vbase->vec_ur_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_sub;
    vbase->vec_ur_neg = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_neg;
    vbase->vec_ur_rev = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_vec_ur_rev;
    vbase->vec_scal_mul_ur = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_vec_scal_mul_ur;
    vbase->vec_conv_ur = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int)) abase_pz_wrapper_vec_conv_ur;
    vbase->vec_reduce = (void (*) (abase_vbase_ptr, void *, void *, unsigned int)) abase_pz_wrapper_vec_reduce;
    vbase->vec_ur_subvec = (void * (*) (abase_vbase_ptr, void *, int)) abase_pz_wrapper_vec_ur_subvec;
    vbase->vec_ur_subvec_const = (const void * (*) (abase_vbase_ptr, const void *, int)) abase_pz_wrapper_vec_ur_subvec_const;
    vbase->vec_ur_coeff_ptr = (void * (*) (abase_vbase_ptr, void *, int)) abase_pz_wrapper_vec_ur_coeff_ptr;
    vbase->vec_ur_coeff_ptr_const = (const void * (*) (abase_vbase_ptr, const void *, int)) abase_pz_wrapper_vec_ur_coeff_ptr_const;
    vbase->vec_elt_stride = (ptrdiff_t (*) (abase_vbase_ptr, int)) abase_pz_wrapper_vec_elt_stride;
    vbase->vec_ur_elt_stride = (ptrdiff_t (*) (abase_vbase_ptr, int)) abase_pz_wrapper_vec_ur_elt_stride;
    vbase->poly_init = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_pz_wrapper_poly_init;
    vbase->poly_clear = (void (*) (abase_vbase_ptr, void *)) abase_pz_wrapper_poly_clear;
    vbase->poly_set = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_poly_set;
    vbase->poly_setmonic = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_poly_setmonic;
    vbase->poly_setcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_poly_setcoef;
    vbase->poly_setcoef_ui = (void (*) (abase_vbase_ptr, void *, unsigned long, unsigned int)) abase_pz_wrapper_poly_setcoef_ui;
    vbase->poly_getcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_pz_wrapper_poly_getcoef;
    vbase->poly_deg = (int (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_poly_deg;
    vbase->poly_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_poly_add;
    vbase->poly_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_poly_sub;
    vbase->poly_add_ui = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long)) abase_pz_wrapper_poly_add_ui;
    vbase->poly_sub_ui = (void (*) (abase_vbase_ptr, void *, const void *, unsigned long)) abase_pz_wrapper_poly_sub_ui;
    vbase->poly_neg = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_poly_neg;
    vbase->poly_scal_mul = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_poly_scal_mul;
    vbase->poly_mul = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_poly_mul;
    vbase->poly_divmod = (void (*) (abase_vbase_ptr, void *, void *, const void *, const void *)) abase_pz_wrapper_poly_divmod;
    vbase->poly_precomp_mod = (void (*) (abase_vbase_ptr, void *, const void *)) abase_pz_wrapper_poly_precomp_mod;
    vbase->poly_mod_pre = (void (*) (abase_vbase_ptr, void *, const void *, const void *, const void *)) abase_pz_wrapper_poly_mod_pre;
    vbase->poly_gcd = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_pz_wrapper_poly_gcd;
    vbase->poly_xgcd = (void (*) (abase_vbase_ptr, void *, void *, void *, const void *, const void *)) abase_pz_wrapper_poly_xgcd;
    vbase->poly_random = (void (*) (abase_vbase_ptr, void *, unsigned int, gmp_randstate_t)) abase_pz_wrapper_poly_random;
    vbase->poly_random2 = (void (*) (abase_vbase_ptr, void *, unsigned int, gmp_randstate_t)) abase_pz_wrapper_poly_random2;
    vbase->poly_cmp = (int (*) (abase_vbase_ptr, const void *, const void *)) abase_pz_wrapper_poly_cmp;
    vbase->poly_asprint = (void (*) (abase_vbase_ptr, char * *, const void *)) abase_pz_wrapper_poly_asprint;
    vbase->poly_fprint = (void (*) (abase_vbase_ptr, FILE *, const void *)) abase_pz_wrapper_poly_fprint;
    vbase->poly_print = (void (*) (abase_vbase_ptr, const void *)) abase_pz_wrapper_poly_print;
    vbase->poly_sscan = (int (*) (abase_vbase_ptr, void *, const char *)) abase_pz_wrapper_poly_sscan;
    vbase->poly_fscan = (int (*) (abase_vbase_ptr, FILE *, void *)) abase_pz_wrapper_poly_fscan;
    /* missing poly_scan */
    vbase->groupsize = (int (*) (abase_vbase_ptr)) abase_pz_wrapper_groupsize;
    vbase->offset = (int (*) (abase_vbase_ptr, int)) abase_pz_wrapper_offset;
    vbase->stride = (int (*) (abase_vbase_ptr)) abase_pz_wrapper_stride;
    vbase->set_ui_at = (void (*) (abase_vbase_ptr, void *, int, unsigned long)) abase_pz_wrapper_set_ui_at;
    vbase->set_ui_all = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_pz_wrapper_set_ui_all;
    vbase->elt_ur_set_ui_at = (void (*) (abase_vbase_ptr, void *, int, unsigned long)) abase_pz_wrapper_elt_ur_set_ui_at;
    vbase->elt_ur_set_ui_all = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_pz_wrapper_elt_ur_set_ui_all;
    vbase->dotprod = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_pz_wrapper_dotprod;
    vbase->mpi_ops_init = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_ops_init;
    vbase->mpi_datatype = (MPI_Datatype (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_datatype;
    vbase->mpi_datatype_ur = (MPI_Datatype (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_datatype_ur;
    vbase->mpi_addition_op = (MPI_Op (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_addition_op;
    vbase->mpi_addition_op_ur = (MPI_Op (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_addition_op_ur;
    vbase->mpi_ops_clear = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_mpi_ops_clear;
    vbase->oo_field_init = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_oo_field_init;
    vbase->oo_field_clear = (void (*) (abase_vbase_ptr)) abase_pz_wrapper_oo_field_clear;
}


/* vim:set ft=cpp: */
