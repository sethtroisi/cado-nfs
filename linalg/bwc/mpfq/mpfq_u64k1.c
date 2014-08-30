/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "mpfq_u64k1.h"

#include "binary-dotprods-backends.h"
#include <inttypes.h>
static int mpfq_u64k1_impl_mpi_attr;     /* for MPI functions */
static MPI_Datatype mpfq_u64k1_impl_mpi_datatype;
static MPI_Datatype mpfq_u64k1_impl_mpi_datatype_ur;
static MPI_Op mpfq_u64k1_impl_mpi_addition_op;
static MPI_Op mpfq_u64k1_impl_mpi_addition_op_ur;
static int mpfq_u64k1_impl_mpi_use_count;   /* several stacked init()/clear() pairs are supported */
/* Active handler: simd_u64k */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_char2 */
/* Options used:{
   family=[ u64k1, u64k2, ],
   k=1,
   tag=u64k1,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     ],
    member_templates_restrict={
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     u64k1=[ u64k1, u64k2, ],
     u64k2=[ u64k1, u64k2, ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=mpfq_vbase,
    global_prefix=mpfq_,
    name=mpfq_vbase,
    substitutions=[
     [ (?^:mpfq_u64k1_elt \*), void *, ],
     [ (?^:mpfq_u64k1_src_elt\b), const void *, ],
     [ (?^:mpfq_u64k1_elt\b), void *, ],
     [ (?^:mpfq_u64k1_dst_elt\b), void *, ],
     [ (?^:mpfq_u64k1_elt_ur \*), void *, ],
     [ (?^:mpfq_u64k1_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_u64k1_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k1_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k1_vec \*), void *, ],
     [ (?^:mpfq_u64k1_src_vec\b), const void *, ],
     [ (?^:mpfq_u64k1_vec\b), void *, ],
     [ (?^:mpfq_u64k1_dst_vec\b), void *, ],
     [ (?^:mpfq_u64k1_vec_ur \*), void *, ],
     [ (?^:mpfq_u64k1_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_u64k1_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k1_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k1_poly \*), void *, ],
     [ (?^:mpfq_u64k1_src_poly\b), const void *, ],
     [ (?^:mpfq_u64k1_poly\b), void *, ],
     [ (?^:mpfq_u64k1_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */


/* Functions operating on the field structure */
/* *simd_u64k::code_for_field_specify */
void mpfq_u64k1_field_specify(mpfq_u64k1_dst_field K MAYBE_UNUSED, unsigned long tag, void * x MAYBE_UNUSED)
{
        if (tag == MPFQ_GROUPSIZE) {
            assert(*(int*)x == 64);
        } else if (tag == MPFQ_PRIME_MPZ) {
            assert(mpz_cmp_ui((mpz_srcptr)x, 2) == 0);
        } else {
            fprintf(stderr, "Unsupported field_specify tag %ld\n", tag);
        }
}


/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */
/* *io::code_for_asprint */
int mpfq_u64k1_asprint(mpfq_u64k1_dst_field K MAYBE_UNUSED, char * * ps, mpfq_u64k1_src_elt x)
{
        /* Hmm, this has never been tested, right ? Attempting a quick fix... */
        const uint64_t * y = x;
        const unsigned int stride = mpfq_u64k1_vec_elt_stride(K,1)/sizeof(uint64_t);
        *ps = mpfq_malloc_check(stride * 16 + 1);
        memset(*ps, ' ', stride * 16);
        int n;
        for(unsigned int i = 0 ; i < stride ; i++) {
            n = snprintf((*ps) + i * 16, 17, "%" PRIx64, y[i]);
            (*ps)[i*16 + n]=',';
        }
        (*ps)[(stride-1) * 16 + n]='\0';
        return (stride-1) * 16 + n;
}

/* *io::code_for_fprint */
int mpfq_u64k1_fprint(mpfq_u64k1_dst_field k, FILE * file, mpfq_u64k1_src_elt x)
{
    char *str;
    int rc;
    mpfq_u64k1_asprint(k,&str,x);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *io::code_for_sscan */
int mpfq_u64k1_sscan(mpfq_u64k1_dst_field k MAYBE_UNUSED, mpfq_u64k1_dst_elt z, const char * str)
{
        char tmp[17];
        uint64_t * y = z;
        const unsigned int stride = mpfq_u64k1_vec_elt_stride(K,1)/sizeof(uint64_t);
        assert(strlen(str) >= 1 * 16);
        int r = 0;
        for(unsigned int i = 0 ; i < stride ; i++) {
            memcpy(tmp, str + i * 16, 16);
            tmp[16]=0;
            if (sscanf(tmp, "%" SCNx64, &(y[i])) == 1) {
                r+=16;
            } else {
                return r;
            }
        }
        return r;
}

/* *io::code_for_fscan */
int mpfq_u64k1_fscan(mpfq_u64k1_dst_field k, FILE * file, mpfq_u64k1_dst_elt z)
{
    char *tmp;
    int allocated, len=0;
    int c, start=0;
    allocated=100;
    tmp = (char *)mpfq_malloc_check(allocated*sizeof(char));
    for(;;) {
        c = fgetc(file);
        if (c==EOF)
            break;
        if (isspace((int)(unsigned char)c)) {
            if (start==0)
                continue;
            else
                break;
        } else {
            if (len==allocated) {
                allocated+=100;
                tmp = (char*)realloc(tmp, allocated*sizeof(char));
            }
            tmp[len]=c;
            len++;
            start=1;
        }
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated*sizeof(char));
    }
    tmp[len]='\0';
    int ret=mpfq_u64k1_sscan(k,z,tmp);
    free(tmp);
    return ret ? len : 0;
}


/* Vector functions */
/* *Mpfq::defaults::vec::alloc::code_for_vec_init, Mpfq::defaults::vec */
void mpfq_u64k1_vec_init(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec * v, unsigned int n)
{
    unsigned int i;
    *v = (mpfq_u64k1_vec) malloc (n*sizeof(mpfq_u64k1_elt));
    for(i = 0; i < n; i++)
        mpfq_u64k1_init(K, (*v) + i);
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_reinit, Mpfq::defaults::vec */
void mpfq_u64k1_vec_reinit(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        unsigned int i;
        *v = (mpfq_u64k1_vec) realloc (*v, m * sizeof(mpfq_u64k1_elt));
        for(i = n; i < m; i+=1)
            mpfq_u64k1_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_u64k1_clear(K, (*v) + i);
        *v = (mpfq_u64k1_vec) realloc (*v, m * sizeof(mpfq_u64k1_elt));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_clear, Mpfq::defaults::vec */
void mpfq_u64k1_vec_clear(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec * v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_u64k1_clear(K, (*v) + i);
    free(*v);
}

/* missing vec_setcoeff_ui */
/* missing vec_scal_mul */
/* missing vec_random2 */
/* *Mpfq::defaults::vec::io::code_for_vec_asprint, Mpfq::defaults::vec */
int mpfq_u64k1_vec_asprint(mpfq_u64k1_dst_field K MAYBE_UNUSED, char * * pstr, mpfq_u64k1_src_vec w, unsigned int n)
{
    if (n == 0) {
        *pstr = (char *)mpfq_malloc_check(4*sizeof(char));
        sprintf(*pstr, "[ ]");
        return strlen(*pstr);
    }
    int alloc = 100;
    int len = 0;
    *pstr = (char *)mpfq_malloc_check(alloc*sizeof(char));
    char *str = *pstr;
    *str++ = '[';
    *str++ = ' ';
    len = 2;
    unsigned int i;
    for(i = 0; i < n; i+=1) {
        if (i) {
            (*pstr)[len++] = ',';
            (*pstr)[len++] = ' ';
        }
        char *tmp;
        mpfq_u64k1_asprint(K, &tmp, w[i]);
        int ltmp = strlen(tmp);
        if (len+ltmp+4 > alloc) {
            alloc = len+ltmp+100 + alloc / 4;
            *pstr = (char *)realloc(*pstr, alloc*sizeof(char));
        }
        strncpy(*pstr+len, tmp, ltmp+4);
        len += ltmp;
        free(tmp);
    }
    (*pstr)[len++] = ' ';
    (*pstr)[len++] = ']';
    (*pstr)[len] = '\0';
    return len;
}

/* *Mpfq::defaults::vec::io::code_for_vec_fprint, Mpfq::defaults::vec */
int mpfq_u64k1_vec_fprint(mpfq_u64k1_dst_field K MAYBE_UNUSED, FILE * file, mpfq_u64k1_src_vec w, unsigned int n)
{
    char *str;
    int rc;
    mpfq_u64k1_vec_asprint(K,&str,w,n);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *Mpfq::defaults::vec::io::code_for_vec_print, Mpfq::defaults::vec */
int mpfq_u64k1_vec_print(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_src_vec w, unsigned int n)
{
    return mpfq_u64k1_vec_fprint(K,stdout,w,n);
}

/* *Mpfq::defaults::vec::io::code_for_vec_sscan, Mpfq::defaults::vec */
int mpfq_u64k1_vec_sscan(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec * w, unsigned int * n, const char * str)
{
    // start with a clean vector
    unsigned int nn;
    int len = 0;
    mpfq_u64k1_vec_reinit(K, w, *n, 0);
    *n = nn = 0;
    while (isspace((int)(unsigned char)str[len]))
        len++;
    if (str[len] != '[')
        return 0;
    len++;
    while (isspace((int)(unsigned char)str[len]))
        len++;
    if (str[len] == ']') {
        len++;
        return len;
    }
    unsigned int i = 0;
    for (;;) {
        if (nn < i+1) {
            mpfq_u64k1_vec_reinit(K, w, nn, i+1);
            *n = nn = i+1;
        }
        int ret = mpfq_u64k1_sscan(K, mpfq_u64k1_vec_coeff_ptr(K, *w, i), str + len);
        if (!ret) {
            *n = 0; /* invalidate data ! */
            return 0;
        }
        i++;
        len += ret;
        while (isspace((int)(unsigned char)str[len]))
            len++;
        if (str[len] == ']') {
            len++;
            break;
        }
        if (str[len] != ',') {
            *n = 0; /* invalidate data ! */
            return 0;
        }
        len++;
        while (isspace((int)(unsigned char)str[len]))
            len++;
    }
    return len;
}

/* *Mpfq::defaults::vec::io::code_for_vec_fscan, Mpfq::defaults::vec */
int mpfq_u64k1_vec_fscan(mpfq_u64k1_dst_field K MAYBE_UNUSED, FILE * file, mpfq_u64k1_vec * w, unsigned int * n)
{
    char *tmp;
    int c;
    int allocated, len=0;
    allocated=100;
    tmp = (char *)mpfq_malloc_check(allocated*sizeof(char));
    int nest = 0, mnest = 0;
    for(;;) {
        c = fgetc(file);
        if (c==EOF) {
            free(tmp);
            return 0;
        }
        if (c == '[') {
            nest++, mnest++;
        }
        if (len==allocated) {
            allocated = len + 10 + allocated / 4;
            tmp = (char*)realloc(tmp, allocated*sizeof(char));
        }
        tmp[len]=c;
        len++;
        if (c == ']') {
            nest--, mnest++;
        }
        if (mnest && nest == 0)
            break;
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated*sizeof(char));
    }
    tmp[len]='\0';
    int ret=mpfq_u64k1_vec_sscan(K,w,n,tmp);
    free(tmp);
    return ret;
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_init, Mpfq::defaults::vec */
void mpfq_u64k1_vec_ur_init(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec_ur * v, unsigned int n)
{
    unsigned int i;
    *v = (mpfq_u64k1_vec_ur) malloc (n*sizeof(mpfq_u64k1_elt_ur));
    for(i = 0; i < n; i+=1)
        mpfq_u64k1_elt_ur_init(K, &( (*v)[i]));
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_reinit, Mpfq::defaults::vec */
void mpfq_u64k1_vec_ur_reinit(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec_ur * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        *v = (mpfq_u64k1_vec_ur) realloc (*v, m * sizeof(mpfq_u64k1_elt_ur));
        unsigned int i;
        for(i = n; i < m; i+=1)
            mpfq_u64k1_elt_ur_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_u64k1_elt_ur_clear(K, (*v) + i);
        *v = (mpfq_u64k1_vec_ur) realloc (*v, m * sizeof(mpfq_u64k1_elt_ur));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_clear, Mpfq::defaults::vec */
void mpfq_u64k1_vec_ur_clear(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_vec_ur * v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_u64k1_elt_ur_clear(K, &( (*v)[i]));
    free(*v);
}

/* missing vec_scal_mul_ur */
/* missing vec_reduce */

/* Polynomial functions */

/* Functions related to SIMD operation */
/* *simd_dotprod::code_for_dotprod */
void mpfq_u64k1_dotprod(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k1_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = xw[0];
    const uint64_t * u0 = xu0[0];
    const uint64_t * u1 = xu1[0];
    dotprod_64K_64(w,u1,u0,n,1);
}


/* Member templates related to SIMD operation */

/* MPI interface */
static void mpfq_u64k1_mpi_op_inner_ur(void *, void *, int *, MPI_Datatype *);
/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
/* Triggered by: mpi_ops_init */
static void mpfq_u64k1_mpi_op_inner_ur(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    mpfq_u64k1_dst_field K;
    MPI_Type_get_attr(*datatype, mpfq_u64k1_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    mpfq_u64k1_vec_ur_add(K, inoutvec, inoutvec, invec, *len);
}

static void mpfq_u64k1_mpi_op_inner(void *, void *, int *, MPI_Datatype *);
/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
/* Triggered by: mpi_ops_init */
static void mpfq_u64k1_mpi_op_inner(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    mpfq_u64k1_dst_field K;
    MPI_Type_get_attr(*datatype, mpfq_u64k1_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    mpfq_u64k1_vec_add(K, inoutvec, inoutvec, invec, *len);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void mpfq_u64k1_mpi_ops_init(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
        if (mpfq_u64k1_impl_mpi_use_count++) return;
    MPI_Type_create_keyval(MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN, &mpfq_u64k1_impl_mpi_attr, NULL);
    MPI_Type_contiguous(mpfq_u64k1_vec_elt_stride(K, 1), MPI_BYTE, &mpfq_u64k1_impl_mpi_datatype);
    MPI_Type_commit(&mpfq_u64k1_impl_mpi_datatype);
    MPI_Type_contiguous(mpfq_u64k1_vec_ur_elt_stride(K, 1), MPI_BYTE, &mpfq_u64k1_impl_mpi_datatype_ur);
    MPI_Type_commit(&mpfq_u64k1_impl_mpi_datatype_ur);
    MPI_Type_set_attr(mpfq_u64k1_impl_mpi_datatype, mpfq_u64k1_impl_mpi_attr, K);
    MPI_Type_set_attr(mpfq_u64k1_impl_mpi_datatype_ur, mpfq_u64k1_impl_mpi_attr, K);
    /* 1 here indicates that our operation is always taken to be
     * commutative */
    MPI_Op_create(&mpfq_u64k1_mpi_op_inner, 1, &mpfq_u64k1_impl_mpi_addition_op);
    MPI_Op_create(&mpfq_u64k1_mpi_op_inner_ur, 1, &mpfq_u64k1_impl_mpi_addition_op_ur);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype */
MPI_Datatype mpfq_u64k1_mpi_datatype(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
    return mpfq_u64k1_impl_mpi_datatype;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype_ur */
MPI_Datatype mpfq_u64k1_mpi_datatype_ur(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
    return mpfq_u64k1_impl_mpi_datatype_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op */
MPI_Op mpfq_u64k1_mpi_addition_op(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
    return mpfq_u64k1_impl_mpi_addition_op;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op_ur */
MPI_Op mpfq_u64k1_mpi_addition_op_ur(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
    return mpfq_u64k1_impl_mpi_addition_op_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_clear */
void mpfq_u64k1_mpi_ops_clear(mpfq_u64k1_dst_field K MAYBE_UNUSED)
{
        if (--mpfq_u64k1_impl_mpi_use_count) return;
    MPI_Op_free(&mpfq_u64k1_impl_mpi_addition_op);
    MPI_Op_free(&mpfq_u64k1_impl_mpi_addition_op_ur);
    MPI_Type_delete_attr(mpfq_u64k1_impl_mpi_datatype, mpfq_u64k1_impl_mpi_attr);
    MPI_Type_delete_attr(mpfq_u64k1_impl_mpi_datatype_ur, mpfq_u64k1_impl_mpi_attr);
    MPI_Type_free(&mpfq_u64k1_impl_mpi_datatype);
    MPI_Type_free(&mpfq_u64k1_impl_mpi_datatype_ur);
    MPI_Type_free_keyval(&mpfq_u64k1_impl_mpi_attr);
}


/* Object-oriented interface */
static void mpfq_u64k1_wrapper_oo_field_clear(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_oo_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_oo_field_clear(vbase);
}

static void mpfq_u64k1_wrapper_oo_field_init(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_oo_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_oo_field_init(vbase);
}

static void mpfq_u64k1_wrapper_mpi_ops_clear(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_mpi_ops_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_mpi_ops_clear(vbase->obj);
}

static MPI_Op mpfq_u64k1_wrapper_mpi_addition_op_ur(mpfq_vbase_ptr);
static MPI_Op mpfq_u64k1_wrapper_mpi_addition_op_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_mpi_addition_op_ur(vbase->obj);
}

static MPI_Op mpfq_u64k1_wrapper_mpi_addition_op(mpfq_vbase_ptr);
static MPI_Op mpfq_u64k1_wrapper_mpi_addition_op(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_mpi_addition_op(vbase->obj);
}

static MPI_Datatype mpfq_u64k1_wrapper_mpi_datatype_ur(mpfq_vbase_ptr);
static MPI_Datatype mpfq_u64k1_wrapper_mpi_datatype_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_mpi_datatype_ur(vbase->obj);
}

static MPI_Datatype mpfq_u64k1_wrapper_mpi_datatype(mpfq_vbase_ptr);
static MPI_Datatype mpfq_u64k1_wrapper_mpi_datatype(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_mpi_datatype(vbase->obj);
}

static void mpfq_u64k1_wrapper_mpi_ops_init(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_mpi_ops_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_mpi_ops_init(vbase->obj);
}

static void mpfq_u64k1_wrapper_dotprod(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_dotprod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec xw MAYBE_UNUSED, mpfq_u64k1_src_vec xu1 MAYBE_UNUSED, mpfq_u64k1_src_vec xu0 MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_dotprod(vbase->obj, xw, xu1, xu0, n);
}

static void mpfq_u64k1_wrapper_elt_ur_set_ui_all(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, unsigned long);
static void mpfq_u64k1_wrapper_elt_ur_set_ui_all(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_set_ui_all(vbase->obj, r, v);
}

static void mpfq_u64k1_wrapper_elt_ur_set_ui_at(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, int, unsigned long);
static void mpfq_u64k1_wrapper_elt_ur_set_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_set_ui_at(vbase->obj, p, k, v);
}

static void mpfq_u64k1_wrapper_set_ui_all(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, unsigned long);
static void mpfq_u64k1_wrapper_set_ui_all(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_u64k1_set_ui_all(vbase->obj, r, v);
}

static void mpfq_u64k1_wrapper_set_ui_at(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, int, unsigned long);
static void mpfq_u64k1_wrapper_set_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_u64k1_set_ui_at(vbase->obj, p, k, v);
}

static int mpfq_u64k1_wrapper_stride(mpfq_vbase_ptr);
static int mpfq_u64k1_wrapper_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_stride(vbase->obj);
}

static int mpfq_u64k1_wrapper_offset(mpfq_vbase_ptr, int);
static int mpfq_u64k1_wrapper_offset(mpfq_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return mpfq_u64k1_offset(vbase->obj, n);
}

static int mpfq_u64k1_wrapper_groupsize(mpfq_vbase_ptr);
static int mpfq_u64k1_wrapper_groupsize(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_groupsize(vbase->obj);
}

static ptrdiff_t mpfq_u64k1_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr, int);
static ptrdiff_t mpfq_u64k1_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_ur_elt_stride(vbase->obj, n);
}

static ptrdiff_t mpfq_u64k1_wrapper_vec_elt_stride(mpfq_vbase_ptr, int);
static ptrdiff_t mpfq_u64k1_wrapper_vec_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_elt_stride(vbase->obj, n);
}

static mpfq_u64k1_src_elt mpfq_u64k1_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr, mpfq_u64k1_src_vec_ur, int);
static mpfq_u64k1_src_elt mpfq_u64k1_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_ur_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_u64k1_dst_elt mpfq_u64k1_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, int);
static mpfq_u64k1_dst_elt mpfq_u64k1_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_ur_coeff_ptr(vbase->obj, v, i);
}

static mpfq_u64k1_src_vec_ur mpfq_u64k1_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr, mpfq_u64k1_src_vec_ur, int);
static mpfq_u64k1_src_vec_ur mpfq_u64k1_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_ur_subvec_const(vbase->obj, v, i);
}

static mpfq_u64k1_dst_vec_ur mpfq_u64k1_wrapper_vec_ur_subvec(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, int);
static mpfq_u64k1_dst_vec_ur mpfq_u64k1_wrapper_vec_ur_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_ur_subvec(vbase->obj, v, i);
}

static void mpfq_u64k1_wrapper_vec_ur_rev(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_vec_ur u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_rev(vbase->obj, w, u, n);
}

static void mpfq_u64k1_wrapper_vec_ur_neg(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_vec_ur u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_neg(vbase->obj, w, u, n);
}

static void mpfq_u64k1_wrapper_vec_ur_sub(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_vec_ur u MAYBE_UNUSED, mpfq_u64k1_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_sub(vbase->obj, w, u, v, n);
}

static void mpfq_u64k1_wrapper_vec_ur_add(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_vec_ur u MAYBE_UNUSED, mpfq_u64k1_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_add(vbase->obj, w, u, v, n);
}

static void mpfq_u64k1_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur x MAYBE_UNUSED, mpfq_u64k1_src_vec_ur w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_u64k1_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_elt_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_elt_ur x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_u64k1_wrapper_vec_ur_set(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur r MAYBE_UNUSED, mpfq_u64k1_src_vec_ur s MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_set(vbase->obj, r, s, n);
}

static void mpfq_u64k1_wrapper_vec_ur_clear(mpfq_vbase_ptr, mpfq_u64k1_vec_ur *, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_clear(vbase->obj, v, n);
}

static void mpfq_u64k1_wrapper_vec_ur_reinit(mpfq_vbase_ptr, mpfq_u64k1_vec_ur *, unsigned int, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_reinit(vbase->obj, v, n, m);
}

static void mpfq_u64k1_wrapper_vec_ur_set_vec(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_set_vec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur w MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_set_vec(vbase->obj, w, u, n);
}

static void mpfq_u64k1_wrapper_vec_ur_set_zero(mpfq_vbase_ptr, mpfq_u64k1_dst_vec_ur, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec_ur r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_set_zero(vbase->obj, r, n);
}

static void mpfq_u64k1_wrapper_vec_ur_init(mpfq_vbase_ptr, mpfq_u64k1_vec_ur *, unsigned int);
static void mpfq_u64k1_wrapper_vec_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_ur_init(vbase->obj, v, n);
}

static int mpfq_u64k1_wrapper_vec_scan(mpfq_vbase_ptr, mpfq_u64k1_vec *, unsigned int *);
static int mpfq_u64k1_wrapper_vec_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_scan(vbase->obj, w, n);
}

static int mpfq_u64k1_wrapper_vec_fscan(mpfq_vbase_ptr, FILE *, mpfq_u64k1_vec *, unsigned int *);
static int mpfq_u64k1_wrapper_vec_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_u64k1_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_fscan(vbase->obj, file, w, n);
}

static int mpfq_u64k1_wrapper_vec_sscan(mpfq_vbase_ptr, mpfq_u64k1_vec *, unsigned int *, const char *);
static int mpfq_u64k1_wrapper_vec_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_sscan(vbase->obj, w, n, str);
}

static int mpfq_u64k1_wrapper_vec_print(mpfq_vbase_ptr, mpfq_u64k1_src_vec, unsigned int);
static int mpfq_u64k1_wrapper_vec_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_print(vbase->obj, w, n);
}

static int mpfq_u64k1_wrapper_vec_fprint(mpfq_vbase_ptr, FILE *, mpfq_u64k1_src_vec, unsigned int);
static int mpfq_u64k1_wrapper_vec_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_u64k1_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_fprint(vbase->obj, file, w, n);
}

static int mpfq_u64k1_wrapper_vec_asprint(mpfq_vbase_ptr, char * *, mpfq_u64k1_src_vec, unsigned int);
static int mpfq_u64k1_wrapper_vec_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_u64k1_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_asprint(vbase->obj, pstr, w, n);
}

static mpfq_u64k1_src_elt mpfq_u64k1_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr, mpfq_u64k1_src_vec, int);
static mpfq_u64k1_src_elt mpfq_u64k1_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_u64k1_dst_elt mpfq_u64k1_wrapper_vec_coeff_ptr(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, int);
static mpfq_u64k1_dst_elt mpfq_u64k1_wrapper_vec_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_coeff_ptr(vbase->obj, v, i);
}

static mpfq_u64k1_src_vec mpfq_u64k1_wrapper_vec_subvec_const(mpfq_vbase_ptr, mpfq_u64k1_src_vec, int);
static mpfq_u64k1_src_vec mpfq_u64k1_wrapper_vec_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_subvec_const(vbase->obj, v, i);
}

static mpfq_u64k1_dst_vec mpfq_u64k1_wrapper_vec_subvec(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, int);
static mpfq_u64k1_dst_vec mpfq_u64k1_wrapper_vec_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec v MAYBE_UNUSED, int i MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_subvec(vbase->obj, v, i);
}

static int mpfq_u64k1_wrapper_vec_is_zero(mpfq_vbase_ptr, mpfq_u64k1_src_vec, unsigned int);
static int mpfq_u64k1_wrapper_vec_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_is_zero(vbase->obj, r, n);
}

static int mpfq_u64k1_wrapper_vec_cmp(mpfq_vbase_ptr, mpfq_u64k1_src_vec, mpfq_u64k1_src_vec, unsigned int);
static int mpfq_u64k1_wrapper_vec_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, mpfq_u64k1_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return mpfq_u64k1_vec_cmp(vbase->obj, u, v, n);
}

static void mpfq_u64k1_wrapper_vec_random(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, unsigned int, gmp_randstate_t);
static void mpfq_u64k1_wrapper_vec_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_u64k1_vec_random(vbase->obj, w, n, state);
}

static void mpfq_u64k1_wrapper_vec_sub(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, mpfq_u64k1_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_sub(vbase->obj, w, u, v, n);
}

static void mpfq_u64k1_wrapper_vec_rev(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_rev(vbase->obj, w, u, n);
}

static void mpfq_u64k1_wrapper_vec_neg(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_neg(vbase->obj, w, u, n);
}

static void mpfq_u64k1_wrapper_vec_add(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, mpfq_u64k1_src_vec u MAYBE_UNUSED, mpfq_u64k1_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_add(vbase->obj, w, u, v, n);
}

static void mpfq_u64k1_wrapper_vec_getcoeff(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt x MAYBE_UNUSED, mpfq_u64k1_src_vec w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_u64k1_vec_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_u64k1_wrapper_vec_setcoeff(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_elt, unsigned int);
static void mpfq_u64k1_wrapper_vec_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec w MAYBE_UNUSED, mpfq_u64k1_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_u64k1_vec_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_u64k1_wrapper_vec_set_zero(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_set_zero(vbase->obj, r, n);
}

static void mpfq_u64k1_wrapper_vec_set(mpfq_vbase_ptr, mpfq_u64k1_dst_vec, mpfq_u64k1_src_vec, unsigned int);
static void mpfq_u64k1_wrapper_vec_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_vec r MAYBE_UNUSED, mpfq_u64k1_src_vec s MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_set(vbase->obj, r, s, n);
}

static void mpfq_u64k1_wrapper_vec_clear(mpfq_vbase_ptr, mpfq_u64k1_vec *, unsigned int);
static void mpfq_u64k1_wrapper_vec_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_clear(vbase->obj, v, n);
}

static void mpfq_u64k1_wrapper_vec_reinit(mpfq_vbase_ptr, mpfq_u64k1_vec *, unsigned int, unsigned int);
static void mpfq_u64k1_wrapper_vec_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    mpfq_u64k1_vec_reinit(vbase->obj, v, n, m);
}

static void mpfq_u64k1_wrapper_vec_init(mpfq_vbase_ptr, mpfq_u64k1_vec *, unsigned int);
static void mpfq_u64k1_wrapper_vec_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_u64k1_vec_init(vbase->obj, v, n);
}

static int mpfq_u64k1_wrapper_scan(mpfq_vbase_ptr, mpfq_u64k1_dst_elt);
static int mpfq_u64k1_wrapper_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt x MAYBE_UNUSED)
{
    return mpfq_u64k1_scan(vbase->obj, x);
}

static int mpfq_u64k1_wrapper_fscan(mpfq_vbase_ptr, FILE *, mpfq_u64k1_dst_elt);
static int mpfq_u64k1_wrapper_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_u64k1_dst_elt z MAYBE_UNUSED)
{
    return mpfq_u64k1_fscan(vbase->obj, file, z);
}

static int mpfq_u64k1_wrapper_sscan(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, const char *);
static int mpfq_u64k1_wrapper_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt z MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_u64k1_sscan(vbase->obj, z, str);
}

static int mpfq_u64k1_wrapper_print(mpfq_vbase_ptr, mpfq_u64k1_src_elt);
static int mpfq_u64k1_wrapper_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_elt x MAYBE_UNUSED)
{
    return mpfq_u64k1_print(vbase->obj, x);
}

static int mpfq_u64k1_wrapper_fprint(mpfq_vbase_ptr, FILE *, mpfq_u64k1_src_elt);
static int mpfq_u64k1_wrapper_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_u64k1_src_elt x MAYBE_UNUSED)
{
    return mpfq_u64k1_fprint(vbase->obj, file, x);
}

static int mpfq_u64k1_wrapper_asprint(mpfq_vbase_ptr, char * *, mpfq_u64k1_src_elt);
static int mpfq_u64k1_wrapper_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * ps MAYBE_UNUSED, mpfq_u64k1_src_elt x MAYBE_UNUSED)
{
    return mpfq_u64k1_asprint(vbase->obj, ps, x);
}

static int mpfq_u64k1_wrapper_is_zero(mpfq_vbase_ptr, mpfq_u64k1_src_elt);
static int mpfq_u64k1_wrapper_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_elt r MAYBE_UNUSED)
{
    return mpfq_u64k1_is_zero(vbase->obj, r);
}

static int mpfq_u64k1_wrapper_cmp(mpfq_vbase_ptr, mpfq_u64k1_src_elt, mpfq_u64k1_src_elt);
static int mpfq_u64k1_wrapper_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_src_elt r MAYBE_UNUSED, mpfq_u64k1_src_elt s MAYBE_UNUSED)
{
    return mpfq_u64k1_cmp(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_elt_ur_sub(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_elt_ur, mpfq_u64k1_src_elt_ur);
static void mpfq_u64k1_wrapper_elt_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s1 MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s2 MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_sub(vbase->obj, r, s1, s2);
}

static void mpfq_u64k1_wrapper_elt_ur_neg(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_elt_ur);
static void mpfq_u64k1_wrapper_elt_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_neg(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_elt_ur_add(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_elt_ur, mpfq_u64k1_src_elt_ur);
static void mpfq_u64k1_wrapper_elt_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s1 MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s2 MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_add(vbase->obj, r, s1, s2);
}

static void mpfq_u64k1_wrapper_elt_ur_set_zero(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur);
static void mpfq_u64k1_wrapper_elt_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_set_zero(vbase->obj, r);
}

static void mpfq_u64k1_wrapper_elt_ur_set_elt(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_elt);
static void mpfq_u64k1_wrapper_elt_ur_set_elt(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED, mpfq_u64k1_src_elt s MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_set_elt(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_elt_ur_set(mpfq_vbase_ptr, mpfq_u64k1_dst_elt_ur, mpfq_u64k1_src_elt_ur);
static void mpfq_u64k1_wrapper_elt_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt_ur r MAYBE_UNUSED, mpfq_u64k1_src_elt_ur s MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_set(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_elt_ur_clear(mpfq_vbase_ptr, mpfq_u64k1_elt_ur *);
static void mpfq_u64k1_wrapper_elt_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_elt_ur * px MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_clear(vbase->obj, px);
}

static void mpfq_u64k1_wrapper_elt_ur_init(mpfq_vbase_ptr, mpfq_u64k1_elt_ur *);
static void mpfq_u64k1_wrapper_elt_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_elt_ur * px MAYBE_UNUSED)
{
    mpfq_u64k1_elt_ur_init(vbase->obj, px);
}

static void mpfq_u64k1_wrapper_neg(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, mpfq_u64k1_src_elt);
static void mpfq_u64k1_wrapper_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, mpfq_u64k1_src_elt s MAYBE_UNUSED)
{
    mpfq_u64k1_neg(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_sub(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, mpfq_u64k1_src_elt, mpfq_u64k1_src_elt);
static void mpfq_u64k1_wrapper_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, mpfq_u64k1_src_elt s1 MAYBE_UNUSED, mpfq_u64k1_src_elt s2 MAYBE_UNUSED)
{
    mpfq_u64k1_sub(vbase->obj, r, s1, s2);
}

static void mpfq_u64k1_wrapper_add(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, mpfq_u64k1_src_elt, mpfq_u64k1_src_elt);
static void mpfq_u64k1_wrapper_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, mpfq_u64k1_src_elt s1 MAYBE_UNUSED, mpfq_u64k1_src_elt s2 MAYBE_UNUSED)
{
    mpfq_u64k1_add(vbase->obj, r, s1, s2);
}

static void mpfq_u64k1_wrapper_random(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, gmp_randstate_t);
static void mpfq_u64k1_wrapper_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_u64k1_random(vbase->obj, r, state);
}

static void mpfq_u64k1_wrapper_set_zero(mpfq_vbase_ptr, mpfq_u64k1_dst_elt);
static void mpfq_u64k1_wrapper_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED)
{
    mpfq_u64k1_set_zero(vbase->obj, r);
}

static void mpfq_u64k1_wrapper_set(mpfq_vbase_ptr, mpfq_u64k1_dst_elt, mpfq_u64k1_src_elt);
static void mpfq_u64k1_wrapper_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_dst_elt r MAYBE_UNUSED, mpfq_u64k1_src_elt s MAYBE_UNUSED)
{
    mpfq_u64k1_set(vbase->obj, r, s);
}

static void mpfq_u64k1_wrapper_clear(mpfq_vbase_ptr, mpfq_u64k1_elt *);
static void mpfq_u64k1_wrapper_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_elt * px MAYBE_UNUSED)
{
    mpfq_u64k1_clear(vbase->obj, px);
}

static void mpfq_u64k1_wrapper_init(mpfq_vbase_ptr, mpfq_u64k1_elt *);
static void mpfq_u64k1_wrapper_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_u64k1_elt * px MAYBE_UNUSED)
{
    mpfq_u64k1_init(vbase->obj, px);
}

static void mpfq_u64k1_wrapper_field_setopt(mpfq_vbase_ptr, unsigned long, void *);
static void mpfq_u64k1_wrapper_field_setopt(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, void * y MAYBE_UNUSED)
{
    mpfq_u64k1_field_setopt(vbase->obj, x, y);
}

static void mpfq_u64k1_wrapper_field_specify(mpfq_vbase_ptr, unsigned long, void *);
static void mpfq_u64k1_wrapper_field_specify(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long tag MAYBE_UNUSED, void * x MAYBE_UNUSED)
{
    mpfq_u64k1_field_specify(vbase->obj, tag, x);
}

static void mpfq_u64k1_wrapper_field_clear(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_field_clear(vbase->obj);
}

static void mpfq_u64k1_wrapper_field_init(mpfq_vbase_ptr);
static void mpfq_u64k1_wrapper_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_u64k1_field_init(vbase->obj);
}

static int mpfq_u64k1_wrapper_field_degree(mpfq_vbase_ptr);
static int mpfq_u64k1_wrapper_field_degree(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_u64k1_field_degree(vbase->obj);
}

static void mpfq_u64k1_wrapper_field_characteristic(mpfq_vbase_ptr, mpz_t);
static void mpfq_u64k1_wrapper_field_characteristic(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpz_t z MAYBE_UNUSED)
{
    mpfq_u64k1_field_characteristic(vbase->obj, z);
}

static unsigned long mpfq_u64k1_wrapper_impl_max_degree();
static unsigned long mpfq_u64k1_wrapper_impl_max_degree()
{
    return mpfq_u64k1_impl_max_degree();
}

static unsigned long mpfq_u64k1_wrapper_impl_max_characteristic_bits();
static unsigned long mpfq_u64k1_wrapper_impl_max_characteristic_bits()
{
    return mpfq_u64k1_impl_max_characteristic_bits();
}

static const char * mpfq_u64k1_wrapper_impl_name();
static const char * mpfq_u64k1_wrapper_impl_name()
{
    return mpfq_u64k1_impl_name();
}

/* Mpfq::engine::oo::oo_field_init */
/* Triggered by: oo */
void mpfq_u64k1_oo_field_init(mpfq_vbase_ptr vbase)
{
    memset(vbase, 0, sizeof(struct mpfq_vbase_s));
    vbase->obj = malloc(sizeof(mpfq_u64k1_field));
    mpfq_u64k1_field_init((mpfq_u64k1_dst_field) vbase->obj);
    vbase->impl_name = (const char * (*) ()) mpfq_u64k1_wrapper_impl_name;
    vbase->impl_max_characteristic_bits = (unsigned long (*) ()) mpfq_u64k1_wrapper_impl_max_characteristic_bits;
    vbase->impl_max_degree = (unsigned long (*) ()) mpfq_u64k1_wrapper_impl_max_degree;
    vbase->field_characteristic = (void (*) (mpfq_vbase_ptr, mpz_t)) mpfq_u64k1_wrapper_field_characteristic;
    /* missing field_characteristic_bits */
    vbase->field_degree = (int (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_field_degree;
    vbase->field_init = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_field_init;
    vbase->field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_field_clear;
    vbase->field_specify = (void (*) (mpfq_vbase_ptr, unsigned long, void *)) mpfq_u64k1_wrapper_field_specify;
    vbase->field_setopt = (void (*) (mpfq_vbase_ptr, unsigned long, void *)) mpfq_u64k1_wrapper_field_setopt;
    vbase->init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_init;
    vbase->clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_clear;
    vbase->set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_wrapper_set;
    /* missing set_ui */
    vbase->set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_set_zero;
    /* missing get_ui */
    /* missing set_mpn */
    /* missing set_mpz */
    /* missing get_mpn */
    /* missing get_mpz */
    vbase->random = (void (*) (mpfq_vbase_ptr, void *, gmp_randstate_t)) mpfq_u64k1_wrapper_random;
    /* missing random2 */
    vbase->add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_u64k1_wrapper_add;
    vbase->sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_u64k1_wrapper_sub;
    vbase->neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_wrapper_neg;
    /* missing mul */
    /* missing sqr */
    /* missing is_sqr */
    /* missing sqrt */
    /* missing pow */
    /* missing powz */
    /* missing frobenius */
    /* missing add_ui */
    /* missing sub_ui */
    /* missing mul_ui */
    /* missing inv */
    vbase->elt_ur_init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_elt_ur_init;
    vbase->elt_ur_clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_elt_ur_clear;
    vbase->elt_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_wrapper_elt_ur_set;
    vbase->elt_ur_set_elt = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_wrapper_elt_ur_set_elt;
    vbase->elt_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_elt_ur_set_zero;
    /* missing elt_ur_set_ui */
    vbase->elt_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_u64k1_wrapper_elt_ur_add;
    vbase->elt_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_wrapper_elt_ur_neg;
    vbase->elt_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_u64k1_wrapper_elt_ur_sub;
    /* missing mul_ur */
    /* missing sqr_ur */
    /* missing reduce */
    vbase->cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *)) mpfq_u64k1_wrapper_cmp;
    /* missing cmp_ui */
    vbase->is_zero = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_u64k1_wrapper_is_zero;
    vbase->asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *)) mpfq_u64k1_wrapper_asprint;
    vbase->fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *)) mpfq_u64k1_wrapper_fprint;
    vbase->print = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_u64k1_wrapper_print;
    vbase->sscan = (int (*) (mpfq_vbase_ptr, void *, const char *)) mpfq_u64k1_wrapper_sscan;
    vbase->fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *)) mpfq_u64k1_wrapper_fscan;
    vbase->scan = (int (*) (mpfq_vbase_ptr, void *)) mpfq_u64k1_wrapper_scan;
    vbase->vec_init = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_init;
    vbase->vec_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned int, unsigned int)) mpfq_u64k1_wrapper_vec_reinit;
    vbase->vec_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_clear;
    vbase->vec_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_set;
    vbase->vec_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_set_zero;
    vbase->vec_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_setcoeff;
    /* missing vec_setcoeff_ui */
    vbase->vec_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_getcoeff;
    vbase->vec_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_add;
    vbase->vec_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_neg;
    vbase->vec_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_rev;
    vbase->vec_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_sub;
    /* missing vec_scal_mul */
    /* missing vec_conv */
    vbase->vec_random = (void (*) (mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t)) mpfq_u64k1_wrapper_vec_random;
    /* missing vec_random2 */
    vbase->vec_cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_cmp;
    vbase->vec_is_zero = (int (*) (mpfq_vbase_ptr, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_is_zero;
    vbase->vec_subvec = (void * (*) (mpfq_vbase_ptr, void *, int)) mpfq_u64k1_wrapper_vec_subvec;
    vbase->vec_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, int)) mpfq_u64k1_wrapper_vec_subvec_const;
    vbase->vec_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, int)) mpfq_u64k1_wrapper_vec_coeff_ptr;
    vbase->vec_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, int)) mpfq_u64k1_wrapper_vec_coeff_ptr_const;
    vbase->vec_asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_asprint;
    vbase->vec_fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_fprint;
    vbase->vec_print = (int (*) (mpfq_vbase_ptr, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_print;
    vbase->vec_sscan = (int (*) (mpfq_vbase_ptr, void *, unsigned int *, const char *)) mpfq_u64k1_wrapper_vec_sscan;
    vbase->vec_fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *, unsigned int *)) mpfq_u64k1_wrapper_vec_fscan;
    vbase->vec_scan = (int (*) (mpfq_vbase_ptr, void *, unsigned int *)) mpfq_u64k1_wrapper_vec_scan;
    vbase->vec_ur_init = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_init;
    vbase->vec_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_set_zero;
    vbase->vec_ur_set_vec = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_set_vec;
    vbase->vec_ur_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned int, unsigned int)) mpfq_u64k1_wrapper_vec_ur_reinit;
    vbase->vec_ur_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_clear;
    vbase->vec_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_set;
    vbase->vec_ur_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_setcoeff;
    vbase->vec_ur_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_getcoeff;
    vbase->vec_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_add;
    vbase->vec_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_sub;
    vbase->vec_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_neg;
    vbase->vec_ur_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_u64k1_wrapper_vec_ur_rev;
    /* missing vec_scal_mul_ur */
    /* missing vec_conv_ur */
    /* missing vec_reduce */
    vbase->vec_ur_subvec = (void * (*) (mpfq_vbase_ptr, void *, int)) mpfq_u64k1_wrapper_vec_ur_subvec;
    vbase->vec_ur_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, int)) mpfq_u64k1_wrapper_vec_ur_subvec_const;
    vbase->vec_ur_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, int)) mpfq_u64k1_wrapper_vec_ur_coeff_ptr;
    vbase->vec_ur_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, int)) mpfq_u64k1_wrapper_vec_ur_coeff_ptr_const;
    vbase->vec_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, int)) mpfq_u64k1_wrapper_vec_elt_stride;
    vbase->vec_ur_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, int)) mpfq_u64k1_wrapper_vec_ur_elt_stride;
    /* missing poly_init */
    /* missing poly_clear */
    /* missing poly_set */
    /* missing poly_setmonic */
    /* missing poly_setcoeff */
    /* missing poly_setcoeff_ui */
    /* missing poly_getcoeff */
    /* missing poly_deg */
    /* missing poly_add */
    /* missing poly_sub */
    /* missing poly_set_ui */
    /* missing poly_add_ui */
    /* missing poly_sub_ui */
    /* missing poly_neg */
    /* missing poly_scal_mul */
    /* missing poly_mul */
    /* missing poly_divmod */
    /* missing poly_precomp_mod */
    /* missing poly_mod_pre */
    /* missing poly_gcd */
    /* missing poly_xgcd */
    /* missing poly_random */
    /* missing poly_random2 */
    /* missing poly_cmp */
    /* missing poly_asprint */
    /* missing poly_fprint */
    /* missing poly_print */
    /* missing poly_sscan */
    /* missing poly_fscan */
    /* missing poly_scan */
    vbase->groupsize = (int (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_groupsize;
    vbase->offset = (int (*) (mpfq_vbase_ptr, int)) mpfq_u64k1_wrapper_offset;
    vbase->stride = (int (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_stride;
    vbase->set_ui_at = (void (*) (mpfq_vbase_ptr, void *, int, unsigned long)) mpfq_u64k1_wrapper_set_ui_at;
    vbase->set_ui_all = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_u64k1_wrapper_set_ui_all;
    vbase->elt_ur_set_ui_at = (void (*) (mpfq_vbase_ptr, void *, int, unsigned long)) mpfq_u64k1_wrapper_elt_ur_set_ui_at;
    vbase->elt_ur_set_ui_all = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_u64k1_wrapper_elt_ur_set_ui_all;
    vbase->dotprod = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_wrapper_dotprod;
    vbase->mpi_ops_init = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_ops_init;
    vbase->mpi_datatype = (MPI_Datatype (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_datatype;
    vbase->mpi_datatype_ur = (MPI_Datatype (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_datatype_ur;
    vbase->mpi_addition_op = (MPI_Op (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_addition_op;
    vbase->mpi_addition_op_ur = (MPI_Op (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_addition_op_ur;
    vbase->mpi_ops_clear = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_mpi_ops_clear;
    vbase->oo_field_init = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_oo_field_init;
    vbase->oo_field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_u64k1_wrapper_oo_field_clear;
}


/* vim:set ft=cpp: */
