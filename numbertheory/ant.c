#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <stdint.h>
#include <sys/resource.h>	/* for getrusage */
#include "macros.h"
#include "mpz_poly.h"

/*{{{ typedefs */

// mpz_t : multiple precision Z-elements ; same for mpq_t
// a matric on mpz_t type
struct mpz_mat_s {
    mpz_t *x;
    unsigned int m,n;
};
typedef struct mpz_mat_s mpz_mat[1];
typedef struct mpz_mat_s * mpz_mat_ptr;
typedef const struct mpz_mat_s * mpz_mat_srcptr;

// same on mpq_t type
struct mpq_mat_s {
    mpq_t *x;
    unsigned int m,n;
};
typedef struct mpq_mat_s mpq_mat[1];
typedef struct mpq_mat_s * mpq_mat_ptr;
typedef const struct mpq_mat_s * mpq_mat_srcptr;
/*}}}*/
/*{{{ entry access*/
mpz_ptr mpz_mat_entry(mpz_mat_ptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpz_srcptr mpz_mat_entry_const(mpz_mat_srcptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpq_ptr mpq_mat_entry(mpq_mat_ptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpq_srcptr mpq_mat_entry_const(mpq_mat_srcptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}
/*}}}*/
/*{{{ init/clear/realloc*/
void mpz_mat_init(mpz_mat_ptr M, unsigned int m, unsigned int n)
{
    M->x = (mpz_t*) ((m*n) ? malloc(m * n * sizeof(mpz_t)) : NULL);

    M->m = m;
    M->n = n;
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_init(mpz_mat_entry(M, i, j));
}

void mpz_mat_clear(mpz_mat_ptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_clear(mpz_mat_entry(M, i, j));
    free(M->x);
}

void mpz_mat_realloc(mpz_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n) return;
    mpz_mat_clear(M);
    mpz_mat_init(M, m, n);
}

void mpq_mat_init(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    M->x = (mpq_t*) ((m*n) ? malloc(m * n * sizeof(mpq_t)) : NULL);
    M->m = m;
    M->n = n;
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_init(mpq_mat_entry(M, i, j));
}

void mpq_mat_clear(mpq_mat_ptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_clear(mpq_mat_entry(M, i, j));
    free(M->x);
}

void mpq_mat_realloc(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n) return;
    mpq_mat_clear(M);
    mpq_mat_init(M, m, n);
}

/*}}}*/
/*{{{ operations on submatrices, and swaps*/
void mpz_mat_submat_swap(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpz_swap(
                    mpz_mat_entry(A0, i0 + i, j0 + j),
                    mpz_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpz_mat_swap(mpz_mat_ptr A, mpz_mat_ptr B)
{
    mpz_mat C;
    memcpy(C, A, sizeof(mpz_mat));
    memcpy(A, B, sizeof(mpz_mat));
    memcpy(B, C, sizeof(mpz_mat));
}

void mpq_mat_submat_swap(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpq_swap(
                    mpq_mat_entry(A0, i0 + i, j0 + j),
                    mpq_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpq_mat_swap(mpq_mat_ptr A, mpq_mat_ptr B)
{
    mpq_mat C;
    memcpy(C, A, sizeof(mpq_mat));
    memcpy(A, B, sizeof(mpq_mat));
    memcpy(B, C, sizeof(mpq_mat));
}

/*}}}*/
/*{{{ set/set_ui*/
void mpz_mat_set(mpz_mat_ptr dst, mpz_mat_srcptr src)
{
    mpz_mat_realloc(dst, src->m, src->n);
    for(unsigned int i = 0 ; i < src->m ; i++)
        for(unsigned int j = 0 ; j < src->n ; j++)
            mpz_set(mpz_mat_entry(dst, i, j), mpz_mat_entry_const(src, i, j));
}

void mpz_mat_set_ui(mpz_mat_ptr M, unsigned int a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_set_ui(mpz_mat_entry(M, i, j), i == j ? a : 0);
}

void mpq_mat_set_ui(mpq_mat_ptr M, unsigned int a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_set_ui(mpq_mat_entry(M, i, j), i == j ? a : 0, 1);
}

void mpq_mat_set(mpq_mat_ptr dst, mpq_mat_srcptr src)
{
    mpq_mat_realloc(dst, src->m, src->n);
    for(unsigned int i = 0 ; i < src->m ; i++)
        for(unsigned int j = 0 ; j < src->n ; j++)
            mpq_set(mpq_mat_entry(dst, i, j), mpq_mat_entry_const(src, i, j));
}

void mpz_mat_urandomm(mpz_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_urandomm(mpz_mat_entry(M, i, j), state, p);
}

void mpq_mat_urandomm(mpq_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_urandomm(mpq_numref(mpq_mat_entry(M, i, j)), state, p);
            mpz_set_ui(mpq_denref(mpq_mat_entry(M, i, j)), 1);
        }
}



/*}}}*/
/*{{{ Joining matrices */
// I didn't write all possible functions there, only those I needed, but it could be completed
void mpz_mat_vertical_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2)
{
	ASSERT_ALWAYS(M1->n == M2->n);
	mpz_mat_realloc(N,M1->m+M2->m,M1->n);
	
	mpz_mat top,bottom;
	mpz_mat_init(top,M1->m,M1->n);
	mpz_mat_init(bottom,M2->m,M2->n);
	mpz_mat_set(top,M1);
	mpz_mat_set(bottom,M2);
	mpz_mat_submat_swap(top,0,0,N,0,0,top->m,top->n);
	mpz_mat_submat_swap(bottom,0,0,N,top->m,0,bottom->m,bottom->n);
	mpz_mat_clear(top);
	mpz_mat_clear(bottom);
	
}

void mpq_mat_vertical_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2)
{
	ASSERT_ALWAYS(M1->n == M2->n);
	mpq_mat_realloc(N,M1->m+M2->m,M1->n);
	
	mpq_mat top,bottom;
	mpq_mat_init(top,M1->m,M1->n);
	mpq_mat_init(bottom,M2->m,M2->n);
	mpq_mat_set(top,M1);
	mpq_mat_set(bottom,M2);
	mpq_mat_submat_swap(top,0,0,N,0,0,top->m,top->n);
	mpq_mat_submat_swap(bottom,0,0,N,top->m,0,bottom->m,bottom->n);
	mpq_mat_clear(top);
	mpq_mat_clear(bottom);
	
}
/*}}}*/
/*{{{ determinant and trace */
// We assume that M is a square Matrix and that its size is not 0
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M)
{
    mpz_set_ui(t, 0);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_add(t, t, mpz_mat_entry_const(M,i,i));
}

// We assume that M is a triangular Matrix and that its size is not 0
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M)
{
    mpz_set_ui(d, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_mul(d, d, mpz_mat_entry_const(M,i,i));
}
/*}}}*/
/*{{{ miscellaneous */
void mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M)
{
    mpz_set_ui(den, 1);
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_lcm(den, den, mpq_denref(mpq_mat_entry_const(M, i, j)));
    mpz_mat_realloc(num, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_ptr dst = mpz_mat_entry(num, i, j);
            mpq_srcptr src = mpq_mat_entry_const(M, i, j);
            mpz_divexact(dst, den, mpq_denref(src));
            mpz_mul(dst, dst, mpq_numref(src));
        }
}

void mpz_mat_to_mpq_mat(mpq_mat_ptr N, mpz_mat_srcptr M)
{
	mpq_mat_realloc(N,M->m,M->n);
	unsigned int i,j;
	for(i = 0 ; i < M->m ; i++) {
		for(j = 0 ; j < M->n ; j++) {
			mpq_set_z(mpq_mat_entry(N,i,j),mpz_mat_entry_const(M,i,j));
		}
	}
}

void mpz_mat_mod_ui(mpz_mat_ptr dst, mpz_mat_srcptr src, unsigned int p)
{
    unsigned int i,j;
    mpz_mat C;
    mpz_mat_init(C,src->m,src->n);
    for (i = 0 ; i < C->m ; i++){
        for (j = 0 ; j < C->n ; j++){
            mpz_fdiv_r_ui(mpz_mat_entry(C,i,j),mpz_mat_entry_const(src,i,j),p);
        }
    }
    mpz_mat_realloc(dst,src->m,src->n);
    mpz_mat_set(dst,C);
    mpz_mat_clear(C);
}
/*}}}*/
/* {{{ row-level operations */
// Return 1 if the k-th line of M is null, 0 else
int mpz_mat_isnull_row(mpz_mat_srcptr M, unsigned int k){
    unsigned int j = 0;
    ASSERT_ALWAYS(k < M->m);
    while((j < M->n) && !mpz_cmp_si(mpz_mat_entry_const(M,k,j),0)){
        j++;
    }
    if(j == M->n){
        return 1;
    }
    return 0;
}

void mpz_mat_swaprows(mpz_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_swap(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i1, j));
    }
}

void mpq_mat_swaprows(mpq_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_swap(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i1, j));
    }
}

/* apply a circular shift on rows [i0...i0+k[ */
void mpz_mat_rotatedownrows(mpz_mat_ptr M, unsigned int i0, unsigned int k)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i0 + k <= M->m);
    if (k <= 1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        for(unsigned int s = k - 1 ; s-- ; ) {
            mpz_swap(mpz_mat_entry(M, i0 + s, j), mpz_mat_entry(M, i0 + s + 1, j));
        }
    }
}

void mpq_mat_rotatedownrows(mpq_mat_ptr M, unsigned int i0, unsigned int k)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i0 + k <= M->m);
    if (k <= 1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        for(unsigned int s = k - 1 ; s-- ; ) {
            mpq_swap(mpq_mat_entry(M, i0 + s, j), mpq_mat_entry(M, i0 + s + 1, j));
        }
    }
}

/* put row perm[k] in row k */
void mpz_mat_permuterows(mpz_mat_ptr M, unsigned int * perm)
{
    mpz_mat Mx;
    mpz_mat_init(Mx, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_swap(mpz_mat_entry(Mx, i, j),
                    mpz_mat_entry(M, perm[i], j));
        }
    }
    mpz_mat_swap(M, Mx);
    mpz_mat_clear(Mx);
}

void mpq_mat_permuterows(mpq_mat_ptr M, unsigned int * perm)
{
    mpq_mat Mx;
    mpq_mat_init(Mx, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpq_swap(mpq_mat_entry(Mx, i, j),
                    mpq_mat_entry(M, perm[i], j));
        }
    }
    mpq_mat_swap(M, Mx);
    mpq_mat_clear(Mx);
}

int permutation_signature(const unsigned int * perm, unsigned int n)
{
    /* can this be done in O(n) only ?? */
    int sign = 1;
    for(unsigned int i = 0 ; i < n ; i++)
        for(unsigned int j = i ; j < n ; j++)
            if (perm[j] < perm[i]) sign*=-1;
    return sign;
}


/* add lambda times row i1 to row i0 */
void mpz_mat_addmulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_addmul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
    }
}

void mpz_mat_addmulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_addmul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}

/* add lambda times row i1 to row i0 */
void mpq_mat_addmulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        /* we have no mpq_addmul, unfortunately */
        mpq_mul(tmp,
                lambda,
                mpq_mat_entry(M, i1, j));
        mpq_add(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                tmp);
    }
    mpq_clear(tmp);
}

/* subtract lambda times row i1 to row i0 */
void mpz_mat_submulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
    }
}

void mpz_mat_submulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}

/* subtract lambda times row i1 to row i0 */
void mpq_mat_submulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        /* we have no mpq_submul, unfortunately */
        mpq_mul(tmp,
                lambda,
                mpq_mat_entry(M, i1, j));
        mpq_sub(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                tmp);
    }
    mpq_clear(tmp);
}


/* add row i1 to row i0 */
void mpz_mat_addrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_add(mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}

/* add row i1 to row i0 */
void mpq_mat_addrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_add(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}


/* subtract row i1 to row i0 */
void mpz_mat_subrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_sub(mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}

/* subtract row i1 to row i0 */
void mpq_mat_subrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_sub(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}

/* multiply row i0 by lambda */
void mpz_mat_mulrow(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
    }
}

void mpz_mat_mulrow_mod(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda, mpz_srcptr p)
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}

/* multiply row i0 by lambda */
void mpq_mat_mulrow(mpq_mat_ptr M, unsigned int i0, mpq_srcptr lambda)
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_mul(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j), lambda);
    }
}

/* this computes an additive combination of n rows into row [didx] of the
 * initial matrix. We require that this destination row be cleared
 * initially.
 */
void mpz_mat_combinerows(mpz_mat_ptr M, unsigned int didx, unsigned int sidx,
        mpz_srcptr * lambda, unsigned int n)
{
    for(unsigned int j = 0 ; j < M->n ; j++) {
        ASSERT_ALWAYS(mpz_cmp_ui(mpz_mat_entry(M, didx, j), 0) == 0);
    }
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_mat_addmulrow(M, didx, sidx + i, lambda[i]);
    }
}

/* }}} */
/*{{{ I/O*/
void mpz_mat_fprint(FILE * stream, mpz_mat_srcptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            gmp_fprintf(stream, " %Zd", mpz_mat_entry_const(M, i, j));
        }
        fprintf(stream, "\n");
    }
}

void mpq_mat_fprint(FILE * stream, mpq_mat_srcptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            gmp_fprintf(stream, " %Qd", mpq_mat_entry_const(M, i, j));
        }
        fprintf(stream, "\n");
    }
}
/*}}}*/
/*{{{ multiplication */
void mpz_mat_multiply(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B)
{
    ASSERT_ALWAYS(A->n == B->m);
    mpz_mat C;
    mpz_mat_init(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j),
                        mpz_mat_entry_const(A, i, k),
                        mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_realloc(D, A->m, B->n);
    mpz_mat_set(D,C);
    mpz_mat_clear(C);
}

void mpz_mat_multiply_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B, unsigned int p)
{
    ASSERT_ALWAYS(A->n == B->m);
    mpz_mat C;
    mpz_mat_init(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j),
                        mpz_mat_entry_const(A, i, k),
                        mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_mod_ui(C,C,p);
    mpz_mat_realloc(D, A->m, B->n);
    mpz_mat_set(D,C);
    mpz_mat_clear(C);
}

void mpz_mat_multiply_by_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_ptr k)
{
    mpz_mat B;
	mpz_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpz_mul(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
    mpz_mat_realloc(C, A->m, A->n);
    mpz_mat_set(C,B);
    mpz_mat_clear(B);
}

void mpz_mat_multiply_by_si(mpz_mat_ptr C, mpz_mat_srcptr A, int k)
{
    mpz_mat B;
	mpz_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpz_mul_si(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
    mpz_mat_realloc(C, A->m, A->n);
    mpz_mat_set(C,B);
    mpz_mat_clear(B);
}

void mpz_mat_multiply_by_ui(mpz_mat_ptr C, mpz_mat_srcptr A, unsigned int k)
{
    mpz_mat B;
	mpz_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpz_mul_ui(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
    mpz_mat_realloc(C, A->m, A->n);
    mpz_mat_set(C,B);
    mpz_mat_clear(B);
}

// Returns A^n for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_power_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned int n){
    unsigned int k;
    mpz_mat C,D;
    ASSERT_ALWAYS(A->n == A->m);
    mpz_mat_init(C,A->m,A->n);
    mpz_mat_init(D,A->m,A->n);
    mpz_mat_set_ui(C,1);
    mpz_mat_fprint(stdout,A);
    for (k = 2 ; k <= n ; k++){
        mpz_mat_multiply(D,C,A);
        mpz_mat_swap(C,D);
    }
    mpz_mat_realloc(B,A->m,A->n);
    mpz_mat_swap(C,B);
    mpz_mat_clear(C);
    mpz_mat_clear(D);
}

// Returns A^n mod pfor n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_power_ui_mod_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned int n, unsigned int p){
    unsigned int k;
    mpz_mat C,D;
    ASSERT_ALWAYS(A->n == A->m);
    mpz_mat_init(C,A->m,A->n);
    mpz_mat_init(D,A->m,A->n);
    mpz_mat_set_ui(C,1);
    for (k = 1 ; k <= n ; k++){
        mpz_mat_multiply_mod_ui(D,C,A,p);
        mpz_mat_swap(C,D);
    }
    mpz_mat_realloc(B,A->m,A->n);
    mpz_mat_swap(C,B);
    mpz_mat_clear(C);
    mpz_mat_clear(D);
}


void mpq_mat_multiply(mpq_mat_ptr D, mpq_mat_srcptr A, mpq_mat_srcptr B)
{
	ASSERT_ALWAYS(A->n == B->m);
    mpq_t z;
    mpq_init(z);
    mpq_mat C;
    mpq_mat_init(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_set_ui(mpq_mat_entry(C, i, j), 0, 1);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpq_mul(z,
                        mpq_mat_entry_const(A, i, k),
                        mpq_mat_entry_const(B, k, j));
                mpq_add(mpq_mat_entry(C, i, j),
                        mpq_mat_entry(C, i, j),
                        z);
            }
        }
    }
    mpq_mat_realloc(D,A->m,B->n);
    mpq_mat_set(D,C);
    mpq_mat_clear(C);
    mpq_clear(z);
}

void mpq_mat_multiply_by_mpq(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_ptr k)
{
    mpq_mat B;
	mpq_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpq_mul(mpq_mat_entry(B, i, j),mpq_mat_entry_const(A,i,j),k);
        }
    }
    mpq_mat_realloc(C, A->m, A->n);
    mpq_mat_set(C,B);
    mpq_mat_clear(B);
}

void mpq_mat_multiply_by_si(mpq_mat_ptr C, mpq_mat_srcptr A, int k)
{
    mpq_mat B;
	mpq_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpz_mul_si(mpq_numref(mpq_mat_entry(B,i,j)),mpq_numref(mpq_mat_entry_const(A,i,j)),k);
			mpz_set(mpq_denref(mpq_mat_entry(B,i,j)),mpq_denref(mpq_mat_entry_const(A,i,j)));
			mpq_canonicalize(mpq_mat_entry(B,i,j));
        }
    }
    mpq_mat_realloc(C, A->m, A->n);
    mpq_mat_set(C,B);
    mpq_mat_clear(B);
}

void mpq_mat_multiply_by_ui(mpq_mat_ptr C, mpq_mat_srcptr A, unsigned int k)
{
	mpq_mat B;
	mpq_mat_init(B,A->m,A->n);
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
			mpz_mul_ui(mpq_numref(mpq_mat_entry(B,i,j)),mpq_numref(mpq_mat_entry_const(A,i,j)),k);
			mpz_set(mpq_denref(mpq_mat_entry(B,i,j)),mpq_denref(mpq_mat_entry_const(A,i,j)));
			mpq_canonicalize(mpq_mat_entry(B,i,j));
        }
    }
    mpq_mat_realloc(C, A->m, A->n);
    mpq_mat_set(C,B);
    mpq_mat_clear(B);
}

/*}}}*/
/* {{{ gaussian reduction over the rationals
 * this is a backend for row gaussian reduction. T receives the
 * transformation matrix. M is modified in place.
 * this includes reordering of the rows.
 *
 * this works reasonably well over Fp, but for rationals we suffer from
 * coefficient blowup.
 */
void mpq_gauss_backend(mpq_mat_ptr M, mpq_mat_ptr T)
{
    unsigned int m = M->m;
    unsigned int n = M->n;
    mpq_mat_realloc(T, m, m);
    mpq_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpq_t tmp1, tmp2;
    mpq_init(tmp1);
    mpq_init(tmp2);
    for(unsigned int j = 0 ; j < n ; j++) {
        unsigned int i1;
        for(i1 = rank ; i1 < M->m ; i1++) {
            if (mpq_cmp_ui(mpq_mat_entry(M, i1, j), 0, 1) != 0) break;
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        mpq_mat_swaprows(M, rank, i1);
        mpq_mat_swaprows(T, rank, i1);
        i1 = rank;
        /* canonicalize this row */
        mpq_inv(tmp1, mpq_mat_entry(M, rank, j));
        mpq_mat_mulrow(M, rank, tmp1);
        mpq_mat_mulrow(T, rank, tmp1);
        for(unsigned int i0 = 0 ; i0 < m ; i0++) {
            if (i0 == rank) continue;
            /* we've normalized row rank, so all it takes is to multiply by
             * M[i0, j]
             */
            mpq_set(tmp2, mpq_mat_entry(M, i0, j));
            mpq_mat_submulrow(M, i0, rank, tmp2);
            mpq_mat_submulrow(T, i0, rank, tmp2);
        }
        rank++;
    }
    mpq_clear(tmp1);
    mpq_clear(tmp2);
}
/* }}} */
/* {{{ gaussian reduction over Z/pZ
 *
 * We sort of return something even in the case where p is not prime.
 */
void mpz_gauss_backend_mod(mpz_mat_ptr M, mpz_mat_ptr T, mpz_srcptr p)
{
    unsigned int m = M->m;
    unsigned int n = M->n;
    mpz_mat_realloc(T, m, m);
    mpz_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpz_t gcd, tmp2, tmp3;
    mpz_init(gcd);
    mpz_init(tmp2);
    mpz_init(tmp3);
    for(unsigned int j = 0 ; j < n ; j++) {
        unsigned int i1 = M->m;
        mpz_set(gcd, p);
        for(unsigned int i = rank ; i < M->m ; i++) {
            mpz_gcd(tmp2, mpz_mat_entry(M, i, j), p);
            if (mpz_cmp(tmp2, gcd) < 0) {
                i1 = i;
                mpz_set(gcd, tmp2);
                if (mpz_cmp_ui(gcd, 1) == 0) break;
            }
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        mpz_mat_swaprows(M, rank, i1);
        mpz_mat_swaprows(T, rank, i1);
        i1 = rank;
        /* canonicalize this row */
        /* gcd is gcd(M[rank, j], p) */
        mpz_divexact(tmp2, mpz_mat_entry(M, rank, j), gcd);
        mpz_divexact(tmp3, p, gcd);
        mpz_invert(tmp2, tmp2, tmp3);
        mpz_mat_mulrow_mod(M, rank, tmp2, p);
        mpz_mat_mulrow_mod(T, rank, tmp2, p);
        ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry(M, rank, j), gcd) == 0);
        for(unsigned int i0 = 0 ; i0 < m ; i0++) {
            if (i0 == rank) continue;
            /* we've normalized row rank, so all it takes is to multiply by
             * M[i0, j]
             */
            mpz_fdiv_q(tmp2, mpz_mat_entry(M, i0, j), gcd);
            /* for i0 > rank, we should have exact quotient */
            mpz_mat_submulrow_mod(M, i0, rank, tmp2, p);
            mpz_mat_submulrow_mod(T, i0, rank, tmp2, p);
            ASSERT_ALWAYS(mpz_cmp_ui(mpz_mat_entry(M, i0, j), 0) == 0);
        }
        rank++;
    }
    mpz_clear(gcd);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
}
/* }}} */
/* {{{ integer heap routines */
/* input: heap[0..n-1[ a heap, and heap[n-1] unused.
 * output: heap[0..n[ a heap, includes x
 */
typedef int (*cmp_t)(void *, unsigned int, unsigned int);

void push_heap(unsigned int * heap, unsigned int x, unsigned int n, cmp_t cmp, void * arg)
{
    unsigned int i = n - 1;
    for(unsigned int i0; i; i = i0) {
        i0 = (i-1) / 2;
        if (cmp(arg, x, heap[i0]) > 0)
            heap[i] = heap[i0];
        else
            break;
    }
    heap[i] = x;
}

/* put heap[0] in x, and keep the heap property on the resulting
 * heap[0..n-1[ ; heap[n-1] is unused on output.
 */
void pop_heap(unsigned int * x, unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    *x = heap[0];
    unsigned int y = heap[n-1];
    /* make heap[0..n-1[ a heap, now */
    /* couldn't there be various strategies here ? */
    unsigned int i = 0;
    for(unsigned int i1 ; 2*i+2 <  n - 1 ; i = i1) {
        unsigned int left = 2*i + 1;
        unsigned int right = 2*i + 2;
        if (cmp(arg, heap[left], heap[right]) > 0)
            i1 = left;
        else
            i1 = right;
        if (cmp(arg, y, heap[i1]) > 0)
            break;
        heap[i] = heap[i1];
    }
    /* there we know that i has no right child, maybe even no left child.
     */
    if (2*i + 1 < n - 1 && cmp(arg, y, heap[2*i+1]) <= 0) {
        heap[i] = heap[2*i+1];
        i=2*i+1;
    }
    heap[i] = y;
}

/* make heap[0..n[ a heap. */
void make_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = 1 ; i <= n ; i++) {
        unsigned int t = heap[i-1];
        push_heap(heap, t, i, cmp, arg);
    }
}

/* must be a heap already */
void sort_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = n ; i ; i--) {
        unsigned int t;
        pop_heap(&t, heap, i, cmp, arg);
        heap[i-1] = t;
    }
}

void sort_heap_inverse(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    sort_heap(heap, n, cmp, arg);
    for(unsigned int i = 0, j = n-1 ; i < j ; i++, j--) {
        int a = heap[i];
        int b = heap[j];
        heap[i] = b;
        heap[j] = a;
    }
}
/* }}} */
/* {{{ Hermite normal form over the integers
 * The algorithm here is randomly cooked from a stupid idea I had. The
 * transformation matrices produced seem to be reasonably good, somewhat
 * better than what magma gives (for overdetermined matrices, of course
 * -- otherwise they're unique). However we do so in an embarrassingly
 * larger amount of time (more than 1000* for 100x100 matrices).
 * Fortunately, our intent is not to run this code on anything larger
 * than 30*30, in which case it's sufficient.
 */
/*{{{ hnf helper algorithm (for column matrices) */
struct mpz_hnf_helper_heap_aux {
    double * dd;
    mpz_mat_ptr A;
};

int mpz_hnf_helper_heap_aux_cmp(struct mpz_hnf_helper_heap_aux * data, unsigned int a, unsigned int b)
{
    /* The heap will have the largest elements on top according to this
     * measure. Thus the comparison on the ->a field must be as per
     * mpz_cmp
     */
    mpz_srcptr xa = mpz_mat_entry_const(data->A, a, 0);
    mpz_srcptr xb = mpz_mat_entry_const(data->A, b, 0);

    assert(mpz_sgn(xa) >= 0);
    assert(mpz_sgn(xb) >= 0);
    int r = mpz_cmp(xa, xb);
    if (r) return r;
    /* among combinations giving the same result, we want the "best" to
     * be the one with *smallest* norm. So instead of the sign of da-db,
     * we want the sign of db-da, here.
     */
    double da = data->dd[a];
    double db = data->dd[b];
    return (db > da) - (da > db);
}

/* this is quite the same as computing the hnf on a column matrix, with
 * the added functionality that entries in a[0..n0[ are reduced with
 * respect to the entries in a[n0..n[. This is equivalent to saying that
 * we apply a transformation matrix which is:
 *       I_n0 R
 * T ==  0    S
 * where S is unimodular. The resulting vector T*A has entries [n0..n[
 * equal to (g,0,...,0), while entries [0..n0[ are reduced modulo g.
 *
 * return +1 or -1 which is the determinant of T.
 */
int mpz_hnf_helper(mpz_mat_ptr T, mpz_mat_ptr dT, mpz_ptr * a, unsigned int n0, unsigned int n)
{
    int signdet=1;
    mpz_mat_realloc(dT, n, n);
    mpz_mat_set_ui(dT, 1);

    if (n == n0) return signdet;

    mpz_t q, r2;
    mpz_mat A;
    mpz_mat_init(A, n, 1);
    mpz_init(q);
    mpz_init(r2);
    unsigned int * heap = malloc(n * sizeof(unsigned int));
    struct mpz_hnf_helper_heap_aux cmpdata[1];
    cmpdata->A = A;
    cmpdata->dd = malloc(n * sizeof(double));
    cmp_t cmp = (cmp_t) &mpz_hnf_helper_heap_aux_cmp;
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_ptr Ai = mpz_mat_entry(A, i, 0);
        mpz_set(Ai, a[i]);
        cmpdata->dd[i] = 0;
        for(unsigned int j = 0 ; j < n ; j++) {
            double d = mpz_get_d(mpz_mat_entry(T, i, j));
            cmpdata->dd[i] += d * d;
        }
        heap[i] = i;
        if (i < n0)
            continue;
        if (mpz_sgn(Ai) == -1) {
            mpz_neg(Ai, Ai);
            signdet = -signdet;
            mpz_set_si(mpz_mat_entry(dT, i, i), -1);
            for(unsigned int j = 0 ; j < n ; j++) {
                mpz_neg(mpz_mat_entry(T, i, j), mpz_mat_entry(T, i, j));
            }
        }
    }
    make_heap(heap + n0, n - n0, cmp, cmpdata);
    unsigned int head;
    for( ; ; ) {
        /* extract the head */
        pop_heap(&head, heap + n0, n - n0, cmp, cmpdata);
        mpz_ptr r0 = mpz_mat_entry(A, head, 0);
        if (mpz_cmp_ui(r0, 0) == 0) {
            /* we're done ! */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        /* reduce A[0..n0[ wrt r0 == A[head] */
        for(unsigned int i = 0 ; i < n0 ; i++) {
            mpz_fdiv_qr(q, r2, mpz_mat_entry(A, i, 0), r0);
            mpz_swap(r2, mpz_mat_entry(A, i, 0));
            mpz_mat_submulrow(T, i, head, q);
            mpz_mat_submulrow(dT, i, head, q);
        }
        if (n0 == n - 1) {
            /* we're actually computing the gcd of one single integer.
             * That makes no sense of course.
             */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        mpz_ptr r1 = mpz_mat_entry(A, heap[n0], 0);
        if (mpz_cmp_ui(r1, 0) == 0) {
            /* we're done ! */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        mpz_fdiv_qr(q, r2, r0, r1);
        mpz_swap(r0, r2);
        mpz_mat_submulrow(T, head, heap[n0], q);
        mpz_mat_submulrow(dT, head, heap[n0], q);
        /* adjust the norm of the row in T */
        cmpdata->dd[head] = 0;
        for(unsigned int j = 0 ; j < n ; j++) {
            double d = mpz_get_d(mpz_mat_entry(T, head, j));
            cmpdata->dd[head] += d * d;
        }
        push_heap(heap + n0, head, n - n0, cmp, cmpdata);
    }
    sort_heap_inverse(heap + n0, n - n0, cmp, cmpdata);
    mpz_mat_permuterows(T, heap);
    mpz_mat_permuterows(dT, heap);
    mpz_mat_permuterows(A, heap);
    signdet*= permutation_signature(heap, n);
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_ptr Ai = mpz_mat_entry(A, i, 0);
        mpz_swap(a[i], Ai);
    }
    free(cmpdata->dd);
    mpz_mat_clear(A);
    mpz_clear(q);
    mpz_clear(r2);
    free(heap);

    return signdet;
}
/* {{{ this computes the row hnf on a column C.
 * The result is always of the form C'=(gcd(C),0,0,...,0). The transform
 * matrix T such that C'=T*C is most interesting.  a is modified in
 * place.
 */
void mpz_gcd_many(mpz_mat_ptr dT, mpz_ptr * a, unsigned int n)
{
    mpz_mat T;
    mpz_mat_init(T, n, n);
    mpz_mat_set_ui(T, 1);
    mpz_hnf_helper(T, dT, a, 0, n);
    mpz_mat_clear(T);
}
/*}}}*/
/*}}}*/

/* return +1 or -1, which is the determinant of the transformation matrix
 * T */
int mpz_hnf_backend(mpz_mat_ptr M, mpz_mat_ptr T)
{
    /* Warning ; for this to work, we need to have extra blank rows in M
     * and T, so we require on input M and T to be overallocated (twice
     * as many rows).
     */
    int signdet = 1;
    unsigned int m = M->m;
    unsigned int n = M->n;
    mpz_mat_realloc(T, m, m);
    mpz_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpz_t quo;
    mpz_init(quo);
    mpz_ptr * colm = malloc(m * sizeof(mpz_t));
    mpz_mat dT, Mx, My;
    mpz_mat_init(dT, 0, 0);
    mpz_mat_init(Mx, 0, 0);
    mpz_mat_init(My, 0, 0);
    for(unsigned int j = 0 ; j < n && rank < m; j++) {
        for(unsigned int i = 0 ; i < m ; i++) {
            colm[i] = mpz_mat_entry(M, i, j);
        }
        signdet *= mpz_hnf_helper(T, dT, colm, rank, m);

        /* apply dT to the submatrix of size m * (n - 1 - j) of M */
        mpz_mat_realloc(Mx, m, n - 1 - j);
        mpz_mat_submat_swap(Mx, 0, 0, M, 0, j + 1, m, n - 1 - j);
        mpz_mat_multiply(My, dT, Mx);
        mpz_mat_submat_swap(M, 0, j + 1, My, 0, 0, m, n - 1 - j);

        mpz_srcptr gcd = colm[rank];
        if (mpz_cmp_ui(gcd, 0) == 0) /* this column is zero -- nothing to do */
            continue;
        rank++;
    }
    free(colm);
    mpz_clear(quo);
    mpz_mat_clear(dT);
    mpz_mat_clear(Mx);
    mpz_mat_clear(My);

    return signdet;
}
/* }}} */

/*{{{ kernel*/
// This is supposed to compute the Kernel of M mod p and to store it in the matrix K. If r is the rank of M, and M is a square matrix n*n, K is a n*(n-r) matrix
void mpz_mat_kernel(mpz_mat_ptr K, mpz_mat_srcptr M, unsigned int p)
{
    mpz_mat T, H;
    mpz_t p2;
    unsigned int r;

    mpz_mat_init(T,M->m,M->m);
    mpz_mat_init(H,M->m,M->n);
    mpz_init_set_ui(p2,p);
    mpz_mat_set(H,M);

    // Storing a reduced matrix of M in H with gauss
    mpz_gauss_backend_mod(H,T,p2);

    // Finding the rank of M and H
    r = H->m;
    while((r>0) && mpz_mat_isnull_row(H,r-1)){
        r--;
    }
    // Kernel is of dimension n-r, and a basis of the kernel is given in the n-r last rows of T
    // We shall keep the convention of magma
    
    // Reallocating K with n-r columns
    mpz_mat_realloc(K, H->m-r, H->n);
    mpz_mat_submat_swap(    K, 0, 0,
                            T, r, 0, 
                            M->m - r, M->m);

    mpz_clear(p2);
    mpz_mat_clear(T);
    mpz_mat_clear(H);
    
}
/* }}} */
/*{{{ inversion*/
void mpq_mat_invert(mpq_mat_ptr dst, mpq_mat_ptr src)
{
	// Inverting B
	mpq_mat inv, id, aux;
	ASSERT_ALWAYS(src->m == src->m);
	mpq_mat_init(inv,src->m,2*src->m);
	mpq_mat_init(id,src->m,src->n);
	mpq_mat_init(aux,src->m,src->n);
	
	mpq_mat_set_ui(id,1);
    mpq_mat_set(aux,src);
    
    mpq_mat_submat_swap(aux,0,0,inv,0,0,src->m,src->n);
    mpq_mat_submat_swap(id,0,0,inv,0,src->n,src->m,src->n);
    
        // B2 now contains B and identity
    mpq_gauss_backend(inv,aux);
        // B2 now contains identity and B^-1
    mpq_mat_submat_swap(dst,0,0,inv,0,src->m,src->m,src->m);
    
    mpq_mat_clear(inv);
    mpq_mat_clear(id);
    mpq_mat_clear(aux);
}
/* }}} */
/*{{{ timer*/
double
seconds (void)
{
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    uint64_t r;
    r = (uint64_t) res->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec;
    return r * 1.0e-6;
}
/*}}}*/
/*{{{ conversion of rows and columns to polynomials*/
void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, unsigned int i)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n);
	unsigned int j;
	for (j = 0 ; j < M->n; j++){
		mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,j));
	}
}

void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, unsigned int j)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->m);
	unsigned int i;
	for (i = 0 ; i < M->m; i++){
		mpz_poly_setcoeff(f,i,mpz_mat_entry_const(M,i,j));
	}
}

void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, unsigned int i)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n);
	mpz_set_si(denom,1);
	unsigned int j;
    for (j = 0 ; j < M->n ; j++) {
        mpz_lcm(denom,denom,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
	for (j = 0 ; j < M->n; j++){
		mpq_t aux;
		mpz_t num;
		mpq_init(aux);
		mpz_init(num);
		mpq_set_z(aux,denom);
		mpq_mul(aux,aux,mpq_mat_entry_const(M,i,j));
		mpz_set(num,mpq_numref(aux));
		mpz_poly_setcoeff(f,j,num);
		mpq_clear(aux);
		mpz_init(num);
	}
}

void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, unsigned int j)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->m);
	mpz_set_si(denom,1);
	unsigned int i;
    for (i = 0 ; i < M->m ; i++) {
        mpz_lcm(denom,denom,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
	for (i = 0 ; i < M->m; i++){
		mpq_t aux;
		mpz_t num;
		mpq_init(aux);
		mpz_init(num);
		mpq_set_z(aux,denom);
		mpq_mul(aux,aux,mpq_mat_entry_const(M,i,j));
		mpz_set(num,mpq_numref(aux));
		mpz_poly_setcoeff(f,i,num);
		mpq_clear(aux);
		mpz_init(num);
	}
}
/*}}}*/


// Prints the polynomial
void print_polynomial(mpz_t* f, int degree)
{
	int i;
	for(i = degree ; i >= 0 ; i--){
		if(mpz_get_ui(f[i]) != 0){
			gmp_printf("%Zd",f[i]);
			if(i > 0){
				printf("*x^%d + ",i);
			}
		}
	}
	printf("\n");
}

// Takes a matrix B containing the generators (w_0, ... w_{n-1}) of an order, and returns the matrix U containing ((w_0)^p, ..., (w_{n-1})^p), reduced mod g and mod p.
void generators_to_power_p(mpq_mat_ptr U, mpq_mat_srcptr B, mpz_poly_ptr g, unsigned int p)
{
    ASSERT_ALWAYS (B->m == B->n) ;

    unsigned int i, j;
    int k;
    mpq_mat_realloc(U,B->m,B->n);

    mpz_t K, tmp, p2; // common denominator
    mpz_init(K);
    mpz_init(tmp);
    mpz_set_si(K,1);
    mpz_poly_t f2;
    mpz_poly_t g2;
    mpz_poly_init(f2,B->n-1);
    mpz_poly_init(g2,B->n-1);
    mpz_init(p2);
    mpz_set_ui(p2,p);

    mpz_t N[B->n]; // will be used to store numerators of w[j] times K 
    for (i = 0 ; i < B->m ; i++){
        mpz_init(N[i]);
    }

    mpz_poly_mod_mpz(g2,g,p2,NULL);

    for (j = 0 ; j < U->m ; j++) {
        mpz_poly_t f;
        mpz_poly_init(f,B->n-1);

        // Putting the LCM of all denominators of coefficients of w[j] in K
        mpz_set_si(K,1);
        for (i = 0 ; i < U-> n ; i++) {
            mpz_lcm(K,K,mpq_denref(mpq_mat_entry_const(B,i,j)));
        }


        // Generating the polynomial
        for (i = 0 ; i < B->n ; i++) {
            mpq_t aux1;
            mpq_t K_rat;
            mpq_init(aux1);
            mpq_init(K_rat);

            mpq_set(aux1,mpq_mat_entry_const(B,i,j));
            mpq_set_z(K_rat,K);
            mpq_mul(aux1,aux1,K_rat);
            mpq_get_num(N[i],aux1);

            mpq_clear(K_rat);
            mpq_clear(aux1);
        }
        mpz_poly_setcoeffs(f,N,B->n-1);
        mpz_poly_set(f2,f);

        // Computing f^p mod g (it returns (K * the corresponding generator)^p
        mpz_poly_power_mod_f(f,f,g,p);

        mpz_pow_ui(K, K, p);

        // Storing w[j] in j-th row of U
        for (k = 0 ; k <= f->deg ; k++) {
            mpz_poly_getcoeff(tmp,k,f);
            mpq_set_num(mpq_mat_entry(U,j,k),tmp);
            mpq_set_den(mpq_mat_entry(U,j,k),K);
            mpq_canonicalize(mpq_mat_entry(U,j,k));
        } 

        mpz_poly_clear(f);
    }

    for (i = 0 ; i < B->m ; i++){
        mpz_clear(N[i]);
    }
    mpz_poly_clear(f2);
    mpz_clear(tmp);
    mpz_clear(p2);
    mpz_poly_clear(g2);
    mpz_clear(K);
}

// Builds the HNF of the matrix containing the vectors of the kernel of the application z -> z^p, and the identity*p
void generators_of_Ip(mpz_mat_ptr I, mpz_mat_srcptr K, unsigned int p){
	mpz_mat J, T0;
	mpz_mat_init(J,K->n,K->n);
	mpz_mat_init(T0,K->n,K->n);
	mpz_mat_realloc(I,K->n,K->n);

	mpz_mat_set_ui(J,1);
	mpz_mat_multiply_by_ui(J,J,p);
	mpz_mat_vertical_join(I,J,K);
	mpz_hnf_backend(I,T0);
	mpz_mat_submat_swap(I,0,0,J,0,0,K->n,K->n);
	mpz_mat_set(I,J);
	
	mpz_mat_clear(T0);
	mpz_mat_clear(J);
}

int main(int argc, char * argv[])/*{{{*/
{
	/*
    unsigned int m = 8;
    unsigned int n = 5;
    if (argc == 3) {
        m = strtoul(argv[1], NULL, 0);
        n = strtoul(argv[2], NULL, 0);
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);

    mpq_mat M;
    mpq_mat T;
    mpz_mat Mz;
    mpz_mat Tz;
    mpz_t p;

    mpq_mat_init(M, m, n);
    mpq_mat_init(T, m, m);

    mpz_mat_init(Mz, m, n);
    mpz_mat_init(Tz, m, m);

    mpz_init(p);

    mpz_set_ui(p, 19);

    if (0) {
	printf("\n\nCas 0.1\n\n");
        mpq_mat_urandomm(M, state, p);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_gauss_backend(M, T);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_mat_fprint(stdout, T);
        printf("\n");
    }

    if (0) {
	printf("\n\nCas 0.2\n\n");
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_gauss_backend_mod(Mz, Tz, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_mat_fprint(stdout, Tz);
        printf("\n");
    }

    if (1) {
	printf("\n\nCas 1\n\n");
        mpz_mat_realloc(Mz, m, n);
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz); printf("\n");
        double t = seconds();
        mpz_hnf_backend(Mz, Tz);
        t = seconds()-t;
        mpz_mat_fprint(stdout, Mz); printf("\n");
        mpz_mat_fprint(stdout, Tz); printf("\n");

        printf("%1.4f\n", t);
    }

    mpz_clear(p);
    mpq_mat_clear(M);
    mpq_mat_clear(T);
    mpz_mat_clear(Mz);
    mpz_mat_clear(Tz);
    gmp_randclear(state);
	*/


	// Here starts my personal work
	// We assume that the polynomial was given in a command line, by giving its coefficient, in reversed order (first the constant coefficient, etc. And the head coefficient in the end);

	/*
	if(argc > 1){

		// Initialisation
		unsigned int degree = argc-2;
		mpz_t f[degree+1];
		mpz_mat mul_alpha, M, N, D, D2;
		mpz_t p;
		mpz_t minus;

		// Storing the coefficients obtained in the command line in the polynomial
		unsigned int i,j, k;
		for(i = 0 ; i <= degree ; i++){
			mpz_init(f[i]);
			mpz_set_str(f[i],argv[i+1],10);
		}
		printf("\nYour polynomial is :\n");
		print_polynomial(f,degree);

		//Initialising each matrix
		mpz_mat_init(mul_alpha,degree,degree);
		mpz_mat_init(M,degree,degree);
		mpz_mat_init(N,degree,degree);
		mpz_mat_init(D,degree,degree);
		mpz_mat_init(D2,2*degree,degree);
		mpz_init(p);
		mpz_init(minus);
		mpz_mat_set_ui(M,1); // M starts as the identity Matrix
		mpz_set_si(minus,-1);
		

		// Filling the coefficients in mul_alpha, the companion matrix of f
		for(i = 0 ; i < degree ; i++){
			// Setting the left part of the companion matrix
			for(j = 0 ; j < degree-1 ; j++){
				if(i == j+1){ mpz_set_ui(mpz_mat_entry(mul_alpha,i,j),1); }
				else{ mpz_set_ui(mpz_mat_entry(mul_alpha,i,j),0); }
			}

			// Computing the coefficients for the column on the right
			mpz_set(p,f[i]);
			mpz_mul(p,p,minus);
			for(k = 1 ; k <= degree-1-i ; k++){
				mpz_mul(p,p,f[degree]);
			}
			mpz_set(mpz_mat_entry(mul_alpha,i,j),p);
		}



		// Filling the coefficients in D, whose determinant must be computed to get the discriminant.
		for(i = 0 ; i < degree ; i++){
			for(j = 0 ; j <= i ; j++){
                            mpz_mat_trace(mpz_mat_entry(D, i-j, j), M);
			}
			mpz_mat_multiply(N,M,mul_alpha);
			mpz_mat_swap(M,N);
		}
		for(j = 1 ; j < degree ; j++){
			for(i = degree-1 ; j <= i ; i--){
				mpz_mat_trace(mpz_mat_entry(D,i,j+(degree-1)-i),M);
			}
			mpz_mat_multiply(N,M,mul_alpha);
			mpz_mat_swap(M,N);
		}


		// Preparing the HNF
		mpz_mat_realloc(M,degree,degree);
                int sign = mpz_hnf_backend(D, M);

                mpz_mat_determinant_triangular(p, D);
                mpz_mul_si(p,p,sign);
		gmp_printf("\nThe discriminant of Z[f_d * alpha] is %Zd\n\n",p);

		mpz_clear(p);
		mpz_clear(minus);
		mpz_mat_clear(mul_alpha);
		mpz_mat_clear(M);
		mpz_mat_clear(N);
		mpz_mat_clear(D);
		mpz_mat_clear(D2);
		for(i = 0 ; i <= degree ; i++){
			mpz_clear(f[i]);
                }
	
	}*/

    // The inputs to this problem are f, one polynomial of degree n, and B, the matrix containing the genereators of one order of the number field obtained with f, as well as p, a prime number

    if (argc != 3) {
        fprintf(stderr, "Usage: ./a.out [filename] [p]\n");
        exit(EXIT_FAILURE);
    }

    unsigned int p = strtoul(argv[2],NULL,0); //19; //atoi(argv[1]);
    FILE * problemfile = fopen(argv[1], "r");

    mpq_mat B, B_inv, T, U;
    mpz_poly_t f, g;
    mpz_mat X, K, I, T0, M;
    mpz_t den;

    printf("Format: [degree] [coeffs] [coeffs of order basis]\n");

    unsigned int n;

    fscanf(problemfile, "%u", &n);

    mpz_poly_init(f,n);
    mpz_poly_init(g,n);
    mpz_mat_init(X,n,n);
    mpz_mat_init(T0,n,n);
    mpz_mat_init(K,n,n);
    mpz_mat_init(M,n,n*n);
    mpz_mat_init(I,n,n);
    mpq_mat_init(B, n, n); // The matrix of generators
    mpq_mat_init(B_inv, n, n); // Its inverse
    mpq_mat_init(T, n, n); // An auxiliary matrix
    mpq_mat_init(U, n, n);
    mpq_mat_set_ui(T, 1);
    mpz_init(den);








    // Filling in f as an example, and storing in g the monic polynomial such as (fd*alpha) is a root of if and only if alpha is a root of f

    for(unsigned int i = 0 ; i <= n ; i++) {
        int c;
        fscanf(problemfile, "%d", &c);
        mpz_poly_setcoeff_si(f,i,c);
    }
    mpz_poly_to_monic(g,f);
    printf("f is : "); mpz_poly_fprintf(stdout,f); printf("\n");
    printf("f^ is : "); mpz_poly_fprintf(stdout,g); printf("\n");

    // Filling in B here as an example ; genereators of the maximal order of the number field of g (tested with magma);
    for(unsigned int i = 0 ; i < n ; i++) {
        int denom;
        fscanf(problemfile, "%d", &denom);
        for(unsigned int j = 0 ; j < n ; j++) {
            int c;
            fscanf(problemfile, "%d", &c);
            mpq_set_si(mpq_mat_entry(B,i,j),c,denom);
        }
    }

    fclose(problemfile);

    printf("generators are :\n"); mpq_mat_fprint(stdout,B); printf("\n");


    // Inverting B
    mpq_mat_invert(B_inv,B);




    // Now building the matrix U, containing all generators to the power of p
    // Generators are polynomials, stored in the matrix B
    generators_to_power_p(U,B,g,p);
    

    // Now multiplying B^-1 and U;
    mpq_mat_multiply(T,B_inv,U);
    mpq_mat_numden(X,den,T);
    mpz_mat_mod_ui(X,X,p);

    printf("Matrix of the application z -> (z^%d mod f mod %d) :\n",p,p); mpz_mat_fprint(stdout,X); printf("\n");

    /* Which power of p ? */
    int k = 1;
    for(unsigned int pk = p ; pk < n ; k++, pk *= p) ;
    mpz_mat_power_ui_mod_ui(X,X,k,p);

    printf("Matrix to the power of %d of the application z -> (z^%d mod f mod %d) :\n",k,p,p); mpz_mat_fprint(stdout,X); printf("\n");


    mpz_mat_kernel(K,X,p);
    printf("F is the application z -> (z^%d mod f mod %d) :\n", p, p);
    printf("A basis of F^%d is :\n",k);
    mpz_mat_fprint(stdout,K); 
    
    
    
    
    
    
    
    
    
    printf("\n\n****** Generating a new order ******\n");
	
	// Getting generators of I_p
	generators_of_Ip(I,K,p);
	printf("Generators of I_p in the basis of the given order O :\n");
	mpz_mat_fprint(stdout,I); printf("\n");
	
	mpq_mat I_rat, I_inv;
	mpq_mat_init(I_rat,n,n);
	mpq_mat_init(I_inv,n,n);
	mpz_mat_to_mpq_mat(I_rat,I);
	mpq_mat_invert(I_inv,I_rat);
	unsigned int i,j;
	
	for (i = 0 ; i < n ; i++){
		for (j = 0 ; j < n ; j++){
			mpz_poly_t gamma, c, aux;
			mpz_t denom_g, denom_c, coeff;
			mpz_mat row, aux_row;
			mpq_mat row_q, c_mat, res;
			mpz_mat_init(row,1,n);
			mpz_mat_init(aux_row,1,n);
			mpq_mat_init(row_q,1,n);
			mpq_mat_init(c_mat,1,n);
			mpq_mat_init(res,n,n);
			mpz_poly_init(gamma,n);
			mpz_poly_init(c,n);
			mpz_poly_init(aux,n);
			mpz_init(denom_g);
			mpz_init(denom_c);
			mpz_init(coeff);
			

			
			 
			printf("i = %d ; j = %d\n",i,j);
			// Storing one generator of B in the polynomial gamma and in denom_g
			mpq_mat_row_to_poly(gamma,denom_g,B,i);
			//Storing one generator of I_p in the polyomial c and in denom_c
			mpz_mat_submat_swap(aux_row,0,0,I,j,0,1,n);
			mpz_mat_set(row,aux_row);
			mpz_mat_submat_swap(aux_row,0,0,I,j,0,1,n);
			mpz_mat_to_mpq_mat(row_q,row);
			mpq_mat_multiply(c_mat,row_q,B);
			mpq_mat_row_to_poly(c,denom_c,c_mat,0);
			
			
			// Computing gamma*c mod g
			mpz_poly_mul_mod_f(aux,gamma,c,g);
			// Storing the result in row_q
			for(int k = 0 ; k <= aux->deg ; k++){
				mpz_poly_getcoeff(coeff,k,aux);
				mpq_set_num(mpq_mat_entry(row_q,0,k),coeff);
				mpz_mul(coeff,denom_c,denom_g);
				mpq_set_den(mpq_mat_entry(row_q,0,k),coeff);
				mpq_canonicalize(mpq_mat_entry(row_q,0,k));
			}
			
			
			//printf("gamma*c in the basis of Ip :\n");
			//mpq_mat_fprint(stdout,row_q);
			
			// Converting row_q (gamma*c mod g) in the basis of I_p (it is supposed to contain only integers)
			mpq_mat_multiply(res,row_q,B_inv);
			mpq_mat_multiply(res,res,I_inv);
			
			mpq_mat_fprint(stdout,res);
			/*
			// Extracting the numerators (remember, integers only) into a mpz_mat
			mpq_mat_numden(row,coeff,res);
			// Computing the same matrix, modulo p (it is supposed to be a vector of n integers, associated to gamma)
			mpz_mat_mod_ui(row,row,p);
			printf("The same, mod %d :\n",p);
			mpz_mat_fprint(stdout,row);
			//mpq_mat_fprint(stdout,res);
			printf("\n\n");
			
			*/
			
			mpq_mat_clear(res);
			mpq_mat_clear(c_mat);
			mpq_mat_clear(row_q);
			mpz_mat_clear(row);
			mpz_mat_clear(aux_row);
			mpz_clear(coeff);
			mpz_clear(denom_g);
			mpz_clear(denom_c);
			mpz_poly_clear(gamma);
			mpz_poly_clear(c);
			mpz_poly_clear(aux);
		}
	}


	mpq_mat_clear(I_rat);
	mpq_mat_clear(I_inv);
	mpz_mat_clear(M);
	mpz_mat_clear(T0);
    mpz_mat_clear(K);
    mpz_mat_clear(X);
    mpz_mat_clear(I);
    mpq_mat_clear(B);
    mpq_mat_clear(B_inv);
    mpq_mat_clear(T);
    mpq_mat_clear(U);
    mpz_clear(den);
    

}
/*}}}*/
