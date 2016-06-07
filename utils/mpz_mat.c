#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "mpz_mat.h"
#include "portability.h"

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
// We assume that M is square
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M)
{
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(t, 0);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_add(t, t, mpz_mat_entry_const(M,i,i));
}

// We assume that M is triangular
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M)
{
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(d, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_mul(d, d, mpz_mat_entry_const(M,i,i));
}

// We assume that M is square
void mpq_mat_trace(mpq_ptr t, mpq_mat_srcptr M)
{
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(t, 0, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpq_add(t, t, mpq_mat_entry_const(M,i,i));
}

// We assume that M is triangular
void mpq_mat_determinant_triangular(mpq_ptr d, mpq_mat_srcptr M)
{
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(d, 1, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpq_mul(d, d, mpq_mat_entry_const(M,i,i));
}
/*}}}*/
/*{{{ miscellaneous */

/* convert to integer matrix divided by lcm of denominator.
 * return 1
 * if den==NULL, return 0 if denominator happens to not be 1 (in which
 * case the matrix return is undefined).
 */

int mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M)
{
    int ret = 1;
    if (den) mpz_set_ui(den, 1);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_srcptr Mij =  mpq_denref(mpq_mat_entry_const(M, i, j));
            if (den) {
                mpz_lcm(den, den, Mij);
            } else {
                ret = 0;
            }
        }
    }
    mpz_mat_realloc(num, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_ptr dst = mpz_mat_entry(num, i, j);
            mpq_srcptr src = mpq_mat_entry_const(M, i, j);
            if (den) {
                mpz_divexact(dst, den, mpq_denref(src));
                mpz_mul(dst, dst, mpq_numref(src));
            } else {
                mpz_set(dst, mpq_numref(src));
            }
        }
    }
    return ret;
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

static int permutation_signature(const unsigned int * perm, unsigned int n)
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

// Returns A^n mod p for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
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

void mpz_mat_add(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B)
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpz_mat C;
    mpz_mat_init(C,A->m,A->n);
    
    for(unsigned int i = 0 ; i < A->m ; i++){
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpz_add(mpz_mat_entry(C,i,j),mpz_mat_entry_const(A,i,j),mpz_mat_entry_const(B,i,j));
        }
    }
    
    mpz_mat_realloc(D,A->m,A->n);
    mpz_mat_set(D,C);
    mpz_mat_clear(C);
}

/* mpz_poly_eval_mpz_mat_mod ? */
/* TODO: Horner ! */
void mpz_mat_in_poly(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f)
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int n = M->n;
    int d = f->deg;
    mpz_mat total, current_power;
    mpz_mat_init(total,n,n);
    mpz_mat_init(current_power,n,n);
    mpz_mat_set_ui(current_power,1);
    
    for (int i = 0 ; i <= d ; i++) {
        mpz_mat current_term;
        mpz_mat_init(current_term,n,n);
        mpz_t aux;
        mpz_init(aux);
        
        mpz_poly_getcoeff(aux,i,f);
        mpz_mat_multiply_by_mpz(current_term,current_power,aux);
        mpz_mat_add(total,total,current_term);
        mpz_mat_multiply(current_power,current_power,M);
        
        mpz_clear(aux);
        mpz_mat_clear(current_term);
        
    }
    
    mpz_mat_set(D,total);
    mpz_mat_clear(current_power);
    mpz_mat_clear(total);
}

/* mpz_poly_eval_mpz_mat_mod_ui ? */
/* TODO: Horner ! */
void mpz_mat_in_poly_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f, unsigned int p)
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int n = M->n;
    int d = f->deg;
    mpz_mat total, current_power;
    mpz_mat_init(total,n,n);
    mpz_mat_init(current_power,n,n);
    mpz_mat_set_ui(current_power,1);
    
    for (int i = 0 ; i <= d ; i++) {
        mpz_mat current_term;
        mpz_mat_init(current_term,n,n);
        mpz_t aux;
        mpz_init(aux);
        
        mpz_poly_getcoeff(aux,i,f);
        mpz_mat_multiply_by_mpz(current_term,current_power,aux);
        mpz_mat_add(total,total,current_term);
        mpz_mat_multiply(current_power,current_power,M);
        
        mpz_clear(aux);
        mpz_mat_clear(current_term);
        
        mpz_mat_mod_ui(total,total,p);
        mpz_mat_mod_ui(current_power,current_power,p);
        
    }
    
    mpz_mat_set(D,total);
    mpz_mat_clear(current_power);
    mpz_mat_clear(total);
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
/* {{{ Gaussian reduction over Z/pZ
 *
 * We sort of return something even in the case where p is not prime.
 *
 * M is the reduced matrix (modification in place of the input), T the
 * transformation matrix (so that T * input-M = output-M)
 *
 * Note that if M is rectangular, this the function we want to call in
 * order to get a row-reduced echelon form of M.
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

static void push_heap(unsigned int * heap, unsigned int x, unsigned int n, cmp_t cmp, void * arg)
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
static void pop_heap(unsigned int * x, unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
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
static void make_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = 1 ; i <= n ; i++) {
        unsigned int t = heap[i-1];
        push_heap(heap, t, i, cmp, arg);
    }
}

/* must be a heap already */
static void sort_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = n ; i ; i--) {
        unsigned int t;
        pop_heap(&t, heap, i, cmp, arg);
        heap[i-1] = t;
    }
}

static void sort_heap_inverse(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
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

static int mpz_hnf_helper_heap_aux_cmp(struct mpz_hnf_helper_heap_aux * data, unsigned int a, unsigned int b)
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
static int mpz_hnf_helper(mpz_mat_ptr T, mpz_mat_ptr dT, mpz_ptr * a, unsigned int n0, unsigned int n)
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
    mpz_mat_realloc(K, H->m-r, H->m);
    mpz_mat_submat_swap(    K, 0, 0,
                            T, r, 0, 
                            M->m - r, M->m);

    mpz_clear(p2);
    mpz_mat_clear(T);
    mpz_mat_clear(H);
    
}
/* }}} */
/*{{{ inversion*/
void mpq_mat_invert(mpq_mat_ptr dst, mpq_mat_srcptr src)
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
/*{{{ equality */
int mpq_mat_eq(mpq_mat_srcptr A, mpq_mat_srcptr B)
{
    if((A->m != B->m) || (A->n != B->n)){
        return 0;
    }
    int test = 1;
    unsigned int i = 0;
    while(test && (i < A->m)){
        unsigned int j = 0;
        while(test && (j < A->n)){
            if(!mpq_equal(mpq_mat_entry_const(A,i,j),mpq_mat_entry_const(B,i,j))){
                test = 0;
            }
            j++;
        }
        i++;
    }
    return test;
}
/* }}} */

