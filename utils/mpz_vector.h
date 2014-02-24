#ifndef MPZ_VECTOR_H_
#define MPZ_VECTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  unsigned int dim; /* dimension of the vector */
  mpz_t *c;         /* its coordinates */
} mpz_vector_struct_t;

typedef mpz_vector_struct_t mpz_vector_t[1];
typedef mpz_vector_struct_t * mpz_vector_ptr;

/* Management of the structure: init, clear, set and swap. */
void mpz_vector_init (mpz_vector_t, unsigned int);
void mpz_vector_clear(mpz_vector_t);
void mpz_vector_swap (mpz_vector_t, mpz_vector_t);
void mpz_vector_set (mpz_vector_t, mpz_vector_t);
void mpz_vector_setcoordinate (mpz_vector_t, unsigned int, mpz_t);
void mpz_vector_setcoordinate_ui (mpz_vector_t, unsigned int, unsigned int);
int mpz_vector_is_coordinate_zero (mpz_vector_t, unsigned int);

/* Implementation of dot product and norm (skew and non-skew version) */
void mpz_vector_dot_product (mpz_t, mpz_vector_t, mpz_vector_t);
void mpz_vector_skew_dot_product (mpz_t, mpz_vector_t, mpz_vector_t, mpz_t);
void mpz_vector_norm (mpz_t, mpz_vector_t);
void mpz_vector_skew_norm (mpz_t, mpz_vector_t, mpz_t);

/* Convert from mpz_vector_t to mpz_poly_t */
void mpz_vector_get_mpz_poly (mpz_poly_t, mpz_vector_t);

/* Operations on vectors */
void mpz_vector_submul (mpz_vector_t r, mpz_t q, mpz_vector_t v);

/* Lagrange algo to reduce lattice of rank 2 */
void mpz_vector_Lagrange (mpz_vector_t, mpz_vector_t,
                          mpz_vector_t, mpz_vector_t, mpz_t);

#ifdef __cplusplus
}
#endif

#endif	/* MPZ_VECTOR_H_ */
