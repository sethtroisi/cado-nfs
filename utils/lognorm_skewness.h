#ifndef LOGNORM_SKEWNESS_H_
#define LOGNORM_SKEWNESS_H_


void mpz_poly_innerproduct_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr);
void mpz_poly_norm_L2 (mpz_t, mpz_t, mpz_poly_srcptr);
void mpz_poly_skew_innerproduct_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr, mpz_t skew);
void mpz_poly_skew_norm_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_t skew);
void mpz_vector_innerproduct_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr);
void mpz_vector_norm_L2 (mpz_t, mpz_t, mpz_vector_srcptr);
void mpz_vector_skew_innerproduct_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr, mpz_t skew);
void mpz_vector_skew_norm_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_t skew);
double double_poly_innerproduct_L2 (double_poly_srcptr, double_poly_srcptr);
double double_poly_norm_L2 (double_poly_srcptr);
double double_poly_skew_innerproduct_L2 (double_poly_srcptr, double_poly_srcptr, double skew);
double double_poly_skew_norm_L2 (double_poly_srcptr, double skew);

void mpz_poly_innerproduct_cir_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr);
void mpz_poly_norm_cir_L2 (mpz_t, mpz_t, mpz_poly_srcptr);
void mpz_poly_skew_innerproduct_cir_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr, mpz_t skew);
void mpz_poly_skew_norm_cir_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_t skew);
void mpz_vector_innerproduct_cir_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr);
void mpz_vector_norm_cir_L2 (mpz_t, mpz_t, mpz_vector_srcptr);
void mpz_vector_skew_innerproduct_cir_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr, mpz_t skew);
void mpz_vector_skew_norm_cir_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_t skew);
double double_poly_innerproduct_cir_L2 (double_poly_srcptr, double_poly_srcptr);
double double_poly_norm_cir_L2 (double_poly_srcptr);
double double_poly_skew_innerproduct_cir_L2 (double_poly_srcptr, double_poly_srcptr, double skew);
double double_poly_skew_norm_cir_L2 (double_poly_srcptr, double skew);

void mpz_poly_innerproduct_sq_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr);
void mpz_poly_norm_sq_L2 (mpz_t, mpz_t, mpz_poly_srcptr);
void mpz_poly_skew_innerproduct_sq_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_poly_srcptr, mpz_t skew);
void mpz_poly_skew_norm_sq_L2 (mpz_t, mpz_t, mpz_poly_srcptr, mpz_t skew);
void mpz_vector_innerproduct_sq_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr);
void mpz_vector_norm_sq_L2 (mpz_t, mpz_t, mpz_vector_srcptr);
void mpz_vector_skew_innerproduct_sq_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_vector_srcptr, mpz_t skew);
void mpz_vector_skew_norm_sq_L2 (mpz_t, mpz_t, mpz_vector_srcptr, mpz_t skew);
double double_poly_innerproduct_sq_L2 (double_poly_srcptr, double_poly_srcptr);
double double_poly_norm_sq_L2 (double_poly_srcptr);
double double_poly_skew_innerproduct_sq_L2 (double_poly_srcptr, double_poly_srcptr, double skew);
double double_poly_skew_norm_sq_L2 (double_poly_srcptr, double skew);

#endif  /* LOGNORM_SKEWNESS_H_ */
