/* header file for common routines to polyselect and kleinjung */

mpz_t* alloc_mpz_array   (int);
mpz_t* realloc_mpz_array (mpz_t *, int, int);
void   clear_mpz_array   (mpz_t *, int);
