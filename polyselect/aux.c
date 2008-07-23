
/************************ arrays of mpz_t ************************************/

/* allocate an array of d coefficients, and initialize it */
mpz_t*
alloc_mpz_array (int d)
{
  mpz_t *f;
  int i;

  f = (mpz_t*) malloc (d * sizeof (mpz_t));
  for (i = 0; i < d; i++)
    mpz_init (f[i]);
  return f;
}

/* reallocate an array having d0 coefficients to d > d0 coefficients */
mpz_t*
realloc_mpz_array (mpz_t *f, int d0, int d)
{
  int i;

  f = (mpz_t*) realloc (f, d * sizeof (mpz_t));
  for (i = d0; i < d; i++)
    mpz_init (f[i]);
  return f;
}

/* free an array of d coefficients */
void
clear_mpz_array (mpz_t *f, int d)
{
  int i;
  
  for (i = 0; i < d; i++)
    mpz_clear (f[i]);
  free (f);
}


