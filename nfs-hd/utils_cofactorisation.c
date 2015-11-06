#include "utils_cofactorisation.h"
#include "utils_norm.h"

#include "ecm/facul.h"
#include "ecm/facul_doit.h"

/*
 * Initialize an array of factors with number, the maximum number of elements.
 *
 * factor: the array of factors.
 * number: the maximum number of elements.
 */
static void factor_init(factor_ptr factor, unsigned int alloc)
{
  ASSERT(alloc > 0);

  factor->alloc = alloc;
  factor->number = 0;
  factor->factorization = (mpz_t * ) malloc(sizeof(mpz_t) * alloc);
  for (unsigned int i = 0; i < alloc; i++) {
    mpz_init(factor->factorization[i]);
  }
}

/*
 * Delete an array of factors.
 *
 * factor: the array of factors.
 */
static void factor_clear(factor_ptr factor)
{
  for (unsigned int i = 0; i < factor->alloc; i++) {
    mpz_clear(factor->factorization[i]);
  }
  free(factor->factorization);
  factor->alloc = 0;
  factor->number = 0;
}

/* obvious! */
static void factor_append(factor_ptr factor, mpz_srcptr z)
{
  if (factor->alloc == factor->number) {
    unsigned int newalloc = factor->alloc + 10;
    factor->factorization = (mpz_t * ) realloc(factor->factorization,
        sizeof(mpz_t) * newalloc);
    for (unsigned int i = factor->alloc; i < newalloc; ++i)
      mpz_init(factor->factorization[i]);
    factor->alloc = newalloc;
  }

  mpz_set(factor->factorization[factor->number], z);
  factor->number++;
}

/*
 * Print an array of factors.
 *
 * factor: an array of factors.
 */
MAYBE_UNUSED static void factor_fprintf(FILE * file, factor_srcptr factor)
{
  fprintf(file, "[");
  if (factor->number != 0) {
    for (unsigned int i = 0; i < factor->number - 1; i++) {
      gmp_fprintf(file, "%Zd, ", factor->factorization[i]);
    }
    gmp_fprintf(file, "%Zd]\n", factor->factorization[factor->number - 1]);
  } else {
    fprintf(file, "]\n");
  }
}

/*
 * TODO: useless.
 * Test if the maximum factor of the array factor is less or equal to B. Return
 *  1 if true, 0 otherwise. If factor is sorted, set sort to 1, 0 otherwise.
 *
 * factor: an array of factors.
 * B: the smoothness bound.
 * sort: 1 if factor is sorted, 0 otherwise.
 */
MAYBE_UNUSED static unsigned int factor_is_smooth(factor_srcptr factor, mpz_t B,
    unsigned int sort)
{
  if (sort) {
    ASSERT(sort == 1);

    if (mpz_cmp(factor->factorization[factor->number - 1], B) > 0) {
      return 0;
    }
  } else {
    ASSERT(sort == 0);

    for (unsigned int i = 0; i < factor->number; i++) {
      if (mpz_cmp(factor->factorization[i], B) > 0) {
        return 0;
      }
    }
  }

  return 1;
}

#ifdef ASSERT_FACTO
/*
 * Return 1 if the factorisation is good, 0 otherwise.
 */
static unsigned int factor_assert(factor_srcptr factor, mpz_srcptr z)
{
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(tmp, 1);
  for (unsigned int i = 0; i < factor->number; i++) {
    mpz_mul(tmp, tmp, factor->factorization[i]);
  }
  unsigned int assert_facto = 0;
  if (!mpz_cmp(tmp, z)) {
    ASSERT(mpz_cmp(tmp, z) == 0);

    assert_facto = 1;
  }
  mpz_clear(tmp);

  return assert_facto;
}
#endif // ASSERT_FACTO

/*
 * Remove all the small factors under a certain bound, and store z_root /
 *  (factors) in z. Return 1 if z_root is entirely factorize.
 */
static int brute_force_factorize_ul(factor_ptr factor, mpz_ptr z,
    mpz_srcptr z_root, unsigned long bound)
{
  prime_info pi;
  prime_info_init (pi);

  mpz_set(z, z_root);
  mpz_t prime_Z;
  mpz_init(prime_Z);
  unsigned long prime = 2;

  for (prime = 2; prime <= bound; prime = getprime_mt(pi)) {
    if (mpz_cmp_ui(z, 1) != 0) {
      mpz_set_ui(prime_Z, prime);
      mpz_t q;
      mpz_init(q);
      mpz_t r;
      mpz_init(r);
      mpz_fdiv_qr(q, r, z, prime_Z);
      while (mpz_cmp_ui(r, 0) == 0) {
        mpz_set(z, q);
        factor_append(factor, prime_Z);
        mpz_fdiv_qr(q, r, z, prime_Z);
      }
      mpz_clear(q);
      mpz_clear(r);
    }
  }

  mpz_clear(prime_Z);
  prime_info_clear (pi);

  if (mpz_cmp_ui(z, 1) == 0) {
    return 1;
  }
  return 0;
}

static int compare_factor(const void * p0, const void * p1)
{
  return(mpz_cmp((mpz_srcptr) p0, (mpz_srcptr) p1));
}

/*
 * To sort by acending value the element of the factor array.
 *
 * factor: array of factors.
 */
static void sort_factor(factor_ptr factor)
{
  qsort(factor->factorization, factor->number,
      sizeof(factor->factorization[0]), compare_factor);
}

/*
 * Print a relation.
 *
 * factor: factorisation of a if the norm is smooth is the corresponding number
 *  field.
 * I:
 * L:
 * a: the polynomial in the original lattice.
 * t: dimension of the lattice.
 * V: number of number fields.
 * size:
 * assert_facto: 1 if the factorisation is correct, 0 otherwise.
 */
static void printf_relation(factor_t * factor, unsigned int * I,
    unsigned int * L, mpz_poly_srcptr a, unsigned int t, unsigned int V,
    unsigned int size, MAYBE_UNUSED unsigned int * assert_facto)
{
  //Print a.
#ifdef PRINT_POLY_RELATION
  printf("# ");
  mpz_poly_fprintf(stdout, a);
#endif // PRINT_POLY_RELATION

  unsigned int index = 0;
  //Print if the factorisation is good.
#ifdef ASSERT_FACTO
  printf("# ");
  for (unsigned int i = 0; i < V - 1; i++) {
    if (index < size) {
      if (i == L[index]) {
        if (I[index]) {
          printf("%u:", assert_facto[index]);
        } else {
          printf(":");
        }
        index++;
      } else {
        printf(":");
      }
    } else {
      printf(":");
    }
  }

  if (index < size) {
    if (V - 1 == L[index]) {
      if (I[index]) {
        printf("%u", assert_facto[index]);
      }
      index++;
    }
  }

  printf("\n");
#endif // ASSERT_FACTO

  //Print the coefficient of a.
  for (int i = 0; i < a->deg; i++) {
    gmp_printf("%Zd,", a->coeff[i]);
  }
  if ((int)t - 1 == a->deg) {
    gmp_printf("%Zd:", a->coeff[a->deg]);
  } else {
    gmp_printf("%Zd,", a->coeff[a->deg]);
    for (int i = a->deg + 1; i < (int)t - 1; i++) {
      printf("0,");
    }
    printf("0:");
  }

  //Print the factorisation in the different number field.
  index = 0;
  for (unsigned int i = 0; i < V - 1; i++) {
    if (index < size) {
      if (i == L[index]) {
        if (I[index]) {
          for (unsigned int j = 0; j < factor[index]->number - 1; j++) {
#ifdef PRINT_DECIMAL
            gmp_printf("%Zd,", factor[index]->factorization[j]);
#else // PRINT_DECIMAL
            printf("%s,",
                mpz_get_str(NULL, 16, factor[index]->factorization[j]));
#endif // PRINT_DECIMAL
          }
#ifdef PRINT_DECIMAL
          gmp_printf("%Zd:", factor[index]->factorization[
                       factor[index]->number - 1]);
#else // PRINT_DECIMAL
          printf("%s:", mpz_get_str(NULL, 16,
                factor[index]->factorization[factor[index]->number - 1]));
#endif // PRINT_DECIMAL
        } else {
          printf(":");
        }
        index++;
      } else {
        printf(":");
      }
    } else {
      printf(":");
    }
  }

  if (index < size) {
    if (V - 1 == L[index]) {
      if (I[index]) {
        for (unsigned int j = 0; j < factor[index]->number - 1; j++) {
#ifdef PRINT_DECIMAL
          gmp_printf("%Zd,", factor[index]->factorization[j]);
#else // PRINT_DECIMAL
          printf("%s,",
                mpz_get_str(NULL, 16, factor[index]->factorization[j]));
#endif // PRINT_DECIMAL
        }
#ifdef PRINT_DECIMAL
        gmp_printf("%Zd", factor[index]->factorization[
                     factor[index]->number - 1]);
#else // PRINT_DECIMAL
        printf("%s", mpz_get_str(NULL, 16,
              factor[index]->factorization[factor[index]->number - 1]));
#endif // PRINT_DECIMAL
      }
      index++;
    }
  }
  printf("\n");
}


typedef struct {
  unsigned long lpb;            // large prime bound = 2^lpb
  unsigned long fbb;            // fbb (the real bound, not its log)
  double BB;                    // square of fbb
  double BBB;                   // cube of fbb
  facul_method_t * methods;     // list of ECMs (null-terminated)
} facul_aux_data;


// FIXME: This has been duplicated from facul.cpp
// (September 2015).
// Maybe merge again, at some point.
// It returns -1 if the factor is not smooth, otherwise the number of
// factors.
// Remark: FACUL_NOT_SMOOTH is just -1.
static int
facul_aux (mpz_t *factors, const struct modset_t m,
    const facul_aux_data *data, int method_start)
{
  int found = 0;
  facul_method_t* methods = data->methods;
  if (methods == NULL)
    return found;

  int i = 0;
  for (i = method_start; methods[i].method != 0; i++)
  {
    struct modset_t fm, cfm;
    int res_fac = 0;

    switch (m.arith) {
      case CHOOSE_UL:
        res_fac = facul_doit_onefm_ul(factors, m.m_ul,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_15UL:
        res_fac = facul_doit_onefm_15ul(factors, m.m_15ul,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_2UL2:
        res_fac = facul_doit_onefm_2ul2 (factors, m.m_2ul2,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      case CHOOSE_MPZ:
        res_fac = facul_doit_onefm_mpz (factors, m.m_mpz,
            methods[i], &fm, &cfm,
            data->lpb, data->BB, data->BBB);
        break;
      default: abort();
    }
    // check our result
    // res_fac contains the number of factors found
    if (res_fac == -1)
    {
      /*
         The cofactor m is not smooth. So, one stops the
         cofactorization.
         */
      found = FACUL_NOT_SMOOTH;
      break;
    }
    if (res_fac == 0)
    {
      /* Zero factor found. If it was the last method for this
         side, then one stops the cofactorization. Otherwise, one
         tries with an other method */
      continue;
    }

    found += res_fac;
    if (res_fac == 2)
      break;

    /*
       res_fac == 1  Only one factor has been found. Hence, our
       factorization is not finished.
       */
    if (fm.arith != CHOOSE_NONE)
    {
      int found2 = facul_aux(factors+res_fac, fm, data, i+1);
      if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
      {
        found = FACUL_NOT_SMOOTH;
        modset_clear(&cfm);
        modset_clear(&fm);
        break;
      }
      else
        found += found2;
      modset_clear(&fm);
    }
    if (cfm.arith != CHOOSE_NONE)
    {
      int found2 = facul_aux(factors+res_fac, cfm, data, i+1);
      if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
      {
        found = FACUL_NOT_SMOOTH;
      }
      else
        found += found2;
      modset_clear(&cfm);
      break;
    }
    break;
  }
  return found;
}


// Returns 1 if norm is smooth, 0 otherwise.
// Factorization is given in factors.
static int call_facul(factor_ptr factors, mpz_srcptr norm_r,
    facul_aux_data * data) {
  unsigned long B = data->fbb;

  mpz_t norm;
  mpz_init(norm);
  mpz_set(norm, norm_r);

  // Trial divide all the small factors.
  int success = brute_force_factorize_ul(factors, norm, norm, B);
  if (success) {
    return 1;
  }

  // TODO: move this test after the next block
  // and do the primality test with fast routines.
  if (mpz_probab_prime_p(norm, 1)) {
    if (mpz_sizeinbase(norm, 2) <= data->lpb) {
      factor_append(factors, norm);
      return 1;
    } else {
      return 0;
    }
  }

  // Prepare stuff for calling facul() machinery.
  // Choose appropriate arithmetic according to size of norm
  // Result is stored in n, which is of type struct modset_t (see facul.h)
  struct modset_t n;
  n.arith = CHOOSE_NONE;
  size_t bits = mpz_sizeinbase(norm, 2);
  if (bits <= MODREDCUL_MAXBITS) {
    ASSERT(mpz_fits_ulong_p(norm));
    modredcul_initmod_ul(n.m_ul, mpz_get_ui(norm));
    n.arith = CHOOSE_UL;
  }
  else if (bits <= MODREDC15UL_MAXBITS)
  {
    unsigned long t[2];
    modintredc15ul_t m;
    size_t written;
    mpz_export(t, &written, -1, sizeof(unsigned long), 0, 0, norm);
    ASSERT_ALWAYS(written <= 2);
    modredc15ul_intset_uls(m, t, written);
    modredc15ul_initmod_int(n.m_15ul, m);
    n.arith = CHOOSE_15UL;
  }
  else if (bits <= MODREDC2UL2_MAXBITS)
  {
    unsigned long t[2];
    modintredc2ul2_t m;
    size_t written;
    mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, norm);
    ASSERT_ALWAYS(written <= 2);
    modredc2ul2_intset_uls (m, t, written);
    modredc2ul2_initmod_int (n.m_2ul2, m);
    n.arith = CHOOSE_2UL2;
  }
  else
  {
    modmpz_initmod_int(n.m_mpz, norm);
    n.arith = CHOOSE_MPZ;
  }
  ASSERT_ALWAYS(n.arith != CHOOSE_NONE);

  // Call the facul machinery.
  // TODO: think about this hard-coded 16...
  mpz_t * fac = (mpz_t *) malloc(sizeof(mpz_t)*16);
  for (int i = 0; i < 16; ++i)
    mpz_init(fac[i]);
  int found = facul_aux(fac, n, data, 0);
  if (found > 0) {
    for (int i = 0; i < found; ++i)
      factor_append(factors, fac[i]);
  }

  for (int i = 0; i < 16; ++i)
    mpz_clear(fac[i]);
  free(fac);
  mpz_clear(norm);
  return found > 0;
}

/*
 * A potential polynomial is found. Factorise it to verify if it gives a
 *  relation.
 *
 * a: the polynomial in the original lattice.
 * f: the polynomials that define the number fields.
 * lpb: the large prime bounds.
 * L:
 * size:
 * t: dimension of the lattice.
 * V: number of number fields.
 */
static void good_polynomial(mpz_poly_srcptr a, mpz_poly_t * f,
    unsigned int * L, unsigned int size, unsigned int t, unsigned int V,
    int main, facul_aux_data *data, unsigned int * nb_rel_found,
    ideal_spq_srcptr special_q, unsigned int q_side,
    MAYBE_UNUSED unsigned int * number_factorisation)
{
  mpz_t res;
  mpz_init(res);

  factor_t * factor = (factor_t * ) malloc(sizeof(factor_t) * size);
  unsigned int * I = (unsigned int * ) malloc(sizeof(unsigned int) * size);

  unsigned int * assert_facto;
#ifdef ASSERT_FACTO
    assert_facto = (unsigned int * ) malloc(sizeof(unsigned int) * size);
    mpz_t res_tmp;
    mpz_init(res_tmp);
#else // ASSERT_FACTO
    assert_facto = NULL;
#endif // ASSERT_FACTO

  unsigned int find = 0;
  if (main == -1) {
    for (unsigned int i = 0; i < size; i++) {
      factor_init(factor[i], 10);
      norm_poly(res, f[L[i]], a);

#ifdef ASSERT_FACTO
      mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

      if (i == q_side) {
        mpz_t spq;
        mpz_init(spq);
        mpz_set_uint64(spq, ideal_spq_get_q(special_q));
        ASSERT(ideal_spq_get_deg_g(special_q) >= 1);
        for (int deg = 0; deg < ideal_spq_get_deg_g(special_q); deg++) {
          mpz_divexact(res, res, spq);
          factor_append(factor[i], spq);
        }
        mpz_clear(spq);
      }

      int is_smooth = call_facul(factor[i], res, &data[L[i]]);

#ifdef ASSERT_FACTO
      assert_facto[i] = factor_assert(factor[i], res_tmp);
      if (assert_facto[i] == 0) {
        number_factorisation[0] = number_factorisation[0] + 1;
      } else {
        number_factorisation[1] = number_factorisation[1] + 1;
      }
#endif // ASSERT_FACTO

      if (is_smooth) {
        find++;
        I[i] = 1;
        sort_factor(factor[i]);
      } else {
        I[i] = 0;
      }
    }
  } else {
    // TODO: this part of code is false for cofactorisation at least.
    ASSERT(main >= 0);
    unsigned int Lmain = 0;
    while (L[Lmain] != (unsigned int)main) {
      Lmain++;
    }

    norm_poly(res, f[main], a);

#ifdef ASSERT_FACTO
      mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

    int is_smooth = call_facul(factor[Lmain], res, &data[main]);

#ifdef ASSERT_FACTO
      assert_facto[Lmain] = factor_assert(factor[Lmain], res);
      if (assert_facto[Lmain] == 0) {
        number_factorisation[0] = number_factorisation[0] + 1;
      } else {
        number_factorisation[1] = number_factorisation[1] + 1;
      }
#endif // ASSERT_FACTO

    if (is_smooth) {
      find = 1;

      for (unsigned int i = 0; i < size; i++) {
        if (i != Lmain) {
          norm_poly(res, f[L[i]], a);

#ifdef ASSERT_FACTO
          mpz_set(res_tmp, res);
#endif // ASSERT_FACTO

          is_smooth = call_facul(factor[i], res, &data[L[i]]);

#ifdef ASSERT_FACTO
          assert_facto[i] = factor_assert(factor[i], res);
          if (assert_facto[i] == 0) {
            number_factorisation[0] = number_factorisation[0] + 1;
          } else {
            number_factorisation[1] = number_factorisation[1] + 1;
          }
#endif // ASSERT_FACTO

          if (is_smooth) {
            find++;
            I[i] = 1;
            sort_factor(factor[i]);
          } else {
            I[i] = 0;
          }
        }
      }
    }
  }

  if (find >= 2) {
    printf_relation(factor, I, L, a, t, V, size, assert_facto);

    (* nb_rel_found)++;

  }

  mpz_clear(res);
  for (unsigned int i = 0; i < size; i++) {
    factor_clear(factor[i]);
  }
  free(factor);
  free(I);

#ifdef ASSERT_FACTO
  mpz_clear(res_tmp);
  free(assert_facto);
#endif
}

/*
 * Return the sum of the element in index. index has size V.
 *
 * index: array of index.
 * V: number of element, number of number fields.
 */
static uint64_t sum_index(uint64_t * index, unsigned int V, int main)
{
  ASSERT(main >= -1);
  ASSERT(V >= 2);

  if (main != -1) {
    return index[main];
  }
  uint64_t sum = 0;
  for (unsigned int i = 0; i < V; i++) {
    sum += index[i];
  }
  return sum;
}

/*
 * Find in all the indices which one has the smallest value.
 *
 *
 *
 */
static unsigned int find_indices(unsigned int ** L,
    uint64_array_t * indices, uint64_t * index, unsigned int V,
    uint64_t max_indices)
{
  * L = (unsigned int * ) malloc(sizeof(unsigned int) * (V));
  unsigned int size = 0;
  uint64_t min = max_indices;
  for (unsigned int i = 0; i < V; i++) {
    if (indices[i]->length != 0 && index[i] < indices[i]->length) {
      if (indices[i]->array[index[i]] < min) {
        min = indices[i]->array[index[i]];
        (*L)[0] = i;
        size = 1;
      } else if (min == indices[i]->array[index[i]]) {
        (*L)[size] = i;
        size++;
      }
    }
  }
  * L = realloc(* L, size * sizeof(unsigned int));
  return size;
}

static unsigned int find_indices_main(unsigned int ** L,
    uint64_array_t * indices, uint64_t * index, unsigned int V, int main)
{
  ASSERT(main >= 0);

  unsigned int size = 1;
  uint64_t target = indices[main]->array[index[main]];
  (*L)[0] = main;
  for (unsigned int i = 0; i < V; i++) {
    if (i != (unsigned int) main && indices[i]->length != 0 && index[i] <
        indices[i]->length) {
      int test = 0;
      while (target > indices[i]->array[index[i]]) {
        index[i] = index[i] + 1;
        if (index[i] == indices[i]->length) {
          test = 1;
          break;
        }
      }
      if (test) {
        ASSERT(test == 1);

        if (target == indices[i]->array[index[i]]) {
          (*L)[size] = main;
          size++;
        }
      }
    }
  }

  * L = realloc(* L, size * sizeof(unsigned int));
  return size;
}

/*
 * For
 */
static void find_relation(uint64_array_t * indices, uint64_t * index,
    uint64_t number_element, mat_Z_srcptr matrix, mpz_poly_t * f,
    sieving_bound_srcptr H, unsigned int V, int main, uint64_t max_indices,
    facul_aux_data *data, unsigned int * nb_rel_found,
    ideal_spq_srcptr special_q, unsigned int q_side,
    MAYBE_UNUSED unsigned int * number_factorisation)
{
  unsigned int * L = (unsigned int *) malloc(V * sizeof(unsigned int));
  unsigned int size = 0;
  if (main == -1) {
    size = find_indices(&L, indices, index, V, max_indices);
  } else {
    size = find_indices_main(&L, indices, index, V, main);
  }

  ASSERT(size >= 1);

  if (size >= 2) {
    mpz_vector_t c;
    mpz_t gcd;
    mpz_init(gcd);
    mpz_vector_init(c, H->t);
    array_index_mpz_vector(c, indices[L[0]]->array[index[L[0]]], H,
        number_element);

    mpz_poly_t a;
    mpz_poly_init(a, 0);
    mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, c);
    mpz_poly_content(gcd, a);

    //a must be irreducible.
    if (mpz_cmp_ui(gcd, 1) == 0 && a->deg > 0 &&
        mpz_cmp_ui(mpz_poly_lc_const(a), 0) > 0) {

      good_polynomial(a, f, L, size, H->t, V, main, data, nb_rel_found,
          special_q, q_side, number_factorisation);
    }

    mpz_poly_clear(a);

    mpz_clear(gcd);
    mpz_vector_clear(c);

    for (unsigned int i = 0; i < size; i++) {
      index[L[i]] = index[L[i]] + 1;
    }
  } else {
    index[L[0]] = index[L[0]] + 1;
  }

  free(L);
}

unsigned int find_relations(uint64_array_t * indices, uint64_t number_element,
    unsigned int * lpb, mat_Z_srcptr matrix, mpz_poly_t * f,
    sieving_bound_srcptr H, unsigned int V, ideal_spq_srcptr special_q,
    unsigned int q_side, int main)
{
  //index[i] is the current index of indices[i].
  uint64_t * index = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  uint64_t length_tot = 0;
  //Maximum of all the indices.
  uint64_t max_indices = 0;
  for (unsigned int i = 0; i < V; i++) {
    index[i] = 0;
    if (i == (unsigned int)main) {
      ASSERT(main >= 0);

      length_tot = indices[main]->length;
    } else if (indices[i]->length != 0 && main == -1) {
      length_tot += indices[i]->length;
      if (max_indices < indices[i]->array[indices[i]->length - 1]) {
        max_indices = indices[i]->array[indices[i]->length - 1];
      }
    }
  }

  //  prepare cofactorization strategy
  facul_aux_data * data;
  data = (facul_aux_data *)malloc(V*sizeof(facul_aux_data));
  ASSERT_ALWAYS(data != NULL);
  for (unsigned int i = 0; i < V; ++i) {
    unsigned int B = 1000; // FIXME: should be a parameter.
    data[i].lpb = (unsigned long) lpb[i];
    data[i].fbb = B;
    data[i].BB = ((double)B) * ((double)B);
    data[i].BBB = ((double)B) * data[i].BB;
    data[i].methods = facul_make_aux_methods(nb_curves95(lpb[i]), 0, 0);
  }

  MAYBE_UNUSED unsigned int * number_factorisation =
    (unsigned int *) malloc(sizeof(unsigned int) * 2);
  //Bad factorisation.
  number_factorisation[0] = 0;
  //Good factorisation.
  number_factorisation[1] = 0;
  unsigned int nb_rel_found = 0;
  if (0 != length_tot) {
    while(sum_index(index, V, main) < length_tot) {
      find_relation(indices, index, number_element, matrix, f, H, V, main,
          max_indices, data, &nb_rel_found, special_q, q_side,
          number_factorisation);
    }
  }
  if (nb_rel_found != 0) {
    printf("# Number of found relations: %u.\n", nb_rel_found);
  } else {
    printf("# No relations\n");
  }

#ifdef ASSERT_FACTO
  printf("# Number of good factorisations: %u.\n", number_factorisation[1]);
  printf("# Number of bad factorisations: %u.\n", number_factorisation[0]);
#endif // ASSERT_FACTO
  free(number_factorisation);

  for (unsigned int i = 0; i < V; ++i)
    facul_clear_aux_methods(data[i].methods);
  free(data);
  free(index);

  return nb_rel_found;
}
