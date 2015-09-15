#include "utils_cofactorisation.h"
#include "utils_mpz.h"
#include "utils_norm.h"

#include "ecm/facul.h"
#include "ecm/facul_doit.h"

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
  printf("# ");
  mpz_poly_fprintf(stdout, a);

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
        printf("%s", mpz_get_str(NULL, 16,
              factor[index]->factorization[factor[index]->number - 1]));
#else // PRINT_DECIMAL
        gmp_printf("%Zd", factor[index]->factorization[
                     factor[index]->number - 1]);
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
static int call_facul(factor_ptr factors, mpz_srcptr norm_r, facul_aux_data * data) {
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
    int main, facul_aux_data *data)
{
  mpz_t res;
  mpz_init(res);
  factor_t * factor = (factor_t * ) malloc(sizeof(factor_t) * size);
  unsigned int * I = (unsigned int * ) malloc(sizeof(unsigned int) * size);

  unsigned int * assert_facto;
#ifdef ASSERT_FACTO
    assert_facto = (unsigned int * ) malloc(sizeof(unsigned int) * size);
#else // ASSERT_FACTO
    assert_facto = NULL;
#endif // ASSERT_FACTO

  unsigned int find = 0;
  if (main == -1) {
    for (unsigned int i = 0; i < size; i++) {
      norm_poly(res, f[L[i]], a);

      int is_smooth = call_facul(factor[i], res, &data[L[i]]);
#ifdef ASSERT_FACTO
      assert_facto[i] = factor_assert(factor[i], res);
#endif

      if (is_smooth) {
        find++;
        I[i] = 1;
        sort_factor(factor[i]);
      } else {
        I[i] = 0;
      }
    }
  } else {
    // TODO: this part of code is not tested.
    ASSERT(main >= 0);
    unsigned int Lmain = 0;
    while (L[Lmain] != (unsigned int)main) {
      Lmain++;
    }

    norm_poly(res, f[main], a);

    int main_is_smooth = call_facul(factor[Lmain], res, &data[main]);

    if (main_is_smooth) {
      find = 1;
#ifdef ASSERT_FACTO
      assert_facto[Lmain] = factor_assert(factor[Lmain], res);
#endif

      for (unsigned int i = 0; i < size; i++) {
        if (i != Lmain) {
          norm_poly(res, f[L[i]], a);

          int is_smooth = call_facul(factor[i], res, &data[L[i]]);
#ifdef ASSERT_FACTO
          assert_facto[i] = factor_assert(factor[i], res);
#endif

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

  }

  mpz_clear(res);
  for (unsigned int i = 0; i < size; i++) {
    factor_clear(factor[i]);
  }
  free(factor);
  free(I);

#ifdef ASSERT_FACTO
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

static unsigned int find_indices_main(unsigned int ** L, uint64_array_t * indices,
    uint64_t * index, unsigned int V, int main)
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
    facul_aux_data *data)
{
  unsigned int * L;
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

#ifdef NUMBER_SURVIVALS
      number_survivals++;
#endif // NUMBER_SURVIVALS

    //a must be irreducible.
    if (mpz_cmp_ui(gcd, 1) == 0 && a->deg > 0 &&
        mpz_cmp_ui(mpz_poly_lc_const(a), 0) > 0) {

#ifdef NUMBER_SURVIVALS
      number_survivals_facto++;
#endif // NUMBER_SURVIVALS

      good_polynomial(a, f, L, size, H->t, V, main, data);
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

void find_relations(uint64_array_t * indices, uint64_t number_element,
    mpz_t * lpb, mat_Z_srcptr matrix, mpz_poly_t * f, sieving_bound_srcptr H,
    unsigned int V, int main)
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
    size_t lpb_bit = mpz_sizeinbase(lpb[i], 2);
    unsigned int B = 1000; // FIXME: should be a parameter.
    data[i].lpb = lpb_bit;
    data[i].fbb = B;
    data[i].BB = ((double)B) * ((double)B);
    data[i].BBB = ((double)B) * data[i].BB;
    data[i].methods = facul_make_aux_methods(nb_curves95(lpb_bit), 0, 0);
  }

  if (0 != length_tot) {
    while(sum_index(index, V, main) < length_tot) {
      find_relation(indices, index, number_element, matrix, f, H, V, main,
          max_indices, data);
    }
  } else {
    printf("# No relations\n");
  }
  for (unsigned int i = 0; i < V; ++i)   
    facul_clear_aux_methods(data[i].methods);
  free(data);
  free(index);
}
