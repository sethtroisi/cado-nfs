#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "macros.h"
#include "sieving_interval.h"
#include "int64_vector.h"
#include <inttypes.h>

void int64_vector_init(int64_vector_ptr v, unsigned int d)
{
  ASSERT(d > 0);
  v->dim = d;
  v->c = (int64_t *) malloc(sizeof(int64_t) * d);
  ASSERT_ALWAYS (v->c != NULL);
  for (unsigned k = 0; k < d; k++) {
    v->c[k] = 0;
  }
}

void int64_vector_clear(int64_vector_ptr v)
{
  free(v->c);
}

void int64_vector_swap(int64_vector_ptr v0, int64_vector_ptr v1)
{
  ASSERT(v0->dim == v1->dim);
  int64_t * tmp = v0->c;
  v0->c = v1->c;
  v1->c = tmp;
}

void int64_vector_set(int64_vector_ptr v, int64_vector_srcptr s)
{
  ASSERT(v->dim == s->dim);
  for (unsigned int i = 0; i < v->dim; i++)
    v->c[i] = s->c[i];
}

void int64_vector_setcoordinate(int64_vector_ptr v, unsigned int i,
                                int64_t z)
{
  ASSERT(i < v->dim);
  v->c[i] = z;
}

int int64_vector_equal(int64_vector_srcptr a, int64_vector_srcptr b)
{
  int r = (a->dim > b->dim) - (b->dim > a->dim);
  if (r) return 1;
  for(int d = a->dim; d >= 0 ; d--) {
    r = a->c[d] - b->c[d];
    if (r) return 1;
  }
  return 0;
}

void int64_vector_fprintf(FILE * file, int64_vector_srcptr v)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < v->dim - 1; i++)
    fprintf(file, "%" PRId64 ", ", v->c[i]);
  fprintf(file, "%" PRId64 "]\n", v->c[v->dim - 1]);
}

void int64_vector_add_one(int64_vector_ptr v, sieving_interval_srcptr H)
{
  int64_vector_add_one_i(v, 0, H);
}

unsigned int int64_vector_add_one_i(int64_vector_ptr v, unsigned int i,
                                    sieving_interval_srcptr H)
{
  ASSERT(v->dim == H->t);
  ASSERT(i < v->dim);

#ifndef NDEBUG
  int64_t tmp = 0;
  int64_t tmp2 = 0;
  for (unsigned int j = 0; j < v->dim; j++) {
    tmp = tmp + v->c[j];
    tmp2 = tmp2 + H->h[j];
  }
  ASSERT(tmp <= tmp2);
#endif

  unsigned int k = i;
  while(k < v->dim) {
    if (v->c[k] == H->h[k] - 1) {
      int64_vector_setcoordinate(v, k, -(int64_t)H->h[k]);
      k++;
    } else {
      break;
    }
  }
  if (k < v->dim) {
    v->c[k]++;
  }

#ifndef NDEBUG
  for (unsigned int j = 0; j < v->dim - 1; j++) {
    tmp = v->c[j];
    ASSERT(tmp >= -(int64_t)H->h[j]);
    ASSERT(tmp < (int64_t)H->h[j]);
  }
  tmp = v->c[v->dim -1];
  ASSERT(tmp >= 0);
  ASSERT(tmp < H->h[v->dim - 1]);
#endif

  return k;
}

void int64_vector_to_mpz_vector(mpz_vector_ptr a, int64_vector_srcptr b)
{
  ASSERT(a->dim == b->dim);

  for (unsigned int i = 0; i < b->dim; i++) {
    mpz_vector_setcoordinate_int64 (a, i, b->c[i]);
  }
}

double int64_vector_norml2sqr(int64_vector_srcptr a)
{
  double norm = 0;
  for (unsigned int i = 0; i < a->dim; i++) {
    norm += (double)(a->c[i] * a->c[i]);
  }
  return norm;
}

double int64_vector_norml2(int64_vector_srcptr a)
{
  return sqrt(int64_vector_norml2sqr(a));
}

int64_t int64_vector_dot_product(int64_vector_srcptr v0, int64_vector_srcptr v1)
{
  ASSERT(v0->dim == v1->dim);

  int64_t dot = 0;
  for (unsigned int i = 0; i < v0->dim; i++) {
    dot += v0->c[i] * v1->c[i];
  }
  return dot;
}
