#include <stdlib.h>
#include <inttypes.h>
#include "list_int64_vector.h"

/*
 * List of int64_vector.
 */

void list_int64_vector_init(list_int64_vector_ptr list)
{
  list->length = 0;
  list->v = (int64_vector_t * ) malloc(sizeof(int64_vector_t) *
      DEFAULT_LENGTH_LIST_INT64_VECTOR);
}

void list_int64_vector_add_int64_vector(list_int64_vector_ptr list,
    int64_vector_srcptr v)
{
  if ((list->length % DEFAULT_LENGTH_LIST_INT64_VECTOR) == 0 && list->length != 0) {
    list->v = realloc(list->v, sizeof(int64_vector_t) * (list->length +
          DEFAULT_LENGTH_LIST_INT64_VECTOR));
  }
  int64_vector_init(list->v[list->length], v->dim);
  int64_vector_set(list->v[list->length], v);
  list->length++;
}

void list_int64_vector_clear(list_int64_vector_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->dim != 0);
    int64_vector_clear(list->v[i]);
  }
  free(list->v);
  list->length = 0;
}

void list_int64_vector_fprintf(FILE * file, list_int64_vector_srcptr list)
{
  fprintf(file, "[");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      int64_vector_fprintf(file, list->v[i]);
    }
  }
  fprintf(file, "]\n");
}

void list_int64_vector_fprintf_comment(FILE * file,
    list_int64_vector_srcptr list)
{
  fprintf(file, "# [");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "# ");
      int64_vector_fprintf(file, list->v[i]);
    }
    fprintf(file, "# ");
  }
  fprintf(file, "]\n");
}

void list_int64_vector_extract_mat_int64(list_int64_vector_ptr list,
    mat_int64_srcptr matrix)
{
  ASSERT(list->length == 0);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, matrix->NumRows);

  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    mat_int64_extract_vector(v_tmp, matrix, i);
    list_int64_vector_add_int64_vector(list, v_tmp);
  }

  int64_vector_clear(v_tmp);
}

int int64_vector_in_list_int64_vector(int64_vector_srcptr vec, list_int64_vector_srcptr list)
{
  ASSERT(list->length > 2);
  ASSERT(vec->dim > 1);
#ifndef NDEBUG
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(vec->dim == list->v[i]->dim);
    for (unsigned int j = 2; j < vec->dim; j++) {
      ASSERT(vec->c[j] == list->v[i]->c[j]);
    }
  }
#endif // NDEBUG

  for (unsigned int i = 0; i < list->length; i++) {
    if (int64_vector_equal(vec, list->v[i]) == 1) {
      return 1;
    }  
  }

  int c = 0;
  unsigned i = 0;
  unsigned j = 0;
  for (i = 0, j = list->length - 1; i < list->length; j =
      i++) {
    if ( ((list->v[i]->c[1] > vec->c[1]) != (list->v[j]->c[1] > vec->c[1])) &&
        (vec->c[0] < (list->v[j]->c[0] - list->v[i]->c[0]) *
         (vec->c[1] - list->v[i]->c[1]) /
         (list->v[j]->c[1] - list->v[i]->c[1]) + list->v[i]->c[0]) ) {
      c = !c;
    }
  }

  return c;
}

static int compare_last(const void * p0, const void * p1)
{
  const int64_vector_t * v0 = (const int64_vector_t * ) p0;
  const int64_vector_t * v1 = (const int64_vector_t * ) p1;
  
  return (int)((*v0)->c[(*v0)->dim - 1] - (*v1)->c[(*v1)->dim - 1]);
}

void list_int64_vector_sort_last(list_int64_vector_ptr list)
{
  qsort(list->v, list->length, sizeof(list->v[0]), compare_last);
}
