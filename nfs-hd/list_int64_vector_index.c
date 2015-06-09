#include "cado.h"
#include "utils.h"
#include "list_int64_vector_index.h"

void list_int64_vector_index_init(list_int64_vector_index_ptr list)
{
  list->length = 0;
  list->v = (int64_vector_index_t * ) malloc(sizeof(int64_vector_index_t) *
      DEFAULT_LENGTH_LIST_INT64_VECTOR_INDEX);
}

void list_int64_vector_index_add_int64_vector_index(list_int64_vector_index_ptr list,
    int64_vector_srcptr v, uint64_t index)
{
  if ((list->length % DEFAULT_LENGTH_LIST_INT64_VECTOR_INDEX) == 0 && list->length != 0) {
    list->v = realloc(list->v, sizeof(int64_vector_index_t) * (list->length +
          DEFAULT_LENGTH_LIST_INT64_VECTOR_INDEX));
  }
  int64_vector_init(list->v[list->length]->vec, v->dim);
  int64_vector_set(list->v[list->length]->vec, v);
  list->v[list->length]->index = index;
  list->length++;
}

void list_int64_vector_index_clear(list_int64_vector_index_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->vec->dim != 0);
    int64_vector_clear(list->v[i]->vec);
  }
  free(list->v);
  list->length = 0;
}

void list_int64_vector_index_fprintf(FILE * file,
    list_int64_vector_index_srcptr list)
{
  fprintf(file, "[");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "%" PRIu64 ": ", list->v[i]->index);
      int64_vector_fprintf(file, list->v[i]->vec);
    }
  }
  fprintf(file, "]\n");
}

void list_int64_vector_index_fprintf_comment(FILE * file,
    list_int64_vector_index_srcptr list)
{
  fprintf(file, "# [");
  if (list->length != 0) {
    fprintf(file, "\n");
    for (unsigned int i = 0; i < list->length; i++) {
      fprintf(file, "# %" PRIu64 ": ", list->v[i]->index);
      int64_vector_fprintf(file, list->v[i]->vec);
    }
    fprintf(file, "# ");
  }
  fprintf(file, "]\n");
}

static int compare_vector_last(const void * p0, const void * p1)
{
  const int64_vector_index_t * v0 = (const int64_vector_index_t * ) p0;
  const int64_vector_index_t * v1 = (const int64_vector_index_t * ) p1;
  
  return (int)((*v0)->vec->c[(*v0)->vec->dim - 1] -
      (*v1)->vec->c[(*v1)->vec->dim - 1]);
}

void list_int64_vector_index_sort_last(list_int64_vector_index_ptr list)
{
  qsort(list->v, list->length, sizeof(list->v[0]), compare_vector_last);
}
