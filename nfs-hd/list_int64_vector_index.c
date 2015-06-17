#include "cado.h"
#include "utils.h"
#include "list_int64_vector_index.h"

void list_int64_vector_index_init(list_int64_vector_index_ptr list,
    unsigned int vector_dim)
{
  memset(list, 0, sizeof(*list));
  list->vector_dim = vector_dim;
}

static void list_int64_vector_index_prepare_write(list_int64_vector_index_ptr list,
    unsigned int index)
{
  if (index >= list->alloc) {
    list->alloc = index + 4 + list->alloc / 4;
    list->v = (int64_vector_index_t * ) realloc(list->v, list->alloc *
        sizeof(int64_vector_index_t));
    for (unsigned int i = list->length; i < list->alloc; i++) {
      int64_vector_init(list->v[i]->vec, list->vector_dim);
      list->v[i]->index = 0;
    }
  }
  if (list->length <= index) {
    list->length = index + 1;
  }
}

void list_int64_vector_index_add_int64_vector_index(
    list_int64_vector_index_ptr list, int64_vector_srcptr v, uint64_t index)
{
  ASSERT(v->dim == list->vector_dim);

  list_int64_vector_index_prepare_write(list, list->length);
  int64_vector_set(list->v[list->length - 1]->vec, v);
  list->v[list->length - 1]->index = index;
}

void list_int64_vector_index_clear(list_int64_vector_index_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->vec->dim != 0);
    int64_vector_clear(list->v[i]->vec);
  }
  free(list->v);
  memset(list, 0, sizeof(*list));
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
