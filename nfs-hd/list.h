#ifndef LIST_H
#define LIST_H 

#include "int64_vector.h"
#include "double_vector.h"

/*
 * List of int64_vector.
 */
typedef struct
{
  int64_vector_t * v;
  unsigned int length;
}  s_list_int64_vector_t;

typedef s_list_int64_vector_t list_int64_vector_t[1];
typedef s_list_int64_vector_t * list_int64_vector_ptr;
typedef const s_list_int64_vector_t * list_int64_vector_srcptr;

/*
 * List of double vector.
 */
typedef struct
{
  double_vector_t * v;
  unsigned int length;
}  s_list_double_vector_t;

typedef s_list_double_vector_t list_double_vector_t[1];
typedef s_list_double_vector_t * list_double_vector_ptr;
typedef const s_list_double_vector_t * list_double_vector_srcptr;

#endif /* LIST_H */
