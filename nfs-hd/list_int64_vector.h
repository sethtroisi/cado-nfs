#ifndef LIST_INT64_VECTOR_H
#define LIST_INT64_VECTOR_H 

#include <stdio.h>
#include <stdint.h>
#include "cado.h"
#include "utils.h"
#include "int64_vector.h"
#include "mat_int64.h"
#include "list.h"

/*
 * Initialise a list of int64_vector.
 * 
 * list: the list.
 */
void list_int64_vector_init(list_int64_vector_ptr list,
    unsigned int vector_dim);

/*
 * Add a vector in the list.
 *
 * list: the list of vector.
 * v: the vector we add.
 */
unsigned int list_int64_vector_add_int64_vector(list_int64_vector_ptr list,
    int64_vector_srcptr v);

/*
 * Delete the list.
 *
 * list: the list.
 */
void list_int64_vector_clear(list_int64_vector_ptr list);

void list_int64_vector_delete_elements(list_int64_vector_ptr list);
/*
 * Write a list in a file.
 *
 * file: the file in which we write.
 * lest: the list.
 */
void list_int64_vector_fprintf(FILE * file, list_int64_vector_srcptr list);

void list_int64_vector_fprintf_comment(FILE * file,
    list_int64_vector_srcptr list);
/*
 * Extract all the vector (in column) of the matrix and put them in a list.
 *
 * list: the list of vector.
 * matrix: the matrix.
 */
void list_int64_vector_extract_mat_int64(list_int64_vector_ptr list,
    mat_int64_srcptr matrix);

/*
 * From PNpoly
 */
int int64_vector_in_polytop_list_int64_vector(int64_vector_srcptr vec, list_int64_vector_srcptr list);

void list_int64_vector_sort_last(list_int64_vector_ptr list);

#endif /* LIST_INT64_VECTOR_H */
