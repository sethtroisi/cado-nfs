#ifndef CLIQUE_REMOVAL_H_
#define CLIQUE_REMOVAL_H_

/* A clique is a connected component of the graph where the nodes are the rows
   and the edges are the columns of weight 2.
   /!\ It is not a clique is the sense of graph theory.
   We will try to use the "connected component" terminology instead.
*/

/********************* comp_t struct (clique) ********************************/

typedef struct {
  uint64_t i; /* smallest row appearing in the connected component */
  float w;   /* Weight of the connected component */
} comp_t;

void comp_print_info_weight_function ();

/******************** uint64_buffer struct ************************************/

/* Double buffer:
 * |---done----|---todo----|---free----|
 * | | | | | | | | | | | | | | | | | | |
 *  ^           ^           ^           ^
 *  begin       next_todo   next_free   end
 */
struct uint64_buffer_s {
  uint64_t *begin, *next_todo, *next_free, *end;
};
typedef struct uint64_buffer_s uint64_buffer_t[1];
typedef struct uint64_buffer_s * uint64_buffer_ptr;
typedef const struct uint64_buffer_s * uint64_buffer_srcptr;

#define UINT64_BUFFER_MIN_SIZE 32

void uint64_buffer_init (uint64_buffer_ptr, size_t);
void uint64_buffer_clear (uint64_buffer_ptr buf);

/********************** main functions ****************************************/
uint64_t compute_one_connected_component (comp_t *, purge_matrix_srcptr,
                                          uint64_buffer_ptr);
void delete_one_connected_component (purge_matrix_ptr, uint64_t,
                                     uint64_buffer_ptr);

void cliques_removal (purge_matrix_ptr, int64_t, unsigned int, int);

#endif /* CLIQUE_REMOVAL_H_ */

