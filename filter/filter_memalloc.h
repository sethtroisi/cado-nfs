#ifndef FILTER_MEMALLOC_H_
#define FILTER_MEMALLOC_H_

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each read relation is expensive, since 
   malloc() allocates some extra information to keep track of every memory 
   blocks. Instead, we allocate memory in big blocks of size BLOCK_SIZE. */


/* memory blocks are allocated of that # of index_t's */
#define BLOCK_SIZE (1<<20)

index_t * my_malloc (unsigned int);
void my_malloc_free_all ();
size_t get_my_malloc_bytes ();

#endif /* FILTER_MEMALLOC_H_ */
