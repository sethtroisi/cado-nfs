#ifndef BW_LVBLOCK_H
#define BW_LVBLOCK_H
/* Long Vector Blocks */

#ifdef	__cplusplus
extern "C" {
#endif

typedef mp_limb_t * bw_vector_block;

#if 0
OLD void _bw_lvblock_alloc_n(bw_vector_block *, coord_t);
OLD void bw_lvblock_set_zero_n(bw_vector_block, coord_t);
OLD #define bw_lvblock_alloc_n(x,n) _bw_lvblock_alloc_n(&(x),n)
OLD size_t bw_lvblock_read_1(bw_vector_block, FILE *);
OLD size_t bw_lvblock_write_1(FILE *, bw_vector_block);
OLD #define bw_lvblock_alloc(x) _bw_lvblock_alloc_n(&(x),ncols)
#define bw_lvblock_complete_reduce(v) bw_lvblock_is_zero(v)
int bw_lvblock_add(bw_vector_block,bw_vector_block);

#endif

void _bw_lvblock_alloc(bw_vector_block *);
#define bw_lvblock_alloc(x) _bw_lvblock_alloc(&(x))

void _bw_lvblock_free(bw_vector_block *);
#define bw_lvblock_free(x) _bw_lvblock_free(&(x))

void bw_lvblock_set_zero_separated(bw_vector_block);
void bw_lvblock_set_zero_full(bw_vector_block);

size_t bw_lvblock_read(bw_vector_block, FILE *);
size_t bw_lvblock_write(FILE *, bw_vector_block);

void bw_lvblock_reduce_separated(bw_vector_block);

void bw_multiproduct_separated(bw_vector_block, bw_vector_block, bw_vector_block);

int bw_lvblock_is_zero_separated(bw_vector_block);

#ifdef	__cplusplus
}
#endif

#endif /* BW_LVBLOCK_H */
