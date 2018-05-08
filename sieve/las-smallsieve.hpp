#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <stdarg.h>
#include "las-smallsieve-types.hpp"
#include "las-types.hpp"
#include "bucket.hpp"

extern void small_sieve_info(const char * what, int side, small_sieve_data_t const & r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t & ssd);
extern void small_sieve_init(small_sieve_data_t & ssd,
                      small_sieve_block_offsets_t & ssd_offsets,
                      std::vector<fb_entry_general>::const_iterator fb_start,
                      std::vector<fb_entry_general>::const_iterator fb_end,
                      std::vector<fb_entry_general>::const_iterator resieve_start,
                      std::vector<fb_entry_general>::const_iterator resieve_end,
                      sieve_info const & si, int side);
extern void small_sieve_start(std::vector<spos_t> & ssdpos, small_sieve_data_t & ssd, unsigned int first_region_index, sieve_info const & si);

extern void small_sieve_start_many(
        std::vector<std::vector<spos_t>> & ssdpos_many,
        small_sieve_data_t & ssd,
        small_sieve_block_offsets_t & ssd_offsets,
        unsigned int first_region_index,
        sieve_info const & si);

extern void small_sieve_copy_start(std::vector<spos_t>& res, std::vector<spos_t> const & base, int bounds[2]);
/*
extern void small_sieve_skip_stride(small_sieve_data_t *ssd, spos_t * ssdpos, unsigned int N, unsigned int interleaving, sieve_info const & si);
*/
extern void sieve_small_bucket_region(unsigned char *S, unsigned int N,
                               small_sieve_data_t const & ssd,
                               std::vector<spos_t> const & ssdpos,
                               sieve_info & si, int side MAYBE_UNUSED,
                               where_am_I & w);

extern void resieve_small_bucket_region (bucket_primes_t *BP,
        unsigned char *S,
        unsigned int N,
        small_sieve_data_t & ssd, std::vector<spos_t> const & ssdpos,
        sieve_info const & si, where_am_I & w MAYBE_UNUSED);

extern void precompute_small_sieve_init_pos(small_sieve_data_t & ssd,
                               std::vector<std::vector<spos_t>> & ssdpos_per_bucket,
                               sieve_info & si,
                               where_am_I & w);



#endif	/* LAS_SMALLSIEVE_HPP_ */
