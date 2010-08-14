#ifndef MF_H_
#define MF_H_

/* A few tools for reading matrices in either ascii or binary format, and
 * extracting info.
 */

#include <stdio.h>
#include <stdint.h>

/* structs of this type are passed, to decide what gets done with inputs.
 */
struct mf_io_file {
    uint64_t size;    // allow more than 4G elements.
    uint64_t alloc;   // allow more than 4G elements.
    uint32_t * p;
    FILE * f;
    int ascii;
};

#ifdef __cplusplus
extern "C" {
#endif

extern void matrix_read_pass(
        struct mf_io_file * m_in,
        struct mf_io_file * m_out,
        struct mf_io_file * rw_out,
        struct mf_io_file * cw_out,
        unsigned int rskip,
        unsigned int cskip,
        int progress);

extern char * build_mat_auxfile(const char * prefix, const char * what, const char * ext);

extern int has_suffix(const char * path, const char * sfx);
int matrix_autodetect_input(struct mf_io_file * mf_in, const char * mfile);

#ifdef __cplusplus
}
#endif

#endif	/* MF_H_ */
