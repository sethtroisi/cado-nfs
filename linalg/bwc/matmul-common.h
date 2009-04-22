#ifndef MATMUL_COMMON_H_
#define MATMUL_COMMON_H_

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include "params.h"
#include "matmul.h"

#ifdef __cplusplus
extern "C" {
#endif

FILE * matmul_common_reload_cache_fopen(size_t, struct matmul_public_s * mm, const char * filename, const char * ext, uint32_t magic);
FILE * matmul_common_save_cache_fopen(size_t, struct matmul_public_s * mm, const char * filename, const char * ext, uint32_t magic);
void matmul_common_init_post(struct matmul_public_s * mm, param_list pl, int suggest);
uint32_t * matmul_common_read_stupid_data(struct matmul_public_s * mm, const char * filename);

extern const char * rowcol[2];  // [0] = "row" [1] = "col"

#ifdef __cplusplus
}
#endif

/* All matmul implementation are peppered with such markers */
#ifndef ASM_COMMENT
#ifdef  __GNUC__
#define ASM_COMMENT(x)  __asm__("#\t" x "\n")
#else
#define ASM_COMMENT(x)  /**/
#endif
#endif

/* I/O with cache files is made easier with these macros ; rather than
 * having to check for errors over and over again... */
#define MATMUL_COMMON_READ_ONE64(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint64_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    final_v__ = storage_v__;                                            \
} while (0)
#define MATMUL_COMMON_READ_ONE32(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint32_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    final_v__ = storage_v__;                                            \
} while (0)
#define MATMUL_COMMON_READ_ONE16(final_v__, file__)  do {               \
    size_t rc;                                                          \
    uint16_t storage_v__;                                               \
    rc = fread(&storage_v__, sizeof(storage_v__), 1, file__);           \
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");   \
    final_v__ = storage_v__;                                            \
} while (0)
#define MATMUL_COMMON_READ_MANY32(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fread(ptr__, sizeof(uint32_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short read from cached matrix file");	\
} while (0)
#define MATMUL_COMMON_READ_MANY16(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fread(ptr__, sizeof(uint16_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short read from cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_ONE64(final_v__, file__)  do {              \
    size_t rc;								\
    uint64_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_ONE32(final_v__, file__)  do {              \
    size_t rc;								\
    uint32_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_ONE16(final_v__, file__)  do {              \
    size_t rc;								\
    uint16_t storage_v__ = final_v__;					\
    rc = fwrite(&storage_v__, sizeof(storage_v__), 1, file__);		\
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_MANY32(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fwrite(ptr__, sizeof(uint32_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short write to cached matrix file");	\
} while (0)
#define MATMUL_COMMON_WRITE_MANY16(ptr__, n__, f__) do {         	\
    size_t rc;								\
    rc = fwrite(ptr__, sizeof(uint16_t), n__, f__);			\
    FATAL_ERROR_CHECK(rc < n__, "Short write to cached matrix file");	\
} while (0)


#endif	/* MATMUL_COMMON_H_ */
