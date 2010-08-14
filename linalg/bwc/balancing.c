#define _POSIX_C_SOURCE 200112L
#define _BSD_SOURCE     /* strdup */
#define _GNU_SOURCE		/* asprintf */
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>

#include "balancing.h"
#include "utils.h"

void balancing_finalize(balancing_ptr bal)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = 0;
    bal->trows = bal->h->nrows;
    bal->tcols = bal->h->ncols;
    if (bal->h->flags & FLAG_PADDING) {
        bal->tcols = bal->trows = MAX(bal->h->nrows, bal->h->ncols);
    }
    if (bal->h->flags & FLAG_ROWPERM) {
        w = cado_crc_lfsr_turn(l, bal->rowperm, bal->trows * sizeof(uint32_t));
    }
    if (bal->h->flags & FLAG_COLPERM) {
        w = cado_crc_lfsr_turn(l, bal->colperm, bal->tcols * sizeof(uint32_t));
    }
    cado_crc_lfsr_clear(l);
    bal->h->checksum = w;
    if (bal->h->flags & FLAG_REPLICATE) {
        // a trick to identify conjugated perms.
        bal->h->checksum &= ~0xff;
    }
}

void balancing_write_inner(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    fprintf(stderr, "Writing balancing data to %s\n", filename);
    pfile = fopen(filename, "w");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    int rc = fwrite(&bal->h, sizeof(struct balancing_header_s), 1, pfile);
    if (bal->h->flags & FLAG_ROWPERM) {
        rc = fwrite(bal->rowperm, sizeof(uint32_t), bal->trows, pfile);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        rc = fwrite(bal->colperm, sizeof(uint32_t), bal->tcols, pfile);
    }
    fclose(pfile);
}

void balancing_write(balancing_ptr bal, const char * mfile, const char * suggest)
{
    if (suggest && strlen(suggest) && suggest[strlen(suggest)-1] != '/') {
        balancing_write_inner(bal, suggest);
    }

    char * dup_prefix=strdup(mfile);
    char * ext[2] = { ".bin", ".txt" };
    for(int j = 0 ; j < 2 ; j++) {
        if (has_suffix(dup_prefix, ext[j])) {
            dup_prefix[strlen(dup_prefix)-strlen(ext[j])]='\0';
            break;
        }
    }

    /* If we've been suggested an output directory, use it ! */
    const char * q = dup_prefix; 
    const char * d = "";
    if (suggest && strlen(suggest)) {
        d=suggest;
        char * last_slash = strrchr(q, '/');
        if (last_slash) {
            q = last_slash + 1;
        }
    }
    char * filename;
    int rc = asprintf(&filename, "%s%s.%dx%d.%08"PRIx32".bin",
            d,q, bal->h->nh, bal->h->nv, bal->h->checksum);

    ASSERT_ALWAYS(rc >= 0);
    free(dup_prefix);
    balancing_write_inner(bal, filename);
    free(filename);
}

void balancing_read_header_inner(balancing_ptr bal, FILE * pfile)
{
    int rc = fread(&bal->h, sizeof(struct balancing_header_s), 1, pfile);
    ASSERT_ALWAYS(rc == 1);
    bal->trows = bal->h->nrows;
    bal->tcols = bal->h->ncols;
    if (bal->h->flags & FLAG_PADDING) {
        bal->tcols = bal->trows = MAX(bal->h->nrows, bal->h->ncols);
    }
}

void balancing_read_header(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    char * derived = derived_filename(filename, "hdr", ".bin");
    ASSERT_ALWAYS(filename);
    pfile = fopen(filename, "r");
    if (pfile == NULL) {
        pfile = fopen(derived, "r");
        if (pfile == NULL) {
            fprintf(stderr, "Cannot read %s nor %s: %s\n",
                    filename, derived, strerror(errno));
            abort();
        }
    }
    balancing_read_header_inner(bal, pfile);
    fclose(pfile);
    free(derived);
}

void balancing_read(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    ASSERT_ALWAYS(filename);
    pfile = fopen(filename, "r");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    balancing_read_header_inner(bal, pfile);
    if (bal->h->flags & FLAG_PADDING) {
        bal->tcols = bal->trows = MAX(bal->h->nrows, bal->h->ncols);
    }
    if (bal->h->flags & FLAG_ROWPERM) {
        bal->rowperm = malloc(bal->trows * sizeof(uint32_t));
        int rc = fread(bal->rowperm, sizeof(uint32_t), bal->trows, pfile);
        ASSERT_ALWAYS(rc == (int) bal->trows);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        bal->colperm = malloc(bal->tcols * sizeof(uint32_t));
        int rc = fread(bal->colperm, sizeof(uint32_t), bal->tcols, pfile);
        ASSERT_ALWAYS(rc == (int) bal->tcols);
    }
    fclose(pfile);
}

void balancing_init(balancing_ptr bal)
{
    memset(bal, 0, sizeof(balancing));
}
void balancing_clear(balancing_ptr bal)
{
    if (bal->colperm) free(bal->colperm);
    if (bal->rowperm) free(bal->rowperm);
    memset(bal, 0, sizeof(balancing));
}
