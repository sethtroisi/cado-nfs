#define _DARWIN_C_SOURCE        /* asprintf */
#define _GNU_SOURCE     /* asprintf */
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "macros.h"
#include "mf.h"
#include "utils.h"

/* This our favorite way of expanding buffer. We always clear to zero */
static inline void expand(struct mf_io_file * m, uint64_t nn)
{
    m->p = realloc(m->p, nn * sizeof(uint32_t));
    if (nn > m->alloc)
        memset(m->p + m->alloc, 0, (nn - m->alloc) * sizeof(uint32_t));
    m->alloc = nn;
}

static void abort_unexpected_eof()
{
    fprintf(stderr, "Unexpected EOF -- (input file truncated ?).\n");
    exit(1);
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

char * build_mat_auxfile(const char * prefix, const char * what, const char * ext)
{
    /* This function gets called only if we have automatically decided on
     * the prefix, which means that normally, the decision has been taken
     * in accordance to the type of the matrix file.
     *
     * We check for a match between the elcted prefix and the matrix
     * filename. If found, that match is chopped.
     */
    char * dup_prefix;
    dup_prefix=strdup(prefix);

    if (has_suffix(dup_prefix, ext)) {
        dup_prefix[strlen(dup_prefix)-strlen(ext)]='\0';
    }
    char * str;
    int rc = asprintf(&str, "%s.%s%s", dup_prefix, what, ext);
    if (rc<0) abort();
    free(dup_prefix);
    return str;
}

extern void matrix_read_pass(
        struct mf_io_file * m_in,
        struct mf_io_file * m_out,
        struct mf_io_file * rw_out,
        struct mf_io_file * cw_out,
        unsigned int rskip,
        unsigned int cskip,
        int progress)
{
    /* This temp buffer is used only if skipping columns and outputting a
     * matrix file -- and even then, we can avoid it if the matrix rows are
     * known to be sorted (although we don't do so at the moment). */
    struct mf_io_file temp[1];
    memset(temp, 0, sizeof(struct mf_io_file));

    /* This gathers several possible ways of reading the matrix. It
     * encompasses the task of both mf_b_fq and mf_a2b_fq. It can also be used
     * to write anf mf_b2a tool, or to supersede the read_easy functions in
     * ../readmat-easy.c
     */
    ASSERT_ALWAYS(m_in);
    ASSERT_ALWAYS(m_in->p == NULL || m_in->f == NULL);
    ASSERT_ALWAYS(m_out == NULL || m_out->p == NULL || m_out->f == NULL);

    int rc;
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;

    /* In some cases we have to seek in the input file. If turns out to
     * be impossible, then avoid. */
    int can_seek = 1;

    unsigned int exp_nr=0;
    unsigned int exp_nc=0;
    int drop_cw_p = 0;
    if (m_in->f && m_in->ascii) {
        rc = fscanf(m_in->f, "%u %u", &exp_nr, &exp_nc);
        if (rc == EOF) abort_unexpected_eof();
        if (cw_out && cw_out->p == NULL)  {
            expand(cw_out, exp_nc - cskip);
            drop_cw_p=1;
        }
        // we don't even set cw_out->size.
        // cw_out->size = exp_nc;
    } else {
        if (cw_out && cw_out->p == NULL)  {
            expand(cw_out, 1000 * 1000);
            drop_cw_p=1;
        }
    }
    unsigned int i;

    uint32_t * ptr = NULL;
    if (!m_in->f) {
        ptr = m_in->p;
    }

    for(i = 0 ; ; i++) {
        uint32_t j;
        /* {{{ read row weight w */
        uint32_t w;
        if (m_in->f) {
            if (m_in->ascii) {
                rc = fscanf(m_in->f, "%"SCNu32, &w);
                if (i >= exp_nr) {
                    if (rc == EOF)
                        break;
                    if (i == exp_nr) {
                        fprintf(stderr,
                                "Warning: reading extra rows at end of file\n");
                    }
                }
                if (rc == EOF) abort_unexpected_eof();
            } else {
                rc = fread(&w, sizeof(uint32_t), 1, m_in->f);
                if (rc == 0 && feof(m_in->f)) break;
            }
        } else {
            if (((uint64_t) (ptr - m_in->p) >= m_in->size))
                break;
            w = *ptr++;
        }
        /* }}} */

        /* reading the row weight to the rw file is deferred to after the
         * point where we've cleaned the columns to be skipped. */

        if ((!m_out && !cw_out) || (i < rskip)) {
            /* {{{ skip fast over. */
            if (m_in->f) {
                if (!m_in->ascii) {
                    if (can_seek) {
                        rc = fseek(m_in->f, w * sizeof(uint32_t), SEEK_CUR);
                        if (rc != 0)
                            can_seek = 0;
                    }
                    if (!can_seek) {
                        for(j = 0 ; j < w ; j++) {
                            uint32_t c;
                            rc = fread(&c, sizeof(uint32_t), 1, m_in->f);
                            if (!rc) abort_unexpected_eof();
                        }
                    }
                } else {
                    /* TODO: do better */
                    for(j = 0 ; j < w ; j++) {
                        uint32_t c;
                        rc = fscanf(m_in->f, "%" SCNu32, &c);
                        if (rc == EOF) abort_unexpected_eof();
                    }
                }
            } else {
                ptr += w;
            }
            /* }}} */
            goto row_done;
        }

        uint32_t * ptr_w = NULL;

        if (m_out) { /* {{{ start row for matrix file */
            if (m_out->f) {
                if (!cskip) {
                    if (m_out->ascii) {
                        fprintf(m_out->f, "%" PRIu32, w);
                    } else {
                        rc = fwrite(&w, sizeof(uint32_t), 1, m_out->f);
                    }
                }
                /* Otherwise we have to defer */
            } else {
                if (m_out->size + w >= m_out->alloc) {
                    expand(m_out, (1 << 20) + (m_out->size+w) * 3/2);
                }
                ptr_w = &(m_out->p[m_out->size]);
                *ptr_w = w;
                m_out->size++;
            }
        } /* }}} */

        uint32_t real_weight = 0;

        if (m_out && cskip && temp->alloc < w) {
            expand(temp, w + w/2);
        }

        for(j = 0 ; j < w ; j++) {
            /* {{{ read this column index */
            uint32_t c;
            if (m_in->f) {
                if (m_in->ascii) {
                    rc = fscanf(m_in->f, "%" SCNu32, &c);
                    if (rc == EOF) abort_unexpected_eof();
                } else {
                    rc = fread(&c, sizeof(uint32_t), 1, m_in->f);
                    if (!rc) abort_unexpected_eof();
                }
            } else {
                c = *ptr++;
            }
            /* }}} */
            if (cskip) {
                if (c < cskip) {
                    continue;
                } else {
                    if (m_out) {
                        temp->p[real_weight] = c - cskip;
                    }
                    real_weight++;
                }
            } else {
                real_weight++;
            }

            if (m_out) { /* {{{ write to matrix file */
                if (m_out->f) {
                    if (!cskip) {
                        if (m_out->ascii) {
                            fprintf(m_out->f, " %"PRIu32, c);
                        } else {
                            rc = fwrite(&c, sizeof(uint32_t), 1, m_out->f);
                        }
                    }
                    // otherwise it's too early, size hasn't been written yet.
                } else {
                    m_out->p[m_out->size++] = c-cskip;
                }
            }/*}}}*/

            if (exp_nc && c >= exp_nc) {
                fprintf(stderr, "Warning: column index %"PRIu32" exceeds header value %u\n", c, exp_nc);
            }
            if (cw_out) {
                if (c-cskip >= cw_out->alloc) {
                    expand(cw_out, c-cskip + cw_out->alloc / 2);
                }
                if (c-cskip >= cw_out->size) {
                    cw_out->size = c-cskip + 1;
                }
                cw_out->p[c-cskip]++;
            }
        }
        w = real_weight;

        if (m_out && cskip) {/* {{{ start row for matrix file */
            if (m_out->f) {
                if (m_out->ascii) {
                    fprintf(m_out->f, "%" PRIu32, w);
                    for(j = 0 ; j < w ; j++) {
                        fprintf(m_out->f, " %"PRIu32, temp->p[j]);
                    }
                } else {
                    rc = fwrite(&w, sizeof(uint32_t), 1, m_out->f);
                    rc = fwrite(temp->p, sizeof(uint32_t), w, m_out->f);
                }
            } else {
                // fix this up. It's the only thing we have to do.
                *ptr_w = w;
            }
        } /* }}} */

        /* row done */
        if (m_out && m_out->f && m_out->ascii)
            fprintf(m_out->f, "\n");

        if (rw_out) { /* {{{ write row weight to rw file */
            if (rw_out->p) rw_out->p[i]=w;
            if (rw_out->f) {
                if (rw_out->ascii) {
                    fprintf(rw_out->f, "%"PRIu32"\n", w);
                } else {
                    rc = fwrite(&w, sizeof(uint32_t), 1, rw_out->f);
                }
            }
        }/*}}}*/

row_done:
        if (progress) {/*{{{*/
            time_t t = time(NULL);
            if (t >= t1) {
                t1 = t + 1;
                int dt = t - t0;
                size_t sz = ftell(m_in->f) >> 20;
                fprintf(stderr, "read %zu MB in %d s (%.1f MB/s)   \r",
                        sz,
                        dt,
                        (double) sz / dt);
            }
        }/*}}}*/
    }
    if (progress)
        fprintf(stderr, "\n");

    /* write cw if requested *//*{{{*/
    if (cw_out) {
        if (cw_out->f) {
            if (cw_out->ascii) {
                uint32_t j;
                for(j = 0 ; j < cw_out->size ; j++) {
                    fprintf(cw_out->f, "%"PRIu32"\n", cw_out->p[j]);
                }
            } else {
                rc = fwrite(cw_out->p, sizeof(uint32_t), cw_out->size, cw_out->f);
            }
        }
    }
    if (drop_cw_p) {
        free(cw_out->p);
        cw_out->p = NULL;
    }/*}}}*/

    if (exp_nr && i != exp_nr) {
        fprintf(stderr,
                "Warning: counted %u rows,"
                " not according with header value %u\n", i, exp_nr);
    }
}

void balancing_compute_checksum(balancing_ptr bal)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = 0;
    if (bal->h->flags & FLAG_ROWPERM) {
        size_t s = bal->h->nrows;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        w = cado_crc_lfsr_turn(l, bal->rowperm, s * sizeof(uint32_t));
    }
    if (bal->h->flags & FLAG_COLPERM) {
        size_t s = bal->h->ncols;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        w = cado_crc_lfsr_turn(l, bal->colperm, s * sizeof(uint32_t));
    }
    cado_crc_lfsr_clear(l);
    bal->h->checksum = w;
    if (bal->h->flags & FLAG_REPLICATE) {
        // a trick to identify conjugated perms.
        bal->h->checksum &= ~0xff;
    }
}

void balancing_write_namefile(balancing_ptr bal, const char * mfile)
{
    char * dup_prefix=strdup(mfile);
    char * ext[2] = { ".bin", ".txt" };
    for(int j = 0 ; j < 2 ; j++) {
        if (has_suffix(dup_prefix, ext[j])) {
            dup_prefix[strlen(dup_prefix)-strlen(ext[j])]='\0';
            break;
        }
    }
    char * filename;
    int rc = asprintf(&filename, "%s.%dx%d.%08"PRIx32".bin",
            dup_prefix, bal->h->nh, bal->h->nv, bal->h->checksum);
    ASSERT_ALWAYS(rc >= 0);
    free(dup_prefix);
    balancing_write(bal, filename);
    free(filename);
}

void balancing_write(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    pfile = fopen(filename, "w");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    int rc = fwrite(&bal->h, sizeof(struct balancing_header_s), 1, pfile);
    if (bal->h->flags & FLAG_ROWPERM) {
        size_t s = bal->h->nrows;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        rc = fwrite(bal->rowperm, sizeof(uint32_t), s, pfile);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        size_t s = bal->h->ncols;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        rc = fwrite(bal->colperm, sizeof(uint32_t), s, pfile);
    }
    fclose(pfile);
}

void balancing_read(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    pfile = fopen(filename, "r");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    int rc = fread(&bal->h, sizeof(struct balancing_header_s), 1, pfile);
    if (bal->h->flags & FLAG_ROWPERM) {
        size_t s = bal->h->nrows;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        bal->rowperm = malloc(s * sizeof(uint32_t));
        rc = fread(bal->rowperm, sizeof(uint32_t), s, pfile);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        size_t s = bal->h->ncols;
        if (bal->h->flags & FLAG_PADDING) {
            s = MAX(bal->h->nrows, bal->h->ncols);
        }
        bal->colperm = malloc(s * sizeof(uint32_t));
        rc = fread(bal->colperm, sizeof(uint32_t), s, pfile);
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
