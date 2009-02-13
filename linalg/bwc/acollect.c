#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include "macros.h"
#include "manu.h"
#include "utils.h"
#include "filenames.h"
#include "bw-common.h"

/* This program is rather standalone. It checks the current directory
 * for files matching the pattern A%u-%u.%u-%u, and concatenates them.
 *
 * It needs to know m, though (auto-detection would be easy).
 */

void usage()
{
    fprintf(stderr, "Usage: acollect [m=<m>] [wdir=<path>] [--remove-old]\n");
    exit(1);
}

struct bw_params bw[1];

int remove_old = 0;

struct afile_s {
    unsigned int n0, n1, j0, j1;
};
typedef struct afile_s afile[1];
typedef struct afile_s * afile_ptr;


typedef int (*sortfunc_t) (const void *, const void *);

int afile_cmp(afile_ptr a, afile_ptr b)
{
    int dj0 = a->j0 - b->j0;
    int dj1 = a->j1 - b->j1;
    int dn0 = a->n0 - b->n0;
    int dn1 = a->n1 - b->n1;

    if (dj0) return dj0;
    if (dj1) return dj1;
    if (dn0) return dn0;
    if (dn1) return dn1;
    return 0;
}

int main(int argc, char * argv[])
{
    int rc = 0;

    param_list pl;
    param_list_init(pl);
    param_list_configure_knob(pl, "--remove-old", &remove_old);
    bw_common_init(bw, pl, argc, argv);
    param_list_clear(pl);


    afile * afiles = NULL;
    int n_afiles = 0;
    int alloc_afiles=1024;
    afiles = realloc(afiles, alloc_afiles * sizeof(afile));
    DIR * dir = opendir(".");
    struct dirent * de;

    for( ; (de = readdir(dir)) != NULL ; ) {
        afile_ptr a = afiles[n_afiles];
        int k;
        int rc = sscanf(de->d_name, A_FILE_PATTERN "%n", &a->n0, &a->n1, &a->j0, &a->j1, &k);
        /* rc is expected to be 4 or 5 depending on our reading of the
         * standard */
        if (rc < 4 || k != (int) strlen(de->d_name)) {
            // fprintf(stderr, "skipped %s\n", de->d_name);
        } else {
            if (a->n1 % CHAR_BIT || a->n0 % CHAR_BIT) {
                fprintf(stderr, "%s has bad boundaries\n",
                        de->d_name);
                exit(1);
            }
            struct stat sbuf[1];
            rc = stat(de->d_name, sbuf);
            if (rc < 0) {
                fprintf(stderr, "stat: %s\n", strerror(errno));
                exit(1);
            }
            ssize_t expected = bw->m * (a->n1-a->n0) / CHAR_BIT * (a->j1 - a->j0);

            if (sbuf->st_size != expected) {
                fprintf(stderr, "%s does not have expected size %zu\n",
                        de->d_name, expected);
                exit(1);
            }

            // fprintf(stderr, "take %s\n", de->d_name);
            n_afiles++;
            if (n_afiles >= alloc_afiles) {
                alloc_afiles += alloc_afiles / 4;
                afiles = realloc(afiles, alloc_afiles * sizeof(afile));
            }
        }
    }
    closedir(dir);
    
    if (n_afiles == 0) {
        free(afiles);
        return 0;
    }

    if (n_afiles == 1) {
        char * tmp;
        rc = asprintf(&tmp, A_FILE_PATTERN,
                afiles[0]->n0,afiles[0]->n1,afiles[0]->j0,afiles[0]->j1);
        printf("%s\n", tmp);
        free(tmp);
        free(afiles);
        return 0;
    }

    qsort(afiles, n_afiles, sizeof(afile), (sortfunc_t) &afile_cmp);

    afile final;
    /* start with unset values */
    final->n0 = UINT_MAX;
    final->n1 = UINT_MAX;
    final->j0 = UINT_MAX;
    final->j1 = UINT_MAX;
    FILE * f = fopen("A.temp", "w");

    int k0 = 0;
    int k1;

    for( ; k0 < n_afiles ; ) {
        unsigned int j0 = afiles[k0]->j0;
        unsigned int j1 = afiles[k0]->j1;
        k1 = k0;
        unsigned int n0 = afiles[k0]->n0;
        unsigned int n1 = afiles[k0]->n1;
        unsigned int n = n0;
        for( ; k1 < n_afiles && afiles[k1]->j0 == j0 ; k1++) {
            if (afiles[k1]->j1 != j1 || afiles[k1]->n0 != n) {
                fprintf(stderr, "Found inconsistent files A%u-%u.%u-%u and A%u-%u.%u-%u\n",
                        afiles[k0]->n0, afiles[k0]->n1, afiles[k0]->j0, afiles[k0]->j1,
                        afiles[k1]->n0, afiles[k1]->n1, afiles[k1]->j0, afiles[k1]->j1);
                exit(1);
            }
            n = afiles[k1]->n1;
        }
        
        /* Files [k0] to [k1] can be collected */

        if (    (final->n0 != UINT_MAX && n0 != final->n0) ||
                (final->n1 != UINT_MAX && n1 != final->n1) ||
                (final->j1 != UINT_MAX && j0 != final->j1))
        {
            fprintf(stderr, "Cannot append A%u-%u.%u-%u to A%u-%u.%u-%u\n",
                    n0,n1,j0,j1,
                    final->n0,final->n1,final->j0,final->j1);
            rc = 1;
            break;
        }

        FILE ** rs = malloc((k1-k0) * sizeof(FILE *));
        int k;
        for(k = k0 ; k < k1 ; k++) {
            char * tmp;
            int rc = asprintf(&tmp, A_FILE_PATTERN,
                    afiles[k]->n0,afiles[k]->n1,afiles[k]->j0,afiles[k]->j1);
            rs[k - k0] = fopen(tmp, "r");
            free(tmp);
            if (rs[k-k0] == NULL) {
                fprintf(stderr, "fopen: %s\n", strerror(errno));
                rc = 2;
                goto bailout;
            }
        }

        char * buf = malloc((n1-n0)/CHAR_BIT);

        final->j1 = j0;

        for(unsigned int j = j0 ; j < j1 ; j++) {
            for(int i = 0 ; i < bw->m ; i++) {
                char * ptr = buf;
                size_t rz;
                size_t sz;
                for(int k = k0 ; k < k1 ; k++) {
                    sz = (afiles[k]->n1 - afiles[k]->n0) / CHAR_BIT;
                    rz = fread(ptr, 1, sz, rs[k-k0]);
                    if (rz < sz) {
                        fprintf(stderr, "fread: short read\n");
                        rc = 2;
                        fflush(f);
                        goto bailout;
                    }
                    ptr += sz;
                }
                sz = (n1-n0)/CHAR_BIT;
                rz = fwrite(buf, 1, sz, f);
                if (rz != sz) {
                    fprintf(stderr, "fwrite: short write\n");
                    rc = 2;
                    fflush(f);
                    goto bailout;
                }
            }
            final->j1++;
        }
        free(buf);

        if (final->j0 == UINT_MAX) { final->j0 = j0; }
        if (final->n0 == UINT_MAX) { final->n0 = n0; }
        if (final->n1 == UINT_MAX) { final->n1 = n1; }

        if (fflush(f) != 0) {
            // we're in trouble
            fprintf(stderr, "fflush(): %s\n", strerror(errno));
            rc = 2;
            break;
        }

        final->j1 = j1;

        for(int k = k0 ; k < k1 ; k++) {
            fclose(rs[k-k0]);
            if (!remove_old) continue;
            char * tmp;
            int rc = asprintf(&tmp, A_FILE_PATTERN,
                    afiles[k]->n0,afiles[k]->n1,afiles[k]->j0,afiles[k]->j1);
            BUG_ON(rc < 0);     // shut up, dammit.
            if (unlink(tmp) < 0) {
                fprintf(stderr, "unlink: %s\n", strerror(errno));
            }
            free(tmp);
        }
        free(rs);

        k0 = k1;
    }
bailout:
    fclose(f);
    if (final->j0 != UINT_MAX) {
        char * tmp;
        int r;
        r = asprintf(&tmp, A_FILE_PATTERN, final->n0,final->n1,final->j0,final->j1);
        r = rename("A.temp", tmp);
        if (r < 0) {
            fprintf(stderr, "rename: %s\n", strerror(errno));
        }
        printf("%s\n",tmp);
        free(tmp);
    }
    free(afiles);
    return rc;
}
