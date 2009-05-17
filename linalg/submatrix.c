#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include "cado_config.h"
#include "macros.h"
#include "readmat.h"
#include "params.h"

const char * filename_in;
const char * filename_out;

sparse_mat_t mat;

void usage()
{
    fprintf(stderr, "Usage: submatrix --in <file> --out <file> <submatrix specs>\n"
            "where the submatrix may be specified as\n"
            "<i0> <j0> <ni> <nj>\n"
            "--anchor <i0>,<j0> --size <ni>,<nj>\n"
            "i0=<i0> i1=<i1> ni=<ni> nj=<nj>\n"
           );
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    argv++,argc--;

    int anchor[2]={-1,-1};
    int size[2]={-1,-1};
    int wild=0;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (wild < 2) {
            anchor[wild++]=atoi(argv[0]);
            argv++,argc--;
            continue;
        }
        if (wild < 4) {
            size[wild++-2]=atoi(argv[0]);
            argv++,argc--;
            continue;
        }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    filename_in = param_list_lookup_string(pl, "in");
    filename_out = param_list_lookup_string(pl, "out");


    param_list_parse_int_and_int(pl, "anchor", anchor, ",");
    param_list_parse_int_and_int(pl, "size", size, ",");
    param_list_parse_int(pl, "i0", &anchor[0]);
    param_list_parse_int(pl, "j0", &anchor[1]);
    param_list_parse_int(pl, "ni", &size[0]);
    param_list_parse_int(pl, "nj", &size[1]);

    if (anchor[0] < 0 || anchor[1] < 0 || size[0] == 0 || size[1] == 0) {
        usage();
        exit(1);
    }

    unsigned int i0 = anchor[0];
    unsigned int j0 = anchor[1];
    unsigned int i1 = i0 + size[0];
    unsigned int j1 = j0 + size[1];

    FILE * f_in;
    if (filename_in) {
        f_in = fopen(filename_in, "r");
        ASSERT_ALWAYS(f_in);
    } else {
        f_in = stdin;
    }

    FILE * f_out;
    if (filename_out) {
        f_out = fopen(filename_out, "w");
        ASSERT_ALWAYS(f_out);
    } else {
        f_out = stdout;
    }

    sparse_mat_init(mat);
    read_matrix_header(f_in, mat);

    if (i1 > mat->nrows || j1 > mat->nrows) {
        fprintf(stderr, "Warning: submatrix extends beyond "
                "actual matrix size.\n"
                "Result is padded with zeros\n");
    }

    fprintf(f_out, "%u %u\n", i1 - i0, j1 - j0);
    unsigned int i = 0;
    for( ; i < i0 ; i++) {
        read_matrix_row(f_in,mat,mat->data,1);
    }
    for( ; i < i1 && i < mat->nrows ; i++) {
        read_matrix_row(f_in,mat,mat->data,1);
        unsigned int l = mat->data[0];
        unsigned int * q = mat->data + 1;
        unsigned int k = 0;
        for( ; k < l && q[k] < j0 ; k++) ;
        unsigned int k0 = k;
        for( ; k < l && q[k] < j1 ; k++) ;
        fprintf(f_out, "%u", k-k0);
        for(k = k0 ; k < l && q[k] < j1 ; k++) {
            fprintf(f_out, " %u", q[k] - j0);
        }
        fprintf(f_out, "\n");
    }
    for( ; i < i1 ; i++) {
        fprintf(f_out, "0\n");
    }

    if (filename_in) fclose(f_in);
    if (filename_out) fclose(f_out);

    sparse_mat_clear(mat);
    param_list_clear(pl);
    return 0;
}
