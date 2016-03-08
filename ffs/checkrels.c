#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <ctype.h>
#include "macros.h"

#include "types.h"
#include "ffspol.h"
#include "params.h"
#include "norm.h"

#define MAX_LINE_SIZE 2048

const char * curr_line;

void fppol_make_monic(fppol_ptr p)
{
    if (!fppol_is_monic(p)) {
        fp_t lc;
        fppol_get_coeff(lc, p, fppol_deg(p));
        fppol_sdiv(p, p, lc);
    }
}
 
const char * read_ffpol_and_advance(fppol_ptr a, const char *str, uint64_t nl)
{
    char s[MAX_LINE_SIZE];
    char *sp = s;
    while (isxdigit(str[0])) {
        sp[0] = str[0];
        sp++;
        str++;
    }
    sp[0] = '\0';

    int ret = fppol_set_str(a, s);

    if (!ret) {
        fprintf(stderr, "Error while reading poly at line %" PRIu64 "\n", nl);
        fprintf(stderr, "Offending line:\n%s", curr_line);
        exit(EXIT_FAILURE);
    }

    return str;
}

const char * advance_char(const char *str, char c, uint64_t nl)
{
    if (str[0] != c) {
        fprintf(stderr, "Error. Expect %c at line %" PRIu64 "\n", c, nl);
        fprintf(stderr, "Offending line:\n%s", curr_line);
        exit(EXIT_FAILURE);
    }
    return str+1;
}


void check_line(const char * line, ffspol_t *ffspol, uint64_t nl)
{
    fppol_t a, b, norm, tmp, factor;
    fppol_init(a);
    fppol_init(b);
    fppol_init(norm);
    fppol_init(factor);
    fppol_init(tmp);

    line = read_ffpol_and_advance(a, line, nl);
    line = advance_char(line, ',', nl);
    line = read_ffpol_and_advance(b, line, nl);
    line = advance_char(line, ':', nl);

    for (int side = 0; side < 2; side++) {
        ffspol_norm(norm, ffspol[side], a, b);
        fppol_set_ti(tmp, 0);
        do {
            line = read_ffpol_and_advance(factor, line, nl);
            fppol_mul(tmp, tmp, factor);
            if ((side == 0) && (line[0] == ':')) {
                line++;
                break;
            }
            if ((side == 1) && (line[0] == '\n'))
                break;
            if (line[0] != ',') {
                fprintf(stderr, "Error. Expect , at line %" PRIu64 "\n", nl);
                fprintf(stderr, "Offending line:\n%s", curr_line);
                exit(EXIT_FAILURE);
            }
            line++;
        } while (1);

        fppol_make_monic(norm);
        fppol_make_monic(tmp);
        if (!fppol_eq(norm, tmp)) {
            fprintf(stderr,
                    "Error: Wrong norm factorization on side %d, line %"
                    PRIu64 "\n", side, nl);
            fprintf(stderr, "Offending line:\n%s", curr_line);
            exit(EXIT_FAILURE);
        }
    }

    fppol_clear(a);
    fppol_clear(b);
    fppol_clear(norm);
    fppol_clear(factor);
    fppol_clear(tmp);
}

void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0             function field polynomial on side 0\n");
    fprintf(stderr, "  pol1             function field polynomial on side 1\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");

    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    int gf = 0;
    char *argv0 = argv[0];
    param_list pl;
    param_list_init(pl);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
    }

    ffspol_t ffspol[2];

    // read function field polynomials
    {
        const char * polstr;
        ffspol_init(ffspol[0]);
        polstr = param_list_lookup_string(pl, "pol0");
        if (polstr != NULL) {
            ffspol_set_str(ffspol[0], polstr);
        } else 
            usage(argv0, "pol0");
        ffspol_init(ffspol[1]);
        polstr = param_list_lookup_string(pl, "pol1");
        if (polstr != NULL) {
            ffspol_set_str(ffspol[1], polstr);
        } else
            usage(argv0, "pol1");
    }

    uint64_t nl = 0, nrels = 0;
    char line[MAX_LINE_SIZE];
    curr_line = line;
    do {
        line[0]='\0';
        char *ret = fgets(line, MAX_LINE_SIZE, stdin);
        nl++;
        if (ret == NULL) {
            if (feof(stdin)) {
                if (line[0] != '\0') {
                    fprintf(stderr, "Incomplete last line\n");
                    return EXIT_FAILURE;
                }
            } else {
                fprintf(stderr,
                        "I/O error while reading line %" PRIu64 "\n", nl);
                return EXIT_FAILURE;
            }
            break;
        }
        int nc = strlen(line);
        ASSERT_ALWAYS(nc >= 1);
        if (line[nc-1] != '\n') {
            fprintf(stderr,
                    "Line too long while reading line %" PRIu64 "\n", nl);
            return EXIT_FAILURE;
        }
        if (line[0] == '#')
            continue;

        check_line(line, ffspol, nl);
        nrels++;
        if ((nrels % 1000000) == 0) {
            fprintf(stderr, "%" PRIu64 " M relations checked.\n",
                    nrels / 1000000);
        }
    } while (1);

    fprintf(stderr, "Total number of relations: %" PRIu64 "\n", nrels);
    param_list_clear(pl);
    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);
    return EXIT_SUCCESS;
}
