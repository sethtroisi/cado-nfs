#ifndef PARAMS_H_
#define PARAMS_H_

#include <stdio.h>
#include <gmp.h>

/* This is by increasing order of priority */
enum parameter_origin { PARAMETER_FROM_FILE, PARAMETER_FROM_CMDLINE, };

struct parameter_s {
    char * key;
    char * value;
    enum parameter_origin origin;
    int parsed;
};
typedef struct parameter_s parameter[1];
typedef struct parameter_s * parameter_ptr;
typedef struct parameter_s const * parameter_srcptr;

struct param_list_s {
    unsigned int alloc;
    unsigned int size;
    parameter * p;
    int consolidated;
};

typedef struct param_list_s param_list[1];
typedef struct param_list_s * param_list_ptr;
typedef struct param_list_s const * param_list_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

// in any case, calls to param_list functions overwrite the previously
// set parameters in the parameter list.

extern void param_list_init(param_list pl);
extern void param_list_clear(param_list pl);


// takes a file, in the Cado-NFS params format, and stores the dictionary
// of parameters to pl.
extern int param_list_read_stream(param_list pl, FILE *f);
extern int param_list_read_file(param_list pl, const char * name);

// sees whether the arguments pointed to by argv[0] and (possibly)
// argv[1] correspond to either -<key> <value>, --<key> <value> or <key>=<value>
// Having key==NULL means that anything will be read.
extern int param_list_update_cmdline(param_list pl, const char * key,
        int * p_argc, char *** p_argv);

// This one allows shorthands. Notice that the alias string has to
// contain the exact for of the wanted alias, which may be either "-x",
// "--x", or "x=" (a terminating = tells the program that the option is
// wanted all in one go).
extern int param_list_update_cmdline_alias(param_list pl, const char * key,
        const char * alias, int * p_argc, char *** p_argv);

extern int param_list_parse_int(param_list, const char *, int *);
extern int param_list_parse_long(param_list, const char *, long *);
extern int param_list_parse_uint(param_list, const char *, unsigned int *);
extern int param_list_parse_ulong(param_list, const char *, unsigned long *);
extern int param_list_parse_double(param_list, const char *, double *);
extern int param_list_parse_string(param_list, const char *, char *, size_t);
extern int param_list_parse_mpz(param_list, const char *, mpz_ptr);

// tells whether everything has been consumed. Otherwise, return the key
// of the first unconsumed argument.
extern int param_list_all_consumed(param_list pl, char ** extraneous);

// warns against unused command-line parameters. This normally indicates
// a user error. parameters ignored from config files are considered
// normal (although note that in some cases, it could be bad as well).
extern int param_list_warn_unused(param_list pl);

// this one is the ``joker'' call.
extern void param_list_add_key(param_list pl,
        const char *, const char *, enum parameter_origin);

// for debugging.
extern void param_list_display(param_list pl, FILE *f);
#ifdef __cplusplus
}
#endif

#endif	/* PARAMS_H_ */
