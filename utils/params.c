#define _BSD_SOURCE     /* strdup */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>              /* isdigit isspace */
#include <limits.h>             /* INT_MIN INT_MAX */
#include <errno.h>

#include "params.h"
#include "macros.h"

void param_list_init(param_list pl)
{
    pl->size = 0;
    pl->alloc = 16;
    pl->p = (parameter *) malloc(pl->alloc * sizeof(parameter));
    pl->consolidated = 1;
}

void param_list_clear(param_list pl)
{
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        free(pl->p[i]->key);
        free(pl->p[i]->value);
    }
    free(pl->p);
    memset(pl, 0, sizeof(pl));
}

static void make_room(param_list pl, unsigned int more)
{
    if (pl->size + more <= pl->alloc) {
        return;
    }

    for( ; pl->size + more > pl->alloc ; ) {
        pl->alloc += pl->alloc >> 1;
    }

    pl->p = (parameter *) realloc(pl->p, pl->alloc * sizeof(parameter));
}

static void param_list_add_key_nodup(param_list pl,
        char * key, char * value, enum parameter_origin o)
{
    make_room(pl, 1);
    pl->p[pl->size]->key = key;
    pl->p[pl->size]->value = value;
    pl->p[pl->size]->origin = o;
    pl->p[pl->size]->parsed = 0;
    pl->size++;
    pl->consolidated = 0;
}

void param_list_add_key(param_list pl,
        const char * key, const char * value, enum parameter_origin o)
{
    param_list_add_key_nodup(pl, strdup(key), strdup(value), o);
}

struct sorting_data {
    parameter_srcptr s;
    int p;
};

typedef int (*sortfunc_t) (const void *, const void *);

int paramcmp(const struct sorting_data * a, const struct sorting_data * b)
{
    int r;
    r = strcmp(a->s->key, b->s->key);
    if (r) return r;
    r = a->s->origin - b->s->origin;
    if (r) return r;
    /* Compare by pointer position so that the sort is stable */
    return a->p - b->p;
}

/* Sort (for searching), and remove duplicates. */
void param_list_consolidate(param_list pl)
{
    if (pl->consolidated) {
        return;
    }
    struct sorting_data * intermediate;
    intermediate = malloc(pl->size * sizeof(struct sorting_data));
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        intermediate[i].p = i;
        intermediate[i].s = pl->p[i];
    }
    qsort(intermediate, pl->size, sizeof(struct sorting_data),
            (sortfunc_t) &paramcmp);

    parameter * np = (parameter *) malloc(pl->alloc * sizeof(parameter));
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        memcpy(np[i], pl->p[intermediate[i].p], sizeof(parameter));
    }
    free(pl->p);
    pl->p = np;
    free(intermediate);

    // now remove duplicates. The sorting has priorities right.
    unsigned int j = 0;
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (i + 1 < pl->size && strcmp(pl->p[i]->key, pl->p[i+1]->key) == 0) {
            free(pl->p[i]->key);
            free(pl->p[i]->value);
        } else {
            // I can't see why there could conceivably be a problem if i == j,
            // but valgrind complains...
            if (i != j) {
                memcpy(pl->p[j], pl->p[i], sizeof(parameter));
            }
            j++;
        }
    }
    pl->size = j;

    pl->consolidated = 1;
}

int param_list_read_stream(param_list pl, FILE *f)
{
    int all_ok=1;
    const int linelen = 512;
    char line[linelen];
    char * newkey;
    char * newvalue;
    while (!feof(f)) {
        if (fgets(line, linelen, f) == NULL)
            break;
        if (line[0] == '#')
            continue;
        // remove possible comment at end of line.
        char * hash;
        if ((hash = strchr(line, '#')) != NULL) {
            *hash = '\0';
        }

        char * p = line;

        // trailing space
        int l = strlen(p);
        for( ; l && isspace(p[l-1]) ; l--);
        p[l] = '\0';

        // leading space.
        for( ; *p && isspace(*p) ; p++, l--);

        // empty ps are ignored.
        if (l == 0)
            continue;

        // look for a left-hand-side.
        l = 0;
        for( ; p[l] && (isalnum(p[l]) || p[l] == '_') ; l++);

        int lhs_length = l;

        if (lhs_length == 0) {
            fprintf(stderr, "Parse error, no usable key for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }

        /* Now we can match (whitespace*)(separator)(whitespace*)(data)
         */
        char * q = p + lhs_length;
        for( ; *q && isspace(*q) ; q++);

        /* match separator, which is one of : = := */
        if (*q == '=') {
            q++;
        } else if (*q == ':') {
            q++;
            if (*q == '=')
                q++;
        } else {
            fprintf(stderr, "Parse error, no separator for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }
        for( ; *q && isspace(*q) ; q++);

        newkey = malloc(lhs_length + 1);
        memcpy(newkey, p, lhs_length);
        newkey[lhs_length]='\0';

        newvalue = strdup(q);

        param_list_add_key_nodup(pl, newkey, newvalue, PARAMETER_FROM_FILE);
    }
    param_list_consolidate(pl);

    return all_ok;
}

int param_list_read_file(param_list pl, const char * name)
{
    FILE * f;
    f = fopen(name, "r");
    if (f == NULL) {
        fprintf(stderr, "Cannot read %s\n", name);
        exit(1);
    }
    int r = param_list_read_stream(pl, f);
    fclose(f);
    return r;
}

int param_list_update_cmdline(param_list pl, const char * key,
        int * p_argc, char *** p_argv)
{
    if (*p_argc == 0)
        return 0;
    const char * a = (*p_argv[0]);
    if (key == NULL) {
        /* A NULL key is a wildcard */
        if (*p_argc >= 2 && a[0] == '-') {
            a++;
            a+= *a == '-';
            int x=0;
            /* parameters _must_ begin by alphabetic characters, or
             * otherwise we suffer to distinguish immediate numerical
             * entries.
             */
            if (!isalpha(a[x]))
                return 0;
            for( ; a[x] && (isalnum(a[x]) || a[x] == '_') ; x++);
            if (a[x] == '\0') {
                param_list_add_key(pl, a, (*p_argv)[1], PARAMETER_FROM_CMDLINE);
                (*p_argv)+=2;
                (*p_argc)-=2;
                return 1;
            }
        } else {
            /* Check for <key>=<value> syntax */
            int x=0;
            for( ; a[x] && (isalnum(a[x]) || a[x] == '_') ; x++);
            if (a[x] == '=' && a[x+1]) {
                char * newkey = malloc(x+1);
                memcpy(newkey, a, x);
                newkey[x]='\0';
                char * newvalue = strdup(a+x+1);
                param_list_add_key_nodup(pl,
                        newkey, newvalue, PARAMETER_FROM_CMDLINE);
                (*p_argv)+=1;
                (*p_argc)-=1;
                return 1;
            }
        }
        return 0;
    }

    if (*p_argc >= 2 && a[0] == '-') {
        a++;
        a+= *a == '-';
        if (strcmp(a, key) == 0) {
            param_list_add_key(pl, key, (*p_argv)[1], PARAMETER_FROM_CMDLINE);
            (*p_argv)+=2;
            (*p_argc)-=2;
            return 1;
        }
    }
    
    if (strncmp(a, key, strlen(key)) == 0) {
        /* Check for <key>=<value> syntax */
        a += strlen(key);
        if (*a == '=' && a[1]) {
            param_list_add_key(pl, key, a+1, PARAMETER_FROM_CMDLINE);
            (*p_argv)+=1;
            (*p_argc)-=1;
            return 1;
        }
    }
    return 0;
}

int param_list_update_cmdline_alias(param_list pl, const char * key,
        const char * alias, int * p_argc, char *** p_argv)
{
    ASSERT_ALWAYS(alias != NULL);
    const char * a = (*p_argv[0]);
    if (alias[strlen(alias)-1] == '=') {
        if (strncmp(a, alias, strlen(alias)) != 0)
            return 0;
        a += strlen(alias);
        if (a[1] == '\0')
            return 0;
        param_list_add_key(pl, key, a+1, PARAMETER_FROM_CMDLINE);
        (*p_argv)+=1;
        (*p_argc)-=1;
        return 1;
    }
    if (strcmp(a, alias) == 0) {
        param_list_add_key(pl, key, (*p_argv)[1], PARAMETER_FROM_CMDLINE);
        (*p_argv)+=2;
        (*p_argc)-=2;
        return 1;
    }
    return 0;
}

int param_strcmp(const char * a, parameter_srcptr b)
{
    return strcmp(a, b->key);
}

static int assoc(param_list pl, const char * key)
{
    void * found;

    param_list_consolidate(pl);
    found = bsearch(key, pl->p, pl->size,
            sizeof(parameter), (sortfunc_t) param_strcmp);
    if (found == NULL)
        return -1;
    parameter * c = (parameter *) found;
    return c-pl->p;
}

int param_list_parse_long(param_list pl, const char * key, long * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    long res;
    res = strtol(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not a long: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return 1;
}

int param_list_parse_int(param_list pl, const char * key, int * r)
{
    long res;
    if (param_list_parse_long(pl, key, &res) == 0)
        return 0;
    if (res > INT_MAX || res < INT_MIN) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within an int: %ld\n",
                key, res);
        exit(1);
    }
    if (r)
        *r  = res;
    return 1;
}

int param_list_parse_intxint(param_list pl, const char * key, int * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    long res[2];
    res[0] = strtol(value, &end, 0);
    if (*end != 'x') {
        fprintf(stderr, "Parse error: parameter for key %s"
                " must match %%dx%%d; got %s\n",
                key, pl->p[v]->value);
        exit(1);
    }
    value = end + 1;
    res[1] = strtol(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s"
                " must match %%dx%%d; got %s\n",
                key, pl->p[v]->value);
        exit(1);
    }
    if (r) {
        r[0] = res[0];
        r[1] = res[1];
    }
    return 1;
}


int param_list_parse_ulong(param_list pl, const char * key, unsigned long * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    unsigned long res;
    res = strtoul(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error:"
                " parameter for key %s is not an ulong: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return 1;
}

int param_list_parse_uint(param_list pl, const char * key, unsigned int * r)
{
    unsigned long res;
    if (param_list_parse_ulong(pl, key, &res) == 0)
        return 0;
    if (res > UINT_MAX) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within an uinsigned int: %ld\n",
                key, res);
        exit(1);
    }
    if (r)
        *r  = res;
    return 1;
}


int param_list_parse_double(param_list pl, const char * key, double * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    double res;
    res = strtod(value, &end);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not an int: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return 1;
}

int param_list_parse_string(param_list pl, const char * key, char * r, size_t n)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    pl->p[v]->parsed=1;
    char * value = pl->p[v]->value;
    if (r && strlen(value) > n-1) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within string buffer"
                " of length %lu\n", key, (unsigned long) n);
        exit(1);
    }
    if (r)
        strncpy(r, value, n);
    return 1;
}

int param_list_parse_mpz(param_list pl, const char * key, mpz_ptr r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    pl->p[v]->parsed=1;
    unsigned int nread;
    int rc;
    char * value = pl->p[v]->value;
    if (r) {
        rc = gmp_sscanf(value, "%Zd%n", r, &nread);
    } else {
        /* scan even when the result is not wanted */
        rc = gmp_sscanf(value, "%*Zd%n", &nread);
    }
    if (rc != 1 || value[nread] != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not an mpz: %s\n",
                key, value);
        exit(1);
    }
    return 1;
}

int param_list_all_consumed(param_list pl, char ** extraneous)
{
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (!pl->p[i]->parsed) {
            pl->p[i]->parsed = 1;
            if (extraneous) {
                *extraneous = pl->p[i]->key;
            }
            return 0;
        }
    }
    return 1;
}

int param_list_warn_unused(param_list pl)
{
    int u = 0;
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (pl->p[i]->origin != PARAMETER_FROM_FILE && !pl->p[i]->parsed) {
            fprintf(stderr, "Warning: unused command-line parameter %s\n",
                    pl->p[i]->key);
            u++;
        }
    }
    return u;
}

void param_list_display(param_list pl, FILE *f)
{
    param_list_consolidate(pl);
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        fprintf(f,"%s=%s\n", pl->p[i]->key, pl->p[i]->value);
    }
}

void param_list_save(param_list pl, const char * filename)
{
    FILE * f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }

    param_list_display(pl, f);
    fclose(f);
}
