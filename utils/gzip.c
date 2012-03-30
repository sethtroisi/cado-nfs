#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "macros.h"
#include "gzip.h"
#include "misc.h"

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

#if 0
const char * suffix = NULL;

const char * copy_suffix_noalloc(const char * name)
{
    const char * p = strrchr(name, '.');
    if (p == NULL)
        p = name + strlen(name);
    return strdup(p);
}

const char * copy_suffix_alloc(const char * name)
{
    return strdup(copy_suffix_noalloc(name));
}
const char * path_remove_suffix(char * name)
{
    char * p = strrchr(name, '.');
    if (p) *p=0;
    return name;
}

#endif

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "gzip -dc %s", "gzip -c --fast > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    /* These two have to be present */
    { "", NULL, NULL },
    { NULL, NULL, NULL },
};

const char * path_basename(const char * path)
{
    const char *p = strrchr(path, '/');
    if (p == NULL) {
        p = path;
    } else {
        p = p + 1;
    }
    return p;
}

int is_supported_compression_format(const char * s)
{
    struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        if (strcmp(r->suffix, s) == 0)
            return 1;
    }
    return 0;
}

/* Return a unix commands list. Exemple :
   cat file_relation1
   gzip -dc file_relation2.gz file_relation3.gz
   bzip2 -dc file_relation4.gz file_relation5.gz
   [empty string]
*/
char **
prempt_open_compressed_rs (char **ficname)
{
  const struct suffix_handler *cp_r = NULL, *r = supported_compression_formats;
  char **cmd;
  size_t s_cmds = 2;
  size_t p_cmds = 0;
  int suffix_choice = 0;

  if (!(cmd = calloc (s_cmds, sizeof(unsigned char *)))) {
    fprintf (stderr, "fopen_compressed_rs: calloc erreur : %s\n", strerror(errno));
    exit (1);
  }
  while (*ficname)
    if (!suffix_choice) {
      if (p_cmds + 1 >= s_cmds) {
	if (!(cmd = realloc (cmd, sizeof(unsigned char *) * (s_cmds<<1)))) {
	  fprintf (stderr, "fopen_compressed_rs: realloc erreur : %s\n", strerror(errno));
	  exit (1);
	}
	memset(&cmd[s_cmds], 0, sizeof(unsigned char *) * s_cmds);
	s_cmds <<= 1;
      }
      if (!(cmd[p_cmds] = malloc(PREMPT_S_CMD))) {
	fprintf (stderr, "fopen_compressed_rs: malloc erreur : %s\n", strerror(errno));
	exit (1);
      }
      for (cp_r = r ; cp_r->suffix ; cp_r++)
	if (has_suffix (*ficname, cp_r->suffix)) break;
      strcpy (cmd[p_cmds], cp_r->pfmt_in ? cp_r->pfmt_in : "cat %s");
      cmd[p_cmds][strlen(cmd[p_cmds]) - 2] = 0; /* "%s" suppress */
      suffix_choice = 1;
      if (strlen (*ficname) + strlen (cmd[p_cmds]) >= PREMPT_S_CMD) {
	fprintf(stderr, "fopen_compressed_rs: PREMPT_S_CMD (%d) too small. Please * 2\n", PREMPT_S_CMD);
	exit (1);
      }
      strcat (cmd[p_cmds], *ficname++);
    }
    else {
      if (has_suffix (*ficname, cp_r->suffix) &&
	  (strlen (*ficname) + strlen (cmd[p_cmds]) + 1 < PREMPT_S_CMD)) {
	strcat (cmd[p_cmds], " ");
	strcat (cmd[p_cmds], *ficname++);
      }
      else {
	suffix_choice = 0;
	p_cmds++;
      }
    }
  return cmd;
}

/* The pipe capacity is 2^16 by default, we can increase it, but it does not
   seem to make a difference, thus we don't change it by default (patch from
   Alain Filbois). */
#define PIPE_CAPACITY 1UL << 20

FILE*
fopen_compressed_r (const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_in) { /* this suffix has an associated deflate command */
            char * command;
            int ret = asprintf(&command, r->pfmt_in, name);
            ASSERT_ALWAYS(ret >= 0);
            f = popen(command, "r");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            /* remove the 'xxx' to change the pipe capacity */
#ifdef F_SETPIPE_SZxxx
            fcntl (fileno (f), F_SETPIPE_SZ, PIPE_CAPACITY);
#endif
            return f;
        } else {
            f = fopen(name, "r");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE*
fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_out) { /* suffix has an associated compression command */
            char * command;
            int ret = asprintf(&command, r->pfmt_out, name);
            ASSERT_ALWAYS(ret >= 0);
            f = popen(command, "w");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            /* remove the 'xxx' to change the pipe capacity */
#ifdef F_SETPIPE_SZxxx
            /* increase the pipe capacity (2^16 by default), thanks to
               Alain Filbois */
            fcntl (fileno (f), F_SETPIPE_SZ, PIPE_CAPACITY);
#endif
            return f;
        } else {
            f = fopen(name, "w");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE * gzip_open(const char * name, const char * mode)
{
    if (strcmp(mode, "r") == 0) {
        return fopen_compressed_r(name, NULL, NULL);
    } else if (strcmp(mode, "w") == 0) {
        return fopen_compressed_w(name, NULL, NULL);
    } else {
        abort();
    }
}

void gzip_close(FILE * f, const char * name)
{
    const char * e = name + strlen(name);
    const struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        const char * ne = e - MIN(strlen(name), strlen(r->suffix));
        if (strcmp(r->suffix, ne) != 0)
            continue;
        if (r->pfmt_out) {
            pclose(f);
            return;
        } else {
            fclose(f);
            return;
        }
    }
}

