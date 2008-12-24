/* popen() and pclose() are POSIX, not C99. So we have to rely on the
 * libc's specifics to enable this function. For glibc, it's the
 * following macro. For other support libraries, it might be something
 * else. See also popen(3) and feature_test_macros(7)
 */
#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gzip.h"

int is_gzip(const char * s)
{
    unsigned int l = strlen(s);
    return l >= 3 && strcmp(s + l - 3, ".gz") == 0;
}

FILE *gzip_open(const char * s, const char * t)
{
    FILE *ofile = NULL;

    if(is_gzip(s)){
	char command[1024];
	if(!strcmp(t, "w")){
            /* first check the file exists and we can write it */
            ofile = fopen (s, "w");
            if (ofile == NULL)
              {
                fprintf (stderr, "Error, cannot write file %s\n", s);
                exit (1);
              }
            fclose (ofile);
	    snprintf(command, sizeof(command), "gzip -c > %s", s);
	    ofile = popen(command, "w");
	}
	else if(!strcmp(t, "r")){
            /* first check the file exists and we can read it */
            ofile = fopen (s, "r");
            if (ofile == NULL)
              {
                fprintf (stderr, "Error, cannot read file %s\n", s);
                exit (1);
              }
            fclose (ofile);
	    snprintf(command, sizeof(command), "gzip -dc %s", s);
	    ofile = popen(command, "r");
	}
	else{
	    fprintf(stderr, "gzip_open: not implemented %s\n", t);
	    return NULL;
	}
    }
    else
	ofile = fopen(s, t);
    return ofile;
}

void gzip_close(FILE *ofile, const char * s)
{
    if(is_gzip(s))
	pclose(ofile);
    else
	fclose(ofile);
}
