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
	    snprintf(command, sizeof(command), "gzip -c > %s", s);
	    ofile = popen(command, "w");
	}
	else if(!strcmp(t, "r")){
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
