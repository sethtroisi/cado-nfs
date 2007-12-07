/* Helper program for rewriting headers after tuning.
   
  Copyright 2007 Emmanuel Thome'.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "replace.h"

void set_hash_define_int(struct hash_define * ptr, const char * name, int x)
{
	ptr->identifier = strdup(name);
	ptr->string = malloc(32);
	snprintf(ptr->string, 32, "%d", x);
}

void set_hash_define(struct hash_define * ptr, const char * name, const char *v)
{
	ptr->identifier = strdup(name);
	ptr->string = strdup(v);
}


typedef int (*sortfunc_t) (const void *, const void *);
int cmp_hash_define(struct hash_define * a, struct hash_define * b)
{
	return strcmp(a->identifier, b->identifier);
}
int search_hash_define(const char * a, struct hash_define * b)
{
	return strcmp(a, b->identifier);
}


void replace(struct hash_define * ptr, unsigned int n, const char * fname)
{

	qsort(ptr, n, sizeof(struct hash_define), (sortfunc_t) & cmp_hash_define);

	FILE * fi, * fo;
	char * fname_out;
	size_t sl = strlen(fname)+10;
	char buffer[1024];
	char id[80];

	fname_out = malloc(sl);
	snprintf(fname_out,sl,"%s.new",fname);

	fi = fopen(fname, "r");
	fo = fopen(fname_out, "w");

	for (;!feof(fi);) {
		const char * x, *y;
		if (fgets(buffer, sizeof(buffer), fi) == NULL) {
			/* EOF, most certainly */
			break;
		}
#define	IS_HDEFINE(X)	(strncmp(X,"#define",7) == 0)
#define	IS_HUNDEF(X)	(strncmp(X,"#undef",6) == 0)
		if (!(IS_HDEFINE(buffer) || IS_HUNDEF(buffer))) {
			fputs(buffer, fo);
			continue;
		}
		x = buffer;
		for(;*x && !isspace(*x);x++);
		for(;*x && isspace(*x);x++);
		for(y=x;*y&&!isspace(*y);y++);
		strncpy(id,x,y-x);
		id[y-x]='\0';
		/* identifier is id */
		struct hash_define * loc;
		loc = bsearch(id, ptr, n, sizeof(struct hash_define), (sortfunc_t) & search_hash_define);
		if (loc == NULL) {
			fputs(buffer, fo);
			continue;
		}
		if (loc->string[0] == '\0') {
			fprintf(stderr,
				"Warning: replacement string %s used twice!\n",
				loc->identifier);
		}
		fprintf(fo,"#define %s\t\t%s\n", id, loc->string);
		loc->string[0]='\0';
	}
	for(unsigned int i = 0 ; i < n ; i++) {
		struct hash_define * loc = ptr + i;
		if (loc->string[0] != '\0') {
			fprintf(stderr, "Warning: replacement string %s not used ; value:\n%s\n", loc->identifier, loc->string);
		}
	}

	fclose(fi);
	fclose(fo);
	rename(fname_out, fname);
}
