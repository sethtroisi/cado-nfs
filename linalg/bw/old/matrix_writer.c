#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "params.h"
#include "types.h"
#include "macros.h"
#include "old-endian.h"
#include "auxfuncs.h"

/* There must be plenty of space so that even a very heavy line fits...
 */
#define	MAX_BUF_SIZE	65536

struct parser {
	char * buf;
	int pos;
	int size;
};

void skipws(struct parser * x) {
	for(;x->pos<x->size && isspace(x->buf[x->pos]);x->pos++);
}
unsigned int numlen(struct parser * x) {
	unsigned int w;
	skipws(x);
	for(w=0;(x->pos+w)<x->size && isdigit(x->buf[x->pos+w]);w++);
	return w;
}
void error(struct parser * x, const char * s) {
	fprintf(stderr,"Parse error: %s at %s\n", s,x->buf+x->pos);
	exit(1);
}
int try_match(struct parser * x, const char * s) {
	skipws(x);
	if (strlen(s)>(x->size-x->pos))
		return 0;
	if (strncmp(x->buf+x->pos,s,strlen(s))==0) {
		x->pos+=strlen(s);
		return 1;
	}
	return 0;
}
/*
int try_match(struct parser * x, char c) {
	skipws(x);
	if (x->pos < x->size && x->buf[x->pos]==c) { x->pos++; return 1; }
	return 0;
}
*/
void want(struct parser * x, const char * s) {
	if (!try_match(x, s)) {
		error(x,"unexpected literal");
	}
}
/*
void want(struct parser * x, char c) {
	if (!try_match(c)) {
		error(x,"unexpected literal");
	}
}
*/
unsigned long want_ui(struct parser * x) {
	unsigned long res;
	char * endptr;
	int w=numlen(x);
	if (w==0) error(x,"want a number");
	res=strtoul(x->buf+x->pos,&endptr,0);
	x->pos=(endptr-x->buf);
	return res;
}
long want_si(struct parser * x) {
	int s=1;
	if (try_match(x, "-")) s=-1;
	return s*want_ui(x);
}
int is_comment(struct parser * x) {
	return try_match(x, "// ");
}

int main(int argc, char * argv[])
{
	char buffer[MAX_BUF_SIZE];
	size_t nr,nc;
	char modulus[200];
	int nscan;
	int i;
	FILE *m, *v;
	struct parser ps;
	char mn[]="matrix.0";
	char vn[]="vector.0";

	if (argc>1 && isdigit(argv[1][0])) {
		mn[strlen(mn)-1]=argv[1][0];
		vn[strlen(mn)-1]=argv[1][0];
		argc--;
		argv++;
	}

	do {
		fgets(buffer,MAX_BUF_SIZE,stdin);
		if (strncmp(buffer,"//",2)!=0) {
			fprintf(stderr,"Bad file format ; wrong header\n");
			exit(1);
		}
#ifdef	__x86_64__
		nscan=sscanf(buffer,"// %ld %*s %ld %*s%*[ ;]%*s %200s",
				&nr,&nc,modulus);
#else
		nscan=sscanf(buffer,"// %d %*s %d %*s%*[ ;]%*s %200s",
				&nr,&nc,modulus);
#endif
	} while (nscan!=3);

	/*
	fprintf(stderr,"%d ROWS\n",nr);
	fprintf(stderr,"%d COLUMNS\n",nc);
	fprintf(stderr,"MODULUS %s\n",modulus);
	*/

	v=fopen(vn,"w");
#ifdef  __x86_64__
	fprintf(v,"WIEDEMANN INSTANCE %ld x %ld OVER Z/%sZ\n",nr,nc,modulus);
#else
	fprintf(v,"WIEDEMANN INSTANCE %d x %d OVER Z/%sZ\n",nr,nc,modulus);
#endif
	for(i=0;i<nr;i++) {
		char null[]={0,0,0,0,0,0,0,0};
		fwrite(null,8,1,v);
	}
	fclose(v);
	m=fopen(mn,"w");
#ifdef  __x86_64__
	fprintf(m,"WIEDEMANN INSTANCE %ld x %ld OVER Z/%sZ\n",nr,nc,modulus);
#else
	fprintf(m,"WIEDEMANN INSTANCE %d x %d OVER Z/%sZ\n",nr,nc,modulus);
#endif
	for(i=0;i<nr;i++) {
		int j, w;
		type32 xx;
		stype32 yy;
		
		fgets(buffer,MAX_BUF_SIZE,stdin);
		ps.buf=buffer;
		ps.pos=0;
		ps.size=strlen(buffer);
		if (is_comment(&ps))
			continue;
		w=want_ui(&ps);
		xx=w;
		DO_BIG_ENDIAN(mswap32(xx));
		fwrite(&xx, sizeof(type32), 1,  m);

		for(j=0;j<w;j++) {
			xx=want_ui(&ps);
			DO_BIG_ENDIAN(mswap32(xx));
			fwrite(&xx, sizeof(type32), 1,  m);
			want(&ps,":");
			yy=want_si(&ps);
			DO_BIG_ENDIAN(mswap32(xx));
			fwrite(&yy, sizeof(type32), 1,  m);
		}
		skipws(&ps);
		if (ps.buf[ps.pos]!=0) {
			fprintf(stderr,
				"Parse error at line: <%s> ; "
				"remaining part <%s>\n",
				ps.buf,ps.buf+ps.pos);
			exit(1);
		}
	}
	fclose(m);
	
	return 0;
}
