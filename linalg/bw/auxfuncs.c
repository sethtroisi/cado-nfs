#include <sys/time.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "manu.h"
#ifdef NEED_MALLOC_H_
#include <malloc.h>
#endif
#include "auxfuncs.h"

#define FILENAME_LENGTH		80
#define MALLOC_LOG_LIMIT	16384

unsigned int ceil_log2(unsigned int d)
{
	unsigned int k;
	for(k=0;d;k++) d>>=1;
	return k;
}

void die(const char * fmt, int signal, ... )
{
	va_list ap;

#ifdef __SUNPRO_C
	va_start(ap,(void)signal);
#else
	va_start(ap,signal);
#endif
	vfprintf(stderr,fmt,ap);
	va_end(ap);
	if (signal>=0)  {
		eternal_sleep();
	} else {
		exit(-signal);
	}
	exit(0);
}

void eternal_sleep(void)
{
	volatile int foobar=0;
	fprintf(stderr,"\nAttach gdb if you want to...\n");
	fflush(stderr);
	for(;foobar==0;sleep(10));
	if (foobar==0) exit(255);
}

void * _my_malloc(size_t s,const char * fn,int lin)
{
#ifdef TRACE_MALLOC
	static int i=0;
	i++;
#endif
#ifdef MALLOC_LOG_LIMIT
	if (s>MALLOC_LOG_LIMIT) {
		int k;
		size_t t=s;
		const char * units[4]={"B", "KB", "MB", "GB"};
		for(k=0;k<3 && ((t&((1L<<10)-1))==0 || (k<2 && t>=(1L<<10)));t>>=10,k++);
		fprintf(stdout,"// Allocating %d%s, file %s, line %d\n",
			(int)t,units[k],fn,lin);
	}
#endif
	if (s!=0)
	{
		void * p = malloc(s);
		if (p!=NULL) return p;
		else
			die("Allocation of %d bytes failed,"
			" not enough memory\n"
			" file %s, line %d\n", 63, s, fn, lin);
        }
	return NULL;
}

int exist(const char * s)
{
	return (access(s,R_OK)==0);
}

static void enforce_limit(int l, rlim_t min, const char * name)
{
	struct rlimit lim;

	if (getrlimit(l,&lim)<0)
		die("getrlimit() failed!\n",31);

	if (lim.rlim_cur == RLIM_INFINITY)
		return;

	if (min==RLIM_INFINITY || lim.rlim_cur < min) {
		die("FATAL: Limit %s is %d while we need at least %d\n",-31,
				name,lim.rlim_cur,(unsigned long) min);
	}
}

void coredump_limit(int allow)
{
	struct rlimit lim;

	int change_limits[]={	RLIMIT_CPU,	/* CPU time in seconds */
				RLIMIT_FSIZE,	/* filesize */
				RLIMIT_DATA,	/* data size */
				RLIMIT_STACK,	/* stack size */
				RLIMIT_CORE,	/* core file size */
#ifdef __SVR4
				RLIMIT_VMEM,
#else
				RLIMIT_RSS,	/* resident set size */
				RLIMIT_AS,	/* address space (vm) limit */
#endif
#ifdef linux
				RLIMIT_NPROC,	/* number of processes */
				RLIMIT_MEMLOCK,	/* locked-in-memory addr space*/
#endif
				RLIMIT_NOFILE	/* number of open files */
	};
	unsigned int i;

	for(i=0;i<sizeof(change_limits)/sizeof(change_limits[0]);i++) {
		if (getrlimit(change_limits[i],&lim)<0)
			die("getrlimit() failed!\n",31);
		lim.rlim_cur=lim.rlim_max;
		if (setrlimit(change_limits[i],&lim) < 0)
			die("setrlimit() failed!\n",31);
	}

	fprintf(stdout,"// All resources unlimited\n");

	if (allow) {
		fprintf(stdout,"// core dump size unlimited\n");
		if (getrlimit(RLIMIT_CORE,&lim) < 0)
			die("getrlimit() failed!\n",31);
		if (lim.rlim_max != RLIM_INFINITY)
			fprintf(stdout,"// (actual limit is %lu bytes)\n",
					(unsigned long) lim.rlim_max);
		lim.rlim_cur=lim.rlim_max;
		if (setrlimit(RLIMIT_CORE,&lim) < 0)
			die("setrlimit() failed!\n",31);
	} else {
		fprintf(stdout,"// core dump size limited to 0kb\n");
		lim.rlim_cur=lim.rlim_max=0;
		if (setrlimit(RLIMIT_CORE,&lim) < 0)
			die("setrlimit() failed!\n",31);
	}
	fflush(stdout);
	fflush(stderr);
}

void check_limit_requirements(void)
{
	enforce_limit(RLIMIT_CPU,3600*10,"cpu time");
	enforce_limit(RLIMIT_FSIZE,1<<30,"file size");
	enforce_limit(RLIMIT_DATA,1<<30,"data size");
	enforce_limit(RLIMIT_STACK,1<<20,"stack size");
#ifdef __SVR4
	enforce_limit(RLIMIT_VMEM,1<<30,"vmem size");
#else
	enforce_limit(RLIMIT_AS,1<<30,"vmem size");
#endif
	enforce_limit(RLIMIT_NOFILE,8,"open files");
}

int mkfname(char *s, const char *fmt, ...)
{
	va_list ap;
#ifndef	HAS_NOT_SNPRINTF
	int res;
#endif

	va_start(ap,fmt);

#ifndef	HAS_NOT_SNPRINTF
	res=vsnprintf(s,FILENAME_LENGTH,fmt,ap);
	if (res>=FILENAME_LENGTH || res==0) {
		va_end(ap);
		return -1;
	}
#else
	if (strlen(fmt)>=FILENAME_LENGTH) {
		strncpy(s,fmt,FILENAME_LENGTH);
		va_end(ap);
		return -1;
	}
	vsprintf(s,fmt,ap);
#endif
	va_end(ap);
	return 0;
}

long fcopy(FILE *dst, FILE *src, size_t kamount)
{
	char buf[BUFSIZ];
	long nr,nw,ncp;
	ssize_t amount=kamount;

	ncp=0;

	for(;!feof(src) && amount>0;) {
		long maxcp=MIN(BUFSIZ,amount);
		for(nr=0;nr < maxcp && !feof(src) && !ferror(src);) {
			nr+=fread(buf+nr,1,maxcp-nr,src);
		}
		if (ferror(src))
			break;
		for(nw=0;nw<nr && !ferror(dst);) {
			nw+=fwrite(buf+nw,1,nr-nw,dst);
		}
		amount-=nw;
		ncp+=nw;
		if (ferror(dst))
			break;
	}
	return ncp;
}

/*
void _BUG(const char * fn,const int lin)
{
	fprintf(stderr,"Bug in %s, line %d ; last errno = %s\n",
			fn,lin,strerror(errno));
	eternal_sleep();
}
*/
	
#ifdef HAS_NOT_SNPRINTF
int vsnprintf(char * s, size_t n, const char * fmt, va_list ap)
{
	if (strlen(fmt)>n)	/* at least handle this case... */
		BUG();
	
	return vsprintf(s,fmt,ap);
}

int snprintf(char * s, size_t n, const char * fmt, ...)
{
	va_list ap;
	int res;
#ifdef __SUNPRO_C
	va_start(ap,(void)fmt);
#else
	va_start(ap,fmt);
#endif
	res=vsnprintf(s,n,fmt,ap);
	va_end(ap);
	return res;
}
#endif	/* HAS_NOT_SNPRINTF */

/* vim:set sw=8: */
