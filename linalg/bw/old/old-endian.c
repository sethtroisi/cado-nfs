#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "old-endian.h"

type32 builtin_endianness;

void reverse(char *a, size_t width, size_t n)
{
    int i;
    char *buf;

    buf = malloc(width);

    for(i = 0 ; i < ((int) n) - i - 1 ; i++)
    {
	memcpy(buf,a+i*width,width);
	memcpy(a+i*width,a+(n-i-1)*width,width);
	memcpy(a+(n-i-1)*width,buf,width);
    }

    free(buf);
}

void check_endianness(FILE * f)
{
    int i;
    type16 test16;

    for(i=0;i<2;i++)
	((char*)&test16)[i]=(char)i;

    for(i=0;i<4;i++)
	((char*)&builtin_endianness)[i]=(char)i;

    if (M_LITTLE_ENDIAN)
    {
	if (test16 == UC16(0x0100)) {
	    if (f) {
		    fprintf(f,"CPU is little endian\n");
	    }
	} else {
	    die("CPU looks like little endian, but 16bits look weird\n",31);
	}
    }

#ifndef KNOWN_LITTLE_ENDIAN
    else if (M_BIG_ENDIAN)
    {
	if (test16 == UC16(0x0001)) {
	    if (f) {
		    fprintf(f,"CPU is big endian\n");
	    }
	} else {
	    die("CPU looks like big endian, but 16bits look weird\n",31);
	}

	mswap32(builtin_endianness);
	if (!M_LITTLE_ENDIAN)
	    die("The mswap32 function does not work correctly\n",31);
	mswap32(builtin_endianness);

	mswap16(test16);
	if (test16 != UC16(0x0100))
	    die("The mswap16 function does not work correctly\n",31);
	mswap16(test16);
    }
    else
	die("Endianness unknown (%08lX). Patch it yourself\n",
		31,builtin_endianness);
#endif
}

#ifndef KNOWN_LITTLE_ENDIAN
void f_WRITE32(type32*d,type32 s)
{
    *d=s;
    DO_BIG_ENDIAN(mswap32(*d));
}

type32 f_READ32(type32* s)
{
    type32 buf;
    buf=*s;
    DO_BIG_ENDIAN(mswap32(*s));
    return buf;
}
#endif
