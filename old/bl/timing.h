

#ifndef TIMING_H_
#define TIMING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#include <sys/time.h>
#include <sys/resource.h>
#define HAS_microseconds
    extern uint64_t microseconds();
    extern double seconds();
    extern uint64_t wallclock_microseconds();
    extern double wallclock_seconds();

#if defined (__i386__) || defined(__x86_64__)
#define	HAS_cputicks
    static inline uint64_t cputicks() {
	uint64_t r;
	__asm__ __volatile__(".align 64\n\t" "xorl %%eax,%%eax\n\t" "cpuid\n\t" "rdtsc\n\t" "movl %%eax,(%0)\n\t" "movl %%edx,4(%0)\n\t" "xorl %%eax,%%eax\n\t" "cpuid\n\t":	/* no output */
			     :"S"(&r)
			     :"eax", "ebx", "ecx", "edx", "memory");
	 return r;
    }
#endif
#ifdef __cplusplus
}
#endif
#endif				/* TIMING_H_ */
