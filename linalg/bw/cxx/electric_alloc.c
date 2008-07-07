
/* undefine to protect underruns instead */
#define PROTECT_OVERRUN

#include <fcntl.h>
#include <sys/mman.h>

void * electric_alloc(size_t s)
{
    /* Use the method of the good old days from electric fence. */
    char *p;
    size_t r = 8192;        /* Any multiple of the page size will do. */
    unsigned int multip = (s+r-1)/r;
    p = (char *)
    mmap(0, (multip + 1) * r, PROT_READ | PROT_WRITE,
                MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#ifdef PROTECT_OVERRUN
    p += (multip + 1) * r;
    mprotect((void*) (p-r), r, PROT_NONE);
    p -= r+s;
#else
    /* protect underrun */
    mprotect(p, r, PROT_NONE);
    p += r;
#endif
    return (void *) p;
}

void electric_free(void * p0, size_t s)
{
    char * p = (char *) p0;
    size_t r = 8192;
    unsigned int multip = (s+r-1)/r;
#ifdef PROTECT_OVERRUN
    p += s;
    mprotect((void*) (p - r), r, PROT_READ | PROT_WRITE);
    p -= (multip + 1) * r;
#else
    p -= r;
    mprotect(p, r, PROT_READ | PROT_WRITE);
#endif
    munmap(p, (multip + 1) * r);
    // actually unmapping crashes, I don't know why.
}

