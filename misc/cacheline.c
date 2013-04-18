#ifndef NBUCKETK
#define NBUCKETK 7
#endif
#define NBUCKET (1<<NBUCKETK)
#define CLSIZE 128
#define ENTRIES_PER_LINE (CLSIZE / sizeof(unsigned int))
#define BUCKET_SIZE ((size_t)1<<(27-NBUCKETK))
/* Pseudo-random sequence */
#define RAND(x) (5*(x)+1)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <emmintrin.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>

/* Since BUCKET_SIZE is a large power of 2, the start of the buckets
   will be a power-of-two size apart in memory, possibly causing a lot
   of cache line collisions. Thus we leave gaps of size bucket_skip
   between the end of a bucket and the start of the next, to make the 
   start address of the buckets not all map to the same few cache lines.
   This did not seem to make much difference in speed.
   Furthermore, the TLB are also set associative, and a large power-of-2
   distance between write pointers causes very frequent mapping of the 
   pages to which the bucket write pointers point to the same TLB set.
   Leaving a page-sized gap fixes this, and makes a noticable difference 
   in speed. */
static const size_t bucket_skip = 4096;
static const size_t fill_iter = BUCKET_SIZE * NBUCKET * 7 / 8;
typedef void (*fillfunc_t)(unsigned int **, unsigned int);


/* Fill buckets, without cacheline buffer. Simply one write to the bucket 
   memory area per update */
void fill(unsigned int **buckets, unsigned int start) {
  unsigned int loc = start;
  size_t i;
  
  for (i = 0; i < fill_iter; i++) {
    const unsigned int bucket = loc >> (32 - NBUCKETK);
    *buckets[bucket]++ = (unsigned int) loc;
    loc = RAND(loc);
  }
}


void fill_cacheline(unsigned int **buckets, unsigned int start) {
  unsigned int loc = start;
  size_t i;
  
  void *cacheline_alloc;
  unsigned int *cacheline;
  // unsigned int cacheline[NBUCKET][ENTRIES_PER_LINE] __attribute__((__aligned__(64)));
  unsigned char nr[NBUCKET];
  
  /* Putting % 128 here makes the code fail half of the time, which is expected,
     but in the cases where it works, it is about 5% faster, which is not expected.
     Oddly, putting %128 here and __aligned__(128) above does NOT produce the speedup */
  // assert ((uintptr_t)&cacheline % 64 == 0);

  i = posix_memalign (&cacheline_alloc, 4096, NBUCKET * CLSIZE + 4096);
  if (i != 0)
    exit(EXIT_FAILURE);
  cacheline = (unsigned int *)cacheline_alloc;
  
  // printf ("cacheline = %p\n", cacheline);
  
  for (i = 0; i < NBUCKET; i++)
    nr[i] = 0;
  
  for (i = 0; i < fill_iter; i++) {
    const unsigned int bucket = loc >> (32 - NBUCKETK);
    cacheline[bucket * ENTRIES_PER_LINE + (nr[bucket]++)] = (unsigned int) loc;
    
    if (nr[bucket] == ENTRIES_PER_LINE) {
      __v2di *cl = (__v2di *) &cacheline[bucket * ENTRIES_PER_LINE];
      __v2di *b = (__v2di *) buckets[bucket];
      __builtin_ia32_movntdq (&(b[0]), cl[0]);
      if (CLSIZE > 16)
        __builtin_ia32_movntdq (&(b[1]), cl[1]);
      if (CLSIZE > 32)
        __builtin_ia32_movntdq (&(b[2]), cl[2]);
      if (CLSIZE > 48)
        __builtin_ia32_movntdq (&(b[3]), cl[3]);
      if (CLSIZE > 64)
        __builtin_ia32_movntdq (&(b[4]), cl[4]);
      if (CLSIZE > 80)
        __builtin_ia32_movntdq (&(b[5]), cl[5]);
      if (CLSIZE > 96)
        __builtin_ia32_movntdq (&(b[6]), cl[6]);
      if (CLSIZE > 112)
        __builtin_ia32_movntdq (&(b[7]), cl[7]);
      buckets[bucket] += ENTRIES_PER_LINE;
      nr[bucket] = 0;
    }
    loc = RAND(loc);
  }
  free (cacheline_alloc);
}

double time()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (ru.ru_utime.tv_sec * 1000000 + ru.ru_utime.tv_usec) / 1.e6;
}

unsigned int * 
bucket_start(unsigned int * const bucket_storage, size_t i)
{
  return bucket_storage + i * (BUCKET_SIZE + bucket_skip);
}

void printpointers(unsigned int ** const buckets, unsigned int * const bucket_storage)
{
  int i;

  printf("Bucket pointers:\n");
  for (i = 0; i < NBUCKET; i++) {
    char *start = (char *) bucket_start(bucket_storage, i);
    char *end = (char *) (buckets[i]);
    printf ("%p - %p, len = %zu\n", start, end, end - start);
    assert ((size_t)(end-start) < BUCKET_SIZE * sizeof(unsigned int));
  }
  printf("\n");
}


int main(int argc, char **argv) {
  int iter = 10;
  unsigned int **buckets;
  unsigned int *bucket_storage;
  int i, j, which, do_which[2] = {1, 1};
  double starttime, endtime;
  fillfunc_t fillfunc[2] =  {&fill, &fill_cacheline};
  const char *which_str[2] = {"Without", "With"};
  
  while (argc > 1) {
    if (strcmp (argv[1], "-cl") == 0) {
      do_which[0] = 0;
      argc--;
      argv++;
      continue;
    }
    
    if (strcmp (argv[1], "-ncl") == 0) {
      do_which[1] = 0;
      argc--;
      argv++;
      continue;
    }

    if (argc > 2 && strcmp (argv[1], "-i") == 0) {
      iter = atoi (argv[2]);
      argc -= 2;
      argv += 2;
      continue;
    }
  }
  
  buckets = (unsigned int **) malloc(NBUCKET * sizeof (unsigned int *));
  {
    void *alloc;
    i = posix_memalign (&alloc, 4096, NBUCKET * (BUCKET_SIZE + bucket_skip) * sizeof(unsigned int));
    if (i != 0)
      exit(EXIT_FAILURE);
    bucket_storage = (unsigned int *) alloc;
  }
  /* Have system allocate everything before we start timing */
  memset (bucket_storage, 0, NBUCKET * (BUCKET_SIZE + bucket_skip) * sizeof(unsigned int));
  
  for (which = 0; which < 2; which++) {
    if (do_which[which]) {
      starttime = time();
      for (j = 0; j < iter; j++) {
        for (i = 0; i < NBUCKET; i++)
          buckets[i] = bucket_start (bucket_storage, i);

        fillfunc[which](buckets, j + 100);
      }
      endtime = time();
      printf ("%s cache line buffers: iter=%d, fill_iter=%lu, %f seconds\n", 
              which_str[which], iter, fill_iter, endtime - starttime);

      // printpointers(buckets, bucket_storage);
    }
  }

  exit (EXIT_SUCCESS);
}
