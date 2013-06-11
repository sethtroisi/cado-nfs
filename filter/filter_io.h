#ifndef FILTER_IO_H_
#define FILTER_IO_H_

#define MAX_FILES 1000000

/*
  1<<NFR = Number of computation threads.
  The best is 2 (NFR = 1).
  if you have 2 cores or less, you could try NFR=1.
  if you have >4 cores, or hyperthreading, AND very slow cores (< 2 Ghz): NFR=2
*/
#define NFR (1)
/* 1<<NNFR = Number of sentences computed of (find root + hashkey of each term)
   computed in one pass. If the number is too big, the buffer between these
   threaded computations and insertRelation (which cannot be parallelized,
   because the hashkey is not bijective) is too big (memory waste) and
   the pipe-line is slow to start;
   If the number is too small, the computation threads are too often in
   nanosleep to keep CPU.
   NNFR=8 to 16, the greatest is the fastest for fast processors & fast
   memory; 14 seems to be the best choice (maximal speed for medium memory use).
*/
#define NNFR (14)

/* Size of relations buffer between parsing & processing.
   CAREFUL! SIZE_BUF_REL must be greater (at least double) than (1<<(NNFR+1)).
   Stocks the sentences precomputed but not insered. 
   About 64K sentences for the optimal.
*/
#define SIZE_BUF_REL MAX((1<<(NNFR+2)),(1<<6))

/* The realistic minimal non-CPU waiting with nanosleep is about
   10 to 40 µs (1<<13 for nanosleep).
   But all the I/O between the threads have been buffered,
   and a thread do a nanosleep only if its buffer is empty.
   So I use here ~2ms (1<<21) to optimize CPU scheduler.
   Max pause is about 4 to 8ms (1<<22, 1<<23); after the program
   is slow down.
*/
static const struct timespec wait_classical = { 0, 1<<21 };

/* What part of relations do we read */
#define NEEDED_ABP 0 //read a,b:p0,...,pk:p'0,...,p'l
#define NEEDED_AB 1  //read a,b
#define NEEDED_H 2 //read a,b:h0,...,hk
#define NEEDED_HMIN 3 //read h0,...,hk only if hi >= min_index
#define NEEDED_ABH 4 //read a,b:h0,...,hk
#define MAX_NEEDED 4

typedef struct {
  int64_t a;
  uint64_t b;
  prime_t *primes; /*if nb<=NB_PRIME_OPT, contains the address of primes_data*/
  prime_t primes_data[NB_PRIMES_OPT];
  weight_t nb;           /* number of primes */
  weight_t nb_alloc;     /* allocated space for primes */
  weight_t nb_above_min_index; /* nb of primes above min_index, must be <=nb */
  index_t num;          /* Relation number */
} buf_rel_t;

typedef struct {
  index_t nrels;       // nb of rels read
  index_t nprimes;     // nb of primes read
  buf_rel_t *buf_data; // buffer for I/O
  unsigned int needed; // ABP, AB, H, ABHMIN
  unsigned int needr;  // non-null if we need to compute r
  index_t min_index;   // store only primes with index >= min_index
  double W;            // Weight (nb of non-zero) of primes >= min_index
  FILE *f_remaining;   // File for writing remaining rels
  FILE *f_deleted;     // File for writing deleted rels
  bit_vector_ptr rel_used; // If not all rels are processed
} buf_arg_t;

typedef struct {
  volatile unsigned int ok;
  unsigned int num, end;
  buf_rel_t *buf_data; // buffer for I/O
} fr_t;

/* For the multithread sync */
volatile unsigned long cpt_rel_a;
volatile unsigned long cpt_rel_b;
static volatile unsigned int end_insertRelation = 0;

static const unsigned char ugly[256] = {
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   255, 255, 255, 255, 255, 255,
 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 
 };

#ifndef HAVE_NANOSLEEP
int nanosleep(const struct timespec *req, struct timespec *rem);
#endif
void prempt_load (prempt_t);
void relation_stream_get_fast (prempt_t, buf_rel_t *, index_t);
int prempt_scan_relations (char **, void* (*)(buf_arg_t *), buf_arg_t *);
void test_and_print_progress_now ();
int is_finish ();
void set_antebuffer_path (char *, const char *);


#endif /* FILTER_IO_H_ */

