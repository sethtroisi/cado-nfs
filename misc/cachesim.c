#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define RAND(x) (5*(x)+1)
#define UNUSED 0xFFFFFFFFFFFFFFFFUL

typedef struct {
  int associativity, lines, linesize, indices, verbose;
  size_t *cache;
  int *lru;
  unsigned long hits, misses;
} cache_t;

typedef struct {
  size_t pagesize, nr_pages;
  size_t *pages;
} pages_t;

int find_lru(const cache_t *cache, const int index)
{
  int i;
  
  assert (index < cache->indices);
  /* Check the lines from the associative set to find the least recently used one.
     The age is increasing, so age == associativity - 1 is the oldest */
  for (i = 0; i < cache->associativity; i++) {
    int line = index + i * cache->indices;
    int age = cache->lru[line];
    assert (0 <= age && age < cache->associativity);
    if (age == cache->associativity - 1)
      return line;
  }
  abort();
}


void update_lru(cache_t *cache, const int line)
{
  const int index = line % cache->indices;
  int i, age;
  
  assert (line < cache->lines);
  age = cache->lru[line];
  assert (age < cache->associativity);
  
  for (i = 0; i < cache->associativity; i++) {
    int update_line = index + i * cache->indices;
    if (cache->lru[update_line] < age)
      cache->lru[update_line]++;
  }
  cache->lru[line] = 0;
}


void check_lru(const cache_t *cache)
{
  int i, j;
  int *counts;
  
  counts = malloc(cache->associativity * sizeof(int));
  
  for (i = 0; i < cache->indices; i++) {
    memset (counts, 0, cache->associativity * sizeof(int));
    for (j = 0; j < cache->associativity; j++) {
      int line = i + j * cache->indices;
      int age = cache->lru[line];
      assert (0 <= age && age < cache->associativity);
      counts[age]++;
    }
    for (j = 0; j < cache->associativity; j++) {
      assert (counts[j] == 1);
    }
  }
  
  free(counts);
}

/* Search the cache tags in the associative set to which address maps */
int is_in_cache(cache_t *cache, size_t line_address)
{
  const int index = line_address % cache->indices;
  int i;
  
  for (i = 0; i < cache->associativity; i++) {
    int line = index + i * cache->indices;
    if (cache->cache[line] == line_address)
      return line;
  }
  return cache->lines;
}

/* Returns 1 if an address was found in cache. If not, the address is added 
   to the cache, replacing another entry */
void access_cache(cache_t *cache, size_t memory_address)
{
  const size_t line_address = memory_address / cache->linesize;
  const int index = line_address % cache->indices;
  size_t oldaddr;
  int line;
  
  line = is_in_cache(cache, line_address);
  if (line < cache->lines) {
    cache->hits++;
    if (cache->verbose) {
      printf ("Access to address %zu (line address %zu) hit in cache line %d\n", 
              memory_address, line_address, line);
    }
    return;
  }

  /* Get the oldest line in this index */
  line = find_lru(cache, index);
  oldaddr = cache->cache[line];
  /* and set it to this address */
  cache->cache[line] = line_address;
  update_lru(cache, line);
  cache->misses++;
  if (cache->verbose) {
    if (oldaddr == UNUSED)
      printf ("Access to address %zu (line address %zu) missed in cache, entered in line %d (was unused)\n", 
              memory_address, line_address, line);
    else
      printf ("Access to address %zu (line address %zu) missed in cache, entered in line %d, replacing line address %zu\n", 
              memory_address, line_address, line, oldaddr);
  }
}

void cache_init(cache_t *cache, const int lines, const int associativity, 
                const int linesize, const int verbose) 
{
  int i;
  assert (lines % associativity == 0);
  
  cache->associativity = associativity;
  cache->lines = lines;
  cache->linesize = linesize;
  cache->indices = lines / associativity;
  cache->cache = malloc(lines * sizeof (size_t));
  cache->lru = malloc(lines * sizeof (int));
  cache->hits = 0;
  cache->misses = 0;
  cache->verbose = verbose;

  /* Init the LRU buffer to a state that satisfies our invariant:
     among the lines in each index, each age in [0, associativity-1] 
     occurs exactly once */
  for (i = 0; i < lines; i++) {
    cache->lru[i] = i / cache->indices;
  }
  
  for (i = 0; i < lines; i++) {
    cache->cache[i] = UNUSED;
  }
}

void cache_printstats(cache_t *cache)
{
  int i;
  
  printf ("Addresses in cache:\n");
  for (i = 0; i < cache->lines; i++) {
    printf ("%zu", cache->cache[i]);
    if (i % 20 == 19)
      printf ("\n");
    else
      printf ("\t");
  }
  
  printf ("\n%lu hits, %lu misses\n", cache->hits, cache->misses);
}

void cache_clear(cache_t *cache) 
{
  free (cache->cache);
  cache->cache = NULL;
  free (cache->lru);
  cache->lru = NULL;
}

void pages_init(pages_t *pages, size_t pagesize, size_t nr_pages)
{
  size_t i;
  pages->pagesize = pagesize;
  pages->nr_pages = nr_pages;
  
  pages->pages = malloc(nr_pages * sizeof (size_t));
  for (i = 0; i < nr_pages; i++)
    pages->pages[i] = i; /* Simple linear layout for now */
}

void pages_clear(pages_t *pages)
{
  free (pages->pages);
  pages->pages = NULL;
}


int 
page_translate(pages_t *pages, size_t address)
{
  size_t vpage, PFN;
  
  vpage = address / pages->pagesize;
  assert (vpage < pages->nr_pages);
  PFN = pages->pages[vpage];
  
  return PFN * pages->pagesize + address % pages->pagesize;
}


void 
simulation(cache_t *cache, pages_t *pages, int argc, char **argv) {
  int i;
  int bucket_start = 123 * pages->pagesize;
  int region_start = 234 * pages->pagesize;
  int update;
  size_t update_stride = 4;
  size_t logical_address, physical_address;
  
  if (argc > 1) {
    update_stride = (size_t) atoi(argv[1]);
  }
  printf ("Using update stride %zu\n", update_stride);
  
  for (i = 0; i < 1000000; i++) {
    /* Read an update from the bucket */
    logical_address = bucket_start + update_stride * i;
    physical_address = page_translate(pages, logical_address);
    access_cache (cache, physical_address);
    
    /* The address where update hits is assumed to be random in the bucket region */
    update = random() % (cache->lines * cache->linesize);
    
    /* Apply that update to the bucket region */
    logical_address = update + region_start;
    physical_address = page_translate(pages, logical_address);
    access_cache (cache, region_start + update);
  }
}


int main(int argc, char **argv) {
  cache_t cache[1];
  pages_t pages[1];
  int verbose = 0;
  int lines, associativity;
  int linesize = 64;
  
  if (argc > 1 && strcmp(argv[1], "-v") == 0) {
    argc--;
    argv++;
    verbose = 1;
  }
  
  if (argc < 4) {
    printf ("%s [-v] cachelines associativity linesize [parameters for simulation]\n", argv[0]);
    exit (EXIT_FAILURE);
  }
  
  lines = atoi (argv[1]); 
  argc--; 
  argv++;
  associativity = atoi (argv[1]); 
  argc--; 
  argv++;
  linesize = atoi (argv[1]); 
  argc--; 
  argv++;
  
  cache_init (cache, lines, associativity, linesize, verbose);
  check_lru(cache);
  pages_init(pages, 4096, 2048);

  simulation(cache, pages, argc, argv);
  
  check_lru(cache);
  cache_printstats(cache);

  cache_clear (cache);
  pages_clear(pages);
  exit (EXIT_FAILURE);
}
