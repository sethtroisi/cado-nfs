/* checks endianness */

#define PRIME_HINT uint16_t

typedef struct {
    uint16_t x;
    PRIME_HINT p;
} bucket_update_t;

int
main ()
{
  bucket_update_t b[1];
  uint32_t *u;

  b->x = 1;
  b->p = 2;
  u = (uint32_t*) b;
  if (u == (b->p << 16) + b->x)
    return 0; /* little endian */
  else
    return 1; /* big endian */
}
