#include "cado.h"
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include "purgedfile.h"
#include "timing.h"
#include "macros.h"
#include "gzip.h"
#include "portability.h"

void
purgedfile_read_firstline (const char *fname, uint64_t *nrows, uint64_t *ncols)
{
  FILE *f_tmp = fopen_maybe_compressed (fname, "rb");
  if (!f_tmp)
  {
    fprintf(stderr, "%s: %s\n", fname, strerror(errno));
    abort();
  }
  int ret = fscanf(f_tmp, "# %" SCNu64 " %" SCNu64 "", nrows, ncols);
  ASSERT_ALWAYS (ret == 2);
  fclose_maybe_compressed(f_tmp, fname);
}
