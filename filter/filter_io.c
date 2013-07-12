#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"

static relation_stream rs;
static char *pmin, *pminlessone;

static char antebuffer[PATH_MAX];         /* "directory/antebuffer" or "cat" */
static char rep_cado[PATH_MAX];           /* directory of cado */

inline int
is_finish ()
{
  return end_insertRelation;
}

void
set_antebuffer_path (char * argv0, const char *path_antebuffer)
{
  set_rep_cado (argv0, rep_cado);
  search_antebuffer (rep_cado, path_antebuffer, antebuffer);
}


void
prempt_load (prempt_t prempt_data)
{
  char **p_files = prempt_data->files, *pprod;
  FILE *f;
  size_t try_load, load;
  char *pmax = &(prempt_data->buf[PREMPT_BUF]);

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

  while (*p_files)
  {
    if (!(f = popen (*p_files, "r")))
    {
      fprintf (stderr, "prempt_load: popen error. %s\n", strerror (errno));
      exit (1);
    }
#if DEBUG >= 1
    char *p, *l;
    /* Write the command that read the files (one file per line) */
    p = strchr(*p_files, '/');
    if (p)
    {
      *p = 0;
      fprintf (stderr, "%s\n", *p_files);
      *p = '/';
      for ( ; ; )
      {
        l = strchr(p, ' ');
        if (l)
        {
          *l = 0;
          fprintf(stderr, "   %-70s\n", p);
          *l = ' ';
          p = strchr(&(l[1]), '/');
          if (!p)
          {
            fprintf(stderr, "%s\n", &(l[1]));
            break;
          }
        }
        else
        {
          fprintf(stderr, "   %-70s\n", p);
          break;
        }
      }
    }
#endif

    pprod = (char *) prempt_data->pprod;
    for ( ; ; )
    {
      if ((pprod != prempt_data->pcons) &&
    (((((PREMPT_BUF + ((size_t) prempt_data->pcons))) - ((size_t) pprod)) &
      (PREMPT_BUF - 1)) <= PREMPT_ONE_READ))
        {
          NANOSLEEP();
        }
      else
      {
        try_load =
         MIN((PREMPT_BUF + ((size_t)pmin)) - ((size_t) pprod), PREMPT_ONE_READ);
        if ((load = fread (pprod, 1, try_load, f)))
        {
          pprod += load;
          if (pprod == pmax)
            pprod = pmin;
          prempt_data->pprod = pprod;
        }
        else if (feof(f))  // we go to the next batch of files
        {
          pclose(f);
          free(*p_files);
          *p_files++ = NULL;
          break;
        }
        else // error
        {
          fprintf (stderr, "prempt_load: load error (%s) from\n%s\n",
                           strerror (errno), *p_files);
          exit (1);
        }
      }
    }
  }
  prempt_data->end = 1;
  pthread_exit (NULL);
}

static inline void
realloc_buffer_primes (buf_rel_t *buf)
{
  if (buf->nb_alloc == NB_PRIMES_OPT)
  {
    buf->nb_alloc += buf->nb_alloc >> 1;
    prime_t *p = buf->primes;
    SMALLOC(buf->primes, buf->nb_alloc, "realloc buffer primes");
    memcpy (buf->primes, p, NB_PRIMES_OPT * sizeof(prime_t));
  }
  else
  {
    buf->nb_alloc += buf->nb_alloc >> 1;
    buf->primes =(prime_t*)realloc(buf->primes, buf->nb_alloc * sizeof(prime_t));
  }
#if DEBUG >= 2
    fprintf (stderr, "realloc_buffer_primes: num=%"PRid" nb_alloc=%u\n",
                      buf->num, buf->nb_alloc);
#endif
}

#define LOAD_ONE(P) { c = *P;  \
       P = ((size_t) (P - pminlessone) & (PREMPT_BUF - 1)) + pmin;}


static inline unsigned char
read_one_prime (p_r_values_t * pr, char **p)
{
  unsigned char v, c;
  LOAD_ONE(*p);
  for (*pr = 0 ; (v = ugly[c]) < 16; )
  {
    *pr = (*pr << 4) + v;
    LOAD_ONE(*p);
  }
  ASSERT_ALWAYS(c == ',' || c == '\n');

  return c;
}

static inline unsigned char
read_a_or_b (uint64_t *n, char **p, unsigned char c, int ab_hexa)
{
  unsigned char v;
  if (ab_hexa)
  {
    for (*n = 0 ; (v = ugly[c]) < 16; )
    {
      *n = (*n << 4) + v;
      LOAD_ONE(*p);
    }
  }
  else
  {
    for (*n = 0 ; (v = ugly[c]) < 10; )
    {
      *n = (*n * 10) + v;
      LOAD_ONE(*p);
    }
  }

  return c;
}

static inline unsigned char
read_b (uint64_t *b, char **p, int ab_hexa)
{
  unsigned char c;
  LOAD_ONE(*p);

  c = read_a_or_b (b, p, c, ab_hexa);
  ASSERT_ALWAYS(c == ':');

  return c;
}

static inline unsigned char
read_a (int64_t *a, char **p, int ab_hexa)
{
  unsigned char c;
  uint64_t n;

  LOAD_ONE(*p);
  if (c == '-')
  {
    *a = -1;
    LOAD_ONE(*p);
  }
  else
    *a = 1;

  c = read_a_or_b (&n, p, c, ab_hexa);

  ASSERT_ALWAYS(c == ',');
  *a *= n;

  return c;
}

static inline unsigned char
skip_ab (char **p)
{
  unsigned char c;
  LOAD_ONE(*p)
  while (c != ':')
    LOAD_ONE(*p)

  return c;
}

#if 0
#ifdef STAT
  __stat_weight += nb_primes_read;
#ifdef STAT_VALUES_COEFF
  {
    unsigned int i;
    for (i = 0; i < mybufrel->nb; i++)
    {
      if (bit_vector_getbit(rel_used, (size_t) mybufrel->num))
      {
        if (abs(mybufrel->primes[i].e) > STAT_VALUES_COEFF_LEN)
          __stat_nbcoeffofvalue[0]++;
        else
          __stat_nbcoeffofvalue[abs(mybufrel->primes[i].e)]++;
      }
    }
  }
#endif
#endif
#endif

static inline void
relation_get_fast_abp (prempt_t prempt_data, buf_rel_t *mybufrel)
{
  char *p;
  unsigned int nb_primes_read;
  int i;
  p_r_values_t pr;
  unsigned char c;
  int side = -1;

  p = (char *) prempt_data->pcons;

#ifndef FOR_FFS
  c = read_a (&(mybufrel->a), &p, 0);
  c = read_b (&(mybufrel->b), &p, 0);
#else
  c = read_a (&(mybufrel->a), &p, 1);
  c = read_b (&(mybufrel->b), &p, 1);
#endif

  nb_primes_read = 0;
  for ( c = 0 ; ; )
  {
    if (c == '\n')
      break;
    if (c == ':')
      side++;

    c = read_one_prime(&pr, &p);

    for (i = nb_primes_read-1; i >= 0; i-- )
      if (mybufrel->primes[i].p == pr)
        mybufrel->primes[i-1].e++;

    if (i == -1)
    {
      if (mybufrel->nb_alloc == nb_primes_read)
        realloc_buffer_primes (mybufrel);

      mybufrel->primes[nb_primes_read++] = (prime_t) {.h=side, .p=pr, .e=1};
    }
  }

  mybufrel->nb = nb_primes_read;
  mybufrel->nb_above_min_index = nb_primes_read + 1;
}

static inline void
relation_get_fast_abh (prempt_t prempt_data, buf_rel_t *mybufrel)
{
  char *p;
  unsigned int nb_primes_read;
  index_t pr;
  unsigned char c;

  p = (char *) prempt_data->pcons;

  c = read_a (&(mybufrel->a), &p, 1);
  c = read_b (&(mybufrel->b), &p, 1);

  nb_primes_read = 0;
  for ( c = 0 ; ; )
  {
    if (c == '\n')
      break;

    c = read_one_prime(&pr, &p);

    if (nb_primes_read > 0 && mybufrel->primes[nb_primes_read-1].h == pr)
        mybufrel->primes[nb_primes_read-1].e++;
    else
    {
      if (mybufrel->nb_alloc == nb_primes_read)
        realloc_buffer_primes (mybufrel);

      mybufrel->primes[nb_primes_read++] = (prime_t) {.h=pr, .p=0, .e=1};
    }
  }

  mybufrel->nb = nb_primes_read;
  mybufrel->nb_above_min_index = nb_primes_read + 1;
}

static inline void
relation_get_fast_ab (prempt_t prempt_data, buf_rel_t *mybufrel)
{
  char *p;
  unsigned char c;

  p = (char *) prempt_data->pcons;

  c = read_a (&(mybufrel->a), &p, 1);
  c = read_b (&(mybufrel->b), &p, 1);

  while (c != '\n')
    LOAD_ONE(p);
}

static inline void
relation_get_fast_abline (prempt_t prempt_data, buf_rel_t *mybufrel)
{
  char *p, *begin;
  unsigned char c;
  unsigned int i = 0;

  p = (char *) prempt_data->pcons;
  begin = p;

#ifndef FOR_FFS
  c = read_a (&(mybufrel->a), &p, 0);
  c = read_b (&(mybufrel->b), &p, 0);
#else
  c = read_a (&(mybufrel->a), &p, 1);
  c = read_b (&(mybufrel->b), &p, 1);
#endif

  do
  {
    LOAD_ONE(begin);
    mybufrel->line[i++] = c;
  } while (c != '\n');

  mybufrel->line[i] = '\0';
}

static inline void
relation_get_fast_line (prempt_t prempt_data, buf_rel_t *mybufrel)
{
  char *p;
  unsigned char c;
  unsigned int i = 0;
  unsigned int nb_primes = 0;

  p = (char *) prempt_data->pcons;

  do
  {
    LOAD_ONE(p);
    mybufrel->line[i++] = c;
    if (c == ':' || c == ',')
      nb_primes++;
      
  } while (c != '\n');

  mybufrel->line[i] = '\0';
  mybufrel->nb = nb_primes;
}

static inline void
relation_get_fast_hmin (prempt_t prempt_data, buf_rel_t *mybufrel, index_t min)
{
  char *p;
  unsigned int nb_primes_read;
  index_t pr;
  unsigned char c;
  weight_t nb_above_min_index = 1; // count the -1 at the end of the relations
                                   // in rel_compact

  p = (char *) prempt_data->pcons;

  c = skip_ab(&p);

  nb_primes_read = 0;
  for ( c = 0 ; ; )
  {
    if (c == '\n')
      break;

    c = read_one_prime(&pr, &p);

    if (nb_primes_read > 0 && mybufrel->primes[nb_primes_read-1].h == pr)
        mybufrel->primes[nb_primes_read-1].e++;
    else
    {
      if (mybufrel->nb_alloc == nb_primes_read)
        realloc_buffer_primes (mybufrel);

      nb_above_min_index += (weight_t) (pr >= min);
      mybufrel->primes[nb_primes_read++] = (prime_t) { .h = pr, .p = 0, .e = 1};
    }
  }

  mybufrel->nb = nb_primes_read;
  mybufrel->nb_above_min_index = nb_above_min_index;
}

info_mat_t
process_rels (char **fic, void* (*callback_fct)(buf_arg_t *),
              void* (*thread_root)(fr_t *), index_t min_index, FILE **outfd,
              bit_vector_ptr rel_used, unsigned int step)
{
  char *pcons, *pcons_old, *pcons_max, *p, **ff;
  pthread_attr_t attr;
  pthread_t thread_load, thread_callback, thread_fr[(1<<NFR)];
  fr_t fr[(1<<NFR)];
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int length_line, i, k;
  int err, needr = (thread_root != NULL);
  char c;
  unsigned int need_del_rels = 0;

  buf_arg_t buf_arg;
  memset (&buf_arg, 0, sizeof(buf_arg_t));
  buf_arg.min_index = min_index;
  buf_arg.fd = outfd;
  buf_arg.rel_used = rel_used;

  if (fic[0] == NULL) // to avoid a seg fault in the case the filelist is empty
  {
    fprintf (stderr, "End of read: 0 relation (filelist was empty)\n");
    return buf_arg.info;
  }

  ASSERT_ALWAYS (step <= MAX_STEP);

  if (step == STEP_PURGE_PASS2)
  {
    ASSERT_ALWAYS(buf_arg.fd != NULL);
    ASSERT_ALWAYS(buf_arg.rel_used != NULL);
    if (buf_arg.fd[1] != NULL)
      need_del_rels = 1;
  }

  buf_arg.rels = (buf_rel_t *) malloc (SIZE_BUF_REL * sizeof (buf_rel_t));
  ASSERT_ALWAYS(buf_arg.rels != NULL);
  fprintf (stderr, "Allocated buffer for rels of %zuMb\n",
                                  (SIZE_BUF_REL * sizeof (buf_rel_t)) >> 20);

#ifdef STAT
  __stat_weight = 0;
#ifdef STAT_VALUES_COEFF
  for (int i = 0; i < STAT_VALUES_COEFF_LEN; i++)
    __stat_nbcoeffofvalue[i] = 0;
#endif
#endif

  MEMSETZERO(fr, (1<<NFR));
  for (i = (1<<NFR); i--; )
    fr[i].buf_data = buf_arg.rels;

  MEMSETZERO(buf_arg.rels, SIZE_BUF_REL);
  for (i = SIZE_BUF_REL; i--; )
  {
    buf_arg.rels[i].primes = buf_arg.rels[i].primes_data;
    buf_arg.rels[i].nb_alloc = NB_PRIMES_OPT;
  }

  if (step == STEP_DUP1 || step == STEP_PURGE_PASS2)
  {
    size_t tmp_alloc = RELATION_MAX_BYTES * sizeof(char);
    for (i = SIZE_BUF_REL; i--; )
    {
      buf_arg.rels[i].line = (char *) malloc (tmp_alloc);
      ASSERT_ALWAYS (buf_arg.rels[i].line != NULL);
    }
  }

  end_insertRelation = 0;
  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  relation_stream_init (rs);
  rs->pipe = 1;
  length_line = 0;
  prempt_data->files = prempt_open_compressed_rs (antebuffer, fic);

  SMALLOC (prempt_data->buf, PREMPT_BUF, "prempt_scan_relations prempt_data");
  pmin = prempt_data->buf;
  pminlessone = pmin - 1;
  prempt_data->pcons = pmin;
  prempt_data->pprod = pmin;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);


  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data)))
  {
    fprintf (stderr, "prempt_scan_relations: pthread_create error 1: %d. %s\n",
                     err, strerror (errno));
    exit (1);
  }

  err = pthread_create (&thread_callback, &attr, (void* (*)(void*))callback_fct,
                        (void*) &buf_arg);
  if (err)
  {
    fprintf (stderr, "prempt_scan_relations: pthread_create error 2: %d. %s\n",
                     err, strerror (errno));
    exit (1);
  }

  if (needr)
  {
    for (i = 0; i < (1<<NFR); i++)
    {
      err = pthread_create (&(thread_fr[i]), &attr,(void*)thread_root, &(fr[i]));
      if (err)
      {
        fprintf (stderr, "prempt_scan_relations: pthread_create error 3: %d. "
                         "%s\n", err, strerror (errno));
        exit (1);
      }
    }
  }

  pcons = (char *) prempt_data->pcons;
  for ( ; ; )
    {
      rs->pos += length_line;
      length_line = 0;
      prempt_data->pcons = pcons;

      while (pcons == prempt_data->pprod)
        if (!prempt_data->end)
          NANOSLEEP();
        else if (pcons == prempt_data->pprod)
          goto end_of_files;

      if (pcons == prempt_data->pprod + sizeof(*pcons))
        NANOSLEEP();

      rs->lnum++;
      if (*pcons != '#')
      {
        while ((((PREMPT_BUF + ((size_t) prempt_data->pprod)) -
               ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES) &&
     !prempt_data->end)
        NANOSLEEP();

        if (pcons > prempt_data->pprod)
        {
          p = &(pcons_max[-1]);
          c = *p;
          *p = '\n';
          pcons_old = pcons;
          while (*pcons++ != '\n');
            *p = c;

          length_line = (pcons - pcons_old);
          if (pcons <= p)
            goto testendline;

          pcons = pmin;
          if (c == '\n')
            goto testendline;
        }

        p = &(((char *) prempt_data->pprod)[-1]);
        c = *p;
        *p = '\n';
        pcons_old = pcons;
        while (*pcons++ != '\n');
        *p = c;
        length_line += (pcons - pcons_old);

      testendline:
        if (c != '\n' && pcons == prempt_data->pprod)
        {
          fprintf (stderr, "prempt_scan_relations: "
                           "the last line has not a carriage return\n");
          exit(1);
        }
        if (length_line > ((unsigned int) RELATION_MAX_BYTES))
        {
          fprintf (stderr, "prempt_scan_relations: relation line size"
                           " (%u) is greater than RELATION_MAX_BYTES (%d)\n",
                           length_line, RELATION_MAX_BYTES);
          exit(1);
        }

        while (cpy_cpt_rel_a == cpt_rel_b + SIZE_BUF_REL)
          NANOSLEEP();

        k = (unsigned int) (cpy_cpt_rel_a & (SIZE_BUF_REL - 1));
        if (cpy_cpt_rel_a + 1 == cpt_rel_b + SIZE_BUF_REL)
          NANOSLEEP();

        buf_arg.rels[k].num = rs->nrels++;

        if (step == STEP_DUP1)
          relation_get_fast_abline (prempt_data, &(buf_arg.rels[k]));
        else if (step == STEP_DUP2_PASS1)
          relation_get_fast_ab (prempt_data, &(buf_arg.rels[k]));
        else if (step == STEP_DUP2_PASS2)
          relation_get_fast_abp (prempt_data, &(buf_arg.rels[k]));
        else if (step == STEP_PURGE_PASS1)
        {
          if (buf_arg.rel_used != NULL &&
              bit_vector_getbit(buf_arg.rel_used, (size_t) buf_arg.rels[k].num))
            relation_get_fast_hmin (prempt_data, &(buf_arg.rels[k]),
                                                     buf_arg.min_index);
        }
        else if (step == STEP_MERGE)
          relation_get_fast_hmin (prempt_data, &(buf_arg.rels[k]),
                                                     buf_arg.min_index);
        else if (step == STEP_PURGE_PASS2)
        {
          if (need_del_rels ||
              bit_vector_getbit(buf_arg.rel_used, (size_t) buf_arg.rels[k].num))
                relation_get_fast_line (prempt_data, &(buf_arg.rels[k]));
        }
        else if (step == STEP_REPLAY)
          relation_get_fast_abh (prempt_data, &(buf_arg.rels[k]));

        /* Delayed find root computation by block of 1<<NNFR */
        if (cpy_cpt_rel_a && !(k & ((1<<NNFR)-1)))
        {
          if (needr)
          {
            i = (k>>NNFR) & ((1<<NFR)-1);
            while (fr[i].ok)
              NANOSLEEP();
            if (k)
            {
              fr[i].num = k - (1<<NNFR);
              fr[i].end = k - 1;
            }
            else
            {
              fr[i].num = SIZE_BUF_REL - (1<<NNFR);
              fr[i].end = SIZE_BUF_REL - 1;
            }
            fr[i].ok = 1;
            if (cpy_cpt_rel_a > (1<<(NFR+NNFR)))
              cpt_rel_a = cpy_cpt_rel_a - (1<<(NFR+NNFR));
          }
          else if (cpy_cpt_rel_a > (1<<(NNFR)))
            cpt_rel_a = cpy_cpt_rel_a - (1<<(NNFR));
        }
        cpy_cpt_rel_a++;
      }
      else
      {
        do
        {
          while (pcons == prempt_data->pprod)
          {
            if (!prempt_data->end)
              NANOSLEEP();
            else if (pcons == prempt_data->pprod)
            {
              fprintf (stderr, "prempt_scan_relations: at the end of files,"
                               " a line without \\n ?\n");
              exit (1);
            }
          }

          p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod
                                             : pcons_max) - 1;
          c = *p;
          *p = '\n';
          pcons_old = pcons;
          while (*pcons++ != '\n');
            *p = c;

          length_line += (pcons - pcons_old);
          err = (pcons > p && c != '\n');
          if (pcons == pcons_max)
            pcons = pmin;
        } while (err);
      }
  }

 end_of_files:
  if (needr)
  {
    if (cpy_cpt_rel_a)
    {
      k = (unsigned int) ((cpy_cpt_rel_a - 1) & (SIZE_BUF_REL - 1));
      if (k & ((1<<NNFR)-1))
      {
        i = ((k>>NNFR)+1) & ((1<<NFR)-1);
        while (fr[i].ok)
          NANOSLEEP();
        fr[i].num = k & ~((1<<NNFR)-1);
        fr[i].end = k;
        fr[i].ok = 1;
      }
    }
    for (i = 0; i < (1<<NFR); i++)
    {
      while (fr[i].ok)
        NANOSLEEP();
      fr[i].ok = 2;
      pthread_join(thread_fr[i], NULL);
    }
  }

  cpt_rel_a = cpy_cpt_rel_a;
  while (cpy_cpt_rel_a != cpt_rel_b)
    NANOSLEEP();

  end_insertRelation = 1;
  pthread_join(thread_callback, NULL);
  /* if (pthread_tryjoin_np (thread_load, NULL)) */
  pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);

  free (prempt_data->buf);
  for (ff = prempt_data->files; *ff; free(*ff++));
  free (prempt_data->files);
  for (i = SIZE_BUF_REL; i-- ; )
    if (buf_arg.rels[i].nb_alloc != NB_PRIMES_OPT )
      SFREE(buf_arg.rels[i].primes);

  if (step == STEP_DUP1 || step == STEP_PURGE_PASS2)
    for (i = SIZE_BUF_REL; i--; )
      free(buf_arg.rels[i].line);

  free(buf_arg.rels);

  relation_stream_trigger_disp_progress(rs);
  fprintf (stderr, "End of read: %"PRid" relations in %.1fs -- %.1f MB/s -- "
                   "%.1f rels/s\n", rs->nrels, rs->dt, rs->mb_s, rs->rels_s);

  buf_arg.info.nrels = rs->nrels;

  relation_stream_clear(rs);

#ifdef STAT
  fprintf (stderr, "STAT: W_total=%"PRIu64"\n", __stat_weight);
#ifdef STAT_VALUES_COEFF
  {
    uint64_t tot = 0;
    int i;
    for (i = 0; i <= STAT_VALUES_COEFF_LEN; i++)
      tot += __stat_nbcoeffofvalue[i];

    fprintf (stderr, "STAT: coeffs with non-zero values: %"PRIu64"\n",
                     __stat_weight);
    for (i = 1; i <= STAT_VALUES_COEFF_LEN; i++)
      fprintf (stderr, "STAT: coeffs of abs value %d: %"PRIu64" (%.2f%%)\n",
               i, __stat_nbcoeffofvalue[i],
               100 * (double) __stat_nbcoeffofvalue[i]/tot);
    fprintf (stderr, "STAT: coeffs of abs value > %d: %"PRIu64" (%.2f%%)\n",
             STAT_VALUES_COEFF_LEN, __stat_nbcoeffofvalue[0],
             100 * (double) __stat_nbcoeffofvalue[0]/tot);
  }
#endif
#endif

  return buf_arg.info;
}


void test_and_print_progress_now ()
{
  if (relation_stream_disp_progress_now_p(rs))
    fprintf(stderr, "Read %"PRid" relations in %.1fs -- %.1f MB/s -- "
                    "%.1f rels/s\n", rs->nrels, rs->dt, rs->mb_s, rs->rels_s);

}
