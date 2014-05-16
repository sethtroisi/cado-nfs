#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

void
stats_init (stats_data_t r, FILE *f, uint8_t max_log_report, const char *verb,
            const char *name, const char *abbrv)
{
  r->last_report = 0;
  r->t0 = wct_seconds();
  r->out = f;
  r->log_report = MIN_LOG_REPORT;
  if (max_log_report >= MIN_LOG_REPORT)
    r->max_log_report = max_log_report;
  else
    r->max_log_report = MIN_LOG_REPORT;
  r->name = name;
  r->abbrv = abbrv;
  r->verb = verb;
}

/* Return 1 if more than 2^r->log_report relations were read since last
   progress report.
   Otherwise return 0 */
int
stats_test_progress (stats_data_t r, uint64_t i)
{
  if ((i >> r->log_report) != r->last_report)
    return 1;
  else
    return 0;
}

/* Print a line of the form: 
      (Done: )(verb) N (name) in Xs -- S(abbrv)/s
*/
void
stats_print_progress (stats_data_t r, uint64_t i, int end)
{
  double t, dt, speed;
  t = wct_seconds();
  dt = t - r->t0;
  speed = dt > 0.01 ? i/dt : 0;
  const char * prefix = (end) ? "Done: " : "";
  fprintf(r->out, "%s%s %" PRIu64 " %s in %.1fs -- %.1f %s/s\n",
                  prefix, r->verb, i, r->name, dt, speed, r->abbrv);
  fflush(r->out);
  if (r->log_report < r->max_log_report)
    r->log_report++;
  r->last_report = (i >> r->log_report);
}

/* Print a line of the form: 
      (Done: )(verb) N (name) in Xs -- S MB/s -- S' (abbrv)/s
*/
void
stats_print_progress_with_MBs (stats_data_t r, uint64_t i, size_t nByte, int end)
{
  double t, dt, speed, mb_s;
  t = wct_seconds();
  dt = t - r->t0;
  speed = dt > 0.01 ? i/dt : 0;
  mb_s = dt > 0.01 ? (nByte/dt * 1.0e-6) : 0;
  const char * prefix = (end) ? "Done: " : "";
  fprintf(r->out, "%s%s %" PRIu64 " %s in %.1fs -- %.1f MB/s -- %.1f %s/s\n",
                  prefix, r->verb, i, r->name, dt, mb_s, speed, r->abbrv);
  fflush(r->out);
  if (r->log_report < r->max_log_report)
    r->log_report++;
  r->last_report = (i >> r->log_report);
}
