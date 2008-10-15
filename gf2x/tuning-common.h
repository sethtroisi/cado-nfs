#ifndef TUNING_COMMON_H_
#define TUNING_COMMON_H_

extern double mulstep;
extern FILE * rp;

#ifndef MAX
#define MAX(a,b)        ((a)<(b) ? (b) : (a))
#endif

#ifndef MIN
#define MIN(a,b)        ((a)<(b) ? (b) : (a))
#endif

#define MINTIME 0.5		/* time resolution */

#define TIME(x, i)				\
  { long j, k = 1;				\
    double s0 = seconds ();			\
    do {					\
      for (j = 0; j < k; j++) (i);		\
      k = 2 * k;				\
      x = seconds () - s0;			\
    } while (x < MINTIME);			\
    (x) = (x) / (k - 1);			\
  }

#ifdef __cplusplus
extern "C" {
#endif

extern void random_wordstring(unsigned long * a, long n);

extern void check(const unsigned long *a, long m,
	   const unsigned long *b, long n,
	   const char * cname, const unsigned long *c,
           const char * dname, const unsigned long *d);

extern void set_tuning_output();
extern int handle_tuning_mulstep(int * p_argc, char *** p_argv);
extern int handle_tuning_outfile(int * p_argc, char *** p_argv);
#ifdef __cplusplus
}
#endif

#endif	/* TUNING_COMMON_H_ */
