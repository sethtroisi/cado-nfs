/* Some data type definition for sieve that aren't used anywhere else */


/* A sieve report, filled in when sieving large factor base primes and
   the new approximate log is below the threshold, or after sieving L1/L2
   primes and scanning the sieve array for reports. 
   The field "a" gets the "a" value of the report, "p" gets the prime that 
   was sieved (or 1 if the report was not filled during sieving) and "l" the 
   approximate log of the remaining cofactor. 
   Let the norm be n = p1 * p2 * ... * pk * P1 * ... * Pn, where p1 ... pk 
   are factor base primes and P1, ..., Pn are large primes. 
   Let c = P1 * ... * Pn. Then the smallest "l" in the sieve reports for a 
   particular "a" value (there can be several reports for one "a"!) satisfies
   |l - log_b(c)| <= SIEVE_PERMISSIBLE_ERROR.
   
   An entry with p == 0 marks the end of the reports array. */
   
typedef struct {
  long a;
  fbprime_t p;
  unsigned char l;
  unsigned char dummy[3];
} sieve_report_t;


typedef struct {
  sieve_report_t *reports;
  size_t alloc;             /* There is space for alloc sieve_report_t's 
			       currently allocated */
  size_t nr;                /* There are currently nr sieve_report_t filled */
} sieve_reportbuffer_t;

#define missinglog_hist_max 32

typedef struct {
  unsigned long
    sum_missingprimes,    /* The sum of missingprimes */
    sum_primes2;          /* The sum of the 2nd largest non-report FB primes */
  
  unsigned int
    guessed_cof_toolarge, /* Number of times the guessed cofactor after 
			   dividing out the missing prime would be > mfb */
    cof_toolarge,         /* Number of times the cofactor after dividing 
			     out all FB primes was > mfb */
    lp_toolarge,          /* Number of times there was a large prime > lpb */
    survivors,            /* Number of surviving candidate relations */
    nr_missingprimes,     /* The number of times we actually found 
			     missingprime */
    nr_primes2,           /* The number of the 2nd largest ... */
    missinglog_hist[missinglog_hist_max], /* The histogram of missinglog */
    missinglog_guessdiscard_hist[missinglog_hist_max]; /* The histogram of
			     missinglog when we discard due to guessed 
			     cofactor > mfb */
} refactor_stats_t;
