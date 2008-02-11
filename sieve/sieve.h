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

