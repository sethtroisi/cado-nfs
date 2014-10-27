#ifndef FM
#define FM

#include <stdio.h>
#include <stdbool.h>

#ifndef BRENT12
#define BRENT12 0
#endif

#ifndef MONTY12
#define MONTY12 1
#endif

#ifndef MONTY16
#define MONTY16 2
#endif

#ifndef PM1_METHOD
#define PM1_METHOD 1
#endif

#ifndef PP1_27_METHOD
#define PP1_27_METHOD 2
#endif

#ifndef PP1_65_METHOD
#define PP1_65_METHOD 3
#endif

#ifndef EC_METHOD
#define EC_METHOD 4
#endif

//todo: clean it in the future!

#define PM1 1
#define PP1 2
#define ECM 3

#define NB_METHOD 4
#define NB_CURVE 3

typedef struct fm {
	unsigned long *method;	//contain: METHOD, CURVE, B1, B2
	double *proba;
	double *time;
	int len_method;		//lenght of the method (default:4)
	int len_proba;		//index of array proba 
	int len_time;		//index of array time
	int len_p_min;
	/*
	   The prime number such that : proba[i] and time[i] equal to 
	   the probability and the time to find a prime number of 
	   len_p_min+i bits with our nmethod.
	 */
} fm_t;

fm_t *fm_create(void);

void fm_free(fm_t * t);

unsigned long *fm_get_method(fm_t * t);

double *fm_get_proba(fm_t * t);

double *fm_get_time(fm_t * t);

int fm_get_len_method(fm_t * t);

int fm_get_len_proba(fm_t * t);

int fm_get_len_time(fm_t * t);

int fm_get_len_p_min(fm_t * t);

void fm_set_method(fm_t * t, unsigned long *value, int len);

void fm_set_proba(fm_t * t, double *value, int len, int len_p_min);

void fm_set_time(fm_t * t, double *value, int len);

void fm_put_zero(fm_t * t);

bool fm_is_zero(fm_t * t);

fm_t *fm_copy(fm_t * t);

int fm_is_equal(fm_t * c1, fm_t * c2);

int fm_print(fm_t *);

int fm_fprint(FILE * output_file, fm_t * t);

#endif				/* FM */
