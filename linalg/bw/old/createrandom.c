#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"

#include <assert.h>
#include <gmp.h>

#include "macros.h"
#include "types.h"
#include "endian.h"
#include "auxfuncs.h"
#include "variables.h"
#include "filenames.h"
#include "relation.h"
#include "banks.h"
#include "bw_scalar.h"

#define BIGGEST_VALUE	256
/* Used if the modulus is large */

void showuse(void)
{
	fprintf(stderr,"Usage : create <bank> <size> <modulus> <density>\n");
	exit(1);
}

int compar(const void * pa, const void * pb)
{
	unsigned int a, b;
	a = (*(unsigned int*)pa);
	b = (*(unsigned int*)pb);
	if (a > b) return 1;
	else if (a < b) return -1;
	else return 0;
}

int main(int argc, char *argv[])
{
	unsigned int	banknum, size, density;
	unsigned int	i,j,k;
	char	      * modulus_str;	
	FILE	      * fm, * fv;
	ext_relation	r;
	unsigned int  * places;
	bw_scalar	lhs;

	if (argc < 5)
		showuse();
	
	banknum = atoi(argv[1]);
	size    = atoi(argv[2]);
	density = atoi(argv[3]);
	modulus_str  = argv[4];

	if (size <=1 || density <=1) {
		fprintf(stderr,"This size/density is not allowed.\n");
		showuse();
	}

	set_all_filenames(banknum);
	if (try_to_write_bank(&fm,&fv) < 0) {
		showuse();
	}

	if (bw_read_modulus_info(modulus_str,1) < 0 || mpz_size(modulus) == 0) {
		showuse();
	}

	printf("ID STRING : WIEDEMANN INSTANCE %u x %u OVER Z/(%s)Z\n",
			size, size, modulus_str);

	printf("Creating matrix file..."); fflush(stdout);
	fprintf(fm, "WIEDEMANN INSTANCE %u x %u OVER Z/(%s)Z\n",
			size, size, modulus_str);

	places = my_malloc(density * sizeof(unsigned int));

	for(i=0;i<size;i++) {
		int local_weight;

		for(j=0;j<density;j++) {
			places[j]=random()%size;
		}
		qsort(places, density, sizeof(unsigned int), &compar);

		for(j=k=local_weight=0;j<density;j++) {
			if (j==k) {
				local_weight++;
				continue;
			}
			if (places[j] == places[k]) {
				places[j]=(unsigned int)-1;
				continue;
			}
			local_weight++;
			k=j;
		}

		ext_relation_alloc(r,local_weight);

#define RSIGN()	((random()&1)?(1):(-1))

		for(j=k=0;j<density;j++) {
			for( ; j < density && places[j]==(unsigned)-1 ; j++);
			if (j == density) break;
			r->SP[k].index = places[j];
			if (mpz_size(modulus) == 1) {
		/*		r->SP[k].expo  = RSIGN() *
					(random() % mpz_get_ui(modulus));
					*/
				r->SP[k].expo  = RSIGN() * (random() % BIGGEST_VALUE);
			} else {
				r->SP[k].expo  = RSIGN() * (random() % BIGGEST_VALUE);
			}
			if (r->SP[k].expo == 0)
				r->SP[k].expo = 1;
			k++;
		}

#undef RSIGN

		ext_relation_write(fm,r);
		ext_relation_free(r);
	}
		
	printf("ok\n"); fflush(stdout);

	fclose(fm);

	printf("Creating vector file..."); fflush(stdout);
	fprintf(fv, "WIEDEMANN INSTANCE %u x %u OVER Z/(%s)Z\n",
			size, size, modulus_str);

	for(i=0;i<size;i++) { 
		bw_scalar_alloc(lhs, BW_SCALAR_SHORT);
		bw_scalar_set_zero(lhs);
		if (mpz_size(modulus) == 1) {
			/* *lhs = (random() % mpz_get_ui(modulus)); */
			*lhs = random() % BIGGEST_VALUE;
		} else {
			*lhs = random() % BIGGEST_VALUE;
		}
		bw_scalar_write(fv,lhs);
		bw_scalar_free(lhs);
	}
	printf("ok\n"); fflush(stdout);

	fclose(fv);

	return 0;
}

/* vim:set sw=8: */

