#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <errno.h>

#include "master_common.h"
#include "gmp-hacks.h"
#include "auxfuncs.h"

const char *a_meta_filename = "A-%02d-%02d";
const char *valu_meta_filename = "val-%02d";
const char *f_meta_filename = "F%02d";
const char *pi_meta_filename = "pi-%d-%d";

void bw_commit_f(bw_nbpoly fpoly, const int * delta)
{
	int i;
	int j;
	int t;
	char filename[FILENAME_LENGTH];
	FILE *f;
	mpz_t blah;

	mpz_init(blah);

	for (j = 0; j < bigdim; j++) {
		sprintf(filename, f_meta_filename, j);
		f = fopen(filename, "w");
		if (f == NULL) {
			perror("writing f");
			eternal_sleep();
		}
		printf("Writing %s\n", filename);
		for (t = 0; t <= delta[j]; t++) {
			for (i = 0; i < n_param; i++) {
				MPZ_SET_MPN(blah,
					    nbmat_scal(nbpoly_coeff (fpoly, t),
						       i, j), bw_allocsize);
				gmp_fprintf(f, "%Zd%c", blah,
					    i == (n_param - 1) ? '\n' : ' ');
			}
			/* Je les écris dans le sens inverse. Pour une
			 * fois... En fait, c'est plus cohérent pour
			 * l'utilisation qui suit, vu que l'évaluation
			 * bête-et-conne n'est pas moins rapide que par
			 * Horner */
		}
		fclose(f);

#if 0
		for (t = 0;; t++) {
			if (!ncol_is_zero(nbmat_col(nbpoly_coeff(fpoly,
								 deg - t),
						    j)))
				break;
		}

		/* t is now the valuation of the (reversed) column. */

		sprintf(filename, valu_meta_filename, j);
		f = fopen(filename, "w");
		if (f == NULL) {
			perror("Writing valuation file");
			eternal_sleep();
		}
		printf("Writing %s\n", filename);
		fprintf(f,
			"COLUMN %d VALUATION %d DEGNOM %d HAPPY %d TIMES\n",
			j, t, deg, chance_list[j]);
		fclose(f);
#endif
	}

	mpz_clear(blah);
}

void print_chance_list(unsigned int t, unsigned int *chance_list)
{
	/*
	   printf("#RESULT T=%d\n", it->t);
	   for(j=0;j<bigdim;j++) {
	   if (chance_list[j]) printf("#RESULT J=%d\n", j);
	   }
	 */

	int j;
	printf("// step %d LOOK [", t);
	for (j = 0; j < bigdim; j++) {
		if (chance_list[j])
			printf(" %d", j);
	}
	printf(" ]\n");
}


int read_data_for_series(bw_mnpoly a)
{
	int i, j, k;
	FILE **a_files;
	int * a_nbys;
	char filename[FILENAME_LENGTH];
	char row[1024];
	mpz_t blah;

	mpz_init(blah);

	a_files = malloc(m_param * n_param * sizeof(FILE *));
	a_nbys = malloc(n_param * sizeof(int));

	for (i = 0; i < m_param; i++) {
		for (j = 0; j < n_param;) {
			int nbys = 0;
			char * ptr;
			int d;

			sprintf(filename, a_meta_filename, i, j);
			a_files[i * n_param + j] = fopen(filename, "r");
			if (a_files[i * n_param + j] == NULL) {
				die("fopen(%s) : %s", errno, filename,
				    strerror(errno));
			}
			printf("Reading file %s", filename);

			/* NOTE : we drop the first coefficient, because
			 * we mean to work with the sequence generated
			 * on x and By, so that we obtain a generator
			 * afterwards.
			 */

			/* We also use this in order to read the number
			 * of coefficients per row */
			fgets(row, sizeof(row), a_files[i * n_param + j]);
			ptr = row;
			for(nbys = 0 ; gmp_sscanf(ptr,"%Zd%n",blah,&d) >= 1 ; ptr += d, nbys++);
			if (nbys == 0) {
				fprintf(stderr, "\nproblem while reading %s, line = %s\n", filename, row);
				abort();
			}
			if (nbys > 1) {
				printf(" [ %d values per row ]", nbys);
			}
			printf("\n");
			a_nbys[j] = nbys;
			j += nbys;
		}
	}

	for (k = 0; k <= total_work; k++) {
		int rc = 0;

		for (i = 0; i < m_param; i++)
			for (j = 0; j < n_param; j+=a_nbys[j]) {
				int l;
				for(l = 0 ; l < a_nbys[j] ; l++) {
					rc +=
						1 == gmp_fscanf(a_files[i * n_param + j],
								"%Zd", blah);

					mpz_fdiv_r(blah, blah, modulus);
					MPN_SET_MPZ(mnmat_scal
							(mnpoly_coeff(a, k), i, j + l),
							bw_allocsize, blah);
				}
			}
		if (rc == 0) {
			break;
		} else if (rc != m_param * n_param) {
			die("problem with A files (not in sync, rc=%d)\n", 1,
			    rc);
		}
	}

	for (i = 0; i < m_param; i++)
		for (j = 0; j < n_param; j+= a_nbys[j]) {
			fclose(a_files[i * n_param + j]);
		}
	free(a_files);
	free(a_nbys);

	printf("Stopped after reading %d coefficients (not counting 1st)\n",
	       k);

	/* XXX k or k-1 */
	return k;
}
void read_mat_file_header(const char *name)
{
	FILE *f = fopen(name, "r");
	if (f == NULL) {
		die("fopen(%s): %s", errno, name, strerror(errno));
	}
	char modulus[512];
	coord_t ncols;

	fscanf(f, "// %d ROWS %d COLUMNS ; MODULUS %s",
	       &nrows, &ncols, modulus);
	bw_read_modulus_info(modulus, 1);

	if (nrows != ncols) {
		die("matrix is not square", 1);
	}


	fclose(f);
}
