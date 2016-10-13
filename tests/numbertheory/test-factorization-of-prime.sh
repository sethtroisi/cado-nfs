#include "ant.hpp"

int main(int argc, char *argv[])
{				/*{{{ */
    /*
       unsigned int m = 8;
       unsigned int n = 5;
       if (argc == 3) {
       m = strtoul(argv[1], NULL, 0);
       n = strtoul(argv[2], NULL, 0);
       }
       gmp_randstate_t state;
       gmp_randinit_default(state);

       mpq_mat M;
       mpq_mat T;
       mpz_mat Mz;
       mpz_mat Tz;
       mpz_t p;

       mpq_mat_init(M, m, n);
       mpq_mat_init(T, m, m);

       mpz_mat_init(Mz, m, n);
       mpz_mat_init(Tz, m, m);

       mpz_init(p);

       mpz_set_ui(p, 19);

       if (0) {
       printf("\n\nCas 0.1\n\n");
       mpq_mat_urandomm(M, state, p);
       mpq_mat_fprint(stdout, M);
       printf("\n");
       mpq_mat_gauss_backend(M, T);
       mpq_mat_fprint(stdout, M);
       printf("\n");
       mpq_mat_fprint(stdout, T);
       printf("\n");
       }

       if (0) {
       printf("\n\nCas 0.2\n\n");
       mpz_mat_urandomm(Mz, state, p);
       mpz_mat_fprint(stdout, Mz);
       printf("\n");
       mpz_mat_gauss_backend_mod(Mz, Tz, p);
       mpz_mat_fprint(stdout, Mz);
       printf("\n");
       mpz_mat_fprint(stdout, Tz);
       printf("\n");
       }

       if (1) {
       printf("\n\nCas 1\n\n");
       mpz_mat_realloc(Mz, m, n);
       mpz_mat_urandomm(Mz, state, p);
       mpz_mat_fprint(stdout, Mz); printf("\n");
       double t = seconds();
       mpz_mat_hnf_backend(Mz, Tz);
       t = seconds()-t;
       mpz_mat_fprint(stdout, Mz); printf("\n");
       mpz_mat_fprint(stdout, Tz); printf("\n");

       printf("%1.4f\n", t);
       }

       mpz_clear(p);
       mpq_mat_clear(M);
       mpq_mat_clear(T);
       mpz_mat_clear(Mz);
       mpz_mat_clear(Tz);
       gmp_randclear(state);
     */

    // The inputs to this problem are f, one polynomial of degree n, and B, the matrix containing the genereators of one order of the number field obtained with f, as well as p, a prime number

    unsigned long seed = clock();

    for( ; argc > 3 ; ) {
        if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "--seed") == 0) {
            seed = atoi(argv[2]);
            argc--,argv++;
            argc--,argv++;
            continue;
        }
        fprintf(stderr, "Usage: ./a.out [options] [filename] [p]\n");
        fprintf(stderr, "Unexpected arg: %s\n", argv[1]);
	exit(EXIT_FAILURE);
    }


    if (argc != 3) {
	fprintf(stderr, "Usage: ./a.out [options] [filename] [p]\n");
	exit(EXIT_FAILURE);
    }

    unsigned int p = strtoul(argv[2], NULL, 0);	//19; //atoi(argv[1]);
    FILE *problemfile = fopen(argv[1], "r");

    if (!problemfile) {
        fprintf(stderr, "%s: %s\n", argv[1], strerror(errno));
        exit(EXIT_FAILURE);
    }

    mpq_mat D;
    mpz_poly f;

    printf("Format: [degree] [coeffs] [coeffs of order basis]\n");

    unsigned int n = 0;
    mpz_poly_init(f, n);
    read_data(&n, f, /*gen,*/ problemfile);	
    fclose(problemfile);
    mpq_mat_init(D, n, n);

    mpz_poly g;
    mpz_poly_init(g, n);
    mpz_poly_to_monic(g, f);
    printf("f  is : ");
    mpz_poly_fprintf(stdout, f);
    //printf("\n");
    //printf("f^ is : ");
    //mpz_poly_fprintf(stdout, g);
    //printf("\n");

    
    
    
    //FILE* f5 = fopen("maxOrderCTime","a+");
    //FILE* f6 = fopen("idealCTime","a+");
    clock_t start, end;
    double cpu_time_used;
     
    start = clock();
    p_maximal_order(D, f, p);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    // Here we print the maximal order in a file
    
    printf("\n");    
    FILE *maxOrderFile = NULL;
    maxOrderFile = fopen("maxOrderC.data","w+");
    mpq_mat_fprint_as_mpz(maxOrderFile, D);
    fclose(maxOrderFile);
    
    
    // Here we print the time for the computation of p-maximal order in a file
    
    FILE *maxOrderTimeFile = NULL;
    maxOrderTimeFile = fopen("maxOrderCTime.data","a+");
    fprintf(maxOrderTimeFile,"%d %f\n",(unsigned int) f->deg, cpu_time_used);
    fclose(maxOrderTimeFile);
    
    
    printf("\n");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
        
    vector<pair<cxx_mpq_mat, int>> ideals;
    
    start = clock();
    factorization_of_prime(ideals, g, p, state);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    sort_matrices(ideals);
    
    // Here we print the time for the computation of p-maximal order in a file
    
    FILE *idealTimeFile = NULL;
    idealTimeFile = fopen("idealCTime.data","a+");
    fprintf(idealTimeFile,"%d %f\n", (unsigned int) f->deg, cpu_time_used);
    fclose(idealTimeFile);
    
    
    // Here we print these ideals with their multiplicities in a file
    
    FILE *idealFile = NULL;
    idealFile = fopen("idealC.data","w+");
    for(unsigned int i = 0 ; i < ideals.size() ; i++){
        //printf("Ideal :\n");
        //mpq_mat_fprint(stdout, ideals[i].first);
        mpq_mat_fprint_as_mpz(idealFile, ideals[i].first);
        //printf("with a multiplicity of %d\n\n",ideals[i].second);
        gmp_fprintf(idealFile, "%d is its multiplicity\n", ideals[i].second);
    }
    fclose(idealFile);
    
    
    //string SBAD = "";
    //string SBADINFO = "";
    //print_comments_for_badideals_above_p(SBAD, SBADINFO, 0, D,g,ideals,p);
    
    gmp_randclear(state);
    
    mpz_poly_clear(g);
    mpz_poly_clear(f);
    mpq_mat_clear(D);



}

/*}}}*/
