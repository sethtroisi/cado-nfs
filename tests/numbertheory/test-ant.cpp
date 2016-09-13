#include "cado.h"
#include "utils.h"
#include "ant.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;
static char ** original_argv;

istream& operator>>(istream& is, cxx_mpz_poly& f)/*{{{*/
{
    int deg;
    if (!(is >> deg) || deg <= 0) {
        is.clear();
        is.setstate(is.rdstate() | ios_base::failbit);
        return is;
    }
    mpz_poly_realloc(f, deg + 1);
    for(int i = 0 ; i <= deg ; i++) {
        if (!(is >> f->coeff[i])) {
            is.clear();
            is.setstate(is.rdstate() | ios_base::failbit);
            return is;
        }
    }
    mpz_poly_cleandeg(f, deg);
    return is;
}
/*}}}*/

static void decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "test", "which test to run");
    param_list_decl_usage(pl, "prime", "prime");
    param_list_decl_usage(pl, "polystr", "polynomial (string)");
    param_list_decl_usage(pl, "polyfile", "polynomial (file)");
    param_list_decl_usage(pl, "out", "output file");
    param_list_decl_usage(pl, "batch", "batch input file with test vectors and expected results");
    param_list_decl_usage(pl, "seed", "seed used for random picks");
}

void usage(param_list_ptr pl, char ** argv, const char * msg = NULL)
{
    param_list_print_usage(pl, argv[0], stderr);
    if (msg) {
        fprintf(stderr, "%s\n", msg);
    }
    exit(EXIT_FAILURE);
}

struct iowrap {
    streambuf * ibuf, * obuf;
private:
    ifstream ifs;
    ofstream ofs;
    istringstream iss;
public:
    iowrap(param_list_ptr pl) {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "polystr")) != NULL) {
            iss.str(tmp);
            ibuf = iss.rdbuf();
        } else if ((tmp = param_list_lookup_string(pl, "polyfile")) != NULL) {
            if (strcmp(tmp, "-") == 0) {
                ibuf = cin.rdbuf();
            } else {
                ifs.open(tmp);
                ibuf = ifs.rdbuf();
            }
        } else {
            usage(pl, original_argv, "Please provide either --polyfile or --polystr");
        }
        if ((tmp = param_list_lookup_string(pl, "out")) != NULL && strcmp(tmp, "-") != 0) {
            ofs.open(tmp);
            obuf = ofs.rdbuf();
        } else {
            obuf = cout.rdbuf();
        }
    }
};

void do_p_maximal_order(param_list_ptr pl) {/*{{{*/
    cxx_mpz_poly f;
    unsigned long p;

    if (!param_list_parse_ulong(pl, "prime", &p)) usage(pl, original_argv, "missing prime argument");

    iowrap io(pl);
    istream in(io.ibuf);
    ostream out(io.obuf);

    if (!(in >> f)) usage(pl, original_argv, "cannot parse polynomial");
    cxx_mpq_mat M = p_maximal_order(f, p);
    cxx_mpz D;
    cxx_mpz_mat A;
    mpq_mat_numden(A, D, M);
    out << "1/" << D << "*\n" << A << endl;
}
/*}}}*/

/* This is over SL_n(Z_p) */
bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A, unsigned long p)
{
    if (M->m != A->m) return false;
    if (M->n != A->n) return false;
    cxx_mpq_mat Mi;
    mpq_mat_inv(Mi, M);
    cxx_mpq_mat AMi;
    mpq_mat_mul(AMi, A, Mi);
    /* check that the p-valuation is zero */
    for(unsigned int i = 0 ; i < AMi->m ; i++) {
        for(unsigned int j = 0 ; j < AMi->n ; j++) {
            mpq_srcptr mij = mpq_mat_entry_const(AMi, i, j);
            if (mpz_divisible_ui_p(mpq_denref(mij), p)) return false;
        }
    }
    return true;
}

bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A)
{
    cxx_mpq_mat Mq(M);
    cxx_mpq_mat Aq(A);
    return sl_equivalent_matrices(Mq, Aq);
}

void do_p_maximal_order_batch(param_list_ptr pl) {/*{{{*/
    cxx_mpz_poly f;
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == NULL)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;
        invalid_argument exc(string("Parse error on input") + s);
        istringstream is0(s);
        if (!(is0 >> f)) usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw exc;
        cxx_mpz d;
        cxx_mpz_mat A(f->deg, f->deg);
        istringstream is1(s);
        unsigned long p;
        if (!(is1 >> p))
            throw exc;
        if (!(is1 >> d))
            throw exc;
        for(unsigned int i = 0 ; i < A->m ; i++) {
            for(unsigned int j = 0 ; j < A->n ; j++) {
                if (!(is1 >> mpz_mat_entry(A, i, j)))
                    throw exc;
            }
        }
        cxx_mpq_mat M;
        mpq_mat_set_mpz_mat_denom(M, A, d);

        cxx_mpq_mat my_M = p_maximal_order(f, p);

        bool ok = sl_equivalent_matrices(M, my_M, p);

        cout << (ok ? "ok" : "NOK") << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << endl;
        test++;
    }
}
/*}}}*/

void do_factorization_of_prime(param_list_ptr pl) {/*{{{*/
    cxx_mpz_poly f;
    unsigned long p;

    if (!param_list_parse_ulong(pl, "prime", &p)) usage(pl, original_argv, "missing prime argument");

    iowrap io(pl);
    istream in(io.ibuf);
    ostream out(io.obuf);

    if (!(in >> f)) usage(pl, original_argv, "cannot parse polynomial");
    cxx_mpq_mat M = p_maximal_order(f, p);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }
    vector<pair<cxx_mpz_mat, int>> F = factorization_of_prime(M, f, p, state);
    gmp_randclear(state);
    out << F.size() << "\n";
    for(unsigned int k = 0 ; k < F.size() ; k++) {
        out << F[k].first << F[k].second << endl;
    }
}
/*}}}*/

void do_factorization_of_prime_batch(param_list_ptr pl) {/*{{{*/
    cxx_mpz_poly f;
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == NULL)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        invalid_argument exc(string("Parse error on input") + s);

        istringstream is0(s);
        if (!(is0 >> f)) usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw exc;

        istringstream is1(s);

        unsigned long p;
        if (!(is1 >> p)) throw exc;

        cxx_mpq_mat M;
        {       /* {{{ read magma's basis of the p-maximal order */
            cxx_mpz_mat A(f->deg, f->deg);
            cxx_mpz d;
            if (!(is1 >> d))
                throw exc;
            for(unsigned int i = 0 ; i < A->m ; i++) {
                for(unsigned int j = 0 ; j < A->n ; j++) {
                    if (!(is1 >> mpz_mat_entry(A, i, j)))
                        throw exc;
                }
            }
            mpq_mat_set_mpz_mat_denom(M, A, d);
        }       // }}}

        vector<pair<cxx_mpz_mat, int> > ideals;
        {       /* {{{ read magma's factorization of p */
            unsigned int nideals;
            if (!(is1 >> nideals))
                throw exc;
            for(unsigned int k = 0 ; k < nideals ; k++) {
                cxx_mpz_mat A(f->deg, f->deg);
                int e;
                for(unsigned int i = 0 ; i < A->m ; i++) {
                    for(unsigned int j = 0 ; j < A->n ; j++) {
                        if (!(is1 >> mpz_mat_entry(A, i, j))) throw exc;
                    }
                }
                if (!(is1 >> e)) throw exc;
                ideals.push_back(make_pair(A, e));
            }
        }       // }}}

        cxx_mpq_mat my_M = p_maximal_order(f, p);
        gmp_randstate_t state;
        gmp_randinit_default(state);

        /* What if we simply ask our code to compute the factorization of
         * p with respect to the order basis which was chosen by magma ?
         * */
        vector<pair<cxx_mpz_mat, int>> my_ideals = factorization_of_prime(M, f, p, state);
        gmp_randclear(state);

        sort(ideals.begin(), ideals.end(), ideal_comparator());
        sort(my_ideals.begin(), my_ideals.end(), ideal_comparator());

        bool ok=(ideals.size() == my_ideals.size());
        for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
            ok = (ideals[k] == my_ideals[k]);
        }
        cout << (ok ? "ok" : "NOK") << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << endl;
        test++;
    }
}
/*}}}*/

int main(int argc, char *argv[]) /*{{{ */
{
    param_list pl;
    param_list_init(pl);

    param_list_configure_alias(pl, "prime", "p");
    param_list_configure_alias(pl, "polystr", "f");
    param_list_configure_alias(pl, "out", "o");

    decl_usage(pl);

    original_argv = argv;

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(pl, original_argv, "unexpected argument");
    }

    const char * tmp = param_list_lookup_string(pl, "test");
    if (!tmp) usage(pl, original_argv, "missing --test argument");

    if (strcmp(tmp, "p-maximal-order") == 0) {
        do_p_maximal_order(pl);
        return 0;
    } else if (strcmp(tmp, "p-maximal-order-batch") == 0) {
        do_p_maximal_order_batch(pl);
        return 0;
    } else if (strcmp(tmp, "factorization-of-prime") == 0) {
        do_factorization_of_prime(pl);
        return 0;
    } else if (strcmp(tmp, "factorization-of-prime-batch") == 0) {
        do_factorization_of_prime_batch(pl);
        return 0;
    } else {
        usage(pl, original_argv, "unknown test");
    }

#if 0
    param_list_clear(pl);
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


#endif
}

/*}}}*/
