#include <iostream>

#include "args.h"
#include "common_arguments.hpp"
#include "prep_arguments.hpp"
#include "arguments.hpp"
#include "config_file.hpp"
#include "select_mpi.h"
#include "fmt.hpp"
#include "random_generation.h"

using namespace std;

common_arguments common;
prep_arguments mine;
config_file_t cf;

unsigned int m, n;
int can_print;

mpz_class p;

const char * params_file = "bw.cfg";

int argparse(int * p_argc, char *** p_argv)
{
    using namespace std;

    ios_base::sync_with_stdio(false);
    cerr.tie(&cout);
    cout.rdbuf()->pubsetbuf(0, 0);
    cerr.rdbuf()->pubsetbuf(0, 0);

    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    can_print = rank == 0 || getenv("CAN_PRINT");

    if (can_print) {
        std::cout << "// This is bw-prep, version " REV << std::endl;
    }
    process_arguments(*p_argc, *p_argv, common, mine, can_print);

    /* {{{ m,n,p are defined on the command line, and saved in the config
     * file later on.  */

    m = mine.m;
    n = mine.n;
    p = mine.p;

    set_config_value(cf, "m", m);
    set_config_value(cf, "n", n);
    set_config_value(cf, "p", p);
    set_config_value(cf, "check", mine.check_interval);

    if (m == 0 || n == 0) {
        cerr << "m and n must be set\n";
        exit(1);
    }

    if (can_print) {
        cout << fmt("// m = %\n") % m;
        cout << fmt("// n = %\n") % n;
        cout << fmt("// p = %\n") % p;
    }
    /*}}} */
    setup_seeding(mine.seed);

    set_config_value(cf, "mpijobs", size);

    return 0;
}

int finish()
{
    if (can_print) {
        cout << "// saving config file" << endl;
        write_config_file(cf, params_file);

        cout << "// done." << endl;
    }

    return 0;
}
