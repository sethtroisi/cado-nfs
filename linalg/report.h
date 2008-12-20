// data structure for reporting actions during the merge; in standard mode
// (say mono proc), this is just a wrap around for a FILE; otherwise,
// it can be used to register things in an array that will be examined and
// flushed from time to time. See the MPI version for more.
typedef struct{
    char type;
    // '0' for the standard stuff
    FILE *outfile;
    // '1' for MPI
    INT **history;
    int mark;
    int bufsize; // says it!
} report_t;

extern void init_rep(report_t *rep, char *outname, sparse_mat_t *mat, int type, int bufsize);
extern void report1(report_t *rep, INT i);
extern void report2(report_t *rep, INT i1, INT i2);
extern void reportn(report_t *rep, INT *ind, int n);

