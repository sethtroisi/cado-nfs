#include <cstdio>
#include <string>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#define DIE_ERRNO_DIAG(tst, func, arg) do {				\
    if ((tst)) {					        	\
        fprintf(stderr, func "(%s): %s\n", arg, strerror(errno));       \
        exit(1);					        	\
    }							        	\
} while (0)

#define LBSIZE  4096

using namespace std;

int main(int argc, char * argv[])
{
    char * perm = NULL;
    char * file = NULL;
    unsigned int nlines = 0;
    vector<long> lsizes;


    argv++,argc--;
    for( ; argc ; ) {
        /*
        if (strcmp(argv[0], "--reverse") == 0) {
            reverse = 1;
        }
        */
        if (!perm) { perm = argv[0]; argv++,argc--; }
        if (!file) { file = argv[0]; argv++,argc--; }
    }


    /* read the file, count the lines */

    long nbytes;
    {
        long o;
        ifstream ffile(file);
        DIE_ERRNO_DIAG(!ffile.is_open(), "open", file);
        for( ;  ; ) {
            string ss;
            o = ffile.tellg();
            lsizes.push_back(o);
            if (!getline(ffile,ss))
                break;
            nlines++;
        }
        nbytes = o;
        ffile.close();
    }

    fprintf(stderr, "Read %u lines, total %ld bytes\n", nlines, nbytes);
    for(unsigned int i = 0 ; i < lsizes.size() - 1 ; i++) {
        lsizes[i] = lsizes[i+1] - lsizes[i];
    }
    lsizes.pop_back();

    vector<unsigned int> p;
    {
        ifstream pfile(perm);
        unsigned int v;
        for( ; pfile >> v ; ) {
            p.push_back(v);
        }
        pfile.close();
    }

    vector<long> loffsets_end(nlines,0);
    {
        for(unsigned int i = 0 ; i < nlines ; i++) {
            assert(p[i] < nlines);
            loffsets_end[p[i]] = lsizes[i];
        }
        long o = 0;
        for(unsigned int i = 0 ; i < nlines ; i++) {
            long no = loffsets_end[i] + o;
            loffsets_end[i] = o;
            o = no;
        }
    }

    char * data = (char *) malloc(nbytes);
    {
        ifstream ffile(file);
        DIE_ERRNO_DIAG(!ffile.is_open(), "open", file);
        for(unsigned int i = 0 ; i < nlines; i++) {
            assert(p[i] < nlines);
            assert(loffsets_end[p[i]] + lsizes[i] <= nbytes);
            ffile.read(data + loffsets_end[p[i]], lsizes[i]);
        }
        ffile.close();
    }
    fwrite(data,1,nbytes,stdout);



    return 0;
}
