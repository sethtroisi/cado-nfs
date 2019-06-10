#include "cado.h"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <pthread.h>
#include <hwloc.h>
#include <mutex>
#include <atomic>
#include <omp.h>
#include "utils.h"
#include "ringbuf.h"

void mf_scan2_decl_usage(cxx_param_list & pl)
{
    param_list_usage_header(pl,
            "This program make one reading pass through a binary matrix, and produces\n"
            "the companion .rw and .cw files.\n"
            "Typical usage:\n"
            "\tmf_scan [<matrix file name> | options...]\n"
            );
    param_list_decl_usage(pl, "withcoeffs", "Handle DLP matrix, with coefficients\n");
    param_list_decl_usage(pl, "mfile", "Input matrix name (free form also accepted)");
    param_list_decl_usage(pl, "rwfile", "Name of the row weight file to write (defaults to auto-determine from matrix name)");
    param_list_decl_usage(pl, "cwfile", "Name of the col weight file to write (defaults to auto-determine from matrix name)");
    param_list_decl_usage(pl, "threads", "Number of threads to use (defaults to auto detect\n");
    param_list_decl_usage(pl, "memory", "Amount of RAM to use for rolling buffer memory");
}

size_t cw_cutoff = 1UL << 16;
std::atomic<std::atomic<uint32_t> *> cw_global_pools[64];
std::mutex mutexes[64];
size_t produced;
std::atomic<size_t> consumed;

#if 0
uint32_t quarter_cutoffs[][4] = {
    /*  0 */ { 1, 1, 1, },
    /*  1 */ { 2, 2, 3, },
    /*  2 */ { 4, 5, 6, },
    /*  3 */ { 9, 11, 13, },
    /*  4 */ { 19, 22, 26, },
    /*  5 */ { 38, 45, 53, },
    /*  6 */ { 76, 90, 107, },
    /*  7 */ { 152, 181, 215, },
    /*  8 */ { 304, 362, 430, },
    /*  9 */ { 608, 724, 861, },
    /* 10 */ { 1217, 1448, 1722, },
    /* 11 */ { 2435, 2896, 3444, },
    /* 12 */ { 4870, 5792, 6888, },
    /* 13 */ { 9741, 11585, 13777, },
    /* 14 */ { 19483, 23170, 27554, },
    /* 15 */ { 38967, 46340, 55108, },
    /* 16 */ { 77935, 92681, 110217, },
    /* 17 */ { 155871, 185363, 220435, },
    /* 18 */ { 311743, 370727, 440871, },
    /* 19 */ { 623487, 741455, 881743, },
    /* 20 */ { 1246974, 1482910, 1763487, },
    /* 21 */ { 2493948, 2965820, 3526975, },
    /* 22 */ { 4987896, 5931641, 7053950, },
    /* 23 */ { 9975792, 11863283, 14107900, },
    /* 24 */ { 19951584, 23726566, 28215801, },
    /* 25 */ { 39903169, 47453132, 56431603, },
    /* 26 */ { 79806338, 94906265, 112863206, },
    /* 27 */ { 159612677, 189812531, 225726412, },
    /* 28 */ { 319225354, 379625062, 451452825, },
    /* 29 */ { 638450708, 759250124, 902905650, },
    /* 30 */ { 1276901416, 1518500249, 1805811301, },
    /* 31 */ { 2553802833, 3037000499, 3611622602, },
};
#endif

inline int get_segment_index(uint32_t c)
{
    unsigned int t = 64 - cado_clz64((uint64_t) c);
#if 1
    return t;
#else
    int j = 0;
    for( ; j < 3 ; j++)
        if (c < quarter_cutoffs[t][j]) break;
    return 4*t+j;
#endif
}
inline uint32_t get_segment_offset(int t)
{
#if 1
    return 1UL << (t-1);
#else
    int t0 = t / 4;
    int j = t % 4;
    uint32_t c0 = 1UL << (t0-1);
    if (j)
        c0 = quarter_cutoffs[t0][j-1];
    return c0;
#endif
}

inline uint32_t get_segment_size(int t)
{
#if 1
    return 1UL << (t-1);
#else
    int t0 = t / 4;
    uint32_t s0 = 1UL << (t0-1);
    int j = t % 4;
    uint32_t c0 = 1UL << (t0-1);
    if (j)
        c0 = quarter_cutoffs[t0][j-1];
    return s0 - (c0 - s0);
#endif
}


ringbuf R;

struct cw_thread {
    std::vector<uint32_t> cw;
    uint32_t colmax=0;
    cw_thread() : cw(cw_cutoff, 0) {};
    void loop() {
        char buffer[1 << 14];
        for(size_t s ; (s = ringbuf_get(R, buffer, sizeof(buffer))) != 0 ; ) {
            consumed += s;
            uint32_t * v = (uint32_t *) buffer;
            ASSERT_ALWAYS(s % sizeof(uint32_t) == 0);
            size_t sv = s / sizeof(uint32_t);
            for(size_t i = 0 ; i < sv ; i++) {
                uint32_t c = v[i];
                colmax = MAX(colmax, c+1);
                if (c < cw_cutoff) {
                    cw[c]++;
                } else {
                    /* Get the bit size */
                    unsigned int t = get_segment_index(c);
                    std::atomic<uint32_t> * x = cw_global_pools[t].load();
                    uint32_t c1 = c-get_segment_offset(t);
                    ASSERT_ALWAYS(c1 < get_segment_size(t));
                    if (!x) {
                        std::lock_guard<std::mutex> dummy(mutexes[t]);
                        x = cw_global_pools[t].load();
                        if (!x) {
                            std::atomic<uint32_t> * xx = new std::atomic<uint32_t>[get_segment_size(t)];
                            for(size_t i = 0 ; i < get_segment_size(t) ; i++) {
                                xx[i]=0;
                            }
                            cw_global_pools[t] = xx;
                            x = xx;
                        }
                    }
                    x[c1]++;
                }
            }
        }
    }
};

int main(int argc, char * argv[])
{
    char * argv0 = argv[0];

    cxx_param_list pl;
    const char * rwfile = NULL;
    const char * cwfile = NULL;
    const char * mfile = NULL;

    unsigned int wild =  0;
    int withcoeffs = 0;

    argv++,argc--;

    mf_scan2_decl_usage(pl);

    param_list_configure_switch(pl, "--withcoeffs", &withcoeffs);

    for(;argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild == 0) {
            mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        cwfile = tmp;
    }

    if (!mfile) {
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    if (strlen(mfile) < 4 || strcmp(mfile + strlen(mfile) - 4, ".bin") != 0) {
        fprintf(stderr, "Warning: matrix file name should end in .bin\n");
    }

    if (!rwfile) {
        char * leakme;
        rwfile = leakme = derived_filename(mfile, "rw", ".bin");
    }

    if (!cwfile) {
        char * leakme;
        cwfile = leakme = derived_filename(mfile, "cw", ".bin");
    }

    if (withcoeffs) abort();    // not implemented yet

    /* Detect hardware */
    hwloc_topology_t topology;
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);
    int depth = hwloc_topology_get_depth(topology);
    int npu = hwloc_get_nbobjs_by_depth(topology, depth-1);
    hwloc_obj_t root = hwloc_get_root_obj(topology);
#if HWLOC_API_VERSION < 0x020000
    uint64_t ram = root->memory.total_memory;
#else
    uint64_t ram = root->total_memory;
#endif

    for(uint64_t x = ram >> 4; x ; x >>= 1) ram |= x;
    ram = ram + 1;

    int threads = npu;
    size_t ringbuf_size = ram / 4;

    param_list_parse_int(pl, "threads", &threads);
    param_list_parse_size_t(pl, "memory", &ringbuf_size);

    ringbuf_init(R, ringbuf_size);

    /* Start with the input */

    FILE * f_in = fopen(mfile, "rb");
    if (f_in == NULL) { perror(mfile); exit(EXIT_FAILURE); }
    FILE * f_rw = fopen(rwfile, "wb");
    if (f_rw == NULL) { perror(rwfile); exit(EXIT_FAILURE); }
    FILE * f_cw = fopen(cwfile, "wb");
    if (f_cw == NULL) { perror(cwfile); exit(EXIT_FAILURE); }

    ASSERT_ALWAYS(threads >= 2);

    int consumers = threads-1;
    cw_thread T[consumers];
    
    double t0 = wct_seconds();
    double last_report = t0;

    omp_set_num_threads(threads);
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        if (t == 0) {
            uint32_t buf[256];
            for( ; ; ) {
                uint32_t row_length;
                int rc = fread32_little(&row_length, 1, f_in);
                if (rc != 1)
                    break;
                rc = fwrite32_little(&row_length, 1, f_rw);
                ASSERT_ALWAYS(rc == 1);
                for( ; row_length ; ) {
                    int s = MIN(row_length, 256);
                    int k = fread32_little(buf, s, f_in);
                    ASSERT_ALWAYS(k == s);
                    ringbuf_put(R, (char *) buf, s * sizeof(uint32_t));
                    produced += s * sizeof(uint32_t);
                    row_length -= s;
                }
                double tt = wct_seconds();
                if (tt > last_report + 1) {
                    char buf1[20];
                    char buf2[20];
                    printf("read %s, parsed %s, in %.1f s\n",
                            size_disp(produced, buf1),
                            size_disp(consumed.load(), buf2),
                            (last_report = tt) - t0);
                }
            }
            ringbuf_mark_done(R);
        } else {
            T[t-1].loop();
        }
    }
#pragma omp barrier
    {
        char buf2[20];
        double tt = wct_seconds();
        printf("parsed %s, in %.1f s\n",
                size_disp(consumed.load(), buf2),
                tt - t0);
    }
    for(int i = 1 ; i < consumers ; i++) {
        for(size_t j = 0 ; j < cw_cutoff ; j++)
            T[0].cw[j] += T[i].cw[j];
        T[0].colmax = MAX(T[0].colmax, T[i].colmax);
    }
    uint32_t colmax = T[0].colmax;
    uint32_t c = 0;
    for( ; c < cw_cutoff && c < colmax ; c++) {
        int rc = fwrite32_little(&T[0].cw[c], 1, f_cw);
        ASSERT_ALWAYS(rc == 1);
    }
    for( ; c < colmax ; ) {
        unsigned int t = get_segment_index(c);
        uint32_t c1 = c-get_segment_offset(t);
        uint32_t max1 = MIN(colmax-get_segment_offset(t), get_segment_size(t));
        uint32_t n1 = max1 - c1;
        std::atomic<uint32_t> * x = cw_global_pools[t].load();
        ASSERT_ALWAYS(x);
        ASSERT_ALWAYS(sizeof(*x) == sizeof(uint32_t));
        int rc = fwrite32_little((uint32_t*) &x[c1], n1, f_cw);
        ASSERT_ALWAYS(rc == (int) n1);
        c += n1;
    }

    ringbuf_clear(R);

    hwloc_topology_destroy(topology);

    fclose(f_rw);
    fclose(f_cw);
    fclose(f_in);
}

