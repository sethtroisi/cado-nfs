#ifndef PARAMS_H_
#define PARAMS_H_

/* Here come all global parameters, passed through the C preprocessor,
 * which alter the code produced */

/* Maximum acceptable length of the ID string present at the beginning of
 * the input files */
#define DATA_IDENT_LEN 80

/* Maximum acceptable length of the name of the modulus. Although this
 * constant is defined to be quite large, there would be no problem in
 * having it restricted to a dozen. In fact, it can be a symbolic name
 * like ``p1''. The program then fetches information about p1 in the
 * file ``input/p1''. Things like known factors go there.
 * The _FMT indicates how tolerant we should be on onput. Its output
 * analogue is simply "%s". */
#define MODULUSNAME_MAX_LENGTH 200
#define MODULUSNAME_FMT "%199s"

/* Maximum length for internal filenames and meta-filenames. This includes
 * the directory component (all paths are relative, so it doesn't grow
 * too high) */
#define FILENAME_LENGTH 80

/* Typical buffer size to hold one line of data from a file */
#define LINEBUF_LENGTH	1024

/* Define this to inhibit all ``assert'' directives, when they have all
 * proven to be true. Supplying OFLAGS=max to the makefile does this too.
 */
/* #define NDEBUG */

/* Define this to put a limit to the size of input data. 0 to disable */
#define SPMLINE_MAXIMUM_LENGTH_ON_INPUT 0

/* This defines the alignment of data in files. It must be equal to the
 * size of the widest type available among all the machines the program
 * is intended to run on */
#define RAW_DATA_ALIGN 64
/* #define RAW_DATA_ALIGN 32 */

/* Define this to compile the checksum-relevant parts of the code. For
 * now, they are unused so it is better to avoid it. This flag is
 * included for consistency reason with the similar source from the
 * DiscreteLog project */
/* #define CHECKSUM_ENABLED */

/* Define this to enable logging of memory allocations above a certain
 * threshold. 0 means : log everything. The threshold is measured in
 * bytes */
#define MALLOC_LOG_LIMIT	16384

/* Take y as random data instead of z */
#define CORRECT_STUPID_WOE


/* Define this is the pthread_getcpuclockid() and clock_gettime() calls
 * are available and functional (the program should be linked with -lrt
 * then).
 *
 * Pointless when we are not concerned with multithreaded code.
 */
#define	HAS_REALTIME_CLOCKS

/* Define if using LinuxThreads */
/* #define	LINUXTHREADS */


#endif /* PARAMS_H_ */
