#ifndef FILENAMES_H_
#define FILENAMES_H_

#include <inttypes.h>

/* Here, all filenames and filename patterns are defined as macros --
 * just for the purpose of having all of them gathered in one place.
 */

#define A_FILE_PATTERN "A%u-%u.%u-%u"

/* different kinds of vector iterates */
#define Y_FILE_BASE "Y"
#define X_FILE_BASE_PATTERN "X"
#define V_FILE_BASE_PATTERN "V%u-%u"
#define S_FILE_BASE_PATTERN "S.sols%u-%u.%u-%u"
#define K_FILE_BASE_PATTERN "K.sols%u-%u"
#define R_FILE_BASE_PATTERN "R.sols%u-%u"       /* short-lived */
#define CHECK_FILE_BASE "C"
#define W_FILE "W"


/* lingen */
#define LINGEN_BOOTSTRAP_FILE   "F_INIT_QUICK"
#define LINGEN_F_FILE   "F"
#define LINGEN_PI_PATTERN       "pi-%u-%u"

/* mksol */
#define F_FILE_SLICE_PATTERN "F%u-%u"
#define F_FILE_SLICE_PATTERN2 "F.sols%u-%u.%u-%u"

#endif	/* FILENAMES_H_ */
