#ifndef FILENAMES_H_
#define FILENAMES_H_

/* Here, all filenames and filename patterns are defined as macros --
 * just for the purpose of having all of them gathered in one place.
 */

#define A_FILE_PATTERN "A%u-%u.%u-%u"

#define TWISTED_EXTENSION ".twisted"

#define X_FILE_BASE "X"
#define X_TWISTED_FILE X_FILE_BASE TWISTED_EXTENSION

/* different kinds of vector iterates */
#define Y_FILE_BASE "Y"
#define V_FILE_BASE_PATTERN "V%u-%u"
#define S_FILE_BASE_PATTERN "S%u-%u"
#define K_FILE_BASE_PATTERN "K"
#define CHECK_FILE_BASE "C"
#define COMMON_VECTOR_ITERATE_PATTERN "%s.%u" TWISTED_EXTENSION

#define W_FILE_BASE "W"

#define W_FILE W_FILE_BASE TWISTED_EXTENSION

/* This shows that the two-level expansion mechanism just can't work in some
 * cases :-( */
#define S_FILE_PATTERN S_FILE_BASE_PATTERN ".%u" TWISTED_EXTENSION
#define K_FILE_PATTERN K_FILE_BASE_PATTERN ".%u" TWISTED_EXTENSION

/* lingen */
#define LINGEN_BOOTSTRAP_FILE   "F_INIT_QUICK"
#define LINGEN_F_FILE   "F"
#define LINGEN_PI_PATTERN       "pi-%u-%u"

/* misc */
#define BW_CONFIG_FILE "bw.cfg"

/* mksol */
#define F_FILE_SLICE_PATTERN "F%u-%u"

#endif	/* FILENAMES_H_ */
