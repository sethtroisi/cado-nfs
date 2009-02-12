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
#define CHECK_FILE_BASE "C"
#define COMMON_VECTOR_ITERATE_PATTERN "%s.%u" TWISTED_EXTENSION

/* misc */
#define BW_CONFIG_FILE "bw.cfg"
#define MATRIX_INFO_FILE_PATTERN "%s.info"


#endif	/* FILENAMES_H_ */
