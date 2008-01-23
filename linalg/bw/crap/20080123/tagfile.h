#ifndef TAGFILE_H_
#define TAGFILE_H_

#ifdef	__cplusplus
extern "C" {
#endif

extern int read_tag_file(void);
extern int write_tag_file(char *);
extern void load_x_vectors(coord_t *, int, int);

#ifdef	__cplusplus
}
#endif

#endif	/* TAGFILE_H_ */
