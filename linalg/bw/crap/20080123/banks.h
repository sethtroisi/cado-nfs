#ifndef BANKS_H_
#define BANKS_H_

#ifdef	__cplusplus
extern "C" {
#endif

extern char modulusname[MODULUSNAME_MAX_LENGTH];

int try_to_write_bank(FILE **, FILE **);
int try_to_load_bank(FILE **, FILE **);

#ifdef	__cplusplus
}
#endif

#endif /* BANKS_H_ */
